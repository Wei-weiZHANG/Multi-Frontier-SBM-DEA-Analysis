# ============================================================================
# Underreporting Sensitivity Analysis: 20% Farms Low-Report Waste_water & Manure
# ============================================================================
# Standalone script: Fixed frontier; random 20% of farms underreport by 10-20%;
# recalc Carbon from modified Manure; compare baseline vs underreported efficiency
# and rank; 100 random draws for distribution; plot LOESS by scale (4 frontiers).
#
# Usage: setwd() to project dir, then source("underreporting_sensitivity.R")
# Or: source("Multi-frontier SBM modules/underreporting_sensitivity.R")
#
# ============================================================================

# ============================================================================
# 0. PACKAGES AND CONFIG
# ============================================================================

required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "lpSolve", "scales", "AER")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Config (override via environment if needed)
xlsx_path   <- Sys.getenv("UNDERREPORT_DATA_PATH", 
  "./swine_farm_data.xlsx")
n_under_rep <- 100L   # number of random "who underreports" draws
pct_farms   <- 0.20   # 20% of farms underreport
u_min       <- 0.10   # underreport factor min (10%)
u_max       <- 0.20   # underreport factor max (20%)
seed        <- 12345L

message("\n================================================================================")
message("UNDERREPORTING SENSITIVITY: Waste_water & Manure low-report 10-20% (20% farms)")
message("================================================================================")
message(sprintf("  Replications: %d  |  Underreporting farms: %.0f%%  |  u in [%.0f%%, %.0f%%]",
                n_under_rep, pct_farms * 100, u_min * 100, u_max * 100))

# Output directory (same base as main analysis or subdir)
results_base <- Sys.getenv("RESULTS_DIR", unset = "results")
out_dir      <- file.path(results_base, "underreporting_sensitivity")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message(sprintf("  Output: %s\n", out_dir))

# ============================================================================
# 1. VARIABLE DEFINITIONS (match multi_frontier_sbm_dea.R)
# ============================================================================

inputs             <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
good_output        <- "Market_pig"
bad_outputs_full   <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission", 
                        "Eutrophication_potential")
constrained_bad_M <- c(3L)
constrained_bad_L <- c(1L, 2L)
constrained_bad_G <- c(4L, 5L)
constrained_bad_A <- c(1L, 2L, 3L, 4L, 5L)

# Weight design (no pre-experiment): exact values from weight configuration table
bad_weights_M <- c(Waste_water = 0.227, Manure = 0.227, Dead_pig = 4.091,
                   Carbon_emission = 0.227, Eutrophication_potential = 0.227)
bad_weights_L <- c(Waste_water = 1.607, Manure = 1.607, Dead_pig = 0.089,
                   Carbon_emission = 0.089, Eutrophication_potential = 0.089)
bad_weights_G <- c(Waste_water = 0.227, Manure = 0.227, Dead_pig = 0.227,
                   Carbon_emission = 4.091, Eutrophication_potential = 4.091)
bad_weights_A <- c(Waste_water = 1, Manure = 1, Dead_pig = 1,
                   Carbon_emission = 1, Eutrophication_potential = 1)
# Ensure order matches bad_outputs_full
bad_weights_M <- bad_weights_M[bad_outputs_full]
bad_weights_L <- bad_weights_L[bad_outputs_full]
bad_weights_G <- bad_weights_G[bad_outputs_full]
bad_weights_A <- bad_weights_A[bad_outputs_full]

# ============================================================================
# 2. SBM AND SCALING (consistent copy from the main analysis)
# ============================================================================

scale_matrix <- function(M) {
  scaling_factors <- apply(M, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
  sweep(M, 2, scaling_factors, "/")
}

sbm_efficiency <- function(X, Y, bad = NULL, constrained_bad_indices = NULL, 
                           bad_weights = NULL, RTS = "vrs", orientation = "output", 
                           show_dmu_progress = TRUE, evaluate_indices = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.null(bad) && !is.matrix(bad)) bad <- as.matrix(bad)
  n <- nrow(X); m <- ncol(X); s_g <- ncol(Y)
  if (!is.null(bad)) { s_b <- ncol(bad) } else { s_b <- 0L; bad <- matrix(0, nrow = n, ncol = 0L) }
  if (is.null(constrained_bad_indices) && s_b > 0) constrained_bad_indices <- seq_len(s_b)
  if (!is.null(bad_weights)) {
    if (length(bad_weights) != s_b) stop("bad_weights length mismatch")
    bad_weights[is.na(bad_weights) | bad_weights < 0] <- 0
  } else { bad_weights <- rep(1, s_b) }
  if (!is.null(evaluate_indices)) {
    evaluate_indices <- sort(unique(as.integer(evaluate_indices)))
    n_eval <- length(evaluate_indices)
  } else { evaluate_indices <- seq_len(n); n_eval <- n }
  eff_scores <- numeric(n)
  show_progress <- show_dmu_progress && (n_eval > 20)
  progress_interval <- if (n_eval > 100) 10 else if (n_eval > 50) 5 else 1
  eval_counter <- 0
  for (k in evaluate_indices) {
    eval_counter <- eval_counter + 1
    if (show_progress && (eval_counter %% progress_interval == 0 || eval_counter == n_eval)) {
      cat(sprintf("\r    SBM: %d/%d DMUs", eval_counter, n_eval))
      if (eval_counter == n_eval) cat("\n")
    }
    x_k <- pmax(X[k, ], 1e-8)
    y_k <- pmax(Y[k, ], 1e-8)
    b_k <- if (s_b > 0) pmax(bad[k, ], 1e-8) else numeric(0)
    num_vars <- 1 + n + m + s_g + s_b
    obj_coef <- numeric(num_vars)
    obj_coef[1] <- 1
    norm_denom <- m + s_b
    for (i in seq_len(m)) obj_coef[1 + n + i] <- -1 / (norm_denom * x_k[i])
    if (s_b > 0) {
      for (i in seq_len(s_b)) {
        w <- if (i %in% constrained_bad_indices && length(bad_weights) >= i) bad_weights[i] else 0
        obj_coef[1 + n + m + s_g + i] <- -w / (norm_denom * b_k[i])
      }
    }
    n_constr <- m + s_g + s_b + 1 + if (RTS == "vrs") 1 else 0
    constr <- matrix(0, nrow = n_constr, ncol = num_vars)
    dir_vec <- character(n_constr); rhs_vec <- numeric(n_constr); row_idx <- 1
    for (i in seq_len(m)) {
      constr[row_idx, 2:(1+n)] <- X[, i]; constr[row_idx, 1 + n + i] <- 1
      constr[row_idx, 1] <- -x_k[i]; dir_vec[row_idx] <- "=="; rhs_vec[row_idx] <- 0; row_idx <- row_idx + 1
    }
    for (j in seq_len(s_g)) {
      constr[row_idx, 2:(1+n)] <- Y[, j]; constr[row_idx, 1 + n + m + j] <- -1
      constr[row_idx, 1] <- -y_k[j]; dir_vec[row_idx] <- "=="; rhs_vec[row_idx] <- 0; row_idx <- row_idx + 1
    }
    if (s_b > 0) {
      cset <- if (length(constrained_bad_indices) > 0) constrained_bad_indices else seq_len(s_b)
      for (j in seq_len(s_b)) {
        constr[row_idx, 2:(1+n)] <- bad[, j]; constr[row_idx, 1 + n + m + s_g + j] <- 1
        constr[row_idx, 1] <- -b_k[j]
        dir_vec[row_idx] <- if (j %in% cset) "==" else "<="
        rhs_vec[row_idx] <- 0; row_idx <- row_idx + 1
      }
    }
    constr[row_idx, 1] <- 1
    for (j in seq_len(s_g)) constr[row_idx, 1 + n + m + j] <- 1 / (s_g * y_k[j])
    dir_vec[row_idx] <- "=="; rhs_vec[row_idx] <- 1; row_idx <- row_idx + 1
    if (RTS == "vrs") {
      constr[row_idx, 2:(1+n)] <- 1; constr[row_idx, 1] <- -1
      dir_vec[row_idx] <- "=="; rhs_vec[row_idx] <- 0
    }
    lp_result <- tryCatch(lpSolve::lp("min", obj_coef, constr, dir_vec, rhs_vec, scale = 1, compute.sens = 0),
                         error = function(e) list(status = 99, objval = NA))
    if (lp_result$status == 0) eff_scores[k] <- pmax(0, pmin(1, lp_result$objval))
    else eff_scores[k] <- NA_real_
  }
  if (!is.null(evaluate_indices) && length(evaluate_indices) < n)
    eff_scores[setdiff(seq_len(n), evaluate_indices)] <- NA_real_
  return(eff_scores)
}

# ============================================================================
# 3. DATA LOAD AND CLEAN
# ============================================================================

if (!file.exists(xlsx_path)) stop("Data file not found: ", xlsx_path)
suppressWarnings(df_raw <- readxl::read_xlsx(xlsx_path, sheet = 3))
all_required <- unique(c(inputs, good_output, bad_outputs_full))
df_clean <- df_raw %>%
  dplyr::select(dplyr::all_of(all_required), 
                dplyr::any_of(c("Farm", "ID", "Report_year", "No.", "No", "Number"))) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_required), ~ suppressWarnings(as.numeric(.x)))) %>%
  tidyr::drop_na(dplyr::all_of(all_required))
for (v in all_required) df_clean[[v]] <- pmax(df_clean[[v]], 1e-8)

n_obs <- nrow(df_clean)
message(sprintf("  Loaded: %d farms\n", n_obs))

# ============================================================================
# 4. BASE SCALING AND FIXED FRONTIER
# ============================================================================

X_base <- as.matrix(df_clean[, inputs])
Y_base <- as.matrix(df_clean[, good_output, drop = FALSE])
Z_base <- as.matrix(df_clean[, bad_outputs_full, drop = FALSE])

scaling_X <- apply(X_base, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
scaling_Y <- apply(Y_base, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
scaling_Z <- apply(Z_base, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))

scale_with_factors <- function(M, factors) sweep(M, 2, factors, "/")

X_base_scaled <- scale_with_factors(X_base, scaling_X)
Y_base_scaled <- scale_with_factors(Y_base, scaling_Y)
Z_base_scaled <- scale_with_factors(Z_base, scaling_Z)

# ============================================================================
# 5. CARBON RECALCULATION (when Manure is modified) – match monte_carlo_sensitivity.R
# ============================================================================

carbon_idx   <- which(bad_outputs_full == "Carbon_emission")
manure_idx   <- which(bad_outputs_full == "Manure")
waste_water_idx <- which(bad_outputs_full == "Waste_water")
energy_idx   <- which(inputs == "Energy_consumption")
EF_Electricity <- 0.5703
EF_E_CH4 <- 1; GWP_CH4 <- 27
EF_M_CH4 <- 5.08; EF_M_N2O <- 0.175; GWP_N2O <- 273
kg_to_tCO2eq <- 1e-3

market_pig_base <- Y_base[, 1]
AP_base <- market_pig_base / 2
manure_base <- Z_base[, manure_idx]
energy_base <- X_base[, energy_idx]
carbon_base <- Z_base[, carbon_idx]

G_E_CH4_base   <- AP_base * EF_E_CH4 * kg_to_tCO2eq * GWP_CH4
G_M_CH4_base   <- AP_base * EF_M_CH4 * kg_to_tCO2eq * GWP_CH4
G_M_N2O_base   <- AP_base * EF_M_N2O * kg_to_tCO2eq * GWP_N2O
energy_MWh_base <- energy_base / 1000
G_Electricity_base <- energy_MWh_base * EF_Electricity

G_Biogas_contribution_prop <- 0.2075
G_Biogas_base_abs <- carbon_base * G_Biogas_contribution_prop
G_Biogas_base <- -G_Biogas_base_abs

valid_idx <- !is.na(G_Biogas_base_abs) & !is.na(manure_base) & 
  G_Biogas_base_abs > 0 & manure_base > 0 & is.finite(G_Biogas_base_abs) & is.finite(manure_base)
n_valid <- sum(valid_idx)
if (n_valid >= 10) {
  lm_biogas <- lm(G_Biogas_base_abs[valid_idx] ~ manure_base[valid_idx])
  G_Biogas_intercept <- coef(lm_biogas)[1]
  G_Biogas_slope <- coef(lm_biogas)[2]
} else {
  G_Biogas_intercept <- mean(G_Biogas_base_abs[valid_idx], na.rm = TRUE)
  G_Biogas_slope <- 0
}
G_Other_base <- pmax(1e-8, carbon_base - (G_E_CH4_base + G_M_CH4_base + G_M_N2O_base + G_Electricity_base + G_Biogas_base))

# Recalculate Carbon for given row indices using new Manure (other components unchanged)
recalc_carbon_for_rows <- function(ix, manure_new, biogas_intercept, biogas_slope) {
  G_Biogas_abs_new <- pmax(1e-8, biogas_intercept + biogas_slope * manure_new)
  G_Biogas_new <- -G_Biogas_abs_new
  carbon_new <- G_E_CH4_base[ix] + G_M_CH4_base[ix] + G_M_N2O_base[ix] +
    G_Electricity_base[ix] + G_Other_base[ix] + G_Biogas_new
  pmax(1e-8, carbon_new)
}

# ============================================================================
# 6. BASELINE EFFICIENCY (all DMUs, fixed frontier = base)
# ============================================================================

frontiers <- list(
  M = list(constrained = constrained_bad_M, weights = as.numeric(bad_weights_M)),
  L = list(constrained = constrained_bad_L, weights = as.numeric(bad_weights_L)),
  G = list(constrained = constrained_bad_G, weights = as.numeric(bad_weights_G)),
  A = list(constrained = constrained_bad_A, weights = as.numeric(bad_weights_A))
)

message("  Computing baseline efficiency (fixed frontier = base data)...")
base_eff <- data.frame(obs_id = seq_len(n_obs), market_pig = df_clean[[good_output]])
for (fid in c("M", "L", "G", "A")) {
  base_eff[[paste0(fid, "_TE_VRS")]] <- sbm_efficiency(
    X_base_scaled, Y_base_scaled, bad = Z_base_scaled,
    constrained_bad_indices = frontiers[[fid]]$constrained,
    bad_weights = frontiers[[fid]]$weights,
    RTS = "vrs", show_dmu_progress = FALSE
  )
}
message("  Baseline done.\n")

# ============================================================================
# 7. UNDERREPORTING REPLICATIONS
# ============================================================================

set.seed(seed)
n_under <- max(1, round(n_obs * pct_farms))

under_eff_list <- vector("list", n_under_rep)
message(sprintf("  Running %d underreporting replications (20%% farms, u~U(%.2f,%.2f))...\n", n_under_rep, u_min, u_max))

for (rep in seq_len(n_under_rep)) {
  if (rep %% 25 == 0 || rep == n_under_rep) cat(sprintf("  Rep %d/%d\n", rep, n_under_rep))
  
  underreport_ix <- sample.int(n_obs, size = n_under, replace = FALSE)
  u_vec <- runif(n_under, min = u_min, max = u_max)
  
  Z_under <- Z_base
  Z_under[underreport_ix, waste_water_idx] <- Z_base[underreport_ix, waste_water_idx] * (1 - u_vec)
  Z_under[underreport_ix, manure_idx]      <- Z_base[underreport_ix, manure_idx]      * (1 - u_vec)
  Z_under[underreport_ix, carbon_idx]   <- recalc_carbon_for_rows(
    underreport_ix, Z_under[underreport_ix, manure_idx], G_Biogas_intercept, G_Biogas_slope)
  
  Z_under_scaled <- scale_with_factors(Z_under, scaling_Z)
  
  X_combined <- rbind(X_base_scaled, X_base_scaled)
  Y_combined <- rbind(Y_base_scaled, Y_base_scaled)
  Z_combined <- rbind(Z_base_scaled, Z_under_scaled)
  perturbed_indices <- (n_obs + 1):(2 * n_obs)
  
  eff_under <- data.frame(obs_id = seq_len(n_obs), market_pig = df_clean[[good_output]], rep = rep)
  for (fid in c("M", "L", "G", "A")) {
    te_all <- sbm_efficiency(X_combined, Y_combined, bad = Z_combined,
                             constrained_bad_indices = frontiers[[fid]]$constrained,
                             bad_weights = frontiers[[fid]]$weights,
                             RTS = "vrs", show_dmu_progress = FALSE,
                             evaluate_indices = perturbed_indices)
    eff_under[[paste0(fid, "_TE_VRS")]] <- te_all[perturbed_indices]
  }
  under_eff_list[[rep]] <- eff_under
}

under_eff_long <- dplyr::bind_rows(under_eff_list)
message("  Underreporting replications done.\n")

# ============================================================================
# 8. LONG-FORMAT DATA FOR PLOT (Baseline vs Underreported, 4 frontiers)
# ============================================================================

scale_col <- "market_pig"
base_long <- base_eff %>%
  dplyr::select(all_of(scale_col), M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS) %>%
  tidyr::pivot_longer(cols = c(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS),
                      names_to = "frontier", values_to = "efficiency") %>%
  dplyr::mutate(frontier = stringr::str_remove(frontier, "_TE_VRS"),
                frontier = dplyr::if_else(frontier == "M", "M", frontier),
                scenario = "Baseline",
                log_scale = log10(pmax(.data[[scale_col]], 1)))
under_long <- under_eff_long %>%
  dplyr::select(all_of(scale_col), rep, M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS) %>%
  tidyr::pivot_longer(cols = c(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS),
                      names_to = "frontier", values_to = "efficiency") %>%
  dplyr::mutate(frontier = stringr::str_remove(frontier, "_TE_VRS"),
                frontier = dplyr::if_else(frontier == "M", "M", frontier),
                scenario = "Underreported",
                log_scale = log10(pmax(.data[[scale_col]], 1)))

plot_data <- dplyr::bind_rows(base_long, under_long) %>%
  dplyr::filter(!is.na(efficiency), !is.na(log_scale), is.finite(10^log_scale), 10^log_scale > 0) %>%
  dplyr::mutate(frontier = dplyr::if_else(frontier == "M", "M", frontier),
                frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE),
                scenario = factor(scenario, levels = c("Baseline", "Underreported")))

# Underreported: mean efficiency per DMU (across 100 reps) for Tobit
under_mean <- under_long %>%
  dplyr::group_by(frontier, market_pig, log_scale) %>%
  dplyr::summarise(efficiency = mean(efficiency, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(log_scale_sq = log_scale^2)

# ============================================================================
# 8b. TOBIT PARAMETERS (Baseline + Underreported) per frontier
# ============================================================================

fit_tobit_safe <- function(data, scenario_name) {
  data <- dplyr::filter(data, !is.na(.data$efficiency), .data$efficiency >= 0, .data$efficiency <= 1,
                        !is.na(.data$log_scale), is.finite(.data$log_scale))
  if (!"log_scale_sq" %in% names(data)) data$log_scale_sq <- data$log_scale^2
  if (nrow(data) < 20) return(list(label = paste0(scenario_name, ": n<20")))
  fit <- tryCatch(
    AER::tobit(efficiency ~ log_scale + log_scale_sq, left = 0, right = 1, data = data),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(label = paste0(scenario_name, ": fit failed")))
  co <- coef(fit)
  vc <- vcov(fit)
  b1 <- co["log_scale"]; b2 <- co["log_scale_sq"]
  se1 <- sqrt(vc["log_scale", "log_scale"])
  se2 <- sqrt(vc["log_scale_sq", "log_scale_sq"])
  tp_scale <- if (abs(b2) > 1e-10) 10^(-b1 / (2 * b2)) else NA_real_
  tp_se <- NA_real_
  if (is.finite(tp_scale) && tp_scale > 100 && tp_scale < 5e5 && !is.na(se1) && !is.na(se2)) {
    boot_tp <- numeric(100)
    for (b in seq_len(100)) {
      idx <- sample.int(nrow(data), replace = TRUE)
      fb <- tryCatch(AER::tobit(efficiency ~ log_scale + log_scale_sq, left = 0, right = 1, data = data[idx, ]),
                     error = function(e) NULL)
      if (!is.null(fb)) {
        cb <- coef(fb)
        if (abs(cb["log_scale_sq"]) > 1e-10)
          boot_tp[b] <- 10^(-cb["log_scale"] / (2 * cb["log_scale_sq"]))
      }
    }
    tp_se <- sd(boot_tp[boot_tp > 100 & boot_tp < 5e5], na.rm = TRUE)
  }
  star1 <- ifelse(2 * (1 - pnorm(abs(b1 / se1))) < 0.05, "*", "")
  star2 <- ifelse(2 * (1 - pnorm(abs(b2 / se2))) < 0.05, "*", "")
  tp_str <- if (is.finite(tp_scale)) sprintf("%.0f", tp_scale) else "—"
  if (is.finite(tp_se) && tp_se > 0) tp_str <- paste0(tp_str, "\u00b1", sprintf("%.0f", tp_se))
  label <- paste0(scenario_name, ":\n",
                  "\u03b2\u2081=", sprintf("%.3f", b1), star1,
                  " \u03b2\u2082=", sprintf("%.3f", b2), star2,
                  "\nTP=", tp_str)
  list(label = label, beta1 = b1, beta2 = b2, tp = tp_scale)
}

# Tobit panels: horizontal layout; right edge aligned with 360000; compact spacing
x_axis_max <- 360000
tobit_annotations <- data.frame()
for (fid in c("M", "L", "G", "A")) {
  x_base  <- x_axis_max * 0.985
  x_under <- 2500
  data_base <- base_long %>% dplyr::filter(frontier == fid) %>%
    dplyr::mutate(log_scale_sq = log_scale^2) %>%
    dplyr::filter(!is.na(efficiency), !is.na(log_scale))
  data_under <- under_mean %>% dplyr::filter(frontier == fid)
  lab_base <- fit_tobit_safe(data_base, "Baseline")
  lab_under <- fit_tobit_safe(data_under, "Underrep.")
  tobit_annotations <- rbind(tobit_annotations,
    data.frame(frontier = factor(fid, levels = c("M", "L", "G", "A"), ordered = TRUE),
               x = x_base,  y = 0.035, label = lab_base$label, hjust = 1),
    data.frame(frontier = factor(fid, levels = c("M", "L", "G", "A"), ordered = TRUE),
               x = x_under, y = 0.035, label = lab_under$label, hjust = 0)
  )
}

# ============================================================================
# 9. LOESS 4-PANEL PLOT (Baseline vs Underreported; Tobit panel)
# ============================================================================

# Design aligned with scale_efficiency_mechanism_analysis.R LOESS_efficiency_scale_4frontiers
SCENARIO_COLORS <- c("Baseline" = "#2166ac", "Underreported" = "#d73027")
LOESS_SPAN <- 0.65

# Draw LOESS + CI in two layers so both ribbons are visible (avoid baseline CI covering underreported CI)
# Baseline first (back), then Underreported (front) so red CI is on top and visible
p_loess <- ggplot(plot_data, aes(x = 10^log_scale, y = efficiency, color = scenario, linetype = scenario, fill = scenario, shape = scenario)) +
  geom_smooth(data = dplyr::filter(plot_data, scenario == "Baseline"),
              method = "loess", se = TRUE, linewidth = 1.2, span = LOESS_SPAN, level = 0.95,
              alpha = 0.2, aes(x = 10^log_scale, y = efficiency, color = scenario, fill = scenario, linetype = scenario),
              inherit.aes = FALSE) +
  geom_smooth(data = dplyr::filter(plot_data, scenario == "Underreported"),
              method = "loess", se = TRUE, linewidth = 1.2, span = LOESS_SPAN, level = 0.95,
              alpha = 0.25, aes(x = 10^log_scale, y = efficiency, color = scenario, fill = scenario, linetype = scenario),
              inherit.aes = FALSE) +
  geom_point(data = dplyr::filter(plot_data, scenario == "Underreported"),
             alpha = 0.45, size = 1.8, show.legend = TRUE) +
  geom_point(data = dplyr::filter(plot_data, scenario == "Baseline"),
             alpha = 0.5, size = 1.9, show.legend = TRUE) +
  scale_x_log10(breaks = c(3000, 10000, 30000, 100000, 360000), labels = scales::comma,
                limits = c(NA, x_axis_max)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(0, 1.25)) +
  scale_color_manual(values = SCENARIO_COLORS, name = "Scenario") +
  scale_fill_manual(values = SCENARIO_COLORS, name = "Scenario") +
  scale_linetype_manual(values = c("Baseline" = "solid", "Underreported" = "dashed"), name = "Scenario") +
  scale_shape_manual(values = c("Baseline" = 16, "Underreported" = 17), name = "Scenario") +
  geom_label(data = tobit_annotations, aes(x = x, y = y, label = label, hjust = hjust),
             inherit.aes = FALSE, vjust = 0, size = 3.8, fontface = "bold",
             color = "black", fill = "white", alpha = 0.85, label.size = 0.5, label.padding = unit(0.3, "lines")) +
  facet_wrap(~frontier, ncol = 2, nrow = 2,
             labeller = labeller(frontier = c(
              "M" = "M-Frontier (Mortality)",
               "L" = "L-Frontier (Local Env)",
               "G" = "G-Frontier (Global)",
               "A" = "A-Frontier (Aggregated)"
             )), drop = FALSE) +
  labs(
    title = "Efficiency vs Scale: Baseline vs Underreported (20% farms, Wastewater & Manure \u221210\u201320%)",
    subtitle = "Fixed frontier; blue solid = baseline; red dashed = underreported (100 replications). Tobit: \u03b2\u2081, \u03b2\u2082, TP=turning point.",
    x = "Production scale (log scale; market pig heads/a)",
    y = "VRS technical efficiency",
    caption = "Points = raw data | Solid/dashed line = LOESS with 95% CI"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "Times", size = 14),
    plot.title = element_text(family = "Times", face = "bold", size = 16, hjust = 0.5, margin = margin(b = 8)),
    plot.subtitle = element_text(family = "Times", size = 13, hjust = 0.5, color = "black", margin = margin(b = 10)),
    axis.text = element_text(family = "Times", size = 13, color = "black"),
    axis.title = element_text(family = "Times", size = 14, face = "bold", color = "black"),
    strip.text = element_text(family = "Times", face = "bold", size = 14, color = "black"),
    plot.caption = element_text(family = "Times", size = 12, hjust = 0.5, color = "black"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    legend.position = "bottom",
    legend.text = element_text(family = "Times", size = 13, color = "black"),
    legend.title = element_text(family = "Times", face = "bold", size = 14, color = "black"),
    plot.margin = margin(5, 4, 5, 5)
  ) +
  guides(
    color = guide_legend(override.aes = list(linetype = c("solid", "dashed"), shape = c(16, 17), size = c(2.1, 1.9)),
                        title = "Scenario", nrow = 1),
    fill = "none", linetype = "none", shape = "none"
  )

# Save plot
fig_path <- file.path(out_dir, "LOESS_underreporting_4frontiers.png")
ggsave(fig_path, p_loess, width = 14, height = 7, dpi = 150)
message(sprintf("  Saved: %s\n", fig_path))

# ============================================================================
# 10. SUMMARY STATS AND CSV EXPORT
# ============================================================================

base_summary <- base_long %>%
  group_by(frontier) %>%
  summarise(mean_eff_baseline = mean(efficiency, na.rm = TRUE),
            sd_baseline = sd(efficiency, na.rm = TRUE), .groups = "drop")
under_summary <- under_long %>%
  group_by(frontier, rep) %>%
  summarise(mean_eff = mean(efficiency, na.rm = TRUE), .groups = "drop") %>%
  group_by(frontier) %>%
  summarise(mean_eff_underreported = mean(mean_eff, na.rm = TRUE),
            sd_underreported = sd(mean_eff, na.rm = TRUE), .groups = "drop")
summary_wide <- base_summary %>%
  left_join(under_summary, by = "frontier") %>%
  mutate(diff_mean = mean_eff_underreported - mean_eff_baseline)

write.csv(summary_wide, file.path(out_dir, "underreporting_summary_by_frontier.csv"), row.names = FALSE)
write.csv(base_eff, file.path(out_dir, "efficiency_baseline.csv"), row.names = FALSE)
readr::write_csv(under_eff_long, file.path(out_dir, "efficiency_underreported_long.csv"))

message("  Summary (mean efficiency by frontier):")
print(summary_wide)
message("\n================================================================================")
message("UNDERREPORTING SENSITIVITY COMPLETE")
message("================================================================================")










