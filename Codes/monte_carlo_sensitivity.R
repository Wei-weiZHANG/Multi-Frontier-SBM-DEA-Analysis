# ============================================================================
# Monte Carlo SBM: Standalone 6-Variable Sensitivity by Scale Group (A-TE-VRS)
# ============================================================================
# Standalone script: Load data (same path as multi_frontier_sbm_dea.R),
# define scale groups (market_pig: <5k, 5-10k, 10-20k, 20-50k, >50k heads/year),
# Frontier = original sample + single-variable perturbed new DMUs (all-6 block not in reference set).
# Single-variable DMUs (>= n per scale x variable) + all-6-perturbed (one per farm); ratio = (eff_single - eff_original)/eff_original, both vs this frontier.
# Stacking: Carbon = electricity then biogas credit; Eutrophication = 0.8*0.8 = 0.64.
#
# Six key variables:
#   1. Electricity (15% reduction): Energy_consumption a%, Carbon_emission c%
#   2. Water use (15% saving): Water_consumption -15%, Water_footprint b%
#   3. FCR (15% reduction): Feed -15%, Eutrophication_potential -20%
#   4. Biogas (15% increase): Carbon_emission d%
#   5. Mortality rate (15% reduction): Dead_pig -15%
#   6. Precise feeding (NUE 0.34→0.5): Eutrophication_potential -20%
#
# Usage: setwd() to project dir, then source("monte_carlo_sensitivity.R")
# ============================================================================

# ============================================================================
# 0. PACKAGES AND CONFIG
# ============================================================================

required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "lpSolve", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Data path: same as main analysis script
xlsx_path <- Sys.getenv("MC_DATA_PATH", "./swine_farm_data.xlsx")
n_new_dmu_per_cell <- 100L
seed <- 12345L
results_base <- Sys.getenv("RESULTS_DIR", unset = "results")
out_dir <- file.path(results_base, "monte_carlo_sensitivity")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Scale groups: market_pig (heads/year)
SCALE_BREAKS <- c(0, 5000, 10000, 20000, 50000, Inf)
SCALE_LABELS <- c("<5k", "5-10k", "10-20k", "20-50k", ">50k")

# Six key variables (display names for plot)
KEY_VARS <- c("Electricity", "Water_use", "FCR", "Biogas", "Mortality_rate", "Precise_feeding")

message("\n================================================================================")
message("MONTE CARLO SBM: 6-Variable A-TE-VRS Sensitivity by Scale Group")
message("================================================================================")
message(sprintf("  New DMUs per (scale group, variable): >= %d  |  Scale groups: %s", n_new_dmu_per_cell, paste(SCALE_LABELS, collapse = ", ")))
message(sprintf("  Key variables: %s", paste(KEY_VARS, collapse = ", ")))
message(sprintf("  Output: %s\n", out_dir))

# ============================================================================
# 1. VARIABLE DEFINITIONS AND WEIGHTS (match main script)
# ============================================================================

inputs <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
good_output <- "Market_pig"
bad_outputs_full <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission",
                       "Water_footprint", "Eutrophication_potential")
# A-frontier only: all six undesirable outputs constrained (equal weights)
constrained_bad_A <- c(1L, 2L, 3L, 4L, 5L, 6L)
bad_weights_A <- c(Waste_water = 1, Manure = 1, Dead_pig = 1,
                  Carbon_emission = 1, Water_footprint = 1, Eutrophication_potential = 1)
bad_weights_A <- bad_weights_A[bad_outputs_full]

# ============================================================================
# 2. SBM AND SCALING
# ============================================================================

scale_matrix <- function(M) {
  scaling_factors <- apply(M, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
  sweep(M, 2, scaling_factors, "/")
}

# reference_indices: when set (e.g. 1:n_obs), the technology set (frontier) uses only these
# DMUs; the evaluated DMU k can be any row index. So new DMUs are evaluated against a fixed frontier.
sbm_efficiency <- function(X, Y, bad = NULL, constrained_bad_indices = NULL,
                           bad_weights = NULL, RTS = "vrs", orientation = "output",
                           show_dmu_progress = TRUE, evaluate_indices = NULL,
                           reference_indices = NULL) {
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
  n_ref <- if (is.null(reference_indices)) n else length(reference_indices)
  if (!is.null(reference_indices)) {
    ref <- as.integer(reference_indices)
    X_ref <- X[ref, , drop = FALSE]
    Y_ref <- Y[ref, , drop = FALSE]
    bad_ref <- if (s_b > 0) bad[ref, , drop = FALSE] else matrix(0, nrow = n_ref, ncol = 0L)
  } else {
    X_ref <- X; Y_ref <- Y; bad_ref <- bad
  }
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
    num_vars <- 1 + n_ref + m + s_g + s_b
    obj_coef <- numeric(num_vars)
    obj_coef[1] <- 1
    norm_denom <- m + s_b
    for (i in seq_len(m)) obj_coef[1 + n_ref + i] <- -1 / (norm_denom * x_k[i])
    if (s_b > 0) {
      for (i in seq_len(s_b)) {
        w <- if (i %in% constrained_bad_indices && length(bad_weights) >= i) bad_weights[i] else 0
        obj_coef[1 + n_ref + m + s_g + i] <- -w / (norm_denom * b_k[i])
      }
    }
    n_constr <- m + s_g + s_b + 1 + if (RTS == "vrs") 1 else 0
    constr <- matrix(0, nrow = n_constr, ncol = num_vars)
    dir_vec <- character(n_constr); rhs_vec <- numeric(n_constr); row_idx <- 1
    for (i in seq_len(m)) {
      constr[row_idx, 2:(1+n_ref)] <- X_ref[, i]; constr[row_idx, 1 + n_ref + i] <- 1
      constr[row_idx, 1] <- -x_k[i]; dir_vec[row_idx] <- "=="; rhs_vec[row_idx] <- 0; row_idx <- row_idx + 1
    }
    for (j in seq_len(s_g)) {
      constr[row_idx, 2:(1+n_ref)] <- Y_ref[, j]; constr[row_idx, 1 + n_ref + m + j] <- -1
      constr[row_idx, 1] <- -y_k[j]; dir_vec[row_idx] <- "=="; rhs_vec[row_idx] <- 0; row_idx <- row_idx + 1
    }
    if (s_b > 0) {
      cset <- if (length(constrained_bad_indices) > 0) constrained_bad_indices else seq_len(s_b)
      for (j in seq_len(s_b)) {
        constr[row_idx, 2:(1+n_ref)] <- bad_ref[, j]; constr[row_idx, 1 + n_ref + m + s_g + j] <- 1
        constr[row_idx, 1] <- -b_k[j]
        dir_vec[row_idx] <- if (j %in% cset) "==" else "<="
        rhs_vec[row_idx] <- 0; row_idx <- row_idx + 1
      }
    }
    constr[row_idx, 1] <- 1
    for (j in seq_len(s_g)) constr[row_idx, 1 + n_ref + m + j] <- 1 / (s_g * y_k[j])
    dir_vec[row_idx] <- "=="; rhs_vec[row_idx] <- 1; row_idx <- row_idx + 1
    if (RTS == "vrs") {
      constr[row_idx, 2:(1+n_ref)] <- 1; constr[row_idx, 1] <- -1
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

scale_with_factors <- function(M, factors) sweep(M, 2, factors, "/")

# ============================================================================
# 3. DATA LOAD, CLEAN, SCALE GROUPS
# ============================================================================

if (!file.exists(xlsx_path)) stop("Data file not found: ", xlsx_path)
suppressWarnings(df_raw <- readxl::read_xlsx(xlsx_path, sheet = 3))
all_required <- unique(c(inputs, good_output, bad_outputs_full))
# Formula parameters from Excel (user-provided, no derived calculation)
formula_cols <- c("Electricity", "EC_noelectricity", "CE_noelectricity", "WF_feed", "CE_nobiogas", "Qbiogas")
df_clean <- df_raw %>%
  dplyr::select(dplyr::all_of(all_required),
                dplyr::any_of(c("Farm", "ID", "Report_year", "No.", "No", "Number", formula_cols))) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_required), ~ suppressWarnings(as.numeric(.x)))) %>%
  tidyr::drop_na(dplyr::all_of(all_required))
for (v in all_required) df_clean[[v]] <- pmax(df_clean[[v]], 1e-8)
# Coerce formula columns to numeric where present
for (v in formula_cols) {
  if (v %in% names(df_clean)) df_clean[[v]] <- suppressWarnings(as.numeric(df_clean[[v]]))
}

df_clean$scale_bin <- cut(
  df_clean[[good_output]],
  breaks = SCALE_BREAKS,
  labels = SCALE_LABELS,
  include.lowest = TRUE,
  right = FALSE
)
df_clean$scale_bin <- factor(df_clean$scale_bin, levels = SCALE_LABELS, ordered = TRUE)
n_obs <- nrow(df_clean)
df_clean$obs_id <- seq_len(n_obs)
message(sprintf("  Loaded: %d farms\n", n_obs))
print(table(df_clean$scale_bin, useNA = "ifany"))

# ============================================================================
# 4. BASE MATRICES AND FORMULA PARAMETERS (from Excel only)
# ============================================================================

X_base <- as.matrix(df_clean[, inputs])
Y_base <- as.matrix(df_clean[, good_output, drop = FALSE])
Z_base <- as.matrix(df_clean[, bad_outputs_full, drop = FALSE])
scaling_X <- apply(X_base, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
scaling_Y <- apply(Y_base, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
scaling_Z <- apply(Z_base, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
X_base_scaled <- scale_with_factors(X_base, scaling_X)
Y_base_scaled <- scale_with_factors(Y_base, scaling_Y)
Z_base_scaled <- scale_with_factors(Z_base, scaling_Z)

# Indices for inputs/bad outputs
carbon_idx <- which(bad_outputs_full == "Carbon_emission")
energy_idx <- which(inputs == "Energy_consumption")
# Formula constants (user-given): a% uses 0.5703, c% uses 0.1229
k_EC_electricity <- 0.5703
k_CE_electricity <- 0.1229
feed_idx <- which(inputs == "Feed")
water_cons_idx <- which(inputs == "Water_consumption")
water_foot_idx <- which(bad_outputs_full == "Water_footprint")
eutro_idx <- which(bad_outputs_full == "Eutrophication_potential")
dead_pig_idx <- which(bad_outputs_full == "Dead_pig")

# Formula parameters: strictly from Excel (no Carbon_emission decomposition)
required_formula_cols <- c("Electricity", "EC_noelectricity", "CE_noelectricity", "WF_feed", "CE_nobiogas", "Qbiogas")
missing_formula <- setdiff(required_formula_cols, names(df_clean))
if (length(missing_formula) > 0)
  stop("Excel must contain formula columns: ", paste(missing_formula, collapse = ", "))
Electricity_vec <- as.numeric(df_clean$Electricity)
EC_noelectricity_vec <- as.numeric(df_clean$EC_noelectricity)
CE_noelectricity_vec <- as.numeric(df_clean$CE_noelectricity)
WF_feed_vec <- as.numeric(df_clean$WF_feed)
CE_nobiogas_vec <- as.numeric(df_clean$CE_nobiogas)
Qbiogas_vec <- as.numeric(df_clean$Qbiogas)
Electricity_vec[is.na(Electricity_vec)] <- 0
EC_noelectricity_vec[is.na(EC_noelectricity_vec)] <- 0
CE_noelectricity_vec[is.na(CE_noelectricity_vec)] <- 0
WF_feed_vec[is.na(WF_feed_vec)] <- 0
CE_nobiogas_vec[is.na(CE_nobiogas_vec)] <- 0
Qbiogas_vec[is.na(Qbiogas_vec)] <- 0

# ============================================================================
# 5. SIX KEY-VARIABLE PERTURBATIONS (return new row: x_row, y_row, z_row)
# ============================================================================

# 1. Electricity: 15% reduction -> Energy_consumption a%, Carbon_emission c%
#    a% = 1 - (Electricity*0.85*1e-3*0.5703 + EC_noelectricity)/Energy_consumption
#    c% = 1 - (Electricity*0.85*1e-3*0.1229 + CE_noelectricity)/Carbon_emission
perturb_electricity <- function(ix) {
  Electricity <- Electricity_vec[ix]
  EC_noelec <- EC_noelectricity_vec[ix]
  CE_noelec <- CE_noelectricity_vec[ix]
  new_EC <- Electricity * 0.85 * 1e-3 * k_EC_electricity + EC_noelec
  new_EC <- pmax(1e-8, new_EC)
  new_CE <- Electricity * 0.85 * 1e-3 * k_CE_electricity + CE_noelec
  new_CE <- pmax(1e-8, new_CE)
  x_row <- X_base[ix, ]; x_row[energy_idx] <- new_EC
  y_row <- Y_base[ix, , drop = FALSE]
  z_row <- Z_base[ix, ]; z_row[carbon_idx] <- new_CE
  list(X = x_row, Y = y_row, Z = z_row)
}

# 2. Water use: 15% saving -> Water_consumption -15%, Water_footprint b%
#    b% = 1 - (0.85*Water_consumption + WF_feed)/Water_footprint  => new_WF = (1-b%)*WF = 0.85*WC + WF_feed
perturb_water <- function(ix) {
  WC <- X_base[ix, water_cons_idx]
  WF_feed <- WF_feed_vec[ix]
  new_WC <- 0.85 * WC
  new_WF <- 0.85 * WC + WF_feed
  new_WF <- pmax(1e-8, new_WF)
  x_row <- X_base[ix, ]; x_row[water_cons_idx] <- pmax(1e-8, new_WC)
  y_row <- Y_base[ix, , drop = FALSE]
  z_row <- Z_base[ix, ]; z_row[water_foot_idx] <- new_WF
  list(X = x_row, Y = y_row, Z = z_row)
}

# 3. FCR: 15% reduction -> Feed -15%, Eutrophication_potential -20%
perturb_fcr <- function(ix) {
  x_row <- X_base[ix, ]; x_row[feed_idx] <- pmax(1e-8, 0.85 * x_row[feed_idx])
  y_row <- Y_base[ix, , drop = FALSE]
  z_row <- Z_base[ix, ]; z_row[eutro_idx] <- pmax(1e-8, 0.8 * z_row[eutro_idx])
  list(X = x_row, Y = y_row, Z = z_row)
}

# 4. Biogas: 15% increase -> Carbon_emission d%
#    d% = 1 - (CE_nobiogas - 1.15*Qbiogas)/Carbon_emission  => new_CE = CE_nobiogas - 1.15*Qbiogas
perturb_biogas <- function(ix) {
  CE_nobiogas <- CE_nobiogas_vec[ix]
  Qbiogas <- Qbiogas_vec[ix]
  new_CE <- CE_nobiogas - 1.15 * Qbiogas
  new_CE <- pmax(1e-8, new_CE)
  x_row <- X_base[ix, ]
  y_row <- Y_base[ix, , drop = FALSE]
  z_row <- Z_base[ix, ]; z_row[carbon_idx] <- new_CE
  list(X = x_row, Y = y_row, Z = z_row)
}

# 5. Mortality rate: 15% reduction -> Dead_pig -15%
perturb_mortality <- function(ix) {
  x_row <- X_base[ix, ]
  y_row <- Y_base[ix, , drop = FALSE]
  z_row <- Z_base[ix, ]; z_row[dead_pig_idx] <- pmax(1e-8, 0.85 * z_row[dead_pig_idx])
  list(X = x_row, Y = y_row, Z = z_row)
}

# 6. Precise feeding: NUE 0.34->0.5 -> Eutrophication_potential -20%
perturb_precise_feeding <- function(ix) {
  x_row <- X_base[ix, ]
  y_row <- Y_base[ix, , drop = FALSE]
  z_row <- Z_base[ix, ]; z_row[eutro_idx] <- pmax(1e-8, 0.8 * z_row[eutro_idx])
  list(X = x_row, Y = y_row, Z = z_row)
}

# All-6 simultaneous perturbation (one per farm). Stacking:
# - Carbon_emission: electricity first (CE_after_elec), then biogas credit: CE_final = CE_after_elec - 0.15*Qbiogas.
# - Eutrophication_potential: FCR -20% and Precise feeding -20% => 0.8*0.8 = 0.64.
perturb_all6 <- function(ix) {
  x_row <- X_base[ix, ]
  y_row <- Y_base[ix, , drop = FALSE]
  z_row <- Z_base[ix, ]

  # 1. Electricity: Energy_consumption, Carbon (first step for stacking)
  E <- Electricity_vec[ix]
  new_EC <- E * 0.85 * 1e-3 * k_EC_electricity + EC_noelectricity_vec[ix]
  new_EC <- pmax(1e-8, new_EC)
  CE_after_elec <- E * 0.85 * 1e-3 * k_CE_electricity + CE_noelectricity_vec[ix]
  CE_after_elec <- pmax(1e-8, CE_after_elec)
  x_row[energy_idx] <- new_EC

  # 2. Water: Water_consumption -15%, Water_footprint
  WC <- x_row[water_cons_idx]
  x_row[water_cons_idx] <- pmax(1e-8, 0.85 * WC)
  z_row[water_foot_idx] <- pmax(1e-8, 0.85 * WC + WF_feed_vec[ix])

  # 3. FCR: Feed -15%, Eutrophication -20% (will multiply by 0.8 again for precise feeding => 0.64)
  x_row[feed_idx] <- pmax(1e-8, 0.85 * x_row[feed_idx])
  z_row[eutro_idx] <- pmax(1e-8, 0.8 * z_row[eutro_idx])

  # 4. Biogas: stack on Carbon (CE_after_elec - additional biogas credit 0.15*Qbiogas)
  z_row[carbon_idx] <- pmax(1e-8, CE_after_elec - 0.15 * Qbiogas_vec[ix])

  # 5. Mortality: Dead_pig -15%
  z_row[dead_pig_idx] <- pmax(1e-8, 0.85 * z_row[dead_pig_idx])

  # 6. Precise feeding: Eutrophication -20% (combined with FCR => 0.8*0.8 = 0.64)
  z_row[eutro_idx] <- pmax(1e-8, 0.8 * z_row[eutro_idx])

  list(X = x_row, Y = y_row, Z = z_row)
}

perturb_funs <- list(
  Electricity = perturb_electricity,
  Water_use = perturb_water,
  FCR = perturb_fcr,
  Biogas = perturb_biogas,
  Mortality_rate = perturb_mortality,
  Precise_feeding = perturb_precise_feeding
)

# ============================================================================
# 6. MONTE CARLO: single-variable DMUs + all-6-perturbed DMUs (one per farm)
#    Frontier = all-6-perturbed DMUs only (no original sample in reference set).
#    Efficiency change ratio: (eff_single - eff_original) / eff_original, both vs this frontier.
# ============================================================================

frontier_A <- list(constrained = constrained_bad_A, weights = as.numeric(bad_weights_A))
set.seed(seed)
new_X_list <- list()
new_Y_list <- list()
new_Z_list <- list()
new_meta <- list(original_ix = integer(0), scale_bin = character(0), variable = character(0))

message("  Phase 1: generating new DMUs (>= ", n_new_dmu_per_cell, " per scale group x variable)...")

for (sg in SCALE_LABELS) {
  idx_in_group <- which(df_clean$scale_bin == sg)
  if (length(idx_in_group) < 1) {
    warning("Scale group ", sg, " has no observations; skipping.")
    next
  }
  for (kv in KEY_VARS) {
    perturb_f <- perturb_funs[[kv]]
    cat(sprintf("\r   %s | %s ...", sg, kv))
    for (s in seq_len(n_new_dmu_per_cell)) {
      ix <- sample(idx_in_group, size = 1)
      new_row <- perturb_f(ix)
      new_X_list[[length(new_X_list) + 1]] <- new_row$X
      new_Y_list[[length(new_Y_list) + 1]] <- c(new_row$Y)
      new_Z_list[[length(new_Z_list) + 1]] <- new_row$Z
      new_meta$original_ix <- c(new_meta$original_ix, ix)
      new_meta$scale_bin <- c(new_meta$scale_bin, sg)
      new_meta$variable <- c(new_meta$variable, kv)
    }
  }
}
cat("\n")
n_new <- length(new_X_list)
message(sprintf("  Generated %d single-variable DMUs.\n", n_new))

X_new <- matrix(unlist(new_X_list), nrow = n_new, byrow = TRUE)
Y_new <- matrix(unlist(new_Y_list), nrow = n_new, byrow = TRUE)
Z_new <- matrix(unlist(new_Z_list), nrow = n_new, byrow = TRUE)

# All-6-perturbed: one DMU per farm (frontier = these only, no original in reference set)
message("  Generating all-6-perturbed DMUs (one per farm, included in frontier)...")
all6_list <- lapply(seq_len(n_obs), perturb_all6)
X_all6 <- matrix(unlist(lapply(all6_list, function(r) r$X)), nrow = n_obs, byrow = TRUE)
Y_all6 <- matrix(unlist(lapply(all6_list, function(r) r$Y)), nrow = n_obs, byrow = TRUE)
Z_all6 <- matrix(unlist(lapply(all6_list, function(r) r$Z)), nrow = n_obs, byrow = TRUE)
message(sprintf("  Generated %d all-6-perturbed DMUs.\n", n_obs))

# Combined: [original | single-variable | all-6-perturbed]. Frontier = original + single-variable new DMUs only.
X_all <- rbind(X_base, X_new, X_all6)
Y_all <- rbind(Y_base, Y_new, Y_all6)
Z_all <- rbind(Z_base, Z_new, Z_all6)
n_total <- n_obs + n_new + n_obs
idx_frontier <- seq_len(n_obs + n_new)  # original sample + single-variable perturbed DMUs (exclude all-6 block)

X_all_scaled <- rbind(
  X_base_scaled,
  scale_with_factors(X_new, scaling_X),
  scale_with_factors(X_all6, scaling_X)
)
Y_all_scaled <- rbind(
  Y_base_scaled,
  scale_with_factors(Y_new, scaling_Y),
  scale_with_factors(Y_all6, scaling_Y)
)
Z_all_scaled <- rbind(
  Z_base_scaled,
  scale_with_factors(Z_new, scaling_Z),
  scale_with_factors(Z_all6, scaling_Z)
)

# Frontier = original + single-variable new DMUs. Ratio = (eff_single - eff_original) / eff_original, both vs this frontier.
message("  Phase 2: computing A-TE-VRS for all DMUs against frontier (original + single-variable, ", length(idx_frontier), " DMUs)...")
eff_all <- sbm_efficiency(X_all_scaled, Y_all_scaled, bad = Z_all_scaled,
                          constrained_bad_indices = frontier_A$constrained,
                          bad_weights = frontier_A$weights,
                          RTS = "vrs", show_dmu_progress = TRUE,
                          evaluate_indices = NULL,
                          reference_indices = idx_frontier)
eff_all <- pmax(0, pmin(1.5, eff_all))
message("  Done.\n")

new_dmu_indices <- (n_obs + 1L):(n_obs + n_new)  # single-variable DMUs only (exclude all-6 block)
eff_original <- eff_all[new_meta$original_ix]
eff_new_dmu <- eff_all[new_dmu_indices]
base_eff_positive <- is.finite(eff_original) & eff_original > 0
ratio_change <- rep(NA_real_, n_new)
ratio_change[base_eff_positive] <- (eff_new_dmu[base_eff_positive] - eff_original[base_eff_positive]) / eff_original[base_eff_positive]

results_df <- data.frame(
  scale_bin = new_meta$scale_bin,
  variable = new_meta$variable,
  sample_id = seq_len(n_new),
  original_ix = new_meta$original_ix,
  base_A_TE_VRS = eff_original,
  new_A_TE_VRS = eff_new_dmu,
  efficiency_change_ratio = ratio_change,
  stringsAsFactors = FALSE
)
results_df <- results_df[is.finite(results_df$efficiency_change_ratio), ]
message(sprintf("  Sensitivity records with valid ratio: %d\n", nrow(results_df)))

# ============================================================================
# 8. AGGREGATE: mean, sd, and 95% CI of efficiency change ratio by scale_bin and variable
# ============================================================================

summary_sensitivity <- results_df %>%
  dplyr::group_by(scale_bin, variable) %>%
  dplyr::summarise(
    mean_ratio = mean(efficiency_change_ratio, na.rm = TRUE),
    sd_ratio = sd(efficiency_change_ratio, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    se_ratio = ifelse(n > 1, sd_ratio / sqrt(n), NA_real_),
    df = pmax(n - 1, 1),
    t_crit = qt(0.975, df),
    ci_half = ifelse(n > 1, t_crit * se_ratio, NA_real_),
    scale_bin = factor(scale_bin, levels = SCALE_LABELS, ordered = TRUE),
    variable = factor(variable, levels = KEY_VARS, ordered = TRUE)
  )
message("  Sensitivity summary (efficiency change ratio by scale group and variable):")
print(summary_sensitivity)

# ============================================================================
# 9. SAVE CSV AND PLOT: bar chart + error bars (publication-ready, 23 x 9 cm)
# ============================================================================

write.csv(results_df, file.path(out_dir, "monte_carlo_6var_sensitivity_raw.csv"), row.names = FALSE)
write.csv(summary_sensitivity, file.path(out_dir, "monte_carlo_6var_sensitivity_summary.csv"), row.names = FALSE)

# Paul Tol "muted" palette
palette_journal <- c(
  Electricity      = "#332288",
  Water_use        = "#88CCEE",
  FCR              = "#44AA77",
  Biogas           = "#AA4499",
  Mortality_rate   = "#DDCC77",
  Precise_feeding  = "#CC6677"
)
palette_journal <- palette_journal[KEY_VARS]
legend_labels <- c(
  Electricity      = "Electricity consumption",
  Water_use        = "Water utilization",
  FCR              = "Feed conversion ratio",
  Biogas           = "Biogas utilization",
  Mortality_rate   = "Mortality ratio",
  Precise_feeding  = "Precise feeding"
)
legend_labels <- legend_labels[KEY_VARS]

p <- ggplot2::ggplot(summary_sensitivity,
                    ggplot2::aes(x = scale_bin, y = mean_ratio, fill = variable)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 1, colour = "grey30", linewidth = 0.25) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.68), width = 0.68) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_ratio - ci_half, ymax = mean_ratio + ci_half),
    position = ggplot2::position_dodge(0.68),
    width = 0.14,
    linewidth = 0.35,
    colour = "grey20",
    na.rm = TRUE
  ) +
  ggplot2::scale_fill_manual(values = palette_journal, labels = legend_labels, name = NULL) +
  ggplot2::scale_y_continuous(
    breaks = seq(-0.15, 0.25, 0.05),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    x = "Scale category (market pig heads/a)",
    y = "Efficiency change ratio"
  ) +
  ggplot2::theme(
    panel.background   = ggplot2::element_rect(fill = "white", colour = NA),
    panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.border       = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.4),
    axis.title         = ggplot2::element_text(size = 9, colour = "black", face = "plain"),
    axis.text          = ggplot2::element_text(size = 8, colour = "black"),
    axis.ticks         = ggplot2::element_line(colour = "grey40", linewidth = 0.25),
    axis.ticks.length  = ggplot2::unit(1.2, "pt"),
    legend.position    = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.title       = ggplot2::element_blank(),
    legend.text        = ggplot2::element_text(size = 8),
    legend.key.size    = ggplot2::unit(4, "mm"),
    legend.spacing.y   = ggplot2::unit(2.2, "mm"),
    legend.background  = ggplot2::element_rect(fill = "white", colour = "grey85", linewidth = 0.25),
    legend.margin      = ggplot2::margin(1.5, 2, 2, 2),
    plot.margin        = ggplot2::margin(5, 5, 5, 5, "mm")
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2))
ggplot2::ggsave(
  file.path(out_dir, "sensitivity_6keyvars_by_scale_bar.png"),
  p, width = 24, height = 6, units = "cm", dpi = 300
)
ggplot2::ggsave(
  file.path(out_dir, "sensitivity_6keyvars_by_scale_bar.pdf"),
  p, width = 23.95 / 2.54, height = 6.52 / 2.54, units = "in", device = "pdf"
)
message(sprintf("  Figure saved: %s (23×9 cm, 300 dpi PNG + PDF)", file.path(out_dir, "sensitivity_6keyvars_by_scale_bar")))

message("\n================================================================================")
message("MONTE CARLO SBM (6-VARIABLE SENSITIVITY) COMPLETE")
message("================================================================================")
message(sprintf("  Scale groups: %s", paste(SCALE_LABELS, collapse = ", ")))
message(sprintf("  Key variables: %s", paste(KEY_VARS, collapse = ", ")))
message(sprintf("  Output dir: %s", out_dir))
message("")








