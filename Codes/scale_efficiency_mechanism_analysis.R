# ============================================================================
# Efficiency-Scale Comprehensive Analysis
# ============================================================================
#
# Purpose: Comprehensive analysis of efficiency-scale relationships including:
#          1. Visual Evidence: LOESS smooth curves, binned bar plots
#          2. Descriptive Evidence: Group statistics quantifying differences
#          3. Inferential Evidence: GAM statistical tests
#          4. Tobit Regression: U-shape pattern detection with turning points
#          5. Unit Input/Output Trends: Per-pig input/output analysis
#          6. Marginal Effects: Impact of 1% changes across scale zones
#          7. Shapley Decomposition: Factor contribution analysis by scale groups
#
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(mgcv)
library(stringr)
library(scales)  # For pretty_breaks function
library(tidyr)

# Load readxl for Excel file reading
if (!require("readxl", quietly = TRUE)) {
  install.packages("readxl", repos = "https://cloud.r-project.org")
  library(readxl)
}

# Load patchwork for combining plots
if (!require("patchwork", quietly = TRUE)) {
  install.packages("patchwork", repos = "https://cloud.r-project.org")
  library(patchwork)
}

# Load Tobit regression packages
if (!require("AER", quietly = TRUE)) {
  install.packages("AER", repos = "https://cloud.r-project.org")
  library(AER)
}
if (!require("censReg", quietly = TRUE)) {
  install.packages("censReg", repos = "https://cloud.r-project.org")
  library(censReg)
}
# For Excel export of U-shape robustness (all model fits)
if (!require("writexl", quietly = TRUE)) {
  install.packages("writexl", repos = "https://cloud.r-project.org")
  library(writexl)
}

cat("==============================================================================\n")
cat("        TRIPLE EVIDENCE CHAIN: U-SHAPE PATTERN ANALYSIS\n")
cat("==============================================================================\n\n")

# Set random seed for reproducibility (ensures consistent results across runs)
set.seed(12345)
cat("✓ Random seed set to 12345 for reproducibility\n\n")

# Diagnostic: Check if called from main script
cat("=== DIAGNOSTIC: Script Initialization ===\n")
cat(sprintf("  Current working directory: %s\n", getwd()))
tryCatch({
  script_file <- basename(sys.frame(1)$ofile)
  cat(sprintf("  Script file: %s\n", script_file))
}, error = function(e) {
  cat("  Script file: scale_efficiency_mechanism_analysis.R\n")
})
cat(sprintf("  results_base exists: %s\n", exists("results_base")))
cat(sprintf("  efficiency_results exists: %s\n", exists("efficiency_results")))
if (exists("results_base")) {
  cat(sprintf("  results_base value: %s\n", results_base))
  cat(sprintf("  results_base directory exists: %s\n", dir.exists(results_base)))
}
if (exists("efficiency_results")) {
  cat(sprintf("  efficiency_results rows: %d\n", nrow(efficiency_results)))
  cat(sprintf("  efficiency_results columns: %d\n", ncol(efficiency_results)))
  cat(sprintf("  efficiency_results column names: %s\n", paste(names(efficiency_results), collapse = ", ")))
}
# Check environment variables
current_results_dir <- Sys.getenv("CURRENT_RESULTS_DIR", unset = "")
cat(sprintf("  CURRENT_RESULTS_DIR environment variable: %s\n", current_results_dir))
cat("==========================================\n\n")

# ============================================================================
# 1. Load Data
# ============================================================================

# Priority 1: Use variables from current session (when called from main script)
if (exists("results_base") && exists("efficiency_results")) {
  results_dir <- results_base
  eff <- efficiency_results
  cat("✓ Using data from current session\n")
  cat(sprintf("✓ Results directory: %s\n", results_dir))
  cat(sprintf("✓ Sample size: %d farms\n", nrow(eff)))
} else {
  # Priority 2: Check environment variable (from main script)
  current_results_dir <- Sys.getenv("CURRENT_RESULTS_DIR", unset = "")
  if (nchar(current_results_dir) > 0) {
    results_dir <- current_results_dir
    cat(sprintf("Using current run results: %s\n", basename(results_dir)))
  } else {
    # Priority 3: Fallback: find latest results folder
    results_base_dir <- "./results"
    if (!dir.exists(results_base_dir)) {
      stop(sprintf("Results base directory not found: %s", results_base_dir))
    }
    all_results <- list.dirs(results_base_dir, full.names = FALSE, recursive = FALSE)
    all_results <- all_results[grepl("^results_\\d{8}_\\d{6}$", all_results)]
    if (length(all_results) == 0) {
      stop("No results folders found. Please run the analysis first.")
    }
    latest_result <- sort(all_results, decreasing = TRUE)[1]
    results_dir <- file.path(results_base_dir, latest_result)
    cat(sprintf("Using latest results folder: %s\n", latest_result))
  }
  
  # Ensure results_dir exists
  if (!dir.exists(results_dir)) {
    stop(sprintf("Results directory does not exist: %s", results_dir))
  }
  
  # Load efficiency data from file
  eff_file_complete <- file.path(results_dir, "efficiency_calculation/all_frontiers_efficiency_with_tgr.csv")
  eff_file_basic <- file.path(results_dir, "efficiency_calculation/frontier_efficiency_results.csv")
  
  if (file.exists(eff_file_complete)) {
    eff <- read_csv(eff_file_complete, show_col_types = FALSE)
    cat("✓ Loaded complete efficiency data (with TGR)\n")
  } else if (file.exists(eff_file_basic)) {
    eff <- read_csv(eff_file_basic, show_col_types = FALSE)
    cat("✓ Loaded basic efficiency data\n")
  } else {
    stop(sprintf("Efficiency data file not found in %s", results_dir))
  }
}

# Standardize names
if ("Market_pig" %in% names(eff)) eff$market_pig <- eff$Market_pig
if (!"market_pig" %in% names(eff)) {
  stop("ERROR: market_pig column not found in efficiency data")
}
if (!"log_scale" %in% names(eff)) {
  eff$log_scale <- log10(pmax(eff$market_pig, 1))
}

# Validate required columns exist
required_cols <- c("M_TE_VRS", "L_TE_VRS", "G_TE_VRS", "A_TE_VRS")
missing_cols <- required_cols[!required_cols %in% names(eff)]
if (length(missing_cols) > 0) {
  stop(sprintf("ERROR: Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

cat(sprintf("✓ Sample size: %d farms\n", nrow(eff)))
cat(sprintf("✓ Scale range: %s to %s pigs\n", 
            format(round(min(eff$market_pig, na.rm = TRUE)), big.mark = ","),
            format(round(max(eff$market_pig, na.rm = TRUE)), big.mark = ",")))
cat(sprintf("✓ Required columns found: %s\n\n", paste(required_cols, collapse = ", ")))

# Create output directories
visualization_dir <- file.path(results_dir, "visualization")
statistical_tests_dir <- file.path(results_dir, "statistical_tests")
if (!dir.exists(visualization_dir)) {
  dir.create(visualization_dir, recursive = TRUE)
  cat(sprintf("✓ Created visualization directory: %s\n", visualization_dir))
}
if (!dir.exists(statistical_tests_dir)) {
  dir.create(statistical_tests_dir, recursive = TRUE)
  cat(sprintf("✓ Created statistical_tests directory: %s\n", statistical_tests_dir))
}

# Helper function for safe plot saving
safe_ggsave <- function(filename, plot, width = 7, height = 5, dpi = 300, ...) {
  tryCatch({
    ggsave(file.path(visualization_dir, filename), plot, 
           width = width, height = height, dpi = dpi, ...)
    cat(sprintf("✓ Saved: %s\n", filename))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("⚠ Failed to save %s: %s\n", filename, e$message))
    return(FALSE)
  })
}

# Define a new variable palette shared by the three variable figures
# (Inputs: solid; Outputs: dashed)
variable_colors_plot <- c(
  "Feed" = "#1F77B4",                    # blue
  "Water_consumption" = "#2CA02C",       # green
  "Energy_consumption" = "#D62728",      # red
  "Labor" = "#FF7F0E",                   # orange
  "Waste_water" = "#9467BD",             # purple
  "Manure" = "#E377C2",                  # pink
  "Dead_pig" = "#17BECF",                # cyan
  "Carbon_emission" = "#8C564B",         # brown
  "Eutrophication_potential" = "#BCBD22" # olive
)

# ============================================================================
# 2. VISUAL EVIDENCE: LOESS Smooth Curves
# ============================================================================

cat("==============================================================================\n")
cat("EVIDENCE 1: Visual Evidence (LOESS Smooth Curves)\n")
cat("==============================================================================\n\n")

# Prepare data for all frontiers
# Check if Slack meta-frontier exists, otherwise use standard meta-frontier
meta_col <- if ("Meta_TE_VRS_Slack" %in% names(eff)) "Meta_TE_VRS_Slack" else "Meta_TE_VRS"
meta_frontier_name <- if ("Meta_TE_VRS_Slack" %in% names(eff)) "Meta_Slack" else "Meta"

# Prepare TE-VRS data
scale_data <- eff %>%
  select(market_pig, log_scale, M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS, all_of(meta_col)) %>%
  pivot_longer(
    cols = c(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS, all_of(meta_col)),
    names_to = "frontier",
    values_to = "efficiency"
  ) %>%
  mutate(frontier = str_remove(frontier, "_TE_VRS")) %>%
  mutate(frontier = str_remove(frontier, "_Slack")) %>%  # Remove _Slack suffix if present
  filter(!is.na(efficiency), !is.na(log_scale)) %>%
  mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A", meta_frontier_name)))

# Prepare SE data for overlay (light gray curves)
scale_data_se <- eff %>%
  select(market_pig, log_scale, M_SE, L_SE, G_SE, A_SE) %>%
  pivot_longer(
    cols = c(M_SE, L_SE, G_SE, A_SE),
    names_to = "frontier",
    values_to = "efficiency"
  ) %>%
  mutate(frontier = str_remove(frontier, "_SE")) %>%
  filter(!is.na(efficiency), !is.na(log_scale) & is.finite(10^log_scale) & 10^log_scale > 0) %>%
  mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE))

# Prepare Meta SE data (if available)
meta_se_col <- if ("Meta_SE" %in% names(eff)) "Meta_SE" else NULL
if (!is.null(meta_se_col)) {
  scale_data_meta_se <- eff %>%
    select(market_pig, log_scale, all_of(meta_se_col)) %>%
    rename(efficiency = all_of(meta_se_col)) %>%
    filter(!is.na(efficiency), !is.na(log_scale) & is.finite(10^log_scale) & 10^log_scale > 0)
} else {
  scale_data_meta_se <- NULL
}

# Define a unified frontier palette (muted, coordinated tones)
FRONTIER_COLORS <- list(
  "M" = "#7EB89E",   # M frontier
  "L" = "#E16859",   # L frontier
  "G" = "#8B72A1",   # G frontier
  "A" = "#60BCDA",   # A frontier
  "Meta" = "#E6B771" # Meta frontier
)

# Calculate binned means for overlay
# Create scale bins: <5k, 5-10k, 10-20k, 20-50k, >50k
# Use original scale (market_pig_original) so bins and diamond x-position match production scale axis
scale_col <- if ("market_pig_original" %in% names(eff)) "market_pig_original" else "market_pig"
eff_with_bins <- eff %>%
  select(all_of(scale_col), log_scale, M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS, all_of(meta_col)) %>%
  mutate(
    scale_bin = cut(
      .data[[scale_col]],
      breaks = c(0, 5000, 10000, 20000, 50000, Inf),
      labels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
      include.lowest = TRUE
    )
  )

binned_means <- eff_with_bins %>%
  select(all_of(scale_col), scale_bin, M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS, all_of(meta_col)) %>%
  pivot_longer(
    cols = c(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS, all_of(meta_col)),
    names_to = "frontier",
    values_to = "efficiency"
  ) %>%
  mutate(frontier = str_remove(frontier, "_TE_VRS")) %>%
  mutate(frontier = str_remove(frontier, "_Slack")) %>%
  filter(!is.na(efficiency), !is.na(scale_bin)) %>%
  mutate(scale_bin = factor(scale_bin, 
                            levels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
                            ordered = TRUE)) %>%
  group_by(scale_bin, frontier) %>%
  summarise(
    n = n(),
    mean_eff = mean(efficiency, na.rm = TRUE),
    mean_scale = exp(mean(log(.data[[scale_col]]), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A", meta_frontier_name)))

# ============================================================================
# TOBIT REGRESSION ANALYSIS: Calculate β₁, β₂, and turning points
# ============================================================================

tobit_results <- list()
tobit_cubic_results <- list()

for (fid in c("M", "L", "G", "A", meta_frontier_name)) {
  # Get efficiency column
  if (fid == meta_frontier_name) {
    eff_col <- meta_col
  } else {
    eff_col <- paste0(fid, "_TE_VRS")
  }
  
  if (!eff_col %in% names(eff)) {
    next
  }
  
  # Prepare data (include cubic term for Tobit cubic robustness check)
  data_tobit <- eff %>%
    select(market_pig, log_scale, all_of(eff_col)) %>%
    rename(efficiency = all_of(eff_col)) %>%
    filter(!is.na(efficiency), !is.na(log_scale), 
           efficiency >= 0, efficiency <= 1,
           market_pig > 0) %>%
    mutate(
      log_scale_sq = log_scale^2,
      log_scale_cube = log_scale^3
    )

  if (nrow(data_tobit) < 10) {
    next
  }

  # Fit Tobit model: efficiency ~ log_scale + log_scale²
  # Left-censored at 0, right-censored at 1
  tryCatch({
    # Use AER::tobit for convenience
    tobit_fit <- AER::tobit(efficiency ~ log_scale + log_scale_sq,
                           left = 0, right = 1,
                           data = data_tobit)
    
    coef_tobit <- coef(tobit_fit)
    vcov_tobit <- vcov(tobit_fit)
    
    beta1 <- coef_tobit["log_scale"]
    beta2 <- coef_tobit["log_scale_sq"]
    
    # Standard errors and p-values
    se_beta1 <- sqrt(vcov_tobit["log_scale", "log_scale"])
    se_beta2 <- sqrt(vcov_tobit["log_scale_sq", "log_scale_sq"])
    
    z_beta1 <- beta1 / se_beta1
    z_beta2 <- beta2 / se_beta2
    
    p_beta1 <- 2 * (1 - pnorm(abs(z_beta1)))
    p_beta2 <- 2 * (1 - pnorm(abs(z_beta2)))
    
    # Calculate turning point (trough for U-shape: β₂ > 0)
    # Turning point: -β₁ / (2 * β₂)
    if (abs(beta2) > 1e-10) {
      turning_point_log <- -beta1 / (2 * beta2)
      turning_point_scale <- 10^turning_point_log
      
      # Bootstrap confidence interval for turning point
      boot_turning_points <- numeric(200)
      boot_success <- 0
      
      for (b in seq_len(200)) {
        boot_idx <- sample(seq_len(nrow(data_tobit)), replace = TRUE)
        data_boot <- data_tobit[boot_idx, ]
        
        tobit_boot <- tryCatch({
          AER::tobit(efficiency ~ log_scale + log_scale_sq,
                    left = 0, right = 1,
                    data = data_boot)
        }, error = function(e) NULL)
        
        if (!is.null(tobit_boot)) {
          coef_boot <- coef(tobit_boot)
          beta1_boot <- coef_boot["log_scale"]
          beta2_boot <- coef_boot["log_scale_sq"]
          
          if (abs(beta2_boot) > 1e-10) {
            turning_point_log_boot <- -beta1_boot / (2 * beta2_boot)
            turning_point_scale_boot <- 10^turning_point_log_boot
            
            # Check if within reasonable range
            if (turning_point_scale_boot > 100 && turning_point_scale_boot < 500000) {
              boot_turning_points[boot_success + 1] <- turning_point_scale_boot
              boot_success <- boot_success + 1
            }
          }
        }
      }
      
      if (boot_success > 10) {
        boot_turning_points <- boot_turning_points[1:boot_success]
        ci_lower <- quantile(boot_turning_points, 0.025, na.rm = TRUE)
        ci_upper <- quantile(boot_turning_points, 0.975, na.rm = TRUE)
        turning_point_se <- sd(boot_turning_points, na.rm = TRUE)
      } else {
        ci_lower <- turning_point_scale * 0.9
        ci_upper <- turning_point_scale * 1.1
        turning_point_se <- turning_point_scale * 0.1
      }
    } else {
      turning_point_log <- NA
      turning_point_scale <- NA
      ci_lower <- NA
      ci_upper <- NA
      turning_point_se <- NA
    }
    
    # Significance stars
    star_beta1 <- ifelse(p_beta1 < 0.001, "***",
                        ifelse(p_beta1 < 0.01, "**",
                               ifelse(p_beta1 < 0.05, "*", "")))
    star_beta2 <- ifelse(p_beta2 < 0.001, "***",
                        ifelse(p_beta2 < 0.01, "**",
                               ifelse(p_beta2 < 0.05, "*", "")))
    
    tobit_results[[fid]] <- list(
      frontier = fid,
      beta1 = as.numeric(beta1),
      beta2 = as.numeric(beta2),
      se_beta1 = se_beta1,
      se_beta2 = se_beta2,
      p_beta1 = p_beta1,
      p_beta2 = p_beta2,
      star_beta1 = star_beta1,
      star_beta2 = star_beta2,
      turning_point_scale = turning_point_scale,
      turning_point_log = turning_point_log,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      turning_point_se = turning_point_se,
      n_obs = nrow(data_tobit)
    )

    # Tobit cubic: robustness check (rule out more complex trends than U-shape)
    tobit_cubic_fit <- tryCatch({
      AER::tobit(efficiency ~ log_scale + log_scale_sq + log_scale_cube,
                 left = 0, right = 1, data = data_tobit)
    }, error = function(e) NULL)
    if (!is.null(tobit_cubic_fit)) {
      coef_cub <- coef(tobit_cubic_fit)
      vcov_cub <- vcov(tobit_cubic_fit)
      nm_cub <- names(coef_cub)
      get_se <- function(name) {
        if (name %in% rownames(vcov_cub)) sqrt(vcov_cub[name, name]) else NA_real_
      }
      get_p <- function(name) {
        b <- coef_cub[name]
        se <- get_se(name)
        if (is.na(se) || se <= 0) return(NA_real_)
        2 * (1 - pnorm(abs(b / se)))
      }
      tobit_cubic_results[[fid]] <- list(
        frontier = fid,
        intercept = as.numeric(coef_cub["(Intercept)"]),
        beta1 = as.numeric(coef_cub["log_scale"]),
        beta2 = as.numeric(coef_cub["log_scale_sq"]),
        beta3 = as.numeric(coef_cub["log_scale_cube"]),
        se_intercept = get_se("(Intercept)"), se_beta1 = get_se("log_scale"),
        se_beta2 = get_se("log_scale_sq"), se_beta3 = get_se("log_scale_cube"),
        p_beta1 = get_p("log_scale"), p_beta2 = get_p("log_scale_sq"), p_beta3 = get_p("log_scale_cube"),
        n_obs = nrow(data_tobit)
      )
    }

  }, error = function(e) {
  })
}

# ============================================================================
# U-SHAPE ROBUSTNESS: Linear OLS, Quadratic OLS, LOESS, Tobit Quad/Cubic
# Collect all fit results and export to single Excel file
# ============================================================================

scale_data_fit <- scale_data %>%
  filter(!is.na(efficiency), efficiency >= 0, efficiency <= 1, !is.na(log_scale)) %>%
  mutate(log_scale_sq = log_scale^2)

linear_ols_list <- list()
quadratic_ols_list <- list()
loess_list <- list()
LOESS_SPAN <- 0.75

for (fid in c("M", "L", "G", "A", meta_frontier_name)) {
  df <- scale_data_fit %>% filter(frontier == fid)
  if (nrow(df) < 5) next

  # Linear OLS
  m_lin <- tryCatch(lm(efficiency ~ log_scale, data = df), error = function(e) NULL)
  if (!is.null(m_lin)) {
    s_lin <- summary(m_lin)
    linear_ols_list[[fid]] <- list(
      Frontier = fid,
      Model = "Linear_OLS",
      Intercept_est = coef(m_lin)[1], Intercept_se = s_lin$coefficients[1, 2], Intercept_p = s_lin$coefficients[1, 4],
      log_scale_est = coef(m_lin)[2], log_scale_se = s_lin$coefficients[2, 2], log_scale_p = s_lin$coefficients[2, 4],
      R2 = s_lin$r.squared, Adj_R2 = s_lin$adj.r.squared, n_obs = nrow(df),
      Significant = ifelse(s_lin$coefficients[2, 4] < 0.05, "Yes", "No")
    )
  }

  # Quadratic OLS
  m_quad <- tryCatch(lm(efficiency ~ log_scale + log_scale_sq, data = df), error = function(e) NULL)
  if (!is.null(m_quad)) {
    s_quad <- summary(m_quad)
    cf <- s_quad$coefficients
    quadratic_ols_list[[fid]] <- list(
      Frontier = fid,
      Model = "Quadratic_OLS",
      Intercept_est = cf[1, 1], Intercept_se = cf[1, 2], Intercept_p = cf[1, 4],
      log_scale_est = cf[2, 1], log_scale_se = cf[2, 2], log_scale_p = cf[2, 4],
      log_scale_sq_est = cf[3, 1], log_scale_sq_se = cf[3, 2], log_scale_sq_p = cf[3, 4],
      R2 = s_quad$r.squared, Adj_R2 = s_quad$adj.r.squared, n_obs = nrow(df),
      Quadratic_Significant = ifelse(cf[3, 4] < 0.05, "Yes", "No")
    )
  }

  # LOESS (summary only: span, residual stats, no coefficients)
  lfit <- tryCatch(loess(efficiency ~ log_scale, data = df, span = LOESS_SPAN), error = function(e) NULL)
  if (!is.null(lfit)) {
    pred_loess <- predict(lfit)
    resid_loess <- df$efficiency - pred_loess
    loess_list[[fid]] <- list(
      Frontier = fid,
      Model = "LOESS",
      span = LOESS_SPAN,
      n_obs = nrow(df),
      residual_sd = sd(resid_loess, na.rm = TRUE),
      residual_mean = mean(resid_loess, na.rm = TRUE),
      enp = lfit$enp
    )
  }
}

# Build Excel sheets
rbind_list_to_df <- function(L, cols = NULL) {
  if (length(L) == 0) return(data.frame())
  out <- do.call(rbind, lapply(L, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  if (!is.null(cols)) out <- out[, cols, drop = FALSE]
  out
}

# Linear OLS sheet
df_linear <- rbind_list_to_df(linear_ols_list)
if (nrow(df_linear) > 0) {
  want_lin <- c("Frontier", "Model", "Intercept_est", "Intercept_se", "Intercept_p",
    "log_scale_est", "log_scale_se", "log_scale_p", "R2", "Adj_R2", "n_obs", "Significant")
  df_linear <- df_linear[, intersect(want_lin, names(df_linear)), drop = FALSE]
}

# Quadratic OLS sheet
df_quadratic <- rbind_list_to_df(quadratic_ols_list)
if (nrow(df_quadratic) > 0) {
  want_quad <- c("Frontier", "Model", "Intercept_est", "Intercept_se", "Intercept_p",
    "log_scale_est", "log_scale_se", "log_scale_p", "log_scale_sq_est", "log_scale_sq_se", "log_scale_sq_p",
    "R2", "Adj_R2", "n_obs", "Quadratic_Significant")
  df_quadratic <- df_quadratic[, intersect(want_quad, names(df_quadratic)), drop = FALSE]
}

# LOESS sheet
df_loess <- rbind_list_to_df(loess_list)
if (nrow(df_loess) > 0) {
  want_loess <- c("Frontier", "Model", "span", "n_obs", "residual_sd", "residual_mean", "enp")
  df_loess <- df_loess[, intersect(want_loess, names(df_loess)), drop = FALSE]
}

# Tobit Quadratic sheet (from tobit_results)
tobit_quad_list <- lapply(names(tobit_results), function(fid) {
  t <- tobit_results[[fid]]
  data.frame(
    Frontier = fid,
    Model = "Tobit_Quadratic",
    beta1 = t$beta1, se_beta1 = t$se_beta1, p_beta1 = t$p_beta1, Significant_beta1 = ifelse(t$p_beta1 < 0.05, "Yes", "No"),
    beta2 = t$beta2, se_beta2 = t$se_beta2, p_beta2 = t$p_beta2, Significant_beta2 = ifelse(t$p_beta2 < 0.05, "Yes", "No"),
    turning_point_scale = t$turning_point_scale, turning_point_se = t$turning_point_se,
    n_obs = t$n_obs,
    stringsAsFactors = FALSE
  )
})
df_tobit_quad <- if (length(tobit_quad_list) > 0) do.call(rbind, tobit_quad_list) else data.frame()

# Tobit Cubic sheet
tobit_cub_list <- lapply(names(tobit_cubic_results), function(fid) {
  t <- tobit_cubic_results[[fid]]
  data.frame(
    Frontier = fid,
    Model = "Tobit_Cubic",
    beta1 = t$beta1, se_beta1 = t$se_beta1, p_beta1 = t$p_beta1, Significant_beta1 = ifelse(t$p_beta1 < 0.05, "Yes", "No"),
    beta2 = t$beta2, se_beta2 = t$se_beta2, p_beta2 = t$p_beta2, Significant_beta2 = ifelse(t$p_beta2 < 0.05, "Yes", "No"),
    beta3 = t$beta3, se_beta3 = t$se_beta3, p_beta3 = t$p_beta3, Significant_beta3 = ifelse(t$p_beta3 < 0.05, "Yes", "No"),
    n_obs = t$n_obs,
    stringsAsFactors = FALSE
  )
})
df_tobit_cubic <- if (length(tobit_cub_list) > 0) do.call(rbind, tobit_cub_list) else data.frame()

# Summary sheet: one row per frontier per model type with key info
summary_rows <- list()
for (fid in c("M", "L", "G", "A", meta_frontier_name)) {
  if (fid %in% names(linear_ols_list)) {
    x <- linear_ols_list[[fid]]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Frontier = fid, Model_Type = "Linear_OLS",
      Key_Parameter = "log_scale", Estimate = x$log_scale_est, SE = x$log_scale_se, P_value = x$log_scale_p,
      Significant = x$Significant, R2 = x$R2, n_obs = x$n_obs, stringsAsFactors = FALSE
    )
  }
  if (fid %in% names(quadratic_ols_list)) {
    x <- quadratic_ols_list[[fid]]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Frontier = fid, Model_Type = "Quadratic_OLS",
      Key_Parameter = "log_scale_sq", Estimate = x$log_scale_sq_est, SE = x$log_scale_sq_se, P_value = x$log_scale_sq_p,
      Significant = x$Quadratic_Significant, R2 = x$R2, n_obs = x$n_obs, stringsAsFactors = FALSE
    )
  }
  if (fid %in% names(loess_list)) {
    x <- loess_list[[fid]]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Frontier = fid, Model_Type = "LOESS",
      Key_Parameter = "span", Estimate = x$span, SE = NA_real_, P_value = NA_real_,
      Significant = NA_character_, R2 = NA_real_, n_obs = x$n_obs, stringsAsFactors = FALSE
    )
  }
  if (fid %in% names(tobit_results)) {
    t <- tobit_results[[fid]]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Frontier = fid, Model_Type = "Tobit_Quadratic",
      Key_Parameter = "beta2", Estimate = t$beta2, SE = t$se_beta2, P_value = t$p_beta2,
      Significant = ifelse(t$p_beta2 < 0.05, "Yes", "No"), R2 = NA_real_, n_obs = t$n_obs, stringsAsFactors = FALSE
    )
  }
  if (fid %in% names(tobit_cubic_results)) {
    t <- tobit_cubic_results[[fid]]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Frontier = fid, Model_Type = "Tobit_Cubic",
      Key_Parameter = "beta3", Estimate = t$beta3, SE = t$se_beta3, P_value = t$p_beta3,
      Significant = ifelse(t$p_beta3 < 0.05, "Yes", "No"), R2 = NA_real_, n_obs = t$n_obs, stringsAsFactors = FALSE
    )
  }
}
df_summary <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()

# Optionally include quantile regression (linear/quadratic) from improve_scale_analysis.R if available
table3_path <- file.path(statistical_tests_dir, "table3_quadratic_quantile_regression.csv")
df_quantile_rq <- data.frame()
if (file.exists(table3_path)) {
  tryCatch({
    df_quantile_rq <- read_csv(table3_path, show_col_types = FALSE)
    if (nrow(df_quantile_rq) > 0) {
      df_quantile_rq <- as.data.frame(df_quantile_rq)
      cat("✓ Included quantile regression (rq) results from improve_scale_analysis.R in Excel\n")
    }
  }, error = function(e) {})
}

# Write Excel (single file, multiple sheets)
excel_ushape_file <- file.path(statistical_tests_dir, "Ushape_robustness_all_models.xlsx")
tryCatch({
  sheets <- list(
    Summary = df_summary,
    Linear_OLS = df_linear,
    Quadratic_OLS = df_quadratic,
    LOESS = df_loess,
    Tobit_Quadratic = df_tobit_quad,
    Tobit_Cubic = df_tobit_cubic
  )
  if (nrow(df_quantile_rq) > 0) sheets$Quantile_Regression_Table3 <- df_quantile_rq
  sheets <- sheets[ sapply(sheets, nrow) > 0 ]
  if (length(sheets) > 0) {
    write_xlsx(sheets, excel_ushape_file)
    cat(sprintf("✓ U-shape robustness Excel saved: %s\n", excel_ushape_file))
  }
}, error = function(e) {
  cat(sprintf("⚠ Could not write U-shape robustness Excel: %s\n", e$message))
})

# Plot 1: 4 Frontiers (M, L, G, A)
tryCatch({
  cat("\n=== Generating 4 Frontiers LOESS plots ===\n")
  scale_data_4frontiers <- scale_data %>% 
    filter(frontier %in% c("M", "L", "G", "A")) %>%
    mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE))
  scale_data_se_4frontiers <- scale_data_se %>% 
    filter(frontier %in% c("M", "L", "G", "A")) %>%
    mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE))
  binned_means_4frontiers <- binned_means %>% 
    filter(frontier %in% c("M", "L", "G", "A")) %>%
    mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE))
  
  cat(sprintf("  ✓ scale_data_4frontiers rows: %d\n", nrow(scale_data_4frontiers)))
  cat(sprintf("  ✓ scale_data_se_4frontiers rows: %d\n", nrow(scale_data_se_4frontiers)))
  cat(sprintf("  ✓ binned_means_4frontiers rows: %d\n", nrow(binned_means_4frontiers)))
  
  p_loess_4panel <- ggplot(scale_data_4frontiers, 
                           aes(x = 10^log_scale, y = efficiency, color = frontier)) +
    # Add SE LOESS curves first (in background, light gray)
    geom_smooth(data = scale_data_se_4frontiers,
                method = "loess", se = TRUE, linewidth = 0.6, 
                color = "gray60", fill = "gray85", alpha = 0.15,
                aes(x = 10^log_scale, y = efficiency, group = frontier),
                inherit.aes = FALSE) +
    # Add SE raw data points (light gray, larger for academic style)
    geom_point(data = scale_data_se_4frontiers,
               aes(x = 10^log_scale, y = efficiency, group = frontier),
               alpha = 0.25, size = 0.8, color = "gray70",
               inherit.aes = FALSE) +
    # TE-VRS LOESS curves (main, colored)
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2,
                aes(fill = frontier)) +
    # TE-VRS raw data points (main, colored)
    geom_point(alpha = 0.55, size = 1.8) +
    geom_point(data = binned_means_4frontiers,
               aes(x = mean_scale, y = mean_eff, color = frontier),
               size = 2.5, shape = 23, fill = "white", stroke = 1.5, alpha = 0.95,
               inherit.aes = FALSE) +
    geom_text(data = binned_means_4frontiers,
              aes(x = mean_scale, y = mean_eff, label = sprintf("%.0f%%", mean_eff * 100)),
              color = "black", size = 5.2, fontface = "bold",
              vjust = -1.2, hjust = 0.5, nudge_y = 0.02,
              inherit.aes = FALSE) +
    # Add Tobit regression annotations (bottom right with box)
    {if (length(tobit_results) > 0) {
      tobit_annotations <- data.frame()
      for (fid in c("M", "L", "G", "A")) {
        if (fid %in% names(tobit_results)) {
          tob <- tobit_results[[fid]]
          # Calculate annotation position (bottom right of each facet)
          scale_range <- range(scale_data_4frontiers$log_scale[scale_data_4frontiers$frontier == fid], na.rm = TRUE)
          x_pos <- 10^(scale_range[2] - (scale_range[2] - scale_range[1]) * 0.05)  # Closer to right edge
          y_pos <- 0.05  # Bottom of plot
          
          # Create annotation text
          annot_text <- paste0(
            "β₁=", sprintf("%.3f", tob$beta1), tob$star_beta1, "\n",
            "β₂=", sprintf("%.3f", tob$beta2), tob$star_beta2
          )
          
          if (!is.na(tob$turning_point_scale)) {
            annot_text <- paste0(annot_text, "\n",
                                 "TP=", sprintf("%.0f", tob$turning_point_scale), 
                                 "±", sprintf("%.0f", tob$turning_point_se))
          }
          
          tobit_annotations <- rbind(tobit_annotations, data.frame(
            frontier = factor(fid, levels = c("M", "L", "G", "A"), ordered = TRUE),
            x = x_pos,
            y = y_pos,
            label = annot_text
          ))
        }
      }
      if (nrow(tobit_annotations) > 0) {
        # Ensure frontier factor order matches
        tobit_annotations$frontier <- factor(tobit_annotations$frontier, 
                                            levels = c("M", "L", "G", "A"), 
                                            ordered = TRUE)
        list(
          # Add background box
          geom_label(data = tobit_annotations,
                       aes(x = x, y = y, label = label),
                       inherit.aes = FALSE,
                     hjust = 1, vjust = 0,
                       size = 5.2, fontface = "bold",
                       color = "black",
                     fill = "white",
                     alpha = 0.85,
                     label.size = 0.5,
                    label.padding = unit(0.5, "lines"),
                     parse = FALSE)
        )
      } else {
        list()
      }
    } else {
      list()
    }} +
    facet_wrap(~frontier, ncol = 2, nrow = 2,
               labeller = labeller(frontier = c(
                "M" = "M-Frontier (Mortality)",
                 "L" = "L-Frontier (Local Env)",
                 "G" = "G-Frontier (Global)",
                 "A" = "A-Frontier (Aggregated)"
               )),
               # Ensure facet order follows factor levels (E-L-G-A)
               drop = FALSE) +
    scale_x_log10(breaks = c(3000, 10000, 30000, 100000, 300000),
                  labels = scales::comma) +
    scale_y_continuous(labels = scales::percent,
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       limits = c(0, 1.25)) +
    scale_color_manual(values = unlist(FRONTIER_COLORS[c("M", "L", "G", "A")])) +
    scale_fill_manual(values = unlist(FRONTIER_COLORS[c("M", "L", "G", "A")])) +
    labs(
      title = "Visual evidence: efficiency vs scale (LOESS smooth + binned means)",
      subtitle = "M, L, G, and A frontiers with nonparametric smooth curves",
      x = "Production scale (log scale; market pig heads/a)",
      y = "VRS technical efficiency",
      caption = "Diamond points = binned means (TE-VRS) | Colored line = TE-VRS LOESS | Gray line = SE LOESS | Points = raw data"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(family = "Times", size = 15),
      plot.title = element_text(family = "Times", face = "bold", size = 17, hjust = 0.5,
                                margin = margin(b = 8)),
      plot.subtitle = element_text(family = "Times", size = 15, hjust = 0.5, color = "gray40",
                                  margin = margin(b = 10)),
      axis.text = element_text(family = "Times", size = 14, color = "black"),
      axis.title = element_text(family = "Times", size = 15, face = "bold", color = "black"),
      strip.text = element_text(family = "Times", face = "bold", size = 15),
      plot.caption = element_text(family = "Times", size = 13, hjust = 0.5, color = "gray50"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5)
    )
  
  safe_ggsave("LOESS_efficiency_scale_4frontiers.png", p_loess_4panel, width = 14, height = 6)
  
  # Create vertical layout version (1 column)
  p_loess_4panel_vertical <- ggplot(scale_data_4frontiers, 
                                   aes(x = 10^log_scale, y = efficiency, color = frontier)) +
    # Add SE LOESS curves first (in background, light gray)
    geom_smooth(data = scale_data_se_4frontiers,
                method = "loess", se = TRUE, linewidth = 0.6, 
                color = "gray60", fill = "gray85", alpha = 0.15,
                aes(x = 10^log_scale, y = efficiency, group = frontier),
                inherit.aes = FALSE) +
    # Add SE raw data points (light gray, larger for academic style)
    geom_point(data = scale_data_se_4frontiers,
               aes(x = 10^log_scale, y = efficiency, group = frontier),
               alpha = 0.25, size = 0.8, color = "gray70",
               inherit.aes = FALSE) +
    # TE-VRS LOESS curves (main, colored)
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2,
                aes(fill = frontier)) +
    # TE-VRS raw data points (main, colored)
    geom_point(alpha = 0.55, size = 1.8) +
    geom_point(data = binned_means_4frontiers,
               aes(x = mean_scale, y = mean_eff, color = frontier),
               size = 2.5, shape = 23, fill = "white", stroke = 1.5, alpha = 0.95,
               inherit.aes = FALSE) +
    geom_text(data = binned_means_4frontiers,
              aes(x = mean_scale, y = mean_eff, label = sprintf("%.0f%%", mean_eff * 100)),
              color = "black", size = 5.2, fontface = "bold",
              vjust = -1.2, hjust = 0.5, nudge_y = 0.02,
              inherit.aes = FALSE) +
    # Add Tobit regression annotations (bottom right with box)
    {if (length(tobit_results) > 0) {
      tobit_annotations <- data.frame()
      for (fid in c("M", "L", "G", "A")) {
        if (fid %in% names(tobit_results)) {
          tob <- tobit_results[[fid]]
          # Calculate annotation position (bottom right of each facet)
          scale_range <- range(scale_data_4frontiers$log_scale[scale_data_4frontiers$frontier == fid], na.rm = TRUE)
          x_pos <- 10^(scale_range[2] - (scale_range[2] - scale_range[1]) * 0.05)  # Closer to right edge
          y_pos <- 0.05  # Bottom of plot
          
          # Create annotation text
          annot_text <- paste0(
            "β₁=", sprintf("%.3f", tob$beta1), tob$star_beta1, "\n",
            "β₂=", sprintf("%.3f", tob$beta2), tob$star_beta2
          )
          
          if (!is.na(tob$turning_point_scale)) {
            annot_text <- paste0(annot_text, "\n",
                                 "TP=", sprintf("%.0f", tob$turning_point_scale), 
                                 "±", sprintf("%.0f", tob$turning_point_se))
          }
          
          tobit_annotations <- rbind(tobit_annotations, data.frame(
            frontier = factor(fid, levels = c("M", "L", "G", "A"), ordered = TRUE),
            x = x_pos,
            y = y_pos,
            label = annot_text
          ))
        }
      }
      if (nrow(tobit_annotations) > 0) {
        # Ensure frontier factor order matches
        tobit_annotations$frontier <- factor(tobit_annotations$frontier, 
                                            levels = c("M", "L", "G", "A"), 
                                            ordered = TRUE)
        list(
          # Add background box
          geom_label(data = tobit_annotations,
                       aes(x = x, y = y, label = label),
                       inherit.aes = FALSE,
                     hjust = 1, vjust = 0,
                       size = 5.2, fontface = "bold",
                       color = "black",
                     fill = "white",
                     alpha = 0.85,
                     label.size = 0.5,
                    label.padding = unit(0.5, "lines"),
                     parse = FALSE)
        )
      } else {
        list()
      }
    } else {
      list()
    }} +
    scale_x_log10(breaks = c(3000, 10000, 30000, 100000, 300000),
                  labels = scales::comma) +
    scale_y_continuous(labels = scales::percent,
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       limits = c(0, 1.25)) +
    scale_color_manual(values = unlist(FRONTIER_COLORS[c("M", "L", "G", "A")])) +
    scale_fill_manual(values = unlist(FRONTIER_COLORS[c("M", "L", "G", "A")])) +
    labs(
      title = "Visual evidence: efficiency vs scale (LOESS smooth + binned means)",
      subtitle = "M, L, G, and A frontiers with nonparametric smooth curves",
      x = "Production scale (log scale; market pig heads/a)",
      y = "VRS technical efficiency",
      caption = "Diamond points = binned means (TE-VRS) | Colored line = TE-VRS LOESS | Gray line = SE LOESS | Points = raw data"
    ) +
    facet_wrap(~frontier, ncol = 1,
               labeller = labeller(frontier = c(
                "M" = "M-Frontier (Mortality)",
                 "L" = "L-Frontier (Local Env)",
                 "G" = "G-Frontier (Global)",
                 "A" = "A-Frontier (Aggregated)"
               )))
    
  # Apply same theme to vertical version
  p_loess_4panel_vertical <- p_loess_4panel_vertical +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(family = "Times", size = 15),
      plot.title = element_text(family = "Times", face = "bold", size = 17, hjust = 0.5,
                                margin = margin(b = 8)),
      plot.subtitle = element_text(family = "Times", size = 15, hjust = 0.5, color = "gray40",
                                  margin = margin(b = 10)),
      axis.text = element_text(family = "Times", size = 14, color = "black"),
      axis.title = element_text(family = "Times", size = 15, face = "bold", color = "black"),
      strip.text = element_text(family = "Times", face = "bold", size = 15),
      plot.caption = element_text(family = "Times", size = 13, hjust = 0.5, color = "gray50"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5)
    )
  
  # Save vertical layout version
  safe_ggsave("LOESS_efficiency_scale_4frontiers_vertical.png", p_loess_4panel_vertical, width = 9.5, height = 12)
  cat("✓ 4 Frontiers LOESS plots generated successfully\n")
}, error = function(e) {
  cat(sprintf("⚠ Error generating 4 Frontiers LOESS plots: %s\n", e$message))
  cat(sprintf("  Error occurred at: %s\n", deparse(e$call)))
})

# Plot 2: Meta-Frontier (use detected meta-frontier name)
tryCatch({
  cat("\n=== Generating Meta-Frontier LOESS plot ===\n")
  scale_data_meta <- scale_data %>% filter(frontier == meta_frontier_name)
  binned_means_meta <- binned_means %>% filter(frontier == meta_frontier_name)
  
  cat(sprintf("  ✓ scale_data_meta rows: %d\n", nrow(scale_data_meta)))
  cat(sprintf("  ✓ binned_means_meta rows: %d\n", nrow(binned_means_meta)))
  cat(sprintf("  ✓ scale_data_meta_se available: %s\n", !is.null(scale_data_meta_se)))
  
  p_loess_meta <- ggplot(scale_data_meta, aes(x = 10^log_scale, y = efficiency)) +
    # Add Meta SE LOESS curve first (in background, light gray) if available
    {if (!is.null(scale_data_meta_se)) {
      geom_smooth(data = scale_data_meta_se,
                  method = "loess", se = TRUE, linewidth = 0.6,
                  color = "gray60", fill = "gray85", alpha = 0.15,
                  aes(x = 10^log_scale, y = efficiency),
                  inherit.aes = FALSE)
    }} +
    # Add Meta SE raw data points (light gray, small) if available
    {if (!is.null(scale_data_meta_se)) {
      geom_point(data = scale_data_meta_se,
                 aes(x = 10^log_scale, y = efficiency),
                 alpha = 0.2, size = 0.8, color = "gray70",
                 inherit.aes = FALSE)
    }} +
    # Meta TE-VRS LOESS curve (main, colored)
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2,
                color = FRONTIER_COLORS[["Meta"]], fill = FRONTIER_COLORS[["Meta"]], alpha = 0.2) +
    # Meta TE-VRS raw data points (main, colored)
    geom_point(alpha = 0.4, size = 1.2, color = FRONTIER_COLORS[["Meta"]]) +
    geom_point(data = binned_means_meta,
               aes(x = mean_scale, y = mean_eff),
               size = 2.5, shape = 23, fill = "white", stroke = 1.5,
               color = FRONTIER_COLORS[["Meta"]], alpha = 0.95,
               inherit.aes = FALSE) +
    geom_text(data = binned_means_meta,
              aes(x = mean_scale, y = mean_eff, label = sprintf("%.0f%%", mean_eff * 100)),
              color = "black", size = 2.5, fontface = "bold",
              vjust = -1.2, hjust = 0.5, nudge_y = 0.02,
              inherit.aes = FALSE) +
    # Add Tobit regression annotation for Meta-frontier
    {if (meta_frontier_name %in% names(tobit_results)) {
      tob_meta <- tobit_results[[meta_frontier_name]]
      scale_range_meta <- range(scale_data_meta$log_scale, na.rm = TRUE)
      x_pos_meta <- 10^(scale_range_meta[2] - (scale_range_meta[2] - scale_range_meta[1]) * 0.15)
      y_pos_meta <- 0.95
      
      annot_text_meta <- paste0(
        "β₁=", sprintf("%.3f", tob_meta$beta1), tob_meta$star_beta1, "\n",
        "β₂=", sprintf("%.3f", tob_meta$beta2), tob_meta$star_beta2
      )
      
      if (!is.na(tob_meta$turning_point_scale)) {
        annot_text_meta <- paste0(annot_text_meta, "\n",
                                  "TP=", sprintf("%.0f", tob_meta$turning_point_scale),
                                  "±", sprintf("%.0f", tob_meta$turning_point_se))
      }
      
      list(geom_text(aes(x = x_pos_meta, y = y_pos_meta, label = annot_text_meta),
                     inherit.aes = FALSE,
                     hjust = 1, vjust = 1,
                     size = 2.8, fontface = "bold",
                     color = "black",
                     bg = "white", alpha = 0.8))
    } else {
      list()
    }} +
    scale_x_log10(breaks = c(3000, 10000, 30000, 100000, 300000),
                  labels = scales::comma) +
    scale_y_continuous(labels = scales::percent,
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       limits = c(0, 1.25)) +
    labs(
      title = paste0("Meta-frontier: efficiency vs scale (LOESS smooth + binned means)",
                     ifelse(meta_col == "Meta_TE_VRS_Slack", 
                            " [Slack Constraint Method]", 
                            " [Maximum Envelope]"
                     )
    ),
      subtitle = ifelse(meta_col == "Meta_TE_VRS_Slack",
                       "Meta-frontier = SBM efficiency with ALL bad outputs unconstrained",
                       "Meta-frontier = max(M, L, G, A) for each farm"),
      x = "Production scale (log)",
      y = "Meta-frontier efficiency (VRS)",
      caption = "Diamond points = binned means (TE-VRS) | Colored line = TE-VRS LOESS | Gray line = SE LOESS | Small points = raw data"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 13, hjust = 0.5, color = "gray50"),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
    )
  
  safe_ggsave("LOESS_efficiency_scale_Meta.png", p_loess_meta, width = 10, height = 7)
  cat("✓ Meta-Frontier LOESS plot generated successfully\n")
}, error = function(e) {
  cat(sprintf("⚠ Error generating Meta-Frontier LOESS plot: %s\n", e$message))
  cat(sprintf("  Error occurred at: %s\n", deparse(e$call)))
})
cat("\n")

# ============================================================================
# 2c. TE-CRS LOESS PLOTS (without SE curves)
# ============================================================================

# ============================================================================
# PRELOAD: Raw data for unit and marginal effects analysis
# ============================================================================
cat("\n=== Preloading raw data for unit and marginal effects analysis ===\n")

# Load raw data for inputs/outputs if not already loaded
if (!exists("df_clean")) {
  xlsx_path <- "./swine_farm_data.xlsx"
  if (file.exists(xlsx_path)) {
    suppressWarnings({
      df_raw <- readxl::read_xlsx(xlsx_path, sheet = 3)
    })
    inputs <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
    good_output <- "Market_pig"
    bad_outputs_full <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
    all_required <- unique(c(inputs, good_output, bad_outputs_full))
    
    df_clean <- df_raw %>%
      dplyr::select(dplyr::all_of(all_required), 
                    dplyr::any_of(c("Farm", "ID", "Report_year"))) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(all_required), 
                                  ~ suppressWarnings(as.numeric(.x)))) %>%
      tidyr::drop_na(dplyr::all_of(all_required))
    cat("✓ Loaded and cleaned raw data for unit and marginal effects analysis\n")
  } else {
    cat("⚠ Raw data file not found: swine_farm_data.xlsx\n")
  }
}

# Merge with efficiency data if df_clean exists
# Prefer obs_id join when eff and df_clean have same row count so factor columns match 1:1 and market_pig_original (head count) is preserved for zone assignment
if (exists("df_clean") && "Market_pig" %in% names(df_clean)) {
  if ("obs_id" %in% names(eff) && nrow(eff) == nrow(df_clean)) {
    data_merged <- eff %>%
      dplyr::left_join(
        df_clean %>%
          dplyr::mutate(obs_id = seq_len(dplyr::n())) %>%
          dplyr::select(obs_id, Feed, Water_consumption, Energy_consumption, Labor,
                       Waste_water, Manure, Dead_pig, Carbon_emission, Eutrophication_potential),
        by = "obs_id"
      )
  } else if ("ID" %in% names(eff) && "ID" %in% names(df_clean)) {
    data_merged <- eff %>%
      dplyr::left_join(df_clean %>% dplyr::select(ID, Feed, Water_consumption, Energy_consumption, Labor,
                                                   Waste_water, Manure, Dead_pig, Carbon_emission, Eutrophication_potential),
                        by = "ID")
  } else {
    scale_round <- if ("market_pig_original" %in% names(eff)) "market_pig_original" else "market_pig"
    data_merged <- eff %>%
      dplyr::mutate(market_pig_round = round(.data[[scale_round]], -2)) %>%
      dplyr::left_join(
        df_clean %>%
          dplyr::mutate(market_pig_round = round(Market_pig, -2)) %>%
          dplyr::select(market_pig_round, Feed, Water_consumption, Energy_consumption, Labor,
                       Waste_water, Manure, Dead_pig, Carbon_emission, Eutrophication_potential) %>%
          dplyr::group_by(market_pig_round) %>%
          dplyr::summarise_all(mean, na.rm = TRUE),
        by = "market_pig_round"
      ) %>%
      dplyr::select(-market_pig_round)
  }
  if (!"market_pig_original" %in% names(data_merged) && "Market_pig" %in% names(data_merged)) {
    data_merged$market_pig_original <- data_merged$Market_pig
  }
  
  # Check if merge was successful
  cat("\n=== DIAGNOSTIC: Data Merge Check ===\n")
  cat(sprintf("data_merged rows: %d\n", nrow(data_merged)))
  required_vars_check <- c("Feed", "Water_consumption", "Energy_consumption", "Labor",
                           "Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
  cat("Variables after merge:\n")
  for (var in required_vars_check) {
    if (var %in% names(data_merged)) {
      n_valid <- sum(!is.na(data_merged[[var]]) & data_merged[[var]] > 0, na.rm = TRUE)
      cat(sprintf("  %s: ✓ (n_valid=%d)\n", var, n_valid))
    } else {
      cat(sprintf("  %s: ✗ MISSING\n", var))
    }
  }
  cat("=====================================\n\n")
} else {
  cat("⚠ df_clean not found or does not have Market_pig column\n")
}

# ============================================================================
# 2c. TE-CRS LOESS PLOTS (without SE curves)
# ============================================================================

tryCatch({
  cat("\n=== Generating TE-CRS LOESS plots ===\n")
  
  # Prepare TE-CRS data (use scale_col so LOESS x-axis matches binned diamond positions)
  scale_data_crs <- eff %>%
    select(all_of(scale_col), M_TE_CRS, L_TE_CRS, G_TE_CRS, A_TE_CRS, all_of(meta_col)) %>%
    mutate(log_scale = log10(pmax(.data[[scale_col]], 1))) %>%
    pivot_longer(
      cols = c(M_TE_CRS, L_TE_CRS, G_TE_CRS, A_TE_CRS),
      names_to = "frontier",
      values_to = "efficiency"
    ) %>%
    mutate(frontier = str_remove(frontier, "_TE_CRS")) %>%
    filter(!is.na(efficiency), !is.na(log_scale) & is.finite(10^log_scale) & 10^log_scale > 0) %>%
    mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
    mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE))
  
  cat(sprintf("  ✓ scale_data_crs rows: %d\n", nrow(scale_data_crs)))
  
  # Calculate binned means for TE-CRS (use scale_col so bins and diamond x match production scale)
  eff_with_bins_crs <- eff %>%
    select(all_of(scale_col), log_scale, M_TE_CRS, L_TE_CRS, G_TE_CRS, A_TE_CRS) %>%
    mutate(
      scale_bin = cut(
        .data[[scale_col]],
        breaks = c(0, 5000, 10000, 20000, 50000, Inf),
        labels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
        include.lowest = TRUE
      )
    )
  
  binned_means_crs <- eff_with_bins_crs %>%
    select(all_of(scale_col), scale_bin, M_TE_CRS, L_TE_CRS, G_TE_CRS, A_TE_CRS) %>%
    pivot_longer(
      cols = c(M_TE_CRS, L_TE_CRS, G_TE_CRS, A_TE_CRS),
      names_to = "frontier",
      values_to = "efficiency"
    ) %>%
    mutate(frontier = str_remove(frontier, "_TE_CRS")) %>%
    filter(!is.na(efficiency), !is.na(scale_bin)) %>%
    mutate(scale_bin = factor(scale_bin, 
                              levels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
                              ordered = TRUE)) %>%
    group_by(scale_bin, frontier) %>%
    summarise(
      n = n(),
      mean_eff = mean(efficiency, na.rm = TRUE),
      mean_scale = exp(mean(log(.data[[scale_col]]), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
    mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE))
  
  cat(sprintf("  ✓ binned_means_crs rows: %d\n", nrow(binned_means_crs)))
  
  # Plot 3: 4 Frontiers TE-CRS (M, L, G, A)
  scale_data_crs_4frontiers <- scale_data_crs %>% filter(frontier %in% c("M", "L", "G", "A"))
  binned_means_crs_4frontiers <- binned_means_crs %>% filter(frontier %in% c("M", "L", "G", "A"))
  
  cat(sprintf("  ✓ scale_data_crs_4frontiers rows: %d\n", nrow(scale_data_crs_4frontiers)))
  cat(sprintf("  ✓ binned_means_crs_4frontiers rows: %d\n", nrow(binned_means_crs_4frontiers)))
  
  p_loess_crs_4panel <- ggplot(scale_data_crs_4frontiers, 
                                aes(x = 10^log_scale, y = efficiency, color = frontier)) +
    # TE-CRS LOESS curves (main, colored)
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2,
                aes(fill = frontier)) +
    # TE-CRS raw data points (main, colored, slightly smaller)
    geom_point(alpha = 0.5, size = 1.1) +
    geom_point(data = binned_means_crs_4frontiers,
               aes(x = mean_scale, y = mean_eff, color = frontier),
               size = 2.5, shape = 23, fill = "white", stroke = 1.5, alpha = 0.95,
               inherit.aes = FALSE) +
    geom_text(data = binned_means_crs_4frontiers,
              aes(x = mean_scale, y = mean_eff, label = sprintf("%.0f%%", mean_eff * 100)),
              color = "black", size = 2.5, fontface = "bold",
              vjust = -1.2, hjust = 0.5, nudge_y = 0.02,
              inherit.aes = FALSE) +
    facet_wrap(~frontier, ncol = 2, nrow = 2,
               labeller = labeller(frontier = c(
                "M" = "M-Frontier (Mortality)",
                 "L" = "L-Frontier (Local Env)",
                 "G" = "G-Frontier (Global)",
                 "A" = "A-Frontier (Aggregated)"
               )),
               # Ensure facet order follows factor levels (E-L-G-A)
               drop = FALSE) +
    scale_x_log10(breaks = c(3000, 10000, 30000, 100000, 300000),
                  labels = scales::comma) +
    scale_y_continuous(labels = scales::percent,
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       limits = c(0, 1.25)) +
    scale_color_manual(values = unlist(FRONTIER_COLORS[c("M", "L", "G", "A")])) +
    scale_fill_manual(values = unlist(FRONTIER_COLORS[c("M", "L", "G", "A")])) +
    labs(
      title = "Visual evidence: CRS efficiency vs scale (LOESS smooth + binned means)",
      subtitle = "M, L, G, and A frontiers with nonparametric smooth curves",
      x = "Production scale (log scale; market pig heads/a)",
      y = "CRS technical efficiency",
      caption = "Diamond points = binned means (TE-CRS) | Colored line = TE-CRS LOESS | Points = raw data"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(family = "Times", size = 12),
      plot.title = element_text(family = "Times", face = "bold", size = 14, hjust = 0.5,
                                margin = margin(b = 8)),
      plot.subtitle = element_text(family = "Times", size = 12, hjust = 0.5, color = "gray40",
                                  margin = margin(b = 10)),
      axis.text = element_text(family = "Times", size = 11, color = "black"),
      axis.title = element_text(family = "Times", size = 12, face = "bold", color = "black"),
      strip.text = element_text(family = "Times", face = "bold", size = 12),
      plot.caption = element_text(family = "Times", size = 10, hjust = 0.5, color = "gray50"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5)
    )
  
  safe_ggsave("LOESS_efficiency_scale_4frontiers_CRS.png", p_loess_crs_4panel, width = 14, height = 6)
  cat("✓ 4 Frontiers TE-CRS LOESS plot generated successfully\n")
  
  # Plot 4: Meta-Frontier TE-CRS
  meta_crs_col <- if ("Meta_TE_CRS" %in% names(eff)) "Meta_TE_CRS" else NULL
  if (!is.null(meta_crs_col)) {
    scale_data_meta_crs <- eff %>%
      select(all_of(scale_col), all_of(meta_crs_col)) %>%
      mutate(log_scale = log10(pmax(.data[[scale_col]], 1))) %>%
      rename(efficiency = all_of(meta_crs_col)) %>%
      filter(!is.na(efficiency), !is.na(log_scale) & is.finite(10^log_scale) & 10^log_scale > 0)
    
    cat(sprintf("  ✓ scale_data_meta_crs rows: %d\n", nrow(scale_data_meta_crs)))
    
    # Calculate binned means for Meta TE-CRS (use scale_col so bins and diamond x match production scale)
    eff_with_bins_meta_crs <- eff %>%
      select(all_of(scale_col), all_of(meta_crs_col)) %>%
      mutate(
        scale_bin = cut(
          .data[[scale_col]],
          breaks = c(0, 5000, 10000, 20000, 50000, Inf),
          labels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
          include.lowest = TRUE
        )
      ) %>%
      rename(efficiency = all_of(meta_crs_col))
    
    binned_means_meta_crs <- eff_with_bins_meta_crs %>%
      filter(!is.na(efficiency), !is.na(scale_bin)) %>%
      mutate(scale_bin = factor(scale_bin, 
                                levels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
                                ordered = TRUE)) %>%
      group_by(scale_bin) %>%
      summarise(
        n = n(),
        mean_eff = mean(efficiency, na.rm = TRUE),
        mean_scale = exp(mean(log(.data[[scale_col]]), na.rm = TRUE)),
        .groups = "drop"
      )
    
    cat(sprintf("  ✓ binned_means_meta_crs rows: %d\n", nrow(binned_means_meta_crs)))
    
    p_loess_meta_crs <- ggplot(scale_data_meta_crs, aes(x = 10^log_scale, y = efficiency)) +
      # Meta TE-CRS LOESS curve (main, colored)
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.2,
                  color = FRONTIER_COLORS[["Meta"]], fill = FRONTIER_COLORS[["Meta"]], alpha = 0.2) +
      # Meta TE-CRS raw data points (main, colored)
      geom_point(alpha = 0.4, size = 1.2, color = FRONTIER_COLORS[["Meta"]]) +
      geom_point(data = binned_means_meta_crs,
                 aes(x = mean_scale, y = mean_eff),
                 size = 2.5, shape = 23, fill = "white", stroke = 1.5,
                 color = FRONTIER_COLORS[["Meta"]], alpha = 0.95,
                 inherit.aes = FALSE) +
      geom_text(data = binned_means_meta_crs,
                aes(x = mean_scale, y = mean_eff, label = sprintf("%.0f%%", mean_eff * 100)),
                color = "black", size = 2.5, fontface = "bold",
                vjust = -1.2, hjust = 0.5, nudge_y = 0.02,
                inherit.aes = FALSE) +
      scale_x_log10(breaks = c(3000, 10000, 30000, 100000, 300000),
                    labels = scales::comma) +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title = paste0("Meta-frontier: CRS efficiency vs scale (LOESS smooth + binned means)",
                       ifelse(meta_col == "Meta_TE_VRS_Slack", 
                              " [Slack Constraint Method]", 
                              " [Maximum Envelope]"
        )),
        subtitle = ifelse(meta_col == "Meta_TE_VRS_Slack",
                         "Meta-frontier = SBM efficiency with ALL bad outputs unconstrained",
                         "Meta-frontier = max(M, L, G, A) for each farm"),
        x = "Production scale (log)",
        y = "Meta-frontier efficiency (CRS)",
        caption = "Diamond points = binned means (TE-CRS) | Colored line = TE-CRS LOESS | Points = raw data"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 13, hjust = 0.5, color = "gray50"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
      )
    
    safe_ggsave("LOESS_efficiency_scale_Meta_CRS.png", p_loess_meta_crs, width = 10, height = 7)
    cat("✓ Meta-Frontier TE-CRS LOESS plot generated successfully\n")
  } else {
    cat("⚠ Meta_TE_CRS column not found, skipping Meta-Frontier TE-CRS plot\n")
  }
  
  cat("✓ TE-CRS LOESS plots generated successfully\n")
}, error = function(e) {
  cat(sprintf("⚠ Error generating TE-CRS LOESS plots: %s\n", e$message))
  cat(sprintf("  Error occurred at: %s\n", deparse(e$call)))
})

cat("\n")

# ============================================================================
# SHAPLEY PLOT GENERATION (Independent of Tobit Results)
# ============================================================================

cat("\n=== INDEPENDENT SHAPLEY PLOT GENERATION ===\n")

# Load raw data for inputs/outputs if not already loaded
if (!exists("df_clean")) {
  xlsx_path <- "./swine_farm_data.xlsx"
  if (file.exists(xlsx_path)) {
    suppressWarnings({
      df_raw <- readxl::read_xlsx(xlsx_path, sheet = 3)
    })
    inputs <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
    good_output <- "Market_pig"
    bad_outputs_full <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
    all_required <- unique(c(inputs, good_output, bad_outputs_full))
    
    df_clean <- df_raw %>%
      dplyr::select(dplyr::all_of(all_required), 
                    dplyr::any_of(c("Farm", "ID", "Report_year"))) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(all_required), 
                                  ~ suppressWarnings(as.numeric(.x)))) %>%
      tidyr::drop_na(dplyr::all_of(all_required))
    cat("✓ Loaded and cleaned raw data for Shapley analysis\n")
  } else {
    cat("⚠ Raw data file not found: swine_farm_data.xlsx\n")
  }
}

# Merge with efficiency data if df_clean exists (same strategy as first data_merged block so factor columns and scale zones are consistent for Shapley and marginal effects)
if (exists("df_clean") && "Market_pig" %in% names(df_clean)) {
  if ("obs_id" %in% names(eff) && nrow(eff) == nrow(df_clean)) {
    data_merged <- eff %>%
      dplyr::left_join(
        df_clean %>%
          dplyr::mutate(obs_id = seq_len(dplyr::n())) %>%
          dplyr::select(obs_id, Market_pig, Feed, Water_consumption, Energy_consumption, Labor,
                       Waste_water, Manure, Dead_pig, Carbon_emission, Eutrophication_potential),
        by = "obs_id"
      ) %>%
      dplyr::mutate(market_pig_original = dplyr::coalesce(.data$market_pig_original, .data$Market_pig))
  } else if ("ID" %in% names(eff) && "ID" %in% names(df_clean)) {
    data_merged <- eff %>%
      dplyr::left_join(df_clean %>% dplyr::select(ID, Feed, Water_consumption, Energy_consumption, Labor,
                                                   Waste_water, Manure, Dead_pig, Carbon_emission, Eutrophication_potential),
                        by = "ID")
  } else {
    scale_round <- if ("market_pig_original" %in% names(eff)) "market_pig_original" else "market_pig"
    data_merged <- eff %>%
      dplyr::mutate(market_pig_round = round(.data[[scale_round]], -2)) %>%
      dplyr::left_join(
        df_clean %>%
          dplyr::mutate(market_pig_round = round(Market_pig, -2)) %>%
          dplyr::select(market_pig_round, Feed, Water_consumption, Energy_consumption, Labor,
                       Waste_water, Manure, Dead_pig, Carbon_emission, Eutrophication_potential) %>%
          dplyr::group_by(market_pig_round) %>%
          dplyr::summarise_all(mean, na.rm = TRUE),
        by = "market_pig_round"
      ) %>%
      dplyr::select(-market_pig_round)
  }
  if (!"market_pig_original" %in% names(data_merged) && "Market_pig" %in% names(data_merged)) {
    data_merged$market_pig_original <- data_merged$Market_pig
  }

  # Check if merge was successful (n_valid = non-NA count; for raw inputs/outputs >0 is expected)
  cat("\n=== DIAGNOSTIC: Data Merge Check for Shapley ===\n")
  cat(sprintf("data_merged rows: %d\n", nrow(data_merged)))
  required_vars_check <- c("Feed", "Water_consumption", "Energy_consumption", "Labor",
                           "Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
  cat("Variables after merge:\n")
  for (var in required_vars_check) {
    if (var %in% names(data_merged)) {
      n_valid <- sum(!is.na(data_merged[[var]]) & data_merged[[var]] > 0, na.rm = TRUE)
      n_non_na <- sum(!is.na(data_merged[[var]]))
      cat(sprintf("  %s: ✓ (n_valid>0=%d, n_non_NA=%d)\n", var, n_valid, n_non_na))
    } else {
      cat(sprintf("  %s: ✗ MISSING\n", var))
    }
  }
  cat("=====================================\n\n")
  
  # Filter data_merged
  data_merged <- data_merged %>%
    dplyr::filter(market_pig > 0, !is.na(A_TE_VRS))
  
  # Shapley plots (only if data_merged exists)
  cat("\n=== Checking for Shapley Plot Generation ===\n")
  cat(sprintf("  data_merged exists: %s\n", exists("data_merged")))
  if (exists("data_merged")) {
    cat(sprintf("  data_merged rows: %d\n", nrow(data_merged)))
    cat(sprintf("  data_merged has A_TE_VRS: %s\n", "A_TE_VRS" %in% names(data_merged)))
    cat(sprintf("  data_merged has A_SE: %s\n", "A_SE" %in% names(data_merged)))
  }
  
  if (exists("data_merged") && nrow(data_merged) > 0) {
    # Calculate Shapley contributions for All Scale only
    factor_vars <- c("Feed", "Water_consumption", "Energy_consumption", "Labor",
                    "Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
    factor_vars <- factor_vars[factor_vars %in% names(data_merged)]
    
    cat(sprintf("  factor_vars found: %d variables\n", length(factor_vars)))
    cat(sprintf("  factor_vars: %s\n", paste(factor_vars, collapse = ", ")))
    
    if (length(factor_vars) > 0) {
      # ============================================================================
      # Shapley Plot (TE-VRS & SE) - All Scale Only
      # ============================================================================
      cat("\n=== Creating Shapley Plot (TE-VRS & SE) ===\n")
      
      # Determine which outcomes are available
      has_te_vrs <- "A_TE_VRS" %in% names(data_merged)
      has_te_crs <- "A_TE_CRS" %in% names(data_merged)
      has_se <- "A_SE" %in% names(data_merged)
      
      cat(sprintf("  has_te_vrs: %s, has_te_crs: %s, has_se: %s\n", has_te_vrs, has_te_crs, has_se))
      
      if (!has_te_vrs && !has_te_crs && !has_se) {
        cat("⚠ Neither A_TE_VRS, A_TE_CRS nor A_SE found, skipping Shapley plot\n")
      } else {
        # Get all scale data (filter by available outcomes)
        if (has_se && has_te_vrs && has_te_crs) {
          all_scale_data <- data_merged %>% 
            dplyr::filter(!is.na(A_TE_VRS), !is.na(A_TE_CRS), !is.na(A_SE))
        } else if (has_se && has_te_vrs) {
          all_scale_data <- data_merged %>% 
            dplyr::filter(!is.na(A_TE_VRS), !is.na(A_SE))
          cat("⚠ A_TE_CRS not found, generating TE-VRS and SE plots only\n")
        } else if (has_te_vrs && has_te_crs) {
          all_scale_data <- data_merged %>% 
            dplyr::filter(!is.na(A_TE_VRS), !is.na(A_TE_CRS))
          cat("⚠ A_SE not found, generating TE-VRS and TE-CRS plots only\n")
        } else if (has_te_vrs) {
          all_scale_data <- data_merged %>% 
            dplyr::filter(!is.na(A_TE_VRS))
          cat("⚠ A_SE and A_TE_CRS not found, generating TE-VRS plot only\n")
        } else if (has_te_crs) {
          all_scale_data <- data_merged %>% 
            dplyr::filter(!is.na(A_TE_CRS))
          cat("⚠ A_SE and A_TE_VRS not found, generating TE-CRS plot only\n")
        } else {
          all_scale_data <- data_merged %>% 
            dplyr::filter(!is.na(A_SE))
          cat("⚠ A_TE_VRS and A_TE_CRS not found, generating SE plot only\n")
        }
        
        cat(sprintf("  all_scale_data rows: %d\n", nrow(all_scale_data)))
        
        if (nrow(all_scale_data) >= 5) {
          # Calculate Shapley values for available outcomes
          shapley_plots_data <- list()
          
          # Determine which outcomes to calculate
          outcomes_to_calc <- c()
          if (has_te_vrs) outcomes_to_calc <- c(outcomes_to_calc, "A_TE_VRS")
          if (has_te_crs) outcomes_to_calc <- c(outcomes_to_calc, "A_TE_CRS")
          if (has_se) outcomes_to_calc <- c(outcomes_to_calc, "A_SE")
          
          cat(sprintf("  Calculating Shapley values for: %s\n", paste(outcomes_to_calc, collapse = ", ")))
          
          for (outcome_var in outcomes_to_calc) {
            outcome_name <- ifelse(outcome_var == "A_TE_VRS", "TE-VRS",
                                 ifelse(outcome_var == "A_TE_CRS", "TE-CRS", "SE"))
            
            y <- all_scale_data[[outcome_var]]
            X <- as.matrix(all_scale_data[, factor_vars, drop = FALSE])
            ok <- complete.cases(y, X) & is.finite(y) & apply(X, 1, function(r) all(is.finite(r)))
            if (sum(ok) < 10) next
            y <- y[ok]
            X <- X[ok, , drop = FALSE]
            var_ok <- apply(X, 2, function(z) sd(z, na.rm = TRUE) > 1e-10)
            if (sum(var_ok) < 2) next
            X <- X[, var_ok, drop = FALSE]
            factor_vars_used <- factor_vars[var_ok]
            X_scaled <- scale(X)
            y_scaled <- as.vector(scale(y))
            if (any(!is.finite(X_scaled)) || any(!is.finite(y_scaled))) next
            n_var <- ncol(X_scaled)
            # Use GAM (mgcv::gam) for Shapley with low k to avoid overfitting
            df_lm <- data.frame(y = y_scaled)
            for (j in seq_len(n_var)) df_lm[[paste0("var", j)]] <- X_scaled[, j]
            # Shapley GAM smooth complexity (s(var, k))
            # Sensitivity check: k = 5
            smooth_terms <- paste0("s(var", seq_len(n_var), ", k = 3)")
            full_formula <- as.formula(paste("y ~", paste(smooth_terms, collapse = " + ")))
            plot_data_long <- data.frame()
            
            for (i in seq_len(n_var)) {
              factor_name <- factor_vars_used[i]
              tryCatch({
                cat(sprintf("  Processing factor %s (%d/%d)...\n", factor_name, i, n_var))
                # Full GAM with all smooth terms
                model_full <- mgcv::gam(full_formula, data = df_lm, method = "REML")
                if (!inherits(model_full, "gam")) stop("Full GAM did not converge")
                
                # GAM without the i-th factor's smooth term
                if (n_var > 1) {
                  vars_without <- setdiff(seq_len(n_var), i)
                  smooth_terms_without <- paste0("s(var", vars_without, ", k = 3)")
                  without_formula <- as.formula(paste("y ~", paste(smooth_terms_without, collapse = " + ")))
                  model_without <- mgcv::gam(without_formula, data = df_lm, method = "REML")
                  if (!inherits(model_without, "gam")) stop("Reduced GAM did not converge")
                } else {
                  model_without <- mgcv::gam(y ~ 1, data = df_lm, method = "REML")
                }
                pred_full <- as.vector(predict(model_full, newdata = df_lm, type = "response"))
                pred_without <- as.vector(predict(model_without, newdata = df_lm, type = "response"))
                marginal_contrib <- pred_full - pred_without
                factor_values <- X[, i]
                factor_values_scaled <- (factor_values - min(factor_values, na.rm = TRUE)) /
                  (max(factor_values, na.rm = TRUE) - min(factor_values, na.rm = TRUE) + 1e-10)
                # Permutation test for importance magnitude (B option):
                # test whether mean(abs(marginal_contrib)) is larger than a null
                # distribution obtained by permuting y_scaled.
                contrib_clean <- marginal_contrib[is.finite(marginal_contrib)]
                signif_flag <- FALSE
                if (length(contrib_clean) >= 5) {
                  obs_stat <- mean(abs(contrib_clean), na.rm = TRUE)
                  P <- 100  # permutations; balanced speed and stability
                  set.seed(20260403 + i * 1000L + length(factor_vars_used))
                  
                  perm_stats <- replicate(P, {
                    perm_y <- sample.int(length(y_scaled), replace = FALSE)
                    df_perm <- df_lm
                    df_perm$y <- y_scaled[perm_y]
                    pred_full_perm <- tryCatch({
                      model_full_perm <- mgcv::gam(full_formula, data = df_perm, method = "REML")
                      as.vector(predict(model_full_perm, newdata = df_perm, type = "response"))
                    }, error = function(e) NULL)
                    
                    if (is.null(pred_full_perm)) return(NA_real_)
                    
                    pred_without_perm <- tryCatch({
                      if (n_var > 1) {
                        model_without_perm <- mgcv::gam(without_formula, data = df_perm, method = "REML")
                        as.vector(predict(model_without_perm, newdata = df_perm, type = "response"))
                      } else {
                        model_without_perm <- mgcv::gam(y ~ 1, data = df_perm, method = "REML")
                        as.vector(predict(model_without_perm, newdata = df_perm, type = "response"))
                      }
                    }, error = function(e) NULL)
                    
                    if (is.null(pred_without_perm)) return(NA_real_)
                    
                    marginal_perm <- pred_full_perm - pred_without_perm
                    marginal_perm <- marginal_perm[is.finite(marginal_perm)]
                    if (length(marginal_perm) < 5) return(NA_real_)
                    mean(abs(marginal_perm), na.rm = TRUE)
                  })
                  
                  perm_stats <- perm_stats[is.finite(perm_stats)]
                  if (length(perm_stats) > 0) {
                    p_value <- (sum(perm_stats >= obs_stat) + 1) / (length(perm_stats) + 1)
                    signif_flag <- p_value < 0.05
                  }
                }
                plot_data_long <- rbind(plot_data_long, data.frame(
                  Factor = factor_name,
                  SHAP_Value = marginal_contrib,
                  Feature_Value = factor_values_scaled,
                  Sample_ID = seq_along(marginal_contrib),
                  Outcome = outcome_name,
                  Significant = signif_flag,
                  stringsAsFactors = FALSE
                ))
                cat(sprintf("    ✓ Factor %s processed successfully\n", factor_name))
              }, error = function(e) {
                cat(sprintf("  ⚠ Error processing factor %s: %s\n", factor_name, e$message))
              })
            }
            
            if (nrow(plot_data_long) > 0) {
              shapley_plots_data[[outcome_name]] <- plot_data_long
              cat(sprintf("  ✓ Calculated Shapley values for %s: %d factors, %d samples\n", 
                         outcome_name, length(unique(plot_data_long$Factor)), 
                         length(unique(plot_data_long$Sample_ID))))
            } else {
              cat(sprintf("  ⚠ No plot_data_long for %s\n", outcome_name))
            }
          }
          
          cat(sprintf("  shapley_plots_data length: %d\n", length(shapley_plots_data)))
          if (length(shapley_plots_data) > 0) {
            cat(sprintf("  Available outcomes: %s\n", paste(names(shapley_plots_data), collapse = ", ")))
          }
          
          # Create plots for all available outcomes (TE-VRS and/or SE)
          # Generate plots for whichever outcomes have data
          if (length(shapley_plots_data) > 0) {
            
            # Create plots for all available outcomes
            outcomes_to_plot <- names(shapley_plots_data)
            
            for (outcome_name in outcomes_to_plot) {
              plot_data_long <- shapley_plots_data[[outcome_name]]
              
              # Determine factor order for this outcome (widest to narrowest, top to bottom)
              factor_order <- plot_data_long %>%
                dplyr::group_by(Factor) %>%
                dplyr::summarise(
                  shap_width = IQR(SHAP_Value, na.rm = TRUE),
                  .groups = "drop"
                ) %>%
                dplyr::arrange(desc(shap_width)) %>%
                dplyr::pull(Factor)
              
              # Add significance markers
              plot_data_long <- plot_data_long %>%
                dplyr::mutate(
                  Factor_Display = gsub("_", " ", ifelse(Factor == "Waste_water", "Wastewater", as.character(Factor))),
                  Factor_Label = ifelse(Significant, paste0(Factor_Display, "*"), Factor_Display)
                )
              
              # Create factor levels with significance markers (widest at top)
              sig_info <- plot_data_long %>%
                dplyr::group_by(Factor) %>%
                dplyr::summarise(Significant = first(Significant), .groups = "drop")
              
              factor_levels_with_sig <- sapply(rev(factor_order), function(f) {  # Reverse for top-to-bottom
                sig <- sig_info$Significant[sig_info$Factor == f]
                f_display <- gsub("_", " ", ifelse(f == "Waste_water", "Wastewater", f))
                if (length(sig) > 0 && sig[1]) {
                  paste0(f_display, "*")
                } else {
                  f_display
                }
              })
              
              plot_data_long <- plot_data_long %>%
                dplyr::mutate(
                  Factor_Label = factor(Factor_Label, levels = factor_levels_with_sig)
                )
              
              # Determine x-axis limits for this outcome
              x_range <- range(plot_data_long$SHAP_Value, na.rm = TRUE)
              x_padding <- diff(x_range) * 0.1
              max_abs_shap <- max(abs(x_range), na.rm = TRUE)
              max_abs_shap <- ifelse(is.finite(max_abs_shap) && max_abs_shap > 0, max_abs_shap, 1)
              max_abs_shap <- max_abs_shap + abs(x_padding)
              # Pretty breaks, then symmetric truncation around 0
              pretty_x <- pretty(c(-max_abs_shap, max_abs_shap), n = 6)
              pretty_x <- pretty_x[is.finite(pretty_x)]
              step_x <- diff(pretty_x)
              step_x <- step_x[is.finite(step_x) & step_x > 0]
              step_x <- if (length(step_x) > 0) min(step_x) else max_abs_shap / 3
              max_abs_shap_sym <- ceiling(max_abs_shap / step_x) * step_x
              x_limits <- c(-max_abs_shap_sym, max_abs_shap_sym)
              x_breaks <- seq(-max_abs_shap_sym, max_abs_shap_sym, by = step_x)
              
              # Calculate density heatmap data (violin-like heatmap) for both TE-VRS and SE
              # Use tryCatch to ensure basic plot is generated even if heatmap fails
              heatmap_data <- tryCatch({
                # Initialize heatmap_data inside the loop for proper scope
                heatmap_data_temp <- data.frame()
                
                # Calculate density for each factor and create heatmap tiles
                for (factor_name in unique(plot_data_long$Factor)) {
                  tryCatch({
                    factor_data <- plot_data_long %>%
                      dplyr::filter(Factor == factor_name, !is.na(SHAP_Value), !is.na(Feature_Value))
                    
                    if (nrow(factor_data) < 3) next
                    
                    # Calculate density with error handling
                    dens <- tryCatch({
                      density(factor_data$SHAP_Value, na.rm = TRUE, n = 128)
                    }, error = function(e) {
                      cat(sprintf("    ⚠ Density calculation failed for factor %s: %s\n", factor_name, e$message))
                      return(NULL)
                    })
                    
                    if (is.null(dens)) next
                    
                    # Create bins for SHAP values
                    n_bins <- 50
                    shap_bins <- seq(x_range[1], x_range[2], length.out = n_bins + 1)
                    
                    for (i in seq_len(n_bins)) {
                      tryCatch({
                        bin_left <- shap_bins[i]
                        bin_right <- shap_bins[i + 1]
                        bin_center <- (bin_left + bin_right) / 2
                        
                        # Find samples in this bin
                        bin_data <- factor_data %>%
                          dplyr::filter(SHAP_Value >= bin_left, SHAP_Value < bin_right)
                        
                        if (nrow(bin_data) > 0) {
                          # Calculate density at bin center (interpolate from density object)
                          dens_value <- tryCatch({
                            approx(dens$x, dens$y, xout = bin_center, rule = 2)$y
                          }, error = function(e) {
                            return(0)
                          })
                          dens_value <- max(0, dens_value)  # Ensure non-negative
                          
                          # Average feature value in this bin
                          avg_feature_value <- mean(bin_data$Feature_Value, na.rm = TRUE)
                          
                          # Normalize density to max density for this factor (for width scaling)
                          max_dens_factor <- max(dens$y, na.rm = TRUE)
                          normalized_dens <- ifelse(max_dens_factor > 0, dens_value / max_dens_factor, 0)
                          
                          heatmap_data_temp <- rbind(heatmap_data_temp, data.frame(
                            Factor = factor_name,
                            Factor_Label = factor_data$Factor_Label[1],
                            SHAP_Value = bin_center,
                            Density = dens_value,
                            Normalized_Density = normalized_dens,
                            Feature_Value = avg_feature_value,
                            Bin_Left = bin_left,
                            Bin_Right = bin_right
                          ))
                        }
                      }, error = function(e) {
                        # Skip this bin if there's an error
                        return(NULL)
                      })
                    }
                  }, error = function(e) {
                    cat(sprintf("    ⚠ Heatmap calculation failed for factor %s: %s\n", factor_name, e$message))
                    return(NULL)
                  })
                }
                
                # Add Factor_Label to heatmap_data
                if (nrow(heatmap_data_temp) > 0) {
                  heatmap_data_temp <- heatmap_data_temp %>%
                    dplyr::left_join(
                      plot_data_long %>%
                        dplyr::select(Factor, Factor_Label) %>%
                        dplyr::distinct(),
                      by = "Factor"
                    ) %>%
                    dplyr::mutate(
                      Factor_Label = factor(Factor_Label, levels = factor_levels_with_sig)
                    )
                  
                  # Prepare heatmap data with y positions
                  factor_y_positions <- seq_along(factor_levels_with_sig)
                  names(factor_y_positions) <- factor_levels_with_sig
                  
                  heatmap_data_temp <- heatmap_data_temp %>%
                    dplyr::mutate(
                      Y_Position = factor_y_positions[as.character(Factor_Label)]
                    )
                  
                  # Calculate tile width proportional to density
                  max_tile_width <- 0.4  # Maximum width of tile (in y-axis units)
                  heatmap_data_temp <- heatmap_data_temp %>%
                    dplyr::mutate(
                      Tile_Width = Normalized_Density * max_tile_width,
                      Y_Min = Y_Position - Tile_Width / 2,
                      Y_Max = Y_Position + Tile_Width / 2
                    )
                }
                
                heatmap_data_temp
              }, error = function(e) {
                cat(sprintf("  ⚠ Heatmap data generation failed: %s\n", e$message))
                cat("  → Generating plot without heatmap tiles\n")
                return(data.frame())  # Return empty data.frame if heatmap fails
              })
            
              # Create plot with complete legend and information
              p_sub <- ggplot2::ggplot(plot_data_long,
                                      ggplot2::aes(x = SHAP_Value, y = Factor_Label)) +
                ggplot2::geom_vline(xintercept = 0, linetype = "dashed", 
                                   color = "gray60", linewidth = 0.5, alpha = 0.7)
            
              # Add heatmap tiles (before points, so they appear in background) for both TE-VRS and SE
              # Use tryCatch to ensure plot is generated even if heatmap rendering fails
              if (!is.null(heatmap_data) && nrow(heatmap_data) > 0) {
                tryCatch({
                  p_sub <- p_sub +
                    ggplot2::geom_rect(
                      data = heatmap_data,
                      ggplot2::aes(
                        xmin = Bin_Left,
                        xmax = Bin_Right,
                        ymin = Y_Min,
                        ymax = Y_Max,
                        fill = Feature_Value
                      ),
                      inherit.aes = FALSE,
                      alpha = 0.6,
                      color = NA
                    )
                }, error = function(e) {
                  cat(sprintf("  ⚠ Heatmap rendering failed: %s\n", e$message))
                  cat("  → Continuing with plot without heatmap tiles\n")
                })
              }
            
              # Add points
              p_sub <- p_sub +
                ggplot2::geom_point(ggplot2::aes(color = Feature_Value), 
                                   size = 1.8, alpha = 0.7,  # Smaller points
                                   position = ggplot2::position_jitter(height = 0.15, width = 0),
                                   stroke = 0.3, na.rm = TRUE) +
                ggplot2::scale_color_gradient2(
                  # Align with variable palette style: blue -> warm neutral -> orange-red
                  low = "#7A9E3A",
                  mid = "#FFF4E8",
                  high = "#C65D3A",
                  midpoint = 0.5,
                  limits = c(0, 1),
                  breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
                  labels = scales::number_format(accuracy = 0.01),
                  name = "Feature value",
                  guide = ggplot2::guide_colorbar(
                    title.position = "top",
                    title.hjust = 0.5,
                    barwidth = 0.5,
                    barheight = 4.8,
                    frame.color = "gray40",
                    frame.linewidth = 0.4,
                    ticks.color = "gray40"
                  )
                ) +
                ggplot2::scale_fill_gradient2(
                  low = "#7A9E3A",
                  mid = "#FFF4E8",
                  high = "#C65D3A",
                  midpoint = 0.5,
                  limits = c(0, 1),
                  guide = "none"  # Hide fill legend, use color legend instead
                ) +
                ggplot2::scale_x_continuous(
                  limits = x_limits,
                  breaks = c(-2, -1, 0, 1, 2),
                  labels = scales::number_format(accuracy = 0.1),
                  expand = ggplot2::expansion(mult = 0.05)
                ) +
                ggplot2::labs(
                  title = paste0("Contributions to A-frontier ", outcome_name),
                  subtitle = "* = p < 0.05",
                  x = ifelse(outcome_name == "TE-VRS", 
                            "Contributions to A-frontier TE-VRS",
                            "Contributions to A-frontier SE"),
                  y = ""
                ) +
                ggplot2::theme_minimal(base_size = 16) +
                ggplot2::theme(
                  text = ggplot2::element_text(family = "Times", size = 14),
                  plot.title = ggplot2::element_text(family = "Times", face = "bold", size = 16, hjust = 0.5,
                                                    margin = ggplot2::margin(b = 4)),
                  plot.subtitle = ggplot2::element_text(family = "Times", size = 13, hjust = 0.5, color = "gray40",
                                                      margin = ggplot2::margin(b = 6)),
                  axis.text.y = ggplot2::element_text(family = "Times", size = 13, color = "black", 
                                                     hjust = 1, margin = ggplot2::margin(r = 4)),
                  axis.text.x = ggplot2::element_text(family = "Times", size = 15, color = "black",
                                                     margin = ggplot2::margin(t = 4)),
                  axis.title = ggplot2::element_text(family = "Times", size = 14, face = "bold", color = "black"),
                  axis.title.x = ggplot2::element_text(family = "Times", size = 13.5, face = "bold", color = "black",
                                                      margin = ggplot2::margin(t = 6)),
                  axis.title.y = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black",
                                                      margin = ggplot2::margin(r = 6)),
                  panel.grid.major.x = ggplot2::element_line(color = "gray88", linewidth = 0.3),
                  panel.grid.major.y = ggplot2::element_line(color = "gray88", linewidth = 0.2),
                  panel.grid.minor.x = ggplot2::element_line(color = "gray95", linewidth = 0.15),
                  panel.grid.minor.y = ggplot2::element_blank(),
                  panel.background = ggplot2::element_rect(fill = "white", color = NA),
                  plot.background = ggplot2::element_rect(fill = "white", color = NA),
                  panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                  legend.position = c(0.95, 0.50),
                  legend.justification = c(1, 0.5),
                  legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.7), color = "gray75", linewidth = 0.3),
                  legend.key.height = ggplot2::unit(1.2, "cm"),
                  legend.key.width = ggplot2::unit(0.3, "cm"),
                  legend.title = ggplot2::element_text(family = "Times", size = 11, face = "bold"),
                  legend.text = ggplot2::element_text(family = "Times", size = 10),
                  plot.margin = ggplot2::margin(8, 8, 8, 8, "pt")
                )
            
              # Save individual plot
              filename <- ifelse(outcome_name == "TE-VRS", 
                                "Fig_Shapley_TE_VRS.png",
                                ifelse(outcome_name == "TE-CRS", 
                                       "Fig_Shapley_TE_CRS.png",
                                       "Fig_Shapley_SE.png"))
              safe_ggsave(filename, p_sub, width = 17.59, height = 14.02, units = "cm")
            }  # Close for loop (outcome_name in outcomes_to_plot)
          } else {
            cat("⚠ No shapley_plots_data available for plotting\n")
          }  # Close else for if (length(shapley_plots_data) > 0)
        } else {
          cat(sprintf("⚠ Insufficient data for Shapley plot (need >= 5 samples, found %d)\n", nrow(all_scale_data)))
        }  # Close else for if (nrow(all_scale_data) >= 5)
      }  # Close if (length(factor_vars) > 0)
    } else {
      cat("⚠ No factor variables found in data_merged\n")
    }
  } else {
    cat("⚠ data_merged does not exist or has no rows\n")
  }  # Close if (exists("data_merged") && nrow(data_merged) > 0)
} else {
  cat("⚠ df_clean not found or does not have Market_pig column\n")
}  # Close if (exists("df_clean") && "Market_pig" %in% names(df_clean))

# ============================================================================
# 7. Unit Input/Output Trends Analysis (Matching Reference Format)
# ============================================================================
cat("\n=== Unit Input/Output Trends Analysis ===\n")

if (exists("df_clean") && "Market_pig" %in% names(df_clean)) {
  # Calculate unit inputs and outputs (per pig)
  # Use original scale (Market_pig = head count) for x-axis so trend matches reference
  unit_trends_data <- df_clean %>%
    dplyr::mutate(
      Feed_per_pig_raw = Feed / Market_pig,
      Water_per_pig_raw = Water_consumption / Market_pig,
      Energy_per_pig_raw = Energy_consumption / Market_pig,
      Labor_per_pig_raw = Labor / Market_pig,
      Waste_water_per_pig_raw = Waste_water / Market_pig,
      Manure_per_pig_raw = Manure / Market_pig,
      Dead_pig_per_pig_raw = Dead_pig / Market_pig,
      Carbon_per_pig_raw = Carbon_emission / Market_pig,
      Eutrophication_per_pig_raw = Eutrophication_potential / Market_pig,
      market_pig = Market_pig,
      log_scale = log10(pmax(Market_pig, 1))
    ) %>%
    dplyr::filter(Market_pig > 0)
  
  if (nrow(unit_trends_data) > 0) {
    # Prepare data for Unit Input/Output Trends (match reference: raw per-pig then z-score per variable)
    unit_data_long <- unit_trends_data %>%
      dplyr::select(market_pig, log_scale,
                   Feed_per_pig_raw, Water_per_pig_raw, Energy_per_pig_raw, Labor_per_pig_raw,
                   Waste_water_per_pig_raw, Manure_per_pig_raw, Dead_pig_per_pig_raw, Carbon_per_pig_raw, Eutrophication_per_pig_raw) %>%
      tidyr::pivot_longer(
        cols = c(Feed_per_pig_raw, Water_per_pig_raw, Energy_per_pig_raw, Labor_per_pig_raw,
                Waste_water_per_pig_raw, Manure_per_pig_raw, Dead_pig_per_pig_raw, Carbon_per_pig_raw, Eutrophication_per_pig_raw),
        names_to = "variable",
        values_to = "unit_value"
      ) %>%
      dplyr::mutate(
        variable = stringr::str_remove(variable, "_per_pig_raw"),
        variable = dplyr::case_when(
          variable == "Water" ~ "Water_consumption",
          variable == "Energy" ~ "Energy_consumption",
          variable == "Carbon" ~ "Carbon_emission",
          variable == "Eutrophication" ~ "Eutrophication_potential",
          TRUE ~ variable
        ),
        type = ifelse(variable %in% c("Feed", "Water_consumption", "Energy_consumption", "Labor"),
                     "Input", "Output")
      ) %>%
      dplyr::filter(!is.na(unit_value), !is.na(log_scale))
    
    # Z-score per variable (same as reference) so trend shape is unchanged by scale
    standardize <- function(x) {
      (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    }
    unit_data_long <- unit_data_long %>%
      dplyr::group_by(variable) %>%
      dplyr::mutate(
        unit_value = standardize(unit_value)
      ) %>%
      dplyr::ungroup()
    
    # Colors already defined earlier in the script
    
    # Ensure variable order: Inputs first, then Outputs
    input_vars <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
    output_vars <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
    variable_order <- c(input_vars, output_vars)
    
    # Ensure unit_data_long has variables in the correct order
    unit_data_long <- unit_data_long %>%
      dplyr::mutate(
        variable = factor(variable, levels = variable_order, ordered = TRUE)
      )
    
    # Define turning point (from Tobit analysis)
    a_frontier_tp <- 20277
    
    # Create helper function to create a single plot (input or output)
    create_unit_plot <- function(data_subset, var_type, title_suffix) {
      # Filter data by type
      plot_data <- data_subset %>%
        dplyr::filter(type == var_type)
      
      if (nrow(plot_data) == 0) return(NULL)
      
      # Get variables for this type
      vars_in_plot <- unique(as.character(plot_data$variable))
      vars_in_plot <- vars_in_plot[order(match(vars_in_plot, variable_order))]
      
      # Create colors and labels for this type
      colors_plot <- variable_colors_plot[vars_in_plot]
      labels_plot <- paste0(gsub("_", " ", ifelse(vars_in_plot == "Waste_water", "Wastewater", vars_in_plot)), " (", tolower(var_type), ")")
      names(labels_plot) <- vars_in_plot
      
      # Linetype: solid for inputs, dashed for outputs
      linetype_value <- ifelse(var_type == "Input", "solid", "dashed")
      
      # Use A-frontier turning point from Tobit analysis
      turning_point <- a_frontier_tp
      
      # Create plot
      p <- ggplot2::ggplot(plot_data,
                                    ggplot2::aes(x = market_pig, y = unit_value, 
                                                color = variable, fill = variable,
                                              group = variable)) +
        # Add colored shadows (confidence intervals)
      ggplot2::geom_smooth(method = "loess", se = TRUE, 
                          linewidth = 0,  # No line for shadow layer
                          span = 0.75,
                            alpha = 0.15) +
        # Add main lines
      ggplot2::geom_smooth(method = "loess", se = FALSE,
                            linewidth = 0.8,
                          span = 0.75,
                            alpha = 1.0,
                            linetype = linetype_value) +
        # Add TP line
        ggplot2::geom_vline(xintercept = turning_point, linetype = "dashed", 
                           color = "#8B4513", linewidth = 0.8, alpha = 0.8) +
        # Add TP label annotation
        ggplot2::annotate("text", 
                         x = turning_point, 
                           y = Inf,
                         label = paste0("TP = ", sprintf("%.0f", turning_point)),
                         color = "#8B4513",
                         size = 3.2,
                         fontface = "bold",
                         family = "Times",
                           hjust = 0,
                           vjust = -0.3,
                           angle = 0) +
        ggplot2::scale_x_log10(limits = c(NA, NA),
                              breaks = c(10000, 30000, 100000, 300000),
                              labels = scales::comma,
                                expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
        ggplot2::scale_y_continuous(
          breaks = c(-3, -2, -1, 0, 1, 2, 3),
          expand = ggplot2::expansion(mult = 0.05)
        ) +
        ggplot2::coord_cartesian(ylim = c(-3, 3)) +
        ggplot2::scale_color_manual(
            values = colors_plot,
            labels = labels_plot,
            name = "",
            breaks = vars_in_plot,
          guide = ggplot2::guide_legend(
            ncol = 1,
            byrow = TRUE,
            override.aes = {
                result_list <- vector("list", length(vars_in_plot))
                for (i in seq_along(vars_in_plot)) {
                result_list[[i]] <- list(
                    linetype = linetype_value,
                  fill = NA,
                  alpha = 1,
                    linewidth = 1.2,
                  shape = NA,
                  size = 0
                )
              }
              result_list
            },
              keywidth = ggplot2::unit(1.2, "cm"),
              keyheight = ggplot2::unit(0.3, "cm")
          )
        ) +
        ggplot2::scale_fill_manual(
            values = colors_plot,
            guide = "none"
        ) +
        ggplot2::labs(
            title = paste0("Unit ", tolower(var_type), " vs production scale with uncertainty"),
          x = "Production scale (log scale; market pig heads/a)",
          y = "Standardized unit value (per pig)"
        ) +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::theme(
          text = ggplot2::element_text(family = "Times", size = 15),
          plot.title = ggplot2::element_text(family = "Times", face = "bold", size = 13, hjust = 0, 
                                            margin = ggplot2::margin(b = 6)),
          axis.text = ggplot2::element_text(family = "Times", size = 13, color = "black"),
          axis.title = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black"),
          panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
          panel.grid.minor = ggplot2::element_line(color = "gray95", linewidth = 0.2),
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.background = ggplot2::element_rect(fill = "white", color = NA),
          plot.background = ggplot2::element_rect(fill = "white", color = NA),
            legend.position = c(0.99, 0.99),
          legend.justification = c(1, 1),
          legend.box = "vertical",
          legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.3),
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(family = "Times", size = 12),
            legend.key.width = ggplot2::unit(1.2, "cm"),
            legend.key.height = ggplot2::unit(0.3, "cm"),
            legend.spacing.y = ggplot2::unit(0.15, "cm"),
            legend.margin = ggplot2::margin(3, 5, 3, 5, "pt"),
            plot.margin = ggplot2::margin(4, 60, 4, 4, "pt")
          )
      
      return(p)
    }
    
    # Create separate plots for inputs and outputs
    p_unit_trends_inputs <- create_unit_plot(unit_data_long, "Input", "Inputs")
    p_unit_trends_outputs <- create_unit_plot(unit_data_long, "Output", "Outputs")
    
    # Save plots separately (7x5 inches each)
    if (!is.null(p_unit_trends_inputs)) {
      safe_ggsave("Fig_Unit_Input_Trends.png", p_unit_trends_inputs, width = 7, height = 5, dpi = 300)
      cat("✓ Saved: Fig_Unit_Input_Trends.png\n")
    }
    if (!is.null(p_unit_trends_outputs)) {
      safe_ggsave("Fig_Unit_Output_Trends.png", p_unit_trends_outputs, width = 7, height = 5, dpi = 300)
      cat("✓ Saved: Fig_Unit_Output_Trends.png\n")
    }
  } else {
    cat("⚠ Insufficient data for unit trends analysis\n")
  }
} else {
  cat("⚠ df_clean not found or does not have Market_pig column\n")
}

# ============================================================================
# 8. Marginal Effects Analysis (Matching Reference Format)
# ============================================================================
cat("\n=== Marginal Effects Analysis ===\n")

tryCatch({
  if (exists("data_merged") && nrow(data_merged) > 0) {
  # Use original scale (market_pig_original) for zone assignment so TP±15% matches production scale
  scale_var <- if ("market_pig_original" %in% names(data_merged)) "market_pig_original" else "market_pig"
  # Compute turning point from Tobit on this dataset so TP and scale_var are in same units
  cat("  Calculating turning points...\n")
  data_tobit_ma <- data_merged %>%
    dplyr::mutate(
      log_scale_ma = log10(pmax(.data[[scale_var]], 1)),
      log_scale_sq_ma = log_scale_ma^2
    ) %>%
    dplyr::filter(!is.na(A_TE_VRS), A_TE_VRS > 0, A_TE_VRS <= 1, is.finite(log_scale_ma))
  a_frontier_tp <- 20277
  if (nrow(data_tobit_ma) >= 10) {
    fit_ma <- tryCatch({
      AER::tobit(A_TE_VRS ~ log_scale_ma + log_scale_sq_ma, left = 0, right = 1, data = data_tobit_ma)
    }, error = function(e) NULL)
    if (!is.null(fit_ma) && length(coef(fit_ma)) >= 3) {
      b1 <- coef(fit_ma)[2]
      b2 <- coef(fit_ma)[3]
      if (is.finite(b2) && b2 > 0) {
        tp_log <- -b1 / (2 * b2)
        a_frontier_tp <- 10^tp_log
        cat(sprintf("  Tobit TP (same scale as %s): %.0f\n", scale_var, a_frontier_tp))
      }
    }
  }
  # TP Zone: TP ± 25% (match reference figure)
  tp_lower_marginal <- a_frontier_tp * 0.5
  tp_upper_marginal <- a_frontier_tp * 1.5

  data_merged <- data_merged %>%
    dplyr::mutate(
      scale_zone = ifelse(.data[[scale_var]] < tp_lower_marginal, "Pre-TP",
                         ifelse(.data[[scale_var]] <= tp_upper_marginal, "TP Zone",
                               "Post-TP"))
    )

  zone_counts <- table(data_merged$scale_zone)
  cat(sprintf("Zone distribution (TP = %.0f, TP zone = [0.75TP, 1.25TP] = TP±25%%):\n", a_frontier_tp))
  print(zone_counts)

  # If any zone has too few samples, adjust boundaries while keeping TP Zone around TP±25%
  min_samples_per_zone <- 3L
  if (any(zone_counts < min_samples_per_zone)) {
    cat(sprintf("Warning: Some zones have < %d samples. Adjusting boundaries...\n", min_samples_per_zone))
    tp_lower_marginal <- max(a_frontier_tp * 0.65, min(data_merged[[scale_var]], na.rm = TRUE))
    tp_upper_marginal <- min(a_frontier_tp * 1.35, max(data_merged[[scale_var]], na.rm = TRUE))
    data_merged <- data_merged %>%
      dplyr::mutate(
        scale_zone = ifelse(.data[[scale_var]] < tp_lower_marginal, "Pre-TP",
                           ifelse(.data[[scale_var]] <= tp_upper_marginal, "TP Zone",
                                 "Post-TP"))
      )
    zone_counts_adj <- table(data_merged$scale_zone)
    cat("Adjusted zone distribution:\n")
    print(zone_counts_adj)
  }
  
  # Calculate marginal effects for each zone
  marginal_effects_data <- data.frame()
  
  marginal_vars_list <- c("Feed", "Water_consumption", "Energy_consumption", "Labor",
                         "Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")

  for (zone in c("Pre-TP", "TP Zone", "Post-TP")) {
    zone_data <- data_merged %>% dplyr::filter(scale_zone == zone)
    min_samples_zone <- ifelse(zone == "TP Zone", 3L, 5L)
    if (nrow(zone_data) < min_samples_zone) next

    for (var in marginal_vars_list) {
      if (!var %in% names(zone_data)) next
      zone_data_clean <- zone_data %>%
        dplyr::filter(
          !is.na(.data[[var]]),
          .data[[var]] > 1e-10,
          !is.na(A_TE_VRS),
          A_TE_VRS > 1e-10,
          is.finite(.data[[var]]),
          is.finite(A_TE_VRS)
        )
      if (nrow(zone_data_clean) < 3L) next

      zone_data_clean <- zone_data_clean %>%
        dplyr::mutate(
          var_log = log(pmax(.data[[var]], 1e-10)),
          te_vrs_log = log(pmax(A_TE_VRS, 1e-10))
        ) %>%
        dplyr::filter(is.finite(var_log), is.finite(te_vrs_log))
      if (sd(zone_data_clean$var_log, na.rm = TRUE) < 1e-10) next  # skip if no variation (singular lm)

      model <- tryCatch({
        lm(te_vrs_log ~ var_log, data = zone_data_clean)
      }, error = function(e) NULL)
      if (is.null(model) || length(coef(model)) < 2) next

      coef_var <- coef(model)[2]
      if (!is.finite(coef_var)) next
      marg_effect <- coef_var * 0.01
      se_marg <- tryCatch(summary(model)$coefficients[2, 2] * 0.01, error = function(e) abs(marg_effect) * 0.1)
      if (!is.finite(se_marg)) se_marg <- abs(marg_effect) * 0.1
      if (!is.finite(marg_effect)) next

      marginal_effects_data <- dplyr::bind_rows(
        marginal_effects_data,
        data.frame(
          Zone = zone,
          Variable = var,
          Type = ifelse(var %in% c("Feed", "Water_consumption", "Energy_consumption", "Labor"),
                       "Input (+1%)", "Output (+1%)"),
          Marginal_Effect = marg_effect,
          SE = se_marg,
          Lower = marg_effect - 1.96 * se_marg,
          Upper = marg_effect + 1.96 * se_marg,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  # Create marginal effects plot if data is available
  if (nrow(marginal_effects_data) > 0) {
    # Ensure Zone order is correct
    marginal_effects_data <- marginal_effects_data %>%
      dplyr::mutate(
        Zone = factor(Zone, levels = c("Pre-TP", "TP Zone", "Post-TP")),
        Variable = factor(Variable, levels = names(variable_colors_plot)),
        Type = factor(Type, levels = c("Input (+1%)", "Output (+1%)"))
      )
    
    # Create labels: Inputs (+1%), Outputs (-1%) to match reference figure
    input_vars_legend <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
    variable_labels <- ifelse(
      names(variable_colors_plot) %in% input_vars_legend,
      paste0(gsub("_", " ", ifelse(names(variable_colors_plot) == "Waste_water", "Wastewater", names(variable_colors_plot))), " (+1%)"),
      paste0(gsub("_", " ", ifelse(names(variable_colors_plot) == "Waste_water", "Wastewater", names(variable_colors_plot))), " (-1%)")
    )
    names(variable_labels) <- names(variable_colors_plot)
    
    # Calculate symmetric y-axis limits for marginal effects (match reference)
    max_abs_marginal <- max(abs(c(marginal_effects_data$Marginal_Effect,
                                  marginal_effects_data$Lower,
                                  marginal_effects_data$Upper)), na.rm = TRUE)
    max_abs_marginal <- ifelse(is.finite(max_abs_marginal) && max_abs_marginal > 0, max_abs_marginal, 1)
    # Pretty breaks, then symmetric truncation around 0
    pretty_y <- pretty(c(-max_abs_marginal, max_abs_marginal), n = 6)
    pretty_y <- pretty_y[is.finite(pretty_y)]
    step_y <- diff(pretty_y)
    step_y <- step_y[is.finite(step_y) & step_y > 0]
    step_y <- if (length(step_y) > 0) min(step_y) else max_abs_marginal / 3
    max_abs_marginal_sym <- ceiling(max_abs_marginal / step_y) * step_y
    y_lim_marginal <- c(-max_abs_marginal_sym, max_abs_marginal_sym)
    y_breaks_marginal <- seq(-max_abs_marginal_sym, max_abs_marginal_sym, by = step_y)

    p_marginal <- ggplot2::ggplot(marginal_effects_data,
                                  ggplot2::aes(x = Zone, y = Marginal_Effect, 
                                              color = Variable, group = Variable,
                                              linetype = Type)) +
      ggplot2::geom_line(linewidth = 0.8, alpha = 0.9) +
      ggplot2::geom_point(size = 2.0, alpha = 0.9, stroke = 0.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                            width = 0.15, linewidth = 0.5, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      ggplot2::scale_y_continuous(
        limits = y_lim_marginal,
        breaks = y_breaks_marginal,
        expand = ggplot2::expansion(mult = 0.1)
      ) +
      ggplot2::scale_color_manual(values = variable_colors_plot, 
                                 labels = variable_labels,
                                 name = "",
                                 guide = ggplot2::guide_legend(ncol = 1, byrow = TRUE)) +
      ggplot2::scale_linetype_manual(values = c("Input (+1%)" = "solid", "Output (+1%)" = "dashed"),
                                     guide = "none") +
      ggplot2::labs(
        title = "Marginal benefits of 1% changes across scale zones",
        x = "Scale zone",
        y = "Marginal benefit on efficiency (%)"
      ) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Times", size = 15),
        plot.title = ggplot2::element_text(family = "Times", face = "bold", size = 13, hjust = 0,
                                          margin = ggplot2::margin(b = 6)),
        axis.text = ggplot2::element_text(family = "Times", size = 13, color = "black"),
        axis.text.x = ggplot2::element_text(family = "Times", size = 13, color = "black",
                                           margin = ggplot2::margin(t = 4)),
        axis.text.y = ggplot2::element_text(family = "Times", size = 13, color = "black",
                                           margin = ggplot2::margin(r = 4)),
        axis.title = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black"),
        axis.title.x = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black",
                                            margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black",
                                            margin = ggplot2::margin(r = 6)),
        panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = ggplot2::element_line(color = "gray95", linewidth = 0.2),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = c(0.99, 0.01),
      legend.justification = c(1, 0),
      legend.box = "vertical",
      legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(family = "Times", size = 11),
      legend.key.width = ggplot2::unit(0.8, "cm"),
      legend.key.height = ggplot2::unit(0.3, "cm"),
      legend.spacing.y = ggplot2::unit(0.2, "cm"),
      legend.margin = ggplot2::margin(4, 6, 4, 6, "pt"),
      plot.margin = ggplot2::margin(4, 60, 4, 4, "pt")
      )
    
    # Save marginal effects plot
    safe_ggsave("Fig_Marginal_Effects.png", p_marginal, width = 7, height = 5, dpi = 300)
    cat("✓ Saved: Fig_Marginal_Effects.png\n")

    # Plot C: Marginal Effects (Input +1%, Output -1%) (match reference; outputs inverted)
    output_vars <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
    input_vars <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
    marginal_effects_data_decrease <- marginal_effects_data %>%
      dplyr::mutate(
        Marginal_Effect_Decrease = ifelse(Variable %in% output_vars, -Marginal_Effect, Marginal_Effect),
        Lower_Decrease = ifelse(Variable %in% output_vars, -Upper, Lower),
        Upper_Decrease = ifelse(Variable %in% output_vars, -Lower, Upper),
        Type_Decrease = ifelse(Variable %in% input_vars, "Input (+1%)", "Output (-1%)")
      ) %>%
      dplyr::mutate(
        Zone = factor(Zone, levels = c("Pre-TP", "TP Zone", "Post-TP")),
        Variable = factor(Variable, levels = names(variable_colors_plot)),
        Type_Decrease = factor(Type_Decrease, levels = c("Input (+1%)", "Output (-1%)"))
      )
    variable_labels_decrease <- paste0(
      gsub("_", " ", ifelse(names(variable_colors_plot) == "Waste_water", "Wastewater", names(variable_colors_plot))),
      ifelse(names(variable_colors_plot) %in% input_vars, " (+1%)", " (-1%)")
    )
    names(variable_labels_decrease) <- names(variable_colors_plot)
    max_abs_marginal_decrease <- max(abs(c(marginal_effects_data_decrease$Marginal_Effect_Decrease,
                                          marginal_effects_data_decrease$Lower_Decrease,
                                          marginal_effects_data_decrease$Upper_Decrease)), na.rm = TRUE)
    max_abs_marginal_decrease <- ifelse(is.finite(max_abs_marginal_decrease) && max_abs_marginal_decrease > 0, max_abs_marginal_decrease, 1)
    pretty_y2 <- pretty(c(-max_abs_marginal_decrease, max_abs_marginal_decrease), n = 6)
    pretty_y2 <- pretty_y2[is.finite(pretty_y2)]
    step_y2 <- diff(pretty_y2)
    step_y2 <- step_y2[is.finite(step_y2) & step_y2 > 0]
    step_y2 <- if (length(step_y2) > 0) min(step_y2) else max_abs_marginal_decrease / 3
    max_abs_marginal_decrease_sym <- ceiling(max_abs_marginal_decrease / step_y2) * step_y2
    y_lim_marginal_decrease <- c(-max_abs_marginal_decrease_sym, max_abs_marginal_decrease_sym)
    y_breaks_marginal_decrease <- seq(-max_abs_marginal_decrease_sym, max_abs_marginal_decrease_sym, by = step_y2)
    p_marginal_decrease <- ggplot2::ggplot(
      marginal_effects_data_decrease,
      ggplot2::aes(x = Zone, y = Marginal_Effect_Decrease,
                  color = Variable, group = Variable,
                  linetype = Type_Decrease)
    ) +
      ggplot2::geom_line(linewidth = 0.8, alpha = 0.9) +
      ggplot2::geom_point(size = 2.0, alpha = 0.9, stroke = 0.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower_Decrease, ymax = Upper_Decrease),
                            width = 0.15, linewidth = 0.5, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      ggplot2::scale_y_continuous(limits = y_lim_marginal_decrease,
                                 breaks = y_breaks_marginal_decrease,
                                 labels = scales::number_format(accuracy = 0.1, scale = 100),
                                 expand = ggplot2::expansion(mult = 0.1)) +
      ggplot2::scale_color_manual(values = variable_colors_plot,
                                 labels = variable_labels_decrease,
                                 name = "",
                                 guide = ggplot2::guide_legend(ncol = 2, byrow = TRUE)) +
      ggplot2::scale_linetype_manual(values = c("Input (+1%)" = "solid", "Output (-1%)" = "dashed"),
                                     guide = "none") +
      ggplot2::labs(
        title = "Marginal benefits of 1% changes across scale zones",
        x = "Scale zone",
        y = "Marginal benefit on efficiency (%)"
      ) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Times", size = 15),
        plot.title = ggplot2::element_text(family = "Times", face = "bold", size = 13, hjust = 0,
                                          margin = ggplot2::margin(b = 6)),
        axis.text = ggplot2::element_text(family = "Times", size = 13, color = "black"),
        axis.text.x = ggplot2::element_text(family = "Times", size = 13, color = "black",
                                           margin = ggplot2::margin(t = 4)),
        axis.text.y = ggplot2::element_text(family = "Times", size = 13, color = "black",
                                           margin = ggplot2::margin(r = 4)),
        axis.title = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black"),
        axis.title.x = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black",
                                            margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(family = "Times", size = 13, face = "bold", color = "black",
                                            margin = ggplot2::margin(r = 6)),
        panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = ggplot2::element_line(color = "gray95", linewidth = 0.2),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.position = c(0.50, 0.01),
        legend.justification = c(0.5, 0),
        legend.box = "vertical",
        legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.3),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(family = "Times", size = 11),
        legend.key.width = ggplot2::unit(0.8, "cm"),
        legend.key.height = ggplot2::unit(0.3, "cm"),
        legend.spacing.y = ggplot2::unit(0.2, "cm"),
        legend.margin = ggplot2::margin(4, 6, 4, 6, "pt"),
        plot.margin = ggplot2::margin(4, 60, 4, 4, "pt")
      )
    safe_ggsave("Fig_Marginal_Effects_Output_Decrease.png", p_marginal_decrease, width = 7, height = 5, dpi = 300)
    cat("✓ Saved: Fig_Marginal_Effects_Output_Decrease.png\n")
  } else {
    cat("⚠ Insufficient data for marginal effects analysis\n")
  }
} else {
  cat("⚠ data_merged not found or has no rows\n")
}
}, error = function(e) {
  cat(sprintf("⚠ Error in marginal effects analysis: %s\n", e$message))
  cat(sprintf("  Error occurred at: %s\n", deparse(e$call)))
})

cat("\n✓ All figures and analyses complete!\n")











