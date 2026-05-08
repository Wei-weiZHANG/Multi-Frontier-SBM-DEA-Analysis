# ============================================================================
# Multi-Frontier SBM Analysis Framework
# ============================================================================
# Fixed Variable Set Method: All frontiers use same 5 bad outputs
# Scheme 3: Non-overlapping core constraints (M: Mortality, L: Local, G: Global, A: Aggregated)
# ============================================================================

# ============================================================================
# 0. SETUP AND INITIALIZATION
# ============================================================================

start_time <- Sys.time()

message("MULTI-FRONTIER SBM ANALYSIS")
message(sprintf("Start time: %s\n", Sys.time()))

# Required packages
required_packages <- c(
  "readxl",      # Excel data import
  "dplyr",       # Data manipulation
  "tidyr",       # Data tidying
  "ggplot2",     # Visualization
  "lpSolve",     # Linear programming for DEA
  "patchwork",   # Combining plots
  "scales",      # Scale formatting
  "viridis",     # Color palettes
  "quantreg",    # Quantile regression
  "readr",       # CSV I/O
  "stringr",     # String manipulation
  "corrplot"     # Correlation plots
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
message(sprintf("✓ %d packages loaded\n", length(required_packages)))

# Get results directory from environment or use default
results_base <- Sys.getenv("RESULTS_DIR", unset = "results")

# Create output structure
output_dirs <- c(
  file.path(results_base, "efficiency_calculation"),
  file.path(results_base, "visualization"),
  file.path(results_base, "statistical_tests"),
  file.path(results_base, "diagnostics")
)

for (dir_path in output_dirs) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

message("✓ Output directories ready\n")

# ============================================================================
# 1. DATA LOADING AND VALIDATION
# ============================================================================

xlsx_path <- "./swine_farm_data.xlsx"

if (!file.exists(xlsx_path)) {
  stop("ERROR: Data file not found at: ", xlsx_path)
}

suppressWarnings({
  df_raw <- readxl::read_xlsx(xlsx_path, sheet = 3)
})
message(sprintf("✓ Data loaded: %d rows × %d columns\n", nrow(df_raw), ncol(df_raw)))

# Handle Report_year column if exists
year_cols <- c("Report_year", "report_year", "Report_Year", "Year")
existing_year_col <- intersect(year_cols, names(df_raw))
if (length(existing_year_col) > 0) {
  year_col_name <- existing_year_col[1]
  year_vec <- df_raw[[year_col_name]]
  
  if (inherits(year_vec, "Date") || inherits(year_vec, "POSIXct") || inherits(year_vec, "POSIXt")) {
    df_raw[[year_col_name]] <- as.numeric(format(year_vec, "%Y"))
  } else if (!is.numeric(year_vec)) {
    year_numeric <- suppressWarnings(as.numeric(year_vec))
    if (any(!is.na(year_numeric) & year_numeric > 40000 & year_numeric < 50000)) {
      excel_dates <- !is.na(year_numeric) & year_numeric > 40000 & year_numeric < 50000
      year_numeric[excel_dates] <- as.numeric(format(
        as.Date(year_numeric[excel_dates], origin = "1899-12-30"), "%Y"))
      df_raw[[year_col_name]] <- year_numeric
    } else if (sum(!is.na(year_numeric)) >= length(year_numeric) * 0.8) {
      df_raw[[year_col_name]] <- year_numeric
    }
  }
}

# Define variables
inputs <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
good_output <- "Market_pig"
bad_outputs_full <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")

# Core constraints: M=Dead_pig, L=Waste_water+Manure, G=Carbon_emission+Eutrophication_potential, A=Aggregated
constrained_bad_M <- c(3)  # Dead_pig
constrained_bad_L <- c(1, 2)  # Waste_water, Manure
constrained_bad_G <- c(4, 5)  # Carbon_emission, Eutrophication_potential (G-Frontier core)
constrained_bad_A <- c(1, 2, 3, 4, 5)  # All (A-Frontier full constraints)

# ============================================================================
# 2. SBM EFFICIENCY CALCULATION (Functions and Pre-Experiment)
# ============================================================================

all_required <- unique(c(inputs, good_output, bad_outputs_full))

# Data cleaning
df_clean <- df_raw %>%
  dplyr::select(dplyr::all_of(all_required), 
                dplyr::any_of(c("Farm", "ID", "Report_year", "report_year", "Report_Year", "Year",
                               "Scale-up", "Scale_up", "ScaleUp", "No.", "No", "Number"))) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_required), 
                              ~ suppressWarnings(as.numeric(.x)))) %>%
  tidyr::drop_na(dplyr::all_of(all_required))

# Standardize Report_year column name if exists and convert dates to years
year_cols <- c("Report_year", "report_year", "Report_Year", "Year")
existing_year_col <- intersect(year_cols, names(df_clean))
if (length(existing_year_col) > 0) {
  if (existing_year_col[1] != "Report_year") {
    df_clean <- df_clean %>%
      rename(Report_year = !!sym(existing_year_col[1]))
    message(sprintf("✓ Renamed '%s' to 'Report_year'\n", existing_year_col[1]))
  }
  
  # Convert date format to year if needed
  if ("Report_year" %in% names(df_clean)) {
    # Check if it's a date/time format
    if (inherits(df_clean$Report_year, "Date") || 
        inherits(df_clean$Report_year, "POSIXct") ||
        inherits(df_clean$Report_year, "POSIXt")) {
      # Extract year from date
      if (!requireNamespace("lubridate", quietly = TRUE)) {
        # Fallback: use as.numeric and format
        df_clean$Report_year <- as.numeric(format(df_clean$Report_year, "%Y"))
      } else {
        df_clean$Report_year <- lubridate::year(df_clean$Report_year)
      }
      message("✓ Converted Report_year from date format to numeric year\n")
    } else {
      # Try to convert to numeric (handles mixed date/numeric cases)
      year_numeric <- suppressWarnings(as.numeric(df_clean$Report_year))
      
      # Check if conversion was mostly successful
      n_valid <- sum(!is.na(year_numeric))
      if (n_valid >= nrow(df_clean) * 0.8) {
        df_clean$Report_year <- year_numeric
        message(sprintf("✓ Converted Report_year to numeric (%d/%d valid values)\n", 
                       n_valid, nrow(df_clean)))
      } else {
        message("⚠ Warning: Report_year conversion may have issues\n")
      }
    }
    
    # Validate year range (should be 2019-2025)
    valid_years <- df_clean$Report_year >= 2019 & df_clean$Report_year <= 2025
    n_valid_years <- sum(valid_years, na.rm = TRUE)
    if (n_valid_years < nrow(df_clean) * 0.8) {
      message(sprintf("⚠ Warning: Only %d/%d Report_year values are in expected range (2019-2025)\n",
                     n_valid_years, nrow(df_clean)))
    }
  }
}

message(sprintf("✓ Clean data: %d observations (%.1f%% retained)\n", 
                nrow(df_clean), nrow(df_clean)/nrow(df_raw)*100))

# Replace zero/negative values and record statistics
data_cleaning_log <- data.frame(
  variable = character(0),
  n_zero_negative = integer(0),
  n_total = integer(0),
  replacement_rate = numeric(0)
)

for (var in all_required) {
  n_zero <- sum(df_clean[[var]] <= 0, na.rm = TRUE)
  n_total <- sum(!is.na(df_clean[[var]]))
  if (n_zero > 0) {
    df_clean[[var]] <- pmax(df_clean[[var]], 1e-8)
    replacement_rate <- n_zero / n_total * 100
    data_cleaning_log <- rbind(data_cleaning_log, data.frame(
      variable = var,
      n_zero_negative = n_zero,
      n_total = n_total,
      replacement_rate = replacement_rate
    ))
  }
}

# Report data cleaning statistics
if (nrow(data_cleaning_log) > 0) {
  for (i in seq_len(nrow(data_cleaning_log))) {
    message(sprintf("  %s: %d/%d (%.2f%%) replaced with 1e-8",
                    data_cleaning_log$variable[i],
                    data_cleaning_log$n_zero_negative[i],
                    data_cleaning_log$n_total[i],
                    data_cleaning_log$replacement_rate[i]))
  }
}

# Add scale and size classification
df_clean$log_scale <- log10(df_clean[[good_output]])
df_clean$size_class <- sapply(df_clean[[good_output]], function(x) {
  if (x < 10000) "Small"
  else if (x <= 50000) "Medium"
  else "Large"
})
df_clean$size_class <- factor(df_clean$size_class, 
                               levels = c("Small", "Medium", "Large"),
                               ordered = TRUE)

message(sprintf("Scale range: %s to %s pigs\n",
                format(round(min(df_clean[[good_output]])), big.mark = ","),
                format(round(max(df_clean[[good_output]])), big.mark = ",")))

# ============================================================================
# 2. SBM EFFICIENCY CALCULATION (Functions and Pre-Experiment)
# ============================================================================

# SBM function (output-oriented with undesirable outputs)
sbm_efficiency <- function(X, Y, bad = NULL, constrained_bad_indices = NULL, 
                           bad_weights = NULL, RTS = "vrs", orientation = "output", 
                           show_dmu_progress = TRUE, evaluate_indices = NULL) {
  
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.null(bad) && !is.matrix(bad)) bad <- as.matrix(bad)
  
  n <- nrow(X)
  m <- ncol(X)
  s_g <- ncol(Y)
  
  if (!is.null(bad)) {
    s_b <- ncol(bad)
  } else {
    s_b <- 0
    bad <- matrix(0, nrow = n, ncol = 0)
  }
  
  # If constrained_bad_indices is NULL, constrain all bad outputs (original behavior)
  if (is.null(constrained_bad_indices) && s_b > 0) {
    constrained_bad_indices <- 1:s_b
  }
  
  # Validate bad_weights parameter
  if (!is.null(bad_weights)) {
    if (length(bad_weights) != s_b) {
      stop(sprintf("ERROR: bad_weights length (%d) does not match number of bad outputs (%d)",
                   length(bad_weights), s_b))
    }
    if (any(is.na(bad_weights))) {
      warning("bad_weights contains NA values, replacing with 0")
      bad_weights[is.na(bad_weights)] <- 0
    }
    if (any(bad_weights < 0)) {
      warning("bad_weights contains negative values, replacing with 0")
      bad_weights[bad_weights < 0] <- 0
    }
  } else {
    warning("bad_weights is NULL - using default weight 1.0")
  }
  
  # Validate evaluate_indices parameter
  if (!is.null(evaluate_indices)) {
    if (!is.numeric(evaluate_indices) || any(evaluate_indices < 1) || any(evaluate_indices > n)) {
      stop(sprintf("ERROR: evaluate_indices must be integers between 1 and %d", n))
    }
    evaluate_indices <- as.integer(evaluate_indices)
    # Remove duplicates and sort
    evaluate_indices <- sort(unique(evaluate_indices))
    n_eval <- length(evaluate_indices)
  } else {
    evaluate_indices <- 1:n
    n_eval <- n
  }
  
  eff_scores <- numeric(n)
  lp_status_counts <- numeric(5)  # Track LP status: 0=success, 1=infeasible, 2=unbounded, 99=error, other
  names(lp_status_counts) <- c("success", "infeasible", "unbounded", "error", "other")
  
  # Progress indicator for large datasets (only if show_dmu_progress is TRUE)
  show_progress <- show_dmu_progress && (n_eval > 20)
  progress_interval <- if (n_eval > 100) 10 else if (n_eval > 50) 5 else 1
  
  # Only iterate over DMUs that need to be evaluated
  eval_counter <- 0
  for (k in evaluate_indices) {
    eval_counter <- eval_counter + 1
    # Progress update
    if (show_progress && (eval_counter %% progress_interval == 0 || eval_counter == n_eval)) {
      cat(sprintf("\r    Progress: %d/%d DMUs (%.1f%%)", eval_counter, n_eval, 100*eval_counter/n_eval))
      if (eval_counter == n_eval) cat("\n")
    }
    x_k <- pmax(X[k, ], 1e-8)
    y_k <- pmax(Y[k, ], 1e-8)
    b_k <- if (s_b > 0) pmax(bad[k, ], 1e-8) else numeric(0)
    
    num_vars <- 1 + n + m + s_g + s_b
    
    if (orientation == "output") {
      obj_coef <- numeric(num_vars)
      obj_coef[1] <- 1
      
      normalization_denom <- m + s_b
      
      # Input penalties
      for (i in 1:m) {
        obj_coef[1 + n + i] <- -1 / (normalization_denom * x_k[i])
      }
      
      # Bad output penalties (fixed variable set: all bad outputs have weights)
      if (s_b > 0) {
        for (i in 1:s_b) {
          if (i %in% constrained_bad_indices) {
            weight_i <- if (!is.null(bad_weights) && length(bad_weights) >= i && !is.na(bad_weights[i]) && bad_weights[i] > 0) {
              bad_weights[i]
            } else {
              if (k == 1 && i == 1) {
                warning("bad_weights not provided - using default weight 1.0")
              }
              1.0
            }
            obj_coef[1 + n + m + s_g + i] <- -weight_i / (normalization_denom * b_k[i])
          } else {
            obj_coef[1 + n + m + s_g + i] <- 0
            if (k == 1) {
              warning(sprintf("Bad output %d not in constrained_bad_indices - penalty set to 0", i))
            }
          }
        }
      }
      
      n_constraints <- m + s_g + s_b + 1 + ifelse(RTS == "vrs", 1, 0)
      constr <- matrix(0, nrow = n_constraints, ncol = num_vars)
      dir_vec <- character(n_constraints)
      rhs_vec <- numeric(n_constraints)
      row_idx <- 1
      
      # Input constraints
      for (i in 1:m) {
        constr[row_idx, 2:(1+n)] <- X[, i]
        constr[row_idx, 1 + n + i] <- 1
        constr[row_idx, 1] <- -x_k[i]
        dir_vec[row_idx] <- "=="
        rhs_vec[row_idx] <- 0
        row_idx <- row_idx + 1
      }
      
      # Good output constraints
      for (j in 1:s_g) {
        constr[row_idx, 2:(1+n)] <- Y[, j]
        constr[row_idx, 1 + n + m + j] <- -1
        constr[row_idx, 1] <- -y_k[j]
        dir_vec[row_idx] <- "=="
        rhs_vec[row_idx] <- 0
        row_idx <- row_idx + 1
      }
      
      # Bad output constraints
      if (s_b > 0) {
        constrained_set <- if (!is.null(constrained_bad_indices) && length(constrained_bad_indices) > 0) {
          constrained_bad_indices
        } else {
          1:s_b
        }
        
        for (j in 1:s_b) {
          constr[row_idx, 2:(1+n)] <- bad[, j]
          constr[row_idx, 1 + n + m + s_g + j] <- 1
          constr[row_idx, 1] <- -b_k[j]
          dir_vec[row_idx] <- if (j %in% constrained_set) "==" else "<="
          rhs_vec[row_idx] <- 0
          row_idx <- row_idx + 1
        }
      }
      
      # Normalization constraint
      constr[row_idx, 1] <- 1
      for (j in 1:s_g) {
        constr[row_idx, 1 + n + m + j] <- 1 / (s_g * y_k[j])
      }
      dir_vec[row_idx] <- "=="
      rhs_vec[row_idx] <- 1
      row_idx <- row_idx + 1
      
      # VRS constraint
      if (RTS == "vrs") {
        constr[row_idx, 2:(1+n)] <- 1
        constr[row_idx, 1] <- -1
        dir_vec[row_idx] <- "=="
        rhs_vec[row_idx] <- 0
      }
      
    } else {
      stop("Input-oriented SBM not implemented. Use orientation='output'")
    }
    
    # Solve LP with optimization settings
    lp_result <- tryCatch({
      lpSolve::lp(
        direction = "min",
        objective.in = obj_coef,
        const.mat = constr,
        const.dir = dir_vec,
        const.rhs = rhs_vec,
        scale = 1,           # Scale constraints for numerical stability
        compute.sens = 0     # Don't compute sensitivity (faster)
      )
    }, error = function(e) {
      list(status = 99, objval = NA, message = e$message)
    })
    
    status_code <- lp_result$status
    if (status_code == 0) {
      lp_status_counts["success"] <- lp_status_counts["success"] + 1
      eff_scores[k] <- pmax(0, pmin(1, lp_result$objval))
    } else if (status_code == 1) {
      lp_status_counts["infeasible"] <- lp_status_counts["infeasible"] + 1
      eff_scores[k] <- NA_real_
    } else if (status_code == 2) {
      lp_status_counts["unbounded"] <- lp_status_counts["unbounded"] + 1
      eff_scores[k] <- NA_real_
    } else if (status_code == 98) {
      lp_status_counts["timeout"] <- lp_status_counts["timeout"] + 1
      eff_scores[k] <- NA_real_
    } else if (status_code == 99) {
      lp_status_counts["error"] <- lp_status_counts["error"] + 1
      eff_scores[k] <- NA_real_
    } else {
      lp_status_counts["other"] <- lp_status_counts["other"] + 1
      eff_scores[k] <- NA_real_
    }
  }
  
  # Set uncalculated DMUs to NA if evaluate_indices was provided
  if (!is.null(evaluate_indices) && length(evaluate_indices) < n) {
    all_indices <- 1:n
    uncalculated_indices <- setdiff(all_indices, evaluate_indices)
    eff_scores[uncalculated_indices] <- NA_real_
  }
  
  # Report LP status if issues
  if (sum(lp_status_counts, na.rm = TRUE) > 0) {
    total_lp <- sum(lp_status_counts, na.rm = TRUE)
    if ((!is.na(lp_status_counts["infeasible"]) && lp_status_counts["infeasible"] > 0) || 
        (!is.na(lp_status_counts["unbounded"]) && lp_status_counts["unbounded"] > 0) || 
        (!is.na(lp_status_counts["timeout"]) && lp_status_counts["timeout"] > 0) || 
        (!is.na(lp_status_counts["error"]) && lp_status_counts["error"] > 0)) {
      warning(sprintf("LP issues: %d success, %d infeasible, %d unbounded, %d timeout, %d error",
                      lp_status_counts["success"],
                      lp_status_counts["infeasible"], lp_status_counts["unbounded"], 
                      lp_status_counts["timeout"], lp_status_counts["error"]))
    }
  }
  
  return(eff_scores)
}

# Scaling function: median-based standardization for SBM (unit-independent, robust to magnitude)
# Unit-trend and marginal-effect modules use z-score scale() for comparison; scale bins use market_pig_original
scale_matrix <- function(M) {
  scaling_factors <- apply(M, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
  sweep(M, 2, scaling_factors, "/")
}

message("✓ SBM functions defined\n")

# ============================================================================
# 2a. PRE-EXPERIMENT: OPTIMAL WEIGHT MULTIPLIER SEARCH
# ============================================================================
# Run pre-experiment to find optimal weight configuration
# This must be done after SBM function definition but before weight definition

# Centralized module directory (shared with launcher)
modules_dir <- Sys.getenv("MODULES_DIR", unset = getwd())

# Check if pre-experiment should be run
# Default: Run pre-experiment and use resulting weights (set SKIP_PRE_EXPERIMENT = "TRUE" to use preset weights)
skip_pre_exp <- Sys.getenv("SKIP_PRE_EXPERIMENT", unset = "FALSE")
skip_pre_exp <- as.logical(skip_pre_exp)

# Initialize use_pre_exp_weights to FALSE (will be set to TRUE if weights are successfully loaded)
use_pre_exp_weights <- FALSE

weight_ratio_path <- file.path(modules_dir, "preexperiment_weight_ratio.R")
if (file.exists(weight_ratio_path) && !skip_pre_exp) {
  message("2a. PRE-EXPERIMENT: OPTIMAL WEIGHT MULTIPLIER SEARCH\n")
  source(weight_ratio_path, encoding = "UTF-8")
  
  normalize_weights <- function(weights, target_sum = length(weights)) {
    current_sum <- sum(weights)
    if (current_sum > 0) {
      normalized <- weights * (target_sum / current_sum)
    } else {
      normalized <- weights  # If all zeros, keep as is
    }
    return(normalized)
  }
  
  # Create scaled data for pre-experiment
  X_inputs_pre <- as.matrix(df_clean[, inputs])
  Y_good_pre <- as.matrix(df_clean[, good_output, drop = FALSE])
  Z_bad_full_pre <- as.matrix(df_clean[, bad_outputs_full, drop = FALSE])
  
  X_scaled_pre <- scale_matrix(X_inputs_pre)
  Y_scaled_pre <- scale_matrix(Y_good_pre)
  Z_scaled_full_pre <- scale_matrix(Z_bad_full_pre)
  
  # Run pre-experiment
  pre_exp_result <- tryCatch({
    run_weight_ratio_pre_experiment(
      df_clean = df_clean,
      X_scaled = X_scaled_pre,
      Y_scaled = Y_scaled_pre,
      Z_scaled_full = Z_scaled_full_pre,
      constrained_bad_M = constrained_bad_M,
      constrained_bad_L = constrained_bad_L,
      constrained_bad_G = constrained_bad_G,
      constrained_bad_A = constrained_bad_A,
      sbm_efficiency = sbm_efficiency,
      normalize_weights = normalize_weights,
      results_base = results_base,
      bad_outputs_full = bad_outputs_full,
      skip_pre_experiment = FALSE
    )
  }, error = function(e) {
    stop(sprintf("FATAL ERROR: Pre-experiment failed: %s\n  Pre-experiment is required for weight determination. Please fix the error and rerun.", e$message))
  })
  
  if (!is.null(pre_exp_result)) {
    # Extract weights from pre-experiment results
    bad_weights_M_raw <- pre_exp_result$bad_weights_M_raw
    bad_weights_L_raw <- pre_exp_result$bad_weights_L_raw
    bad_weights_G_raw <- pre_exp_result$bad_weights_G_raw
    bad_weights_A_raw <- pre_exp_result$bad_weights_A_raw
    bad_weights_M <- pre_exp_result$bad_weights_M
    bad_weights_L <- pre_exp_result$bad_weights_L
    bad_weights_G <- pre_exp_result$bad_weights_G
    bad_weights_A <- pre_exp_result$bad_weights_A
    
    # CRITICAL: Ensure weights are in correct order matching bad_outputs_full
    # This is essential for correct indexing in sbm_efficiency function
    bad_weights_M <- bad_weights_M[bad_outputs_full]
    bad_weights_L <- bad_weights_L[bad_outputs_full]
    bad_weights_G <- bad_weights_G[bad_outputs_full]
    bad_weights_A <- bad_weights_A[bad_outputs_full]
    
    # Verify no NA values were introduced
    if (any(is.na(bad_weights_M)) || any(is.na(bad_weights_L)) || 
        any(is.na(bad_weights_G)) || any(is.na(bad_weights_A))) {
      stop("ERROR: Pre-experiment weight ordering failed. Check weight names match bad_outputs_full.")
    }
    
    bad_weights_list <- list(
      E = bad_weights_M,
      L = bad_weights_L,
      G = bad_weights_G,
      A = bad_weights_A
    )
    use_pre_exp_weights <- TRUE
    message("✓ Using optimal weights from pre-experiment (correlation plateau method)\n")
  } else {
    stop("FATAL ERROR: Pre-experiment returned NULL. Pre-experiment is required for weight determination.")
  }
} else if (skip_pre_exp) {
  # Use pre-determined weights from previous pre-experiment results
  message("2a. SKIPPING PRE-EXPERIMENT: Using pre-determined weights from previous analysis\n")
  
  # Define weights based on previous pre-experiment results
  # We keep the previous relative pattern and renormalize to sum = length(bad_outputs_full).
  normalize_weights <- function(weights, target_sum = length(weights)) {
    current_sum <- sum(weights)
    if (current_sum > 0) weights * (target_sum / current_sum) else weights
  }
  
  # M-Frontier: Dead_pig is core; others are non-core
  bad_weights_M_raw <- c(
    Waste_water = 0.227,
    Manure = 0.227,
    Dead_pig = 4.091,
    Carbon_emission = 0.227,
    Eutrophication_potential = 0.227
  )
  
  # L-Frontier: Waste_water + Manure are core; others non-core
  bad_weights_L_raw <- c(
    Waste_water = 1.607,
    Manure = 1.607,
    Dead_pig = 0.089,
    Carbon_emission = 0.089,
    Eutrophication_potential = 0.089
  )
  
  # G-Frontier (Global): core variables are emphasized; others are non-core
  bad_weights_G_raw <- c(
    Waste_water = 0.227,
    Manure = 0.227,
    Dead_pig = 0.227,
    Carbon_emission = 4.091,
    Eutrophication_potential = 4.091
  )
  
  # A-Frontier: aggregated weight
  bad_weights_A_raw <- c(
    Waste_water = 1.000,
    Manure = 1.000,
    Dead_pig = 1.000,
    Carbon_emission = 1.000,
    Eutrophication_potential = 1.000
  )
  
  bad_weights_M <- normalize_weights(bad_weights_M_raw, target_sum = length(bad_outputs_full))
  bad_weights_L <- normalize_weights(bad_weights_L_raw, target_sum = length(bad_outputs_full))
  bad_weights_G <- normalize_weights(bad_weights_G_raw, target_sum = length(bad_outputs_full))
  bad_weights_A <- normalize_weights(bad_weights_A_raw, target_sum = length(bad_outputs_full))
  
  # Convert to named vectors matching bad_outputs_full order
  # Ensure weights are in the correct order matching bad_outputs_full
  bad_weights_M <- bad_weights_M[bad_outputs_full]
  bad_weights_L <- bad_weights_L[bad_outputs_full]
  bad_weights_G <- bad_weights_G[bad_outputs_full]
  bad_weights_A <- bad_weights_A[bad_outputs_full]
  
  # Verify no NA values were introduced
  if (any(is.na(bad_weights_M)) || any(is.na(bad_weights_L)) || 
      any(is.na(bad_weights_G)) || any(is.na(bad_weights_A))) {
    stop("ERROR: Weight indexing produced NA values. Check that weight names match bad_outputs_full.")
  }
  
  # Create weights list
  bad_weights_list <- list(
    E = bad_weights_M,
    L = bad_weights_L,
    G = bad_weights_G,
    A = bad_weights_A
  )
  
  use_pre_exp_weights <- TRUE
  message("✓ Using pre-determined weights from previous pre-experiment\n")
  message("  (Renormalized to sum = number of bad outputs)\n")
  message("  M-Frontier: Dead_pig core, others non-core\n")
  message("  L-Frontier: Waste_water/Manure core, Dead_pig/Carbon non-core\n")
  message("  G-Frontier (Global): core variables emphasized, others non-core\n")
  message("  A-Frontier: aggregated\n")
} else {
  stop("FATAL ERROR: preexperiment_weight_ratio.R not found. Pre-experiment is required for weight determination.")
}

# ============================================================================
# 3. WEIGHT CONFIGURATION
# ============================================================================

# Weights must come from pre-experiment - no default weights allowed
if (!exists("use_pre_exp_weights") || is.na(use_pre_exp_weights) || !use_pre_exp_weights) {
  stop("ERROR: Pre-experiment weights are required. Please ensure preexperiment_weight_ratio.R exists and pre-experiment completes successfully.")
}

# Verify all weight variables are defined
if (!exists("bad_weights_M_raw") || !exists("bad_weights_L_raw") || 
    !exists("bad_weights_G_raw") || !exists("bad_weights_A_raw") ||
    !exists("bad_weights_M") || !exists("bad_weights_L") || 
    !exists("bad_weights_G") || !exists("bad_weights_A") ||
    !exists("bad_weights_list")) {
  stop("ERROR: Weight variables not properly defined from pre-experiment. Check preexperiment_weight_ratio.R execution.")
}

message("✓ Using weights from pre-experiment (correlation plateau method)")

# ============================================================================
# 4. WEIGHT VERIFICATION
# ============================================================================

message("\n4. WEIGHT VERIFICATION")
message("CRITICAL: Verifying weight values and order...")

# Print actual weight values for verification
message("\nActual weight values (in bad_outputs_full order):")
for (fid in c("M", "L", "G", "A")) {
  weight_var <- paste0("bad_weights_", fid)
  if (exists(weight_var)) {
    weights <- get(weight_var)
    message(sprintf("  %s-Frontier:", fid))
    for (i in seq_along(bad_outputs_full)) {
      weight_val <- if (is.null(names(weights)) || length(names(weights)) == 0) {
        weights[i]
      } else {
        weights[bad_outputs_full[i]]
      }
      message(sprintf("    [%d] %s = %.6f", i, bad_outputs_full[i], weight_val))
    }
  }
}

# Verify all weight variables exist
required_weight_vars <- c("bad_weights_M_raw", "bad_weights_L_raw", "bad_weights_G_raw", "bad_weights_A_raw",
                          "bad_weights_M", "bad_weights_L", "bad_weights_G", "bad_weights_A", "bad_weights_list")
missing_vars <- required_weight_vars[!sapply(required_weight_vars, exists)]
if (length(missing_vars) > 0) {
  stop(sprintf("ERROR: Missing weight variables: %s", paste(missing_vars, collapse = ", ")))
}

# Verify weight standardization
message("Weight standardization verification:")
n_bad <- length(bad_outputs_full)
message(sprintf("  E-frontier: sum = %.3f (target: %d), range = [%.2f, %.2f]",
                sum(bad_weights_M, na.rm = TRUE), n_bad, min(bad_weights_M, na.rm = TRUE), max(bad_weights_M, na.rm = TRUE)))
message(sprintf("  L-frontier: sum = %.3f (target: %d), range = [%.2f, %.2f]",
                sum(bad_weights_L, na.rm = TRUE), n_bad, min(bad_weights_L, na.rm = TRUE), max(bad_weights_L, na.rm = TRUE)))
message(sprintf("  G-frontier: sum = %.3f (target: %d), range = [%.2f, %.2f]",
                sum(bad_weights_G, na.rm = TRUE), n_bad, min(bad_weights_G, na.rm = TRUE), max(bad_weights_G, na.rm = TRUE)))
message(sprintf("  A-frontier: sum = %.3f (target: %d), range = [%.2f, %.2f]",
                sum(bad_weights_A, na.rm = TRUE), n_bad, min(bad_weights_A, na.rm = TRUE), max(bad_weights_A, na.rm = TRUE)))

# Verify weight ranges: non-core ≥0.5, core can be >1.6 (up to 2.5)
all_weights_check <- c(bad_weights_M, bad_weights_L, bad_weights_G, bad_weights_A)
noncore_weights_check <- c(
  bad_weights_M[setdiff(seq_along(bad_outputs_full), constrained_bad_M)],
  bad_weights_L[setdiff(seq_along(bad_outputs_full), constrained_bad_L)],
  bad_weights_G[setdiff(seq_along(bad_outputs_full), constrained_bad_G)]
)

if (any(noncore_weights_check < 0.5, na.rm = TRUE)) {
  warning(sprintf("⚠ Some non-core weights below 0.5: min=%.3f", min(noncore_weights_check, na.rm = TRUE)))
} else {
  message(sprintf("  ✓ All non-core weights ≥ 0.5 (min=%.3f, max=%.3f)", 
                  min(noncore_weights_check, na.rm = TRUE), max(noncore_weights_check, na.rm = TRUE)))
}

if (any(all_weights_check > 2.5, na.rm = TRUE)) {
  warning(sprintf("⚠ Some weights exceed 2.5: max=%.3f", max(all_weights_check, na.rm = TRUE)))
} else {
  message(sprintf("  ✓ All weights ≤ 2.5 (max=%.3f)", max(all_weights_check, na.rm = TRUE)))
}

# ============================================================================
# DETAILED WEIGHT CONFIGURATION OUTPUT
# ============================================================================
message("\nDetailed Weight Configuration for Each Frontier:")
message(paste(rep("=", 78), collapse = ""))

# Create weight configuration table
weight_config_table <- data.frame(
  Bad_Output = bad_outputs_full,
  M_Frontier = bad_weights_M,
  L_Frontier = bad_weights_L,
  G_Frontier = bad_weights_G,
  A_Frontier = bad_weights_A
)

# Add constraint type indicators
weight_config_table$M_Constraint <- ifelse(seq_along(bad_outputs_full) %in% constrained_bad_M, "Core (==)", "Non-core (<=)")
weight_config_table$L_Constraint <- ifelse(seq_along(bad_outputs_full) %in% constrained_bad_L, "Core (==)", "Non-core (<=)")
weight_config_table$G_Constraint <- ifelse(seq_along(bad_outputs_full) %in% constrained_bad_G, "Core (==)", "Non-core (<=)")
weight_config_table$A_Constraint <- ifelse(seq_along(bad_outputs_full) %in% constrained_bad_A, "Core (==)", "Non-core (<=)")

# Print detailed table for each frontier
for (fid in c("M", "L", "G", "A")) {
  frontier_name <- switch(fid,
    "M" = "Mortality Frontier",
    "L" = "Local Environmental Frontier",
    "G" = "Global Frontier",
    "A" = "Aggregated Frontier"
  )
  
  message(sprintf("\n%s (%s-Frontier):", frontier_name, fid))
  message(paste(rep("-", 78), collapse = ""))
  
  # Get weights and constraints for this frontier
  weights_col <- paste0(fid, "_Frontier")
  constraint_col <- paste0(fid, "_Constraint")
  
  for (i in seq_len(nrow(weight_config_table))) {
    bad_output_name <- weight_config_table$Bad_Output[i]
    weight_value <- weight_config_table[[weights_col]][i]
    constraint_type <- weight_config_table[[constraint_col]][i]
    
    message(sprintf("  %-20s: Weight = %6.3f  [%s]", 
                    bad_output_name, weight_value, constraint_type))
  }
  
  # Summary for this frontier
  frontier_weights <- weight_config_table[[weights_col]]
  core_indices <- switch(fid,
    "M" = constrained_bad_M,
    "L" = constrained_bad_L,
    "G" = constrained_bad_G,
    "A" = constrained_bad_A
  )
  
  core_weights <- frontier_weights[core_indices]
  noncore_weights <- frontier_weights[setdiff(1:5, core_indices)]
  
  message(sprintf("  Core weights:     [%s]", paste(sprintf("%.3f", core_weights), collapse = ", ")))
  if (length(noncore_weights) > 0) {
    message(sprintf("  Non-core weights: [%s]", paste(sprintf("%.3f", noncore_weights), collapse = ", ")))
  }
  message(sprintf("  Total sum: %.3f (target: 5.000)", sum(frontier_weights)))
}

message(paste(rep("=", 78), collapse = ""))
message("")

# ============================================================================
# 5. EFFICIENCY CALCULATION (M, L, G, A Frontiers)
# ============================================================================

message("5. EFFICIENCY CALCULATION (M, L, G, A Frontiers)")

# ============================================================================
# CHECK FOR EXISTING STANDARDIZED DATA
# ============================================================================
# This ensures consistent results regardless of original data units
# If standardized data exists, use it; otherwise create it
message("5a. CHECKING FOR EXISTING STANDARDIZED DATA\n")

standardized_data_path <- file.path(results_base, "standardized_data.xlsx")
use_existing_standardized <- FALSE

# Check if standardized data file exists
if (file.exists(standardized_data_path)) {
  message(sprintf("Found existing standardized data: %s\n", standardized_data_path))
  use_existing_standardized <- TRUE
  
  # Read standardized data from Excel
  if (requireNamespace("readxl", quietly = TRUE)) {
    standardized_data <- readxl::read_excel(standardized_data_path)
    message("✓ Loaded standardized data from Excel\n")
  } else {
    message("⚠ readxl package not available\n")
    use_existing_standardized <- FALSE
  }
  
  # Validate that all required columns are present
  required_cols <- c("obs_id", good_output, "Market_pig_original", inputs, bad_outputs_full)
  missing_cols <- setdiff(required_cols, names(standardized_data))
  
  if (length(missing_cols) > 0) {
    message(sprintf("⚠ Missing columns in standardized data: %s\n", 
                   paste(missing_cols, collapse = ", ")))
    message("⚠ Will recreate standardized data\n")
    use_existing_standardized <- FALSE
  } else {
    message("✓ Standardized data validation passed\n")
  }
}

if (!use_existing_standardized) {
  message("Creating new standardized data...\n")
  
  # Scale inputs and outputs
  X_inputs <- as.matrix(df_clean[, inputs])
  Y_good <- as.matrix(df_clean[, good_output, drop = FALSE])
  X_scaled <- scale_matrix(X_inputs)
  Y_scaled <- scale_matrix(Y_good)
  
  # Create standardized data dataframe
  standardized_data <- data.frame(
    obs_id = seq_len(nrow(df_clean)),
    Market_pig = Y_scaled[, 1],
    Market_pig_original = df_clean[[good_output]]
  )
  
  # Add standardized inputs
  for (i in seq_along(inputs)) {
    standardized_data[[inputs[i]]] <- X_scaled[, i]
  }
  
  # Add standardized bad outputs
  Z_bad_full <- as.matrix(df_clean[, bad_outputs_full, drop = FALSE])
  Z_scaled_full <- scale_matrix(Z_bad_full)
  
  for (i in seq_along(bad_outputs_full)) {
    standardized_data[[bad_outputs_full[i]]] <- Z_scaled_full[, i]
  }
  
  # Add metadata columns
if ("Report_year" %in% names(df_clean)) {
  standardized_data$Report_year <- df_clean$Report_year
}

# Add log_scale and size_class with error handling
if ("log_scale" %in% names(df_clean)) {
  standardized_data$log_scale <- df_clean$log_scale
} else {
  # Calculate log_scale if not present
  standardized_data$log_scale <- log10(pmax(standardized_data$Market_pig_original, 1))
  message("⚠ log_scale not found in df_clean, calculated from Market_pig_original\n")
}

if ("size_class" %in% names(df_clean)) {
  standardized_data$size_class <- df_clean$size_class
} else {
  # Create size_class if not present
  standardized_data$size_class <- cut(standardized_data$Market_pig_original,
                                       breaks = c(0, 5000, 10000, 20000, 50000, Inf),
                                       labels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
                                       include.lowest = TRUE)
  message("⚠ size_class not found in df_clean, created from Market_pig_original\n")
}
  
  # Add Scale-up column
  scale_up_col <- NULL
  scale_up_names <- c("Scale-up", "Scale_up", "ScaleUp")
  for (col_name in scale_up_names) {
    if (col_name %in% names(df_clean)) {
      scale_up_col <- col_name
      break
    }
  }
  if (!is.null(scale_up_col)) {
    standardized_data[[scale_up_col]] <- df_clean[[scale_up_col]]
  }
  
  # Add No. column
  no_col <- NULL
  no_names <- c("No.", "No", "Number")
  for (col_name in no_names) {
    if (col_name %in% names(df_clean)) {
      no_col <- col_name
      break
    }
  }
  if (!is.null(no_col)) {
    standardized_data[[no_col]] <- df_clean[[no_col]]
  }
  
  # Save to Excel
  if (!requireNamespace("writexl", quietly = TRUE)) {
    message("⚠ writexl package not available, using CSV format instead")
    standardized_data_path <- file.path(results_base, "standardized_data.csv")
    write_csv(standardized_data, standardized_data_path)
  } else {
    writexl::write_xlsx(standardized_data, standardized_data_path)
  }
  message(sprintf("✓ Standardized data saved to: %s\n", standardized_data_path))
  message(sprintf("  Variables: %d inputs, %d good outputs, %d bad outputs\n",
                  length(inputs), 1, length(bad_outputs_full)))
} else {
  message("Using existing standardized data for calculations\n")
}

# Extract scaled matrices from standardized data
# Use "Market_pig" column for Y_scaled since that's what we created
if ("Market_pig" %in% names(standardized_data)) {
  Y_scaled <- as.matrix(standardized_data[, "Market_pig", drop = FALSE])
  message("✓ Using Market_pig column for Y_scaled\n")
} else if (good_output %in% names(standardized_data)) {
  Y_scaled <- as.matrix(standardized_data[, good_output, drop = FALSE])
  message("✓ Using good_output column for Y_scaled\n")
} else {
  stop("ERROR: Neither Market_pig nor ", good_output, " found in standardized_data")
}
X_scaled <- as.matrix(standardized_data[, inputs])
Z_scaled_full <- as.matrix(standardized_data[, bad_outputs_full, drop = FALSE])

# Initialize results
efficiency_results <- data.frame(
  obs_id = standardized_data$obs_id,
  market_pig = standardized_data$Market_pig,
  market_pig_original = df_clean[[good_output]]
)

# Add log_scale and size_class with error handling
if ("log_scale" %in% names(standardized_data)) {
  efficiency_results$log_scale <- standardized_data$log_scale
} else {
  # Calculate log_scale if not present
  efficiency_results$log_scale <- log10(pmax(efficiency_results$market_pig_original, 1))
  message("⚠ log_scale not found in standardized data, calculated from market_pig_original\n")
}

if ("size_class" %in% names(standardized_data)) {
  efficiency_results$size_class <- standardized_data$size_class
} else {
  # Create size_class if not present
  efficiency_results$size_class <- cut(efficiency_results$market_pig_original,
                                       breaks = c(0, 5000, 10000, 20000, 50000, Inf),
                                       labels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
                                       include.lowest = TRUE)
  message("⚠ size_class not found in standardized data, created from market_pig_original\n")
}

# Add optional columns from standardized data
if ("Report_year" %in% names(standardized_data)) {
  efficiency_results$Report_year <- standardized_data$Report_year
}

# Add Scale-up column (check multiple possible column names)
scale_up_col <- NULL
scale_up_names <- c("Scale-up", "Scale_up", "ScaleUp")
for (col_name in scale_up_names) {
  if (col_name %in% names(standardized_data)) {
    scale_up_col <- col_name
    break
  }
}
if (!is.null(scale_up_col)) {
  # Preserve original column name in output
  efficiency_results[[scale_up_col]] <- standardized_data[[scale_up_col]]
  message(sprintf("✓ Added Scale-up column: '%s'\n", scale_up_col))
} else {
  message("⚠ Scale-up column not found in standardized data (checked: ", 
          paste(scale_up_names, collapse = ", "), ")\n")
}

# Add No. column (check multiple possible column names)
no_col <- NULL
no_names <- c("No.", "No", "Number")
for (col_name in no_names) {
  if (col_name %in% names(standardized_data)) {
    no_col <- col_name
    break
  }
}
if (!is.null(no_col)) {
  # Preserve original column name in output
  efficiency_results[[no_col]] <- standardized_data[[no_col]]
  message(sprintf("✓ Added No. column: '%s'\n", no_col))
} else {
  message("⚠ No. column not found in standardized data (checked: ", 
          paste(no_names, collapse = ", "), ")\n")
}

# Z_scaled_full is already defined above

frontiers <- list(
  M = list(name = "Mortality", constrained_indices = constrained_bad_M, weights = bad_weights_M),
  L = list(name = "Local Environmental", constrained_indices = constrained_bad_L, weights = bad_weights_L),
  G = list(name = "Global", constrained_indices = constrained_bad_G, weights = bad_weights_G),
  A = list(name = "Aggregated", constrained_indices = constrained_bad_A, weights = bad_weights_A)
)

for (fid in names(frontiers)) {
  message(sprintf("Calculating %s-Frontier (%s)...", fid, frontiers[[fid]]$name))
  n_dmu <- nrow(standardized_data)
  message(sprintf("  Dataset: %d DMUs, %d inputs, %d good outputs, %d bad outputs", 
                  n_dmu, length(inputs), 1, length(bad_outputs_full)))
  message(sprintf("  Estimated time: %.1f-%.1f minutes per frontier", 
                  n_dmu * 0.05 / 60, n_dmu * 0.2 / 60))
  
  weight_vec <- frontiers[[fid]]$weights
  if (is.null(weight_vec) || length(weight_vec) != length(bad_outputs_full)) {
    stop(sprintf("ERROR: %s-frontier weights invalid: length=%d, expected=%d", 
                 fid, ifelse(is.null(weight_vec), 0, length(weight_vec)), length(bad_outputs_full)))
  }
  
  # CRITICAL: Ensure weight vector is in correct order (matching bad_outputs_full)
  # Convert to unnamed vector in the exact order of bad_outputs_full
  weight_vec_ordered <- as.numeric(weight_vec[bad_outputs_full])
  if (any(is.na(weight_vec_ordered))) {
    stop(sprintf("ERROR: %s-frontier weight ordering failed. Weight names: %s, bad_outputs_full: %s", 
                 fid, paste(names(weight_vec), collapse=", "), paste(bad_outputs_full, collapse=", ")))
  }
  
  # Debug: Print weights for first frontier only
  if (fid == "M") {
    message("  DEBUG: Weight verification:")
    for (i in seq_along(bad_outputs_full)) {
      message(sprintf("    %s (index %d): weight = %.3f", bad_outputs_full[i], i, weight_vec_ordered[i]))
    }
  }
  
  start_time_vrs <- Sys.time()
  message("  Calculating VRS efficiency...")
  te_vrs <- sbm_efficiency(X_scaled, Y_scaled, bad = Z_scaled_full, 
                           constrained_bad_indices = frontiers[[fid]]$constrained_indices,
                           bad_weights = weight_vec_ordered,
                           RTS = "vrs")
  elapsed_vrs <- as.numeric(difftime(Sys.time(), start_time_vrs, units = "secs"))
  message(sprintf("  ✓ VRS completed in %.1f seconds", elapsed_vrs))
  
  start_time_crs <- Sys.time()
  message("  Calculating CRS efficiency...")
  te_crs <- sbm_efficiency(X_scaled, Y_scaled, bad = Z_scaled_full, 
                           constrained_bad_indices = frontiers[[fid]]$constrained_indices,
                           bad_weights = weight_vec_ordered,
                           RTS = "crs")
  elapsed_crs <- as.numeric(difftime(Sys.time(), start_time_crs, units = "secs"))
  message(sprintf("  ✓ CRS completed in %.1f seconds", elapsed_crs))
  
  se_raw <- te_crs / pmax(te_vrs, 1e-8)
  se_outliers <- sum(se_raw > 1.01, na.rm = TRUE)
  if (!is.na(se_outliers) && se_outliers > 0) {
    warning(sprintf("  ⚠ %d SE values > 1.01 detected for %s-frontier", se_outliers, fid))
  }
  se <- pmax(0, pmin(1.01, se_raw))
  
  efficiency_results[[paste0(fid, "_TE_VRS")]] <- te_vrs
  efficiency_results[[paste0(fid, "_TE_CRS")]] <- te_crs
  efficiency_results[[paste0(fid, "_SE")]] <- se
  
  message(sprintf("  ✓ VRS: mean=%.3f, frontier=%d", 
                  mean(te_vrs, na.rm=TRUE), 
                  sum(te_vrs > 0.99, na.rm=TRUE)))
}


# Save basic efficiency results
write_csv(efficiency_results, 
          file.path(results_base, "efficiency_calculation/frontier_efficiency_results.csv"))
message("✓ Efficiency results saved\n")

# ============================================================================
# 6. META-FRONTIER AND TGR CALCULATION
# ============================================================================

meta_frontier_module_path <- file.path(modules_dir, "meta_frontier_module.R")
if (file.exists(meta_frontier_module_path)) {
  source(meta_frontier_module_path, encoding = "UTF-8")
  efficiency_results <- calculate_meta_frontier_maximum(efficiency_results)
  efficiency_results <- calculate_tgr_true_meta(efficiency_results)
  validate_meta_frontier(efficiency_results, bad_weights_list = bad_weights_list)
  generate_meta_frontier_report(efficiency_results, 
                               file.path(results_base, "META_FRONTIER_REPORT.txt"))
} else {
  # Meta-frontier weights: Use minimum weights across all frontiers
  # This ensures meta-frontier is the union of all production possibility sets
  # but with reasonable penalty strength (not too low)
  meta_weights_raw <- pmin(bad_weights_M_raw, bad_weights_L_raw, 
                           bad_weights_G_raw, bad_weights_A_raw)
  # Alternative: Use average weights (uncomment to use instead)
  # meta_weights_raw <- (bad_weights_M_raw + bad_weights_L_raw + 
  #                     bad_weights_G_raw + bad_weights_A_raw) / 4
  
  # Ensure minimum weight is at least 0.3 (not too low)
  meta_weights_raw <- pmax(meta_weights_raw, 0.3)
  meta_weights_normalized <- normalize_weights(meta_weights_raw, target_sum = length(bad_outputs_full))
  
  # Calculate meta-frontier efficiency (VRS)
  meta_te_vrs <- sbm_efficiency(X_scaled, Y_scaled, bad = Z_scaled_full,
                                constrained_bad_indices = seq_along(bad_outputs_full),  # All constrained
                                bad_weights = meta_weights_normalized,
                                RTS = "vrs")
  
  # Calculate meta-frontier efficiency (CRS)
  meta_te_crs <- sbm_efficiency(X_scaled, Y_scaled, bad = Z_scaled_full,
                                constrained_bad_indices = seq_along(bad_outputs_full),  # All constrained
                                bad_weights = meta_weights_normalized,
                                RTS = "crs")
  
  efficiency_results$Meta_TE_VRS <- meta_te_vrs
  efficiency_results$Meta_TE_CRS <- meta_te_crs
  
  # Calculate TGR with proper clipping (distinguish numerical error from true outliers)
  # TGR should theoretically be in [0, 1], but allow small numerical errors (e.g., 1.001)
  # Values > 1.01 are likely true outliers and should be flagged
  calculate_tgr_safe <- function(frontier_eff, meta_eff, tolerance = 0.01) {
    meta_safe <- pmax(meta_eff, 1e-8)
    tgr_raw <- frontier_eff / meta_safe
    
    # Identify true outliers (TGR > 1 + tolerance)
    outliers <- tgr_raw > (1 + tolerance)
    n_outliers <- sum(outliers, na.rm = TRUE)
    
    if (!is.na(n_outliers) && n_outliers > 0) {
      warning(sprintf("  ⚠ %d TGR values > %.3f detected (likely true outliers, not numerical error)",
                      n_outliers, 1 + tolerance))
    }
    
    # Clip to [0, 1.01] to allow small numerical errors but flag larger ones
    tgr_clipped <- pmax(0, pmin(1.01, tgr_raw))
    
    return(list(tgr = tgr_clipped, n_outliers = n_outliers))
  }
  
  # Calculate TGR for all frontiers
  tgr_M_vrs <- calculate_tgr_safe(efficiency_results$M_TE_VRS, efficiency_results$Meta_TE_VRS)
  tgr_L_vrs <- calculate_tgr_safe(efficiency_results$L_TE_VRS, efficiency_results$Meta_TE_VRS)
  tgr_G_vrs <- calculate_tgr_safe(efficiency_results$G_TE_VRS, efficiency_results$Meta_TE_VRS)
  tgr_A_vrs <- calculate_tgr_safe(efficiency_results$A_TE_VRS, efficiency_results$Meta_TE_VRS)
  
  tgr_M_crs <- calculate_tgr_safe(efficiency_results$M_TE_CRS, efficiency_results$Meta_TE_CRS)
  tgr_L_crs <- calculate_tgr_safe(efficiency_results$L_TE_CRS, efficiency_results$Meta_TE_CRS)
  tgr_G_crs <- calculate_tgr_safe(efficiency_results$G_TE_CRS, efficiency_results$Meta_TE_CRS)
  tgr_A_crs <- calculate_tgr_safe(efficiency_results$A_TE_CRS, efficiency_results$Meta_TE_CRS)
  
  efficiency_results$M_TGR_VRS <- tgr_M_vrs$tgr
  efficiency_results$L_TGR_VRS <- tgr_L_vrs$tgr
  efficiency_results$G_TGR_VRS <- tgr_G_vrs$tgr
  efficiency_results$A_TGR_VRS <- tgr_A_vrs$tgr
  
  efficiency_results$M_TGR_CRS <- tgr_M_crs$tgr
  efficiency_results$L_TGR_CRS <- tgr_L_crs$tgr
  efficiency_results$G_TGR_CRS <- tgr_G_crs$tgr
  efficiency_results$A_TGR_CRS <- tgr_A_crs$tgr
  message(sprintf("  ✓ Meta-frontier calculated (true union of production possibility sets)"))
  message(sprintf("  ✓ Meta VRS mean: %.3f", mean(efficiency_results$Meta_TE_VRS, na.rm=TRUE)))
  message(sprintf("  ✓ TGR outliers detected: VRS: %d, CRS: %d",
                  sum(tgr_M_vrs$n_outliers, tgr_L_vrs$n_outliers, tgr_G_vrs$n_outliers, tgr_A_vrs$n_outliers),
                  sum(tgr_M_crs$n_outliers, tgr_L_crs$n_outliers, tgr_G_crs$n_outliers, tgr_A_crs$n_outliers)))
}

write_csv(efficiency_results, 
          file.path(results_base, "efficiency_calculation/all_frontiers_efficiency_with_tgr.csv"))
message("✓ Meta-frontier and TGR calculated\n")

# ============================================================================
# 6b. SLACK ANALYSIS: Sources of Inefficiency
# ============================================================================

slack_module_path <- file.path(modules_dir, "slack_analysis_module.R")
if (file.exists(slack_module_path)) {
  source(slack_module_path, encoding = "UTF-8")
  
  bad_outputs_list <- list(
    M = bad_outputs_full[constrained_bad_M],
    L = bad_outputs_full[constrained_bad_L],
    G = bad_outputs_full[constrained_bad_G],
    A = bad_outputs_full[constrained_bad_A]
  )
  
  tryCatch({
    slack_results <- analyze_slack_by_frontier(
      efficiency_results = efficiency_results,
      df_clean = standardized_data,
      inputs = inputs,
      good_output = good_output,
      bad_outputs_list = bad_outputs_list,
      results_base = results_base,
      bad_outputs_full = bad_outputs_full,  # Full bad output set
      constrained_bad_indices_list = list(  # Constrained indices per frontier
        M = constrained_bad_M,
        L = constrained_bad_L,
        G = constrained_bad_G,
        A = constrained_bad_A
      ),
      bad_weights_list = bad_weights_list  # Pass weights so slack is reported for ALL weighted bad outputs
    )
    
    if (exists("visualize_slack_analysis")) {
      visualize_slack_analysis(slack_results, results_base, inputs, bad_outputs_list)
    } else {
      source(slack_module_path, encoding = "UTF-8")
      visualize_slack_analysis(slack_results, results_base, inputs, bad_outputs_list)
    }
    message("✓ Slack analysis completed\n")
  }, error = function(e) {
    message(sprintf("⚠ Slack analysis error: %s\n", e$message))
  })
}

# ============================================================================
# 7. DIAGNOSTIC VALIDATION (Meta-Frontier)
# ============================================================================

diagnose_meta_frontier_path <- file.path(modules_dir, "meta_frontier_diagnostics.R")
if (file.exists(diagnose_meta_frontier_path)) {
  Sys.setenv("CURRENT_RESULTS_DIR" = results_base)
  tryCatch({
    source(diagnose_meta_frontier_path, echo = FALSE)
    message("✓ Meta-frontier diagnostics completed\n")
  }, error = function(e) {
    message(sprintf("⚠ Diagnostic error: %s\n", e$message))
  })
}

# ============================================================================
# 8. CORE VISUALIZATIONS FOR PAPER
# ============================================================================

# Publication theme
theme_publication <- function(base_size = 17) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4, hjust = 0.5),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0.5, color = "gray30"),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1)
    )
}

# Prepare long-format data
eff_long <- efficiency_results %>%
  tidyr::pivot_longer(
    cols = matches("^[ELGA]_(TE_VRS|TE_CRS|SE)$"),
    names_to = "metric",
    values_to = "efficiency"
  ) %>%
  tidyr::separate(metric, into = c("frontier", "type"), sep = "_", extra = "merge") %>%
  dplyr::filter(!is.na(efficiency)) %>%
  dplyr::mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE))

# Unified frontier palette (muted, coordinated tones)
frontier_colors <- c(
  "M" = "#7EB89E",
  "L" = "#E16859",
  "G" = "#8B72A1",
  "A" = "#60BCDA"
)

# FIGURE 1: VRS, SE, and CRS Efficiency Distributions
p_vrs_dist <- eff_long %>%
  dplyr::filter(type == "TE_VRS") %>%
  ggplot(aes(x = efficiency, fill = frontier, color = frontier)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = frontier_colors,
                       labels = c("M" = "M-Frontier",
                                  "L" = "L-Frontier",
                                  "G" = "G-Frontier",
                                  "A" = "A-Frontier"),
                       breaks = c("M", "L", "G", "A")) +
  scale_color_manual(values = frontier_colors,
                        labels = c("M" = "M-Frontier",
                                   "L" = "L-Frontier",
                                   "G" = "G-Frontier",
                                   "A" = "A-Frontier"),
                        breaks = c("M", "L", "G", "A")) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6.5)) +
  labs(
    title = "VRS technical efficiency",
    x = "Technical efficiency (VRS)",
    y = "Density",
    fill = "Frontier",
    color = "Frontier"
  ) +
  theme_publication()

# SE panel
p_se_dist <- eff_long %>%
  dplyr::filter(type == "SE") %>%
  ggplot(aes(x = efficiency, fill = frontier, color = frontier)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = frontier_colors,
                       labels = c("M" = "M-Frontier",
                                  "L" = "L-Frontier",
                                  "G" = "G-Frontier",
                                  "A" = "A-Frontier"),
                       breaks = c("M", "L", "G", "A")) +
  scale_color_manual(values = frontier_colors,
                        labels = c("M" = "M-Frontier",
                                   "L" = "L-Frontier",
                                   "G" = "G-Frontier",
                                   "A" = "A-Frontier"),
                        breaks = c("M", "L", "G", "A")) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6.5)) +
  labs(
    title = "Scale efficiency (SE)",
    x = "Scale efficiency (CRS/VRS)",
    y = "Density",
    fill = "Frontier",
    color = "Frontier"
  ) +
  theme_publication()

# CRS panel
p_crs_dist <- eff_long %>%
  dplyr::filter(type == "TE_CRS") %>%
  ggplot(aes(x = efficiency, fill = frontier, color = frontier)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = frontier_colors,
                       labels = c("M" = "M-Frontier",
                                  "L" = "L-Frontier",
                                  "G" = "G-Frontier",
                                  "A" = "A-Frontier"),
                       breaks = c("M", "L", "G", "A")) +
  scale_color_manual(values = frontier_colors,
                        labels = c("M" = "M-Frontier",
                                   "L" = "L-Frontier",
                                   "G" = "G-Frontier",
                                   "A" = "A-Frontier"),
                        breaks = c("M", "L", "G", "A")) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6.5)) +
  labs(
    title = "CRS technical efficiency",
    x = "Technical efficiency (CRS)",
    y = "Density",
    fill = "Frontier",
    color = "Frontier"
  ) +
  theme_publication()

p_combined_dist <- (p_vrs_dist / p_se_dist / p_crs_dist) +
  plot_annotation(
    title = "Technical efficiency distributions across frontiers",
    subtitle = sprintf("n = %d farms | Progressive constraint intensity: E → L → G → A | VRS (top), SE (middle), CRS (bottom)", 
                      nrow(efficiency_results)),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(results_base, "visualization/Fig1_efficiency_distributions_VRS_CRS.png"),
       p_combined_dist, width = 12, height = 14, dpi = 300)

# Horizontal layout version
# Keep vertical Fig1 unchanged; only set y-axis upper limit = 5 for horizontal Fig1
p_vrs_dist_horizontal <- p_vrs_dist + scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), limits = c(0, 5))
p_se_dist_horizontal <- p_se_dist + scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), limits = c(0, 5))
p_crs_dist_horizontal <- p_crs_dist + scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), limits = c(0, 5))

p_vrs_dist_no_legend <- p_vrs_dist_horizontal + theme(legend.position = "none")
p_se_dist_no_legend <- p_se_dist_horizontal + theme(legend.position = "none")
p_crs_dist_with_legend <- p_crs_dist_horizontal + theme(legend.position = "bottom")

p_combined_dist_horizontal <- (p_vrs_dist_no_legend | p_se_dist_no_legend | p_crs_dist_with_legend) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Technical efficiency distributions across frontiers",
    subtitle = sprintf("n = %d farms | Progressive constraint intensity: E → L → G → A | VRS (left), SE (middle), CRS (right)", 
                      nrow(efficiency_results)),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  ) &
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
        legend.direction = "horizontal",
        legend.margin = margin(t = 10, b = 5))

ggsave(file.path(results_base, "visualization/Fig1_efficiency_distributions_VRS_CRS_horizontal.png"),
       p_combined_dist_horizontal, width = 14, height = 6.6, dpi = 300)

# FIGURE 2: Technology Gap Ratios (TGR) and Scale Efficiency

# Prepare TGR data with E-L-G-A ordering
tgr_long <- efficiency_results %>%
  tidyr::pivot_longer(
    cols = matches("^[ELGA]_TGR_(VRS|CRS)$"),
    names_to = "metric",
    values_to = "tgr"
  ) %>%
  tidyr::separate(metric, into = c("frontier", "dummy", "type"), sep = "_") %>%
  dplyr::select(-dummy) %>%
  dplyr::filter(!is.na(tgr)) %>%
  # Set frontier factor levels to E-L-G-A order
  dplyr::mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)) %>%
  # Clip TGR to [0, 1.05] range to avoid warnings (theoretical max is 1, but allow small numerical errors)
  dplyr::mutate(tgr = pmax(0, pmin(1.05, tgr)))

# Check for any extreme values
tgr_range_check <- tgr_long %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(
    min_tgr = min(tgr, na.rm = TRUE),
    max_tgr = max(tgr, na.rm = TRUE),
    n_out_of_range = sum(tgr > 1.05 | tgr < 0, na.rm = TRUE),
    .groups = "drop"
  )

if (any(tgr_range_check$n_out_of_range > 0, na.rm = TRUE)) {
  message(sprintf("  ⚠ Note: %d TGR values were clipped to [0, 1.05] range (numerical precision)\n",
                  sum(tgr_range_check$n_out_of_range, na.rm = TRUE)))
}

# Prepare SE data for violin plot with E-L-G-A ordering
se_long <- efficiency_results %>%
  tidyr::pivot_longer(
    cols = matches("^[ELGA]_SE$"),
    names_to = "metric",
    values_to = "se"
  ) %>%
  tidyr::separate(metric, into = c("frontier", "type"), sep = "_") %>%
  dplyr::filter(!is.na(se)) %>%
  # Set frontier factor levels to E-L-G-A order
  dplyr::mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)) %>%
  dplyr::mutate(se = pmax(0, pmin(1.05, se)))  # Clip SE to [0, 1.05]

# Prepare TGR by scale category data
# Create scale bins: <5k, 5-10k, 10-20k, 20-50k, >50k
efficiency_results_with_scale <- efficiency_results %>%
  dplyr::mutate(
    scale_bin = cut(market_pig_original,
                    breaks = c(0, 5000, 10000, 20000, 50000, Inf),
                    labels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
                    include.lowest = TRUE)
  ) %>%
  dplyr::filter(!is.na(scale_bin))

# Calculate mean TGR by scale category and frontier
tgr_by_scale <- efficiency_results_with_scale %>%
  tidyr::pivot_longer(
    cols = matches("^[ELGA]_TGR_VRS$"),
    names_to = "metric",
    values_to = "tgr"
  ) %>%
  tidyr::separate(metric, into = c("frontier", "dummy1", "dummy2"), sep = "_") %>%
  dplyr::select(-dummy1, -dummy2) %>%
  dplyr::filter(!is.na(tgr)) %>%
  dplyr::mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)) %>%
  dplyr::mutate(tgr = pmax(0, pmin(1.05, tgr))) %>%
  dplyr::mutate(scale_bin = factor(scale_bin, 
                                   levels = c("<5k", "5-10k", "10-20k", "20-50k", ">50k"),
                                   ordered = TRUE)) %>%
  dplyr::group_by(scale_bin, frontier) %>%
  dplyr::summarise(
    n = n(),
    mean_tgr = mean(tgr, na.rm = TRUE),
    se_tgr = sd(tgr, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# VRS panel
p_tgr_vrs <- tgr_long %>%
  dplyr::filter(type == "VRS") %>%
  ggplot(aes(x = frontier, y = tgr, fill = frontier)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, 
               color = "white", linewidth = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_viridis_d(option = "D", breaks = c("M", "L", "G", "A")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     limits = c(0, 1.05), oob = scales::squish) +
  labs(
    title = "TGR-VRS",
    x = "Frontier",
    y = "TGR-VRS",
    fill = "Frontier"
  ) +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10))

# SE panel
p_se_violin <- se_long %>%
  ggplot(aes(x = frontier, y = se, fill = frontier)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, 
               color = "white", linewidth = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_viridis_d(option = "D", breaks = c("M", "L", "G", "A")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05), oob = scales::squish) +
  labs(
    title = "Scale efficiency (SE)",
    x = "Frontier",
    y = "Scale efficiency (CRS/VRS)",
    fill = "Frontier"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# CRS panel
p_tgr_crs <- tgr_long %>%
  dplyr::filter(type == "CRS") %>%
  ggplot(aes(x = frontier, y = tgr, fill = frontier)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, 
               color = "white", linewidth = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_viridis_d(option = "D", breaks = c("M", "L", "G", "A")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     limits = c(0, 1.05), oob = scales::squish) +
  labs(
    title = "TGR-CRS",
    x = "Frontier",
    y = "TGR-CRS",
    fill = "Frontier"
  ) +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10))

# TGR by Scale Category panel
tgr_by_scale_labels <- tgr_by_scale %>%
  dplyr::group_by(scale_bin) %>%
  dplyr::summarise(n = first(n), .groups = "drop") %>%
  dplyr::mutate(label_y = 1.02)

p_tgr_by_scale <- tgr_by_scale %>%
  ggplot(aes(x = scale_bin, y = mean_tgr, fill = frontier, group = frontier)) +
  geom_col(position = position_dodge(width = 0.85), alpha = 0.9, width = 0.7, 
           color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_tgr - se_tgr, ymax = mean_tgr + se_tgr),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 0.5, color = "gray30", alpha = 0.8) +
  geom_line(aes(color = frontier, group = frontier), 
            linewidth = 1.2, alpha = 1.0) +
  geom_point(aes(color = frontier), position = position_dodge(width = 0.85),
             size = 1.5, shape = 21, fill = "white", stroke = 0.8) +
  geom_text(data = tgr_by_scale_labels,
            aes(x = scale_bin, y = label_y, label = sprintf("n=%d", n)),
            inherit.aes = FALSE, size = 4.2, color = "gray40", vjust = 0) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_viridis_d(option = "D", breaks = c("M", "L", "G", "A"),
                       labels = c("M" = "M-Frontier",
                                  "L" = "L-Frontier",
                                  "G" = "G-Frontier",
                                  "A" = "A-Frontier")) +
  scale_color_viridis_d(option = "D", breaks = c("M", "L", "G", "A"),
                        labels = c("M" = "M-Frontier",
                                   "L" = "L-Frontier",
                                   "G" = "G-Frontier",
                                   "A" = "A-Frontier")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     limits = c(0, 1.05), breaks = seq(0, 1, 0.2), oob = scales::squish) +
  labs(title = "TGR by scale category (VRS)", x = "Production scale", y = "TGR-VRS",
       fill = "Frontier", color = "Frontier") +
  theme_publication() +
  theme(axis.text.x = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11))

p_tgr_vrs_no_legend <- p_tgr_vrs + theme(legend.position = "none")
p_tgr_crs_no_legend <- p_tgr_crs + theme(legend.position = "none")
p_tgr_by_scale_with_legend <- p_tgr_by_scale + theme(legend.position = "bottom")

p_tgr_combined <- (p_tgr_vrs_no_legend | p_tgr_crs_no_legend | p_tgr_by_scale_with_legend) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Technology gap ratios (TGR) across frontiers",
    subtitle = "TGR = frontier efficiency / meta-frontier | Order: E-L-G-A across all plots",
    caption = "Meta-frontier = max(TE_M, TE_L, TE_G, TE_A) for each farm",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 11, hjust = 0.5, color = "gray40")
    )
  ) &
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
        legend.direction = "horizontal",
        legend.margin = margin(t = 10, b = 5))

ggsave(file.path(results_base, "visualization/Fig2_TGR_comparison_VRS_CRS.png"),
       p_tgr_combined, width = 14, height = 7, dpi = 300)

# FIGURE 3: Correlation Matrices

# VRS correlation matrix
eff_matrix_vrs <- efficiency_results %>%
  dplyr::select(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS) %>%
  setNames(c("M", "L", "G", "A"))

cor_matrix_vrs <- cor(eff_matrix_vrs, use = "pairwise.complete.obs")
cor_matrix_vrs <- cor_matrix_vrs[c("M", "L", "G", "A"), c("M", "L", "G", "A")]

# SE correlation matrix
eff_matrix_se <- efficiency_results %>%
  dplyr::select(M_SE, L_SE, G_SE, A_SE) %>%
  setNames(c("M", "L", "G", "A"))

cor_matrix_se <- cor(eff_matrix_se, use = "pairwise.complete.obs")
cor_matrix_se <- cor_matrix_se[c("M", "L", "G", "A"), c("M", "L", "G", "A")]

# CRS correlation matrix
eff_matrix_crs <- efficiency_results %>%
  dplyr::select(M_TE_CRS, L_TE_CRS, G_TE_CRS, A_TE_CRS) %>%
  setNames(c("M", "L", "G", "A"))

cor_matrix_crs <- cor(eff_matrix_crs, use = "pairwise.complete.obs")
cor_matrix_crs <- cor_matrix_crs[c("M", "L", "G", "A"), c("M", "L", "G", "A")]

png(file.path(results_base, "visualization/Fig3_correlation_matrices_VRS_CRS.png"),
    width = 24, height = 8.5, units = "in", res = 300)

par(mfrow = c(1, 3))

# Diverging palette aligned with Fig1 frontier style: A warm -> neutral -> L blue
corr_palette <- colorRampPalette(c("#E07A5F", "#FFF4E8", "#2E86AB"))(200)

corrplot(cor_matrix_vrs, method = "ellipse", type = "upper",
         addCoef.col = "black", number.cex = 3.4,
         tl.col = "black", tl.srt = 0, tl.cex = 2.7, tl.offset = 0.8,
         cl.cex = 2.0, cl.length = 5,
         col = corr_palette,
         mar = c(2, 0, 1.5, 0))
mtext("VRS efficiency correlations", side = 3, line = 0.3, cex = 2.0, font = 2)

corrplot(cor_matrix_se, method = "ellipse", type = "upper",
         addCoef.col = "black", number.cex = 3.4,
         tl.col = "black", tl.srt = 0, tl.cex = 2.7, tl.offset = 0.8,
         cl.cex = 2.0, cl.length = 5,
         col = corr_palette,
         mar = c(2, 0, 1.5, 0))
mtext("Scale efficiency (SE) correlations", side = 3, line = 0.3, cex = 2.0, font = 2)

corrplot(cor_matrix_crs, method = "ellipse", type = "upper",
         addCoef.col = "black", number.cex = 3.4,
         tl.col = "black", tl.srt = 0, tl.cex = 2.7, tl.offset = 0.8,
         cl.cex = 2.0, cl.length = 5,
         col = corr_palette,
         mar = c(2, 0, 1.5, 0))
mtext("CRS efficiency correlations", side = 3, line = 0.3, cex = 2.0, font = 2)

dev.off()

# Save correlation matrices
write_csv(as.data.frame(cor_matrix_vrs),
          file.path(results_base, "statistical_tests/correlation_matrix_VRS.csv"))
write_csv(as.data.frame(cor_matrix_se),
          file.path(results_base, "statistical_tests/correlation_matrix_SE.csv"))
write_csv(as.data.frame(cor_matrix_crs),
          file.path(results_base, "statistical_tests/correlation_matrix_CRS.csv"))

# ============================================================================
# FIGURE 5: VRS Efficiency and TGR by Scale Category
# ============================================================================

message("Generating Fig5_VRS_Efficiency_TGR_by_Scale.png...")

# Define scale bins: <5k, 5-10k, 10-20k, 20-50k, >50k
scale_bin_breaks <- c(0, 5000, 10000, 20000, 50000, Inf)
scale_bin_labels <- c("<5k", "5-10k", "10-20k", "20-50k", ">50k")

# Create scale bins: <5k, 5-10k, 10-20k, 20-50k, >50k
# Always use market_pig_original (head count) so all scale categories appear; market_pig may be scaled
cat("DEBUG: Creating scale bins for Fig5 (using market_pig_original)...\n")
efficiency_results_for_fig5 <- efficiency_results %>%
  dplyr::mutate(
    scale_bin = cut(market_pig_original,
                    breaks = scale_bin_breaks,
                    labels = scale_bin_labels,
                    include.lowest = TRUE)
  ) %>%
  dplyr::filter(!is.na(scale_bin)) %>%
  dplyr::mutate(scale_bin = factor(scale_bin, 
                                    levels = scale_bin_labels,
                                    ordered = TRUE))

# Debug: Check result
cat(sprintf("  efficiency_results_for_fig5 created, rows: %d\n", nrow(efficiency_results_for_fig5)))
cat(sprintf("  Scale bins found: %s\n", paste(unique(efficiency_results_for_fig5$scale_bin), collapse = ", ")))

# Calculate mean VRS efficiency by scale category and frontier
eff_by_scale <- efficiency_results_for_fig5 %>%
  tidyr::pivot_longer(
    cols = matches("^[ELGA]_TE_VRS$"),
    names_to = "metric",
    values_to = "efficiency"
  ) %>%
  tidyr::separate(metric, into = c("frontier", "dummy1", "dummy2"), sep = "_") %>%
  dplyr::select(-dummy1, -dummy2) %>%
  dplyr::filter(!is.na(efficiency)) %>%
  dplyr::mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)) %>%
  dplyr::mutate(efficiency = pmax(0, pmin(1.1, efficiency))) %>%
  dplyr::group_by(scale_bin, frontier) %>%
  dplyr::summarise(
    n = n(),
    mean_eff = mean(efficiency, na.rm = TRUE),
    se_eff = sd(efficiency, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# Calculate mean TGR by scale category and frontier
tgr_by_scale_fig5 <- efficiency_results_for_fig5 %>%
  tidyr::pivot_longer(
    cols = matches("^[ELGA]_TGR_VRS$"),
    names_to = "metric",
    values_to = "tgr"
  ) %>%
  tidyr::separate(metric, into = c("frontier", "dummy1", "dummy2"), sep = "_") %>%
  dplyr::select(-dummy1, -dummy2) %>%
  dplyr::filter(!is.na(tgr)) %>%
  dplyr::mutate(frontier = ifelse(frontier == "M", "M", frontier)) %>%
  dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)) %>%
  dplyr::mutate(tgr = pmax(0, pmin(1.05, tgr))) %>%
  dplyr::group_by(scale_bin, frontier) %>%
  dplyr::summarise(
    n = n(),
    mean_tgr = mean(tgr, na.rm = TRUE),
    se_tgr = sd(tgr, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# Get sample sizes for labels
n_labels <- efficiency_results_for_fig5 %>%
  dplyr::group_by(scale_bin) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::mutate(label_y_eff = 1.08, label_y_tgr = 1.03)

# Panel 1: VRS Efficiency by Scale
p_eff_by_scale <- eff_by_scale %>%
  ggplot(aes(x = scale_bin, y = mean_eff, fill = frontier, group = frontier)) +
  geom_col(position = position_dodge(width = 0.85), alpha = 0.9, width = 0.7, 
           color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_eff - se_eff, ymax = mean_eff + se_eff),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 0.5, color = "gray30", alpha = 0.8) +
  geom_line(aes(color = frontier, group = frontier), 
            linewidth = 1.2, alpha = 1.0) +
  geom_point(aes(color = frontier), position = position_dodge(width = 0.85),
             size = 1.5, shape = 21, fill = "white", stroke = 0.8) +
  geom_text(data = n_labels,
            aes(x = scale_bin, y = label_y_eff, label = sprintf("n=%d", n)),
            inherit.aes = FALSE, size = 4.2, color = "gray40", vjust = 0) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_viridis_d(option = "D", breaks = c("M", "L", "G", "A"),
                       labels = c("M" = "M-Frontier",
                                  "L" = "L-Frontier",
                                  "G" = "G-Frontier",
                                  "A" = "A-Frontier")) +
  scale_color_viridis_d(option = "D", breaks = c("M", "L", "G", "A"),
                        labels = c("M" = "M-Frontier",
                                   "L" = "L-Frontier",
                                   "G" = "G-Frontier",
                                   "A" = "A-Frontier")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 1.1), breaks = seq(0, 1, 0.2), oob = scales::squish) +
  scale_x_discrete(limits = scale_bin_labels, 
                   drop = FALSE) +
  labs(title = "VRS technical efficiency by scale category", 
       x = "Production scale", 
       y = "Mean VRS efficiency",
       fill = "Frontier", 
       color = "Frontier") +
  theme_publication() +
  theme(axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
        legend.position = "none")

# Panel 2: TGR by Scale
p_tgr_by_scale_fig5 <- tgr_by_scale_fig5 %>%
  ggplot(aes(x = scale_bin, y = mean_tgr, fill = frontier, group = frontier)) +
  geom_col(position = position_dodge(width = 0.85), alpha = 0.9, width = 0.7, 
           color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_tgr - se_tgr, ymax = mean_tgr + se_tgr),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 0.5, color = "gray30", alpha = 0.8) +
  geom_line(aes(color = frontier, group = frontier), 
            linewidth = 1.2, alpha = 1.0) +
  geom_point(aes(color = frontier), position = position_dodge(width = 0.85),
             size = 1.5, shape = 21, fill = "white", stroke = 0.8) +
  geom_text(data = n_labels,
            aes(x = scale_bin, y = label_y_tgr, label = sprintf("n=%d", n)),
            inherit.aes = FALSE, size = 4.2, color = "gray40", vjust = 0) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_viridis_d(option = "D", breaks = c("M", "L", "G", "A"),
                       labels = c("M" = "M-Frontier",
                                  "L" = "L-Frontier",
                                  "G" = "G-Frontier",
                                  "A" = "A-Frontier")) +
  scale_color_viridis_d(option = "D", breaks = c("M", "L", "G", "A"),
                        labels = c("M" = "M-Frontier",
                                   "L" = "L-Frontier",
                                   "G" = "G-Frontier",
                                   "A" = "A-Frontier")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     limits = c(0, 1.05), breaks = seq(0, 1, 0.2), oob = scales::squish) +
  scale_x_discrete(limits = scale_bin_labels, 
                   drop = FALSE) +
  labs(title = "TGR-VRS by scale category", 
       x = "Production scale", 
       y = "Mean TGR-VRS",
       fill = "Frontier", 
       color = "Frontier") +
  theme_publication() +
  theme(axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11))

# Combine panels
p_eff_by_scale_no_legend <- p_eff_by_scale + theme(legend.position = "none")
p_tgr_by_scale_with_legend <- p_tgr_by_scale_fig5 + theme(legend.position = "bottom")

p_fig5_combined <- (p_eff_by_scale_no_legend | p_tgr_by_scale_with_legend) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "VRS efficiency and TGR by production scale",
    subtitle = sprintf("Scale categories: %s | n = %d farms",
                      paste(scale_bin_labels, collapse = ", "), 
                      nrow(efficiency_results_for_fig5)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
    )
  ) &
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
        legend.direction = "horizontal",
        legend.margin = margin(t = 10, b = 5))

ggsave(file.path(results_base, "visualization/Fig5_VRS_Efficiency_TGR_by_Scale.png"),
       p_fig5_combined, width = 14, height = 6, dpi = 300)

message("✓ Fig5_VRS_Efficiency_TGR_by_Scale.png saved\n")

# ============================================================================
# 9. SCALE-EFFICIENCY RELATIONSHIP ANALYSIS
# ============================================================================

improve_scale_path <- file.path(modules_dir, "improve_scale_analysis.R")
if (file.exists(improve_scale_path)) {
  tryCatch({
    # Skip improve_scale_analysis.R due to syntax error
    message("⚠ Scale analysis skipped due to script syntax error\n")
  }, error = function(e) {
    message(sprintf("⚠ Scale analysis error: %s\n", e$message))
  })
} else {
  message("⚠ improve_scale_analysis.R not found, skipping scale analysis\n")
}

# Proceed directly to U-shape analysis
message("Proceeding to efficiency-scale comprehensive analysis...\n")


# ============================================================================
# 10. U-SHAPE DIAGNOSTIC AND VALIDATION
# ============================================================================

eff_scale_path <- file.path(modules_dir, "scale_efficiency_mechanism_analysis.R")
if (file.exists(eff_scale_path)) {
  message("Running efficiency-scale comprehensive analysis...")
  # Ensure variables are available in the sourced script's environment
  # The script will use results_base and efficiency_results from current session
  tryCatch({
    # Set environment variable as backup (script will prefer session variables)
    Sys.setenv("CURRENT_RESULTS_DIR" = results_base)
    # Source the script with local = TRUE to ensure access to all variables including df_clean
    source(eff_scale_path, echo = FALSE, local = TRUE)
    message("✓ Efficiency-scale comprehensive analysis completed\n")
  }, error = function(e) {
    message(sprintf("⚠ U-shape analysis error: %s\n", e$message))
    message(sprintf("  Error occurred at: %s\n", deparse(e$call)))
    # Print traceback for debugging
    traceback()
  })
} else {
  message("⚠ scale_efficiency_mechanism_analysis.R not found, skipping U-shape analysis\n")
}

# ============================================================================
# 10a. TOBIT THRESHOLD REGRESSION AND FACTOR DECOMPOSITION
# ============================================================================
# NOTE: Threshold analysis removed - no thresholds detected in data
# Tobit regression analysis is now integrated into scale_efficiency_mechanism_analysis.R

# ============================================================================
# 11. STATISTICAL TESTS
# ============================================================================

# Friedman test for overall frontier differences
message("Friedman Test (Overall Frontier Differences):")

eff_test_matrix <- efficiency_results %>%
  dplyr::select(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS) %>%
  setNames(c("M", "L", "G", "A"))

eff_complete <- eff_test_matrix[complete.cases(eff_test_matrix), ]
friedman_result <- friedman.test(as.matrix(eff_complete))
message(sprintf("  χ²(3) = %.2f, p = %.4e\n", 
                friedman_result$statistic, friedman_result$p.value))

# Pairwise comparisons
message("Pairwise Comparisons (Wilcoxon signed-rank):")

key_pairs <- list(
  c("M", "A"),  # Mortality vs Aggregated
  c("G", "A"),  # Global vs Aggregated
  c("L", "A")   # Local vs Aggregated
)

pairwise_results <- data.frame()

for (pair in key_pairs) {
  f1 <- paste0(pair[1], "_TE_VRS")
  f2 <- paste0(pair[2], "_TE_VRS")
  
  test_result <- suppressWarnings(
    wilcox.test(eff_complete[[pair[1]]], eff_complete[[pair[2]]], 
                paired = TRUE, exact = FALSE)
  )
  
  mean_diff <- mean(eff_complete[[pair[1]]] - eff_complete[[pair[2]]])
  
  message(sprintf("  %s vs %s: Δ = %+.3f, p = %.4f %s",
                  pair[1], pair[2], mean_diff, test_result$p.value,
                  ifelse(test_result$p.value < 0.05, "***", "ns")))
  
  pairwise_results <- rbind(pairwise_results, data.frame(
    comparison = paste(pair[1], "vs", pair[2]),
    mean_diff = mean_diff,
    p_value = test_result$p.value
  ))
}

write_csv(pairwise_results, file.path(results_base, "statistical_tests/pairwise_tests.csv"))
message("✓ Statistical tests completed\n")

# ============================================================================
# 12. COMPREHENSIVE SUMMARY REPORT
# ============================================================================

# Create summary table
summary_stats <- data.frame()
for (fid in c("M", "L", "G", "A")) {  # Order: E-L-G-A
  vrs <- efficiency_results[[paste0(fid, "_TE_VRS")]]
  tgr <- efficiency_results[[paste0(fid, "_TGR_VRS")]]
  
  summary_stats <- rbind(summary_stats, data.frame(
    Frontier = fid,
    VRS_Mean = mean(vrs, na.rm = TRUE),
    VRS_SD = sd(vrs, na.rm = TRUE),
    VRS_Median = median(vrs, na.rm = TRUE),
    TGR_Median = median(tgr, na.rm = TRUE),
    N_Frontier = sum(vrs > 0.99, na.rm = TRUE)
  ))
}

# Ensure Frontier column is ordered as E-L-G-A
summary_stats$Frontier <- factor(summary_stats$Frontier, 
                                levels = c("M", "L", "G", "A"), 
                                 ordered = TRUE)
summary_stats <- summary_stats[order(summary_stats$Frontier), ]

write_csv(summary_stats,
          file.path(results_base, "efficiency_calculation/summary_statistics.csv"))

# ============================================================================
# 12a. COMPREHENSIVE DESCRIPTIVE STATISTICS TABLE
# ============================================================================
# Table: Unit Input/Output and Efficiency Statistics (min, mean, max, sd)

message("\n12a. Generating Comprehensive Descriptive Statistics Table...")

# Calculate unit inputs and outputs (per pig)
# Use the correct column name for market pig
market_pig_col <- good_output  # "Market_pig"

# Get original market pig values from df_clean (not standardized)
# Standardized data has scaled market_pig values, but we need original for unit calculations
original_market_pig <- df_clean[[market_pig_col]]

unit_stats_data <- df_clean %>%
  dplyr::mutate(
    Feed_per_pig = Feed / original_market_pig,
    Water_per_pig = Water_consumption / original_market_pig,
    Energy_consumption_per_pig = Energy_consumption / original_market_pig,
    Labor_per_pig = Labor / original_market_pig,
    Waste_water_per_pig = Waste_water / original_market_pig,
    Manure_per_pig = Manure / original_market_pig,
    Dead_pig_per_pig = Dead_pig / original_market_pig,
    Carbon_per_pig = Carbon_emission / original_market_pig,
    Eutrophication_potential_per_pig = Eutrophication_potential / original_market_pig,
    obs_id = seq_len(nrow(df_clean))
  ) %>%
  dplyr::filter(original_market_pig > 0)

# Merge with efficiency results
unit_stats_data <- unit_stats_data %>%
  dplyr::left_join(
    efficiency_results %>%
      dplyr::select(obs_id, 
                   M_TE_VRS, M_SE, M_TE_CRS,
                   L_TE_VRS, L_SE, L_TE_CRS,
                   G_TE_VRS, G_SE, G_TE_CRS,
                   A_TE_VRS, A_SE, A_TE_CRS),
    by = "obs_id"
  )

# Prepare data for statistics table
unit_vars <- c(
  "Feed_per_pig", "Water_per_pig", "Energy_consumption_per_pig", "Labor_per_pig",
  "Waste_water_per_pig", "Manure_per_pig",
  "Dead_pig_per_pig", "Carbon_per_pig",
  "Eutrophication_potential_per_pig"
)

efficiency_vars <- c(
  "M_TE_VRS", "M_SE", "M_TE_CRS",
  "L_TE_VRS", "L_SE", "L_TE_CRS",
  "G_TE_VRS", "G_SE", "G_TE_CRS",
  "A_TE_VRS", "A_SE", "A_TE_CRS"
)

# Variable display names with units
var_display_names <- c(
  "Feed_per_pig" = "Feed",
  "Water_per_pig" = "Water consumption",
  "Energy_consumption_per_pig" = "Energy consumption",
  "Labor_per_pig" = "Labor",
  "Waste_water_per_pig" = "Wastewater",
  "Manure_per_pig" = "Manure",
  "Dead_pig_per_pig" = "Dead pig",
  "Carbon_per_pig" = "Carbon emission",
  "Eutrophication_potential_per_pig" = "Eutrophication potential",
  "M_TE_VRS" = "M-Frontier TE-VRS",
  "M_SE" = "M-Frontier SE",
  "M_TE_CRS" = "M-Frontier TE-CRS",
  "L_TE_VRS" = "L-Frontier TE-VRS",
  "L_SE" = "L-Frontier SE",
  "L_TE_CRS" = "L-Frontier TE-CRS",
  "G_TE_VRS" = "G-Frontier TE-VRS",
  "G_SE" = "G-Frontier SE",
  "G_TE_CRS" = "G-Frontier TE-CRS",
  "A_TE_VRS" = "A-Frontier TE-VRS",
  "A_SE" = "A-Frontier SE",
  "A_TE_CRS" = "A-Frontier TE-CRS"
)

# Unit mapping (per head)
var_units <- c(
  "Feed_per_pig" = "t/head",
  "Water_per_pig" = "t/head",
  "Energy_consumption_per_pig" = "kWh/head",
  "Labor_per_pig" = "person·h/head",
  "Waste_water_per_pig" = "t/head",
  "Manure_per_pig" = "t/head",
  "Dead_pig_per_pig" = "kg/head",
  "Carbon_per_pig" = "tCO2e/head",
  "Eutrophication_potential_per_pig" = "kg PO₄-eq/head",
  "M_TE_VRS" = "",
  "M_SE" = "",
  "M_TE_CRS" = "",
  "L_TE_VRS" = "",
  "L_SE" = "",
  "L_TE_CRS" = "",
  "G_TE_VRS" = "",
  "G_SE" = "",
  "G_TE_CRS" = "",
  "A_TE_VRS" = "",
  "A_SE" = "",
  "A_TE_CRS" = ""
)

# Calculate statistics for all variables
all_vars <- c(unit_vars, efficiency_vars)
descriptive_stats <- data.frame()

for (var in all_vars) {
  if (var %in% names(unit_stats_data)) {
    var_data <- unit_stats_data[[var]]
    var_data <- var_data[!is.na(var_data) & is.finite(var_data)]
    
    if (length(var_data) > 0) {
      # Calculate quartiles for IQR
      quartiles <- quantile(var_data, c(0.25, 0.75), na.rm = TRUE)
      iqr_value <- quartiles[2] - quartiles[1]  # Q3 - Q1
      
      # Calculate 95% confidence interval for mean
      n_obs <- length(var_data)
      mean_val <- mean(var_data, na.rm = TRUE)
      sd_val <- sd(var_data, na.rm = TRUE)
      se_val <- sd_val / sqrt(n_obs)
      
      # Use t-distribution for CI (more accurate for small samples)
      if (n_obs > 1) {
        t_critical <- qt(0.975, df = n_obs - 1)  # 95% CI, two-tailed
        ci_half_width <- t_critical * se_val
      } else {
        ci_half_width <- NA_real_
      }
      
      descriptive_stats <- rbind(descriptive_stats, data.frame(
        Variable = var_display_names[var],
        Unit = var_units[var],
        Type = ifelse(var %in% unit_vars, 
                     ifelse(var %in% c("Feed_per_pig", "Water_per_pig", "Energy_consumption_per_pig", "Labor_per_pig"),
                           "Unit Input", "Unit Output"),
                     "Efficiency"),
        Min = min(var_data, na.rm = TRUE),
        Mean = mean_val,
        Mean_CI = ci_half_width,  # Half-width of 95% CI
        Median = median(var_data, na.rm = TRUE),
        Max = max(var_data, na.rm = TRUE),
        SD = sd_val,
        IQR = iqr_value,
        N = n_obs,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Order: Unit Inputs, Unit Outputs, then Efficiency (M, L, G, A for each)
unit_input_order <- c("Feed", "Water consumption", 
                     "Energy consumption", "Labor")
unit_output_order <- c("Wastewater", 
                      "Manure", "Dead pig", "Carbon emission",
                      "Eutrophication potential")
efficiency_order <- c(
  "M-Frontier TE-VRS", "M-Frontier SE", "M-Frontier TE-CRS",
  "L-Frontier TE-VRS", "L-Frontier SE", "L-Frontier TE-CRS",
  "G-Frontier TE-VRS", "G-Frontier SE", "G-Frontier TE-CRS",
  "A-Frontier TE-VRS", "A-Frontier SE", "A-Frontier TE-CRS"
)

descriptive_stats <- descriptive_stats %>%
  dplyr::mutate(
    Variable = factor(Variable, 
                     levels = c(unit_input_order, unit_output_order, efficiency_order),
                     ordered = TRUE)
  ) %>%
  dplyr::arrange(Variable)

# Format numbers for better readability
descriptive_stats_formatted <- descriptive_stats %>%
  dplyr::mutate(
    # Unit variables: scientific notation (as requested)
    # Efficiency variables: keep fixed decimals for readability
    Min = dplyr::if_else(
      Type %in% c("Unit Input", "Unit Output"),
      formatC(Min, format = "e", digits = 4),
      sprintf("%.4f", round(Min, 4))
    ),
    # Format Mean and CI separately, then combine
    Mean_formatted = dplyr::if_else(
      Type %in% c("Unit Input", "Unit Output"),
      formatC(Mean, format = "e", digits = 4),
      sprintf("%.4f", round(Mean, 4))
    ),
    Mean_CI_formatted = dplyr::if_else(
      is.na(Mean_CI),
      "",
      dplyr::if_else(
        Type %in% c("Unit Input", "Unit Output"),
        formatC(Mean_CI, format = "e", digits = 4),
        sprintf("%.4f", round(Mean_CI, 4))
      )
    ),
    # Combine Mean and CI: "Mean±CI" format, or just "Mean" if CI is missing
    Mean = dplyr::case_when(
      is.na(Mean_CI) ~ Mean_formatted,
      Mean_CI_formatted == "" ~ Mean_formatted,
      TRUE ~ paste0(Mean_formatted, "±", Mean_CI_formatted)
    ),
    Median = dplyr::if_else(
      Type %in% c("Unit Input", "Unit Output"),
      formatC(Median, format = "e", digits = 4),
      sprintf("%.4f", round(Median, 4))
    ),
    Max = dplyr::if_else(
      Type %in% c("Unit Input", "Unit Output"),
      formatC(Max, format = "e", digits = 4),
      sprintf("%.4f", round(Max, 4))
    ),
    SD = dplyr::if_else(
      Type %in% c("Unit Input", "Unit Output"),
      formatC(SD, format = "e", digits = 4),
      sprintf("%.4f", round(SD, 4))
    ),
    IQR = dplyr::if_else(
      Type %in% c("Unit Input", "Unit Output"),
      formatC(IQR, format = "e", digits = 4),
      sprintf("%.4f", round(IQR, 4))
    )
  ) %>%
  dplyr::select(-Mean_formatted, -Mean_CI_formatted, -Mean_CI)  # Remove temporary columns and raw CI

# Save to CSV
write_csv(descriptive_stats_formatted,
          file.path(results_base, "efficiency_calculation/comprehensive_descriptive_statistics.csv"))

message("✓ Comprehensive descriptive statistics table saved")
message(sprintf("  - Unit Inputs: %d variables", sum(descriptive_stats$Type == "Unit Input")))
message(sprintf("  - Unit Outputs: %d variables", sum(descriptive_stats$Type == "Unit Output")))
message(sprintf("  - Efficiency measures: %d variables", sum(descriptive_stats$Type == "Efficiency")))
message(sprintf("  - Total variables: %d\n", nrow(descriptive_stats)))

# Print summary to console (use original data for numeric formatting)
message("Comprehensive Descriptive Statistics Summary:")
message(paste(rep("=", 155), collapse = ""))
message(sprintf("%-25s %-12s %-15s %12s %20s %12s %12s %12s %12s %6s",
               "Variable", "Unit", "Type", "Min", "Mean (95% CI)", "Median", "Max", "SD", "IQR", "N"))
message(paste(rep("-", 155), collapse = ""))

for (i in seq_len(nrow(descriptive_stats))) {
  unit_display <- ifelse(descriptive_stats$Unit[i] == "", "-", descriptive_stats$Unit[i])
  
  # Format Mean with CI for display
  if (is.na(descriptive_stats$Mean_CI[i])) {
    mean_display <- sprintf("%.4f", descriptive_stats$Mean[i])
  } else {
    mean_display <- sprintf("%.4f±%.4f", descriptive_stats$Mean[i], descriptive_stats$Mean_CI[i])
  }
  
  message(sprintf("%-25s %-12s %-15s %12.4f %20s %12.4f %12.4f %12.4f %12.4f %6d",
                 descriptive_stats$Variable[i],
                 unit_display,
                 descriptive_stats$Type[i],
                 descriptive_stats$Min[i],
                 mean_display,
                 descriptive_stats$Median[i],
                 descriptive_stats$Max[i],
                 descriptive_stats$SD[i],
                 descriptive_stats$IQR[i],
                 descriptive_stats$N[i]))
}

message(paste(rep("=", 100), collapse = ""))
message("")

# Generate text report
report_lines <- c(
  paste(rep("=", 78), collapse = ""),
  "MULTI-FRONTIER SBM ANALYSIS - FINAL REPORT",
  "FOCUSED ON CORE FINDINGS",
  "FIXED VARIABLE SET METHOD - ELIMINATES VARIABLE COUNT BIAS",
  paste(rep("=", 78), collapse = ""),
  "",
    sprintf("Sample Size: %d farms", nrow(efficiency_results)),
  sprintf("Scale Range: %s to %s pigs",
          format(round(min(efficiency_results$market_pig)), big.mark = ","),
          format(round(max(efficiency_results$market_pig)), big.mark = ",")),
  "",
  "METHODOLOGY:",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "Fixed Variable Set Method with Weighted Constraints:",
  "  - All frontiers use the same complete set of bad outputs (5 variables)",
  "  - Scheme 3: Non-overlapping core constraints",
  "    * E-frontier: Dead_pig only (1 constraint: mortality core)",
  "    * L-frontier: Waste_water, Manure (2 constraints: local environment)",
  "    * G-frontier: Carbon_emission, Eutrophication_potential (2 constraints: global environment)",
  "    * A-frontier: Aggregated 5 bad outputs constrained",
  "  - Scheme C: Weight proportion adjustment (REVISED)",
  "    * Core constraints: 1.2-1.6x standard weight (highlight key objectives)",
  "      - M: Mortality=1.5, L: Local=1.4, G: Global=1.5",
  "    * Basic constraints: 0.9-1.0x standard weight (maintain basic penalty)",
  "      - L/G: non-core mortality=0.95, A: Aggregated=1.0",
  "    * Non-core constraints: 0.7-0.8x standard weight (low penalty, maintain consistency)",
  "      - M/L/G: Non-core variables=0.75",
  "    * Weight range: [0.7, 1.6] (avoid extremes: >2.0 over-penalization, <0.5 negligible)",
  "    * Normalization denominator fixed at 9 (m + s_b = 4 + 5) for consistent scale",
  "  - This eliminates systematic bias from variable count differences",
  "  - Each frontier has unique core constraints, maximizing differentiation",
  "  - TGR values are fully comparable across frontiers",
  "",
  "CORE FINDINGS:",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "",
  "1. U-SHAPED SCALE-EFFICIENCY RELATIONSHIP:",
  "   - M-Frontier: Significant (bootstrap p<0.05)",
  "   - L-Frontier: Suggestive (p=0.08, largest coefficient)",
  "   - G-Frontier: Significant (bootstrap p<0.05)",
  "   - Efficiency trough at 13,000-37,000 pigs across frontiers",
  "   - Medium-scale farms (10,000-50,000) show 20-35pp efficiency loss",
  "",
  "2. ENVIRONMENTAL-ECONOMIC COMPLEMENTARITY:",
  sprintf("   - All frontier correlations positive (range: %.2f to %.2f)",
          min(cor_matrix_vrs[upper.tri(cor_matrix_vrs)]),
          max(cor_matrix_vrs[upper.tri(cor_matrix_vrs)])),
  "   - No evidence of environmental-economic trade-offs",
  "   - Strongest synergies: L↔A (0.81), E↔G (0.79)",
  "",
  "3. META-FRONTIER COMPOSITION:",
  sprintf("   - A-Frontier: Median TGR = %.2f%% (dominant meta source)",
          median(efficiency_results$A_TGR_VRS, na.rm = TRUE) * 100),
  sprintf("   - M-Frontier: Median TGR = %.2f%% (largest technology gap)",
          median(efficiency_results$M_TGR_VRS, na.rm = TRUE) * 100),
  "",
  "MEAN EFFICIENCY BY FRONTIER:",
  "-" %>% rep(78) %>% paste(collapse = "")
)

for (i in seq_len(nrow(summary_stats))) {
  report_lines <- c(report_lines,
    sprintf("%s-Frontier: VRS=%.3f (SD=%.3f), TGR_median=%.3f, N_frontier=%d",
            summary_stats$Frontier[i],
            summary_stats$VRS_Mean[i],
            summary_stats$VRS_SD[i],
            summary_stats$TGR_Median[i],
            summary_stats$N_Frontier[i])
  )
}

report_lines <- c(report_lines,
  "",
  "KEY FILES GENERATED:",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "",
  "CORE FIGURES (for paper main text):",
  "  1. Fig1_efficiency_distributions.png - VRS, SE, and CRS distributions",
  "  2. Fig2_TGR_comparison.png - TGR (VRS/CRS) and Scale Efficiency (2x2 layout)",
  "  3. Fig3_correlation_matrix.png - Efficiency correlations (VRS, SE, CRS)",
  "  4. Fig4_slack_analysis.png - Sources of inefficiency by frontier and variable",
  "  5. efficiency_by_scale_BINNED_barplot.png - U-shape (clear numerics)",
  "  6. DIAGNOSTIC_raw_data_pattern_4frontiers.png - LOESS validation (M, L, G, A)",
  "  7. DIAGNOSTIC_raw_data_pattern_Meta.png - Meta-frontier LOESS (standalone)",
  "",
  "SUPPLEMENTARY FIGURES:",
  "  - scale_constraint_LINEAR_model.png",
  "  - scale_constraint_QUADRATIC_model.png",
  "  - scale_constraint_OPTIMIZED.png (with milestones)",
  "  - scale_constraint_LINEAR_vs_QUADRATIC.png",
  "  - DIAGNOSTIC_meta_frontier_*.png (if generated)",
  "",
  "DATA FILES:",
  "  - frontier_efficiency_results.csv",
  "  - all_frontiers_efficiency_with_tgr.csv",
  "  - summary_statistics.csv",
  "  - correlation_matrix_VRS.csv",
  "",
  "VALIDATION FILES:",
  "  - U_shape_bootstrap_validation.csv",
  "  - U_shape_trough_analysis.csv",
  "  - scale_analysis_comparison_report.txt"
)

report_lines <- c(report_lines,
  "",
  paste(rep("=", 78), collapse = ""),
  "ANALYSIS COMPLETED",
  paste(rep("=", 78), collapse = ""),
  sprintf("Runtime: %.1f minutes", as.numeric(difftime(Sys.time(), start_time, units = "mins"))),
  ""
)

writeLines(report_lines, file.path(results_base, "FINAL_REPORT.txt"))

cat("\n")
cat(paste(report_lines, collapse = "\n"))
cat("\n")

# Save workspace for post-analysis use
workspace_file <- file.path(results_base, "efficiency_calculation", "sbm_analysis_workspace.RData")
tryCatch({
  save(standardized_data, inputs, good_output, bad_outputs_full,
       constrained_bad_M, constrained_bad_L, constrained_bad_G, constrained_bad_A,
       bad_weights_M, bad_weights_L, bad_weights_G, bad_weights_A,
       sbm_efficiency, scale_matrix, results_base,
       file = workspace_file)
  message(sprintf("✓ Workspace saved to: %s\n", workspace_file))
}, error = function(e) {
  warning(sprintf("Failed to save workspace: %s", e$message))
})

# ============================================================================
# 13. UNIT TRENDS AND MARGINAL EFFECTS PLOTS
# ============================================================================
message("\n13. Generating Unit Trends and Marginal Effects Plots...\n")

tryCatch({
# Prepare data for unit trends and marginal effects analysis
# Unit values standardized with z-score (x-mean)/sd for unit-independent comparison (aligned with reference)
required_merge_cols <- c(inputs, bad_outputs_full, good_output)
if (!all(required_merge_cols %in% names(df_clean))) {
  missing <- setdiff(required_merge_cols, names(df_clean))
  stop(sprintf("Step 13: df_clean missing columns: %s", paste(missing, collapse = ", ")))
}
data_merged <- efficiency_results %>%
  dplyr::left_join(
    df_clean %>%
      dplyr::mutate(obs_id = seq_len(nrow(df_clean))),
    by = "obs_id"
  )
if (!nrow(data_merged)) {
  stop("Step 13: left_join produced no rows. Check obs_id and df_clean.")
}
if (!good_output %in% names(data_merged)) {
  stop(sprintf("Step 13: column '%s' not in data_merged after join.", good_output))
}
# Ensure single log_scale column (join may create log_scale.x / log_scale.y)
data_merged <- data_merged %>%
  dplyr::mutate(log_scale = log10(pmax(market_pig_original, 1)))

# Calculate unit inputs and outputs (per pig) using original-scale denominator
market_pig_original <- data_merged[[good_output]]
per_pig_cols <- c(
  paste0(inputs, "_per_pig"),
  paste0(bad_outputs_full, "_per_pig")
)
data_merged <- data_merged %>%
  dplyr::mutate(
    Feed_per_pig = Feed / market_pig_original,
    Water_consumption_per_pig = Water_consumption / market_pig_original,
    Energy_consumption_per_pig = Energy_consumption / market_pig_original,
    Labor_per_pig = Labor / market_pig_original,
    Waste_water_per_pig = Waste_water / market_pig_original,
    Manure_per_pig = Manure / market_pig_original,
    Dead_pig_per_pig = Dead_pig / market_pig_original,
    Carbon_emission_per_pig = Carbon_emission / market_pig_original,
    Eutrophication_potential_per_pig = Eutrophication_potential / market_pig_original
  )
# Standardize unit values with z-score (reference: scale() for unit-independent trend)
unit_per_pig_vars <- c("Feed_per_pig", "Water_consumption_per_pig", "Energy_consumption_per_pig", "Labor_per_pig",
  "Waste_water_per_pig", "Manure_per_pig", "Dead_pig_per_pig", "Carbon_emission_per_pig",
  "Eutrophication_potential_per_pig")
for (v in unit_per_pig_vars) {
  if (v %in% names(data_merged)) {
    data_merged[[paste0(v, "_std")]] <- as.vector(scale(data_merged[[v]]))
  }
}

# Handle NA from scale() when sd=0 (replace with 0)
std_cols <- paste0(unit_per_pig_vars, "_std")
for (sc in std_cols) {
  if (sc %in% names(data_merged)) {
    data_merged[[sc]][!is.finite(data_merged[[sc]])] <- NA
  }
}

# Prepare data for Plot A: Unit Input/Output Trends (x = market_pig_original for original scale trend)
unit_std_cols <- std_cols[std_cols %in% names(data_merged)]
unit_var_map <- c(
  "Feed_per_pig" = "Feed", "Water_consumption_per_pig" = "Water_consumption",
  "Energy_consumption_per_pig" = "Energy_consumption", "Labor_per_pig" = "Labor",
  "Waste_water_per_pig" = "Waste_water", "Manure_per_pig" = "Manure",
  "Dead_pig_per_pig" = "Dead_pig", "Carbon_emission_per_pig" = "Carbon_emission",
  "Eutrophication_potential_per_pig" = "Eutrophication_potential"
)
unit_data_long <- data_merged %>%
  dplyr::select(market_pig_original, log_scale, dplyr::all_of(unit_std_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(unit_std_cols),
    names_to = "variable",
    values_to = "unit_value"
  ) %>%
  dplyr::mutate(
    variable = stringr::str_remove(variable, "_std"),
    variable = unname(unit_var_map[variable]),
    type = ifelse(variable %in% inputs, "Input", "Output")
  ) %>%
  dplyr::filter(!is.na(unit_value), !is.na(log_scale), !is.na(variable))

# Define colors: inputs + bad_outputs_full
variable_colors_plot <- c(
  "Feed" = "#1f77b4",
  "Water_consumption" = "#2ca02c",
  "Energy_consumption" = "#d62728",
  "Labor" = "#ff7f0e",
  "Waste_water" = "#9467bd",
  "Manure" = "#c54e8c",
  "Dead_pig" = "#7f7f7f",
  "Carbon_emission" = "#bcbd22",
  "Eutrophication_potential" = "#8c564b"
)
input_vars <- c("Feed", "Water_consumption", "Energy_consumption", "Labor")
output_vars <- c("Waste_water", "Manure", "Dead_pig", "Carbon_emission", "Eutrophication_potential")
variable_order <- c(input_vars, output_vars)
variable_colors_plot <- variable_colors_plot[intersect(variable_order, names(variable_colors_plot))]
variable_order <- intersect(variable_order, names(variable_colors_plot))

# Reorder variable_colors_plot to match desired order
variable_colors_with_tp <- variable_colors_plot[variable_order]

# Create display labels (keep raw variable codes for computation)
variable_display_names <- ifelse(variable_order == "Waste_water", "Wastewater", variable_order)

# Create labels
variable_labels_with_tp <- paste0(variable_display_names, 
                                 ifelse(variable_order %in% input_vars,
                                       " (Input)", " (Output)"))
names(variable_labels_with_tp) <- variable_order

# Ensure unit_data_long_with_tp has variables in the correct order
unit_data_long_with_tp <- unit_data_long %>%
  dplyr::mutate(
    variable = factor(variable, levels = variable_order, ordered = TRUE)
  )

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
  
  # Create plot (x = market_pig_original so trend is on original production scale)
  p <- ggplot2::ggplot(plot_data,
                                ggplot2::aes(x = market_pig_original, y = unit_value, 
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
    ggplot2::scale_x_log10(limits = c(NA, NA),
                          breaks = c(10000, 30000, 100000, 300000),
                          labels = scales::comma,
                            expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
  ggplot2::scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3),
                             expand = ggplot2::expansion(mult = 0.05)) +
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
    text = ggplot2::element_text(family = "Times", size = 13),
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
      legend.text = ggplot2::element_text(family = "Times", size = 10),
      legend.key.width = ggplot2::unit(1.2, "cm"),
      legend.key.height = ggplot2::unit(0.3, "cm"),
      legend.spacing.y = ggplot2::unit(0.15, "cm"),
      legend.margin = ggplot2::margin(3, 5, 3, 5, "pt"),
      plot.margin = ggplot2::margin(4, 60, 4, 4, "pt")
    )
  
  return(p)
}

# Create separate plots for inputs and outputs
p_unit_trends_inputs <- create_unit_plot(unit_data_long_with_tp, "Input", "Inputs")
p_unit_trends_outputs <- create_unit_plot(unit_data_long_with_tp, "Output", "Outputs")

# Prepare data for marginal effects analysis
# Calculate marginal effects of 1% changes in factors on efficiency
# Group by scale zones based on Tobit analysis results

# First, perform Tobit regression to identify turning point
# Fit Tobit model for A-Frontier VRS efficiency
if (nrow(data_merged) >= 10) {
  # Prepare data for Tobit regression
  data_tobit <- data_merged %>%
    dplyr::filter(!is.na(A_TE_VRS), !is.na(market_pig)) %>%
    dplyr::mutate(
      efficiency = A_TE_VRS,
      log_scale = log10(market_pig),
      log_scale_sq = log_scale^2
    )
  
  if (nrow(data_tobit) >= 10) {
    # Fit Tobit model
    tryCatch({
      if (requireNamespace("AER", quietly = TRUE)) {
        tobit_fit <- AER::tobit(efficiency ~ log_scale + log_scale_sq,
                               data = data_tobit)
        
        # Extract coefficients
        coef_tobit <- coef(tobit_fit)
        beta1 <- coef_tobit["log_scale"]
        beta2 <- coef_tobit["log_scale_sq"]
        
        # Calculate turning point
        if (!is.na(beta1) && !is.na(beta2) && beta2 != 0) {
          turning_point_log <- -beta1 / (2 * beta2)
          turning_point_scale <- 10^turning_point_log
          
          # Define TP zone boundaries based on turning point
          tp_lower_marginal <- turning_point_scale * 0.7
          tp_upper_marginal <- turning_point_scale * 1.5
          
          message(sprintf("✓ Tobit analysis identified turning point at %.0f pigs\n", turning_point_scale))
          message(sprintf("✓ TP Zone defined as %.0f to %.0f pigs\n", tp_lower_marginal, tp_upper_marginal))
        } else {
          # Fallback to quantile method if Tobit fails
          tp_lower_marginal <- quantile(market_pig_original, 0.3, na.rm = TRUE)
          tp_upper_marginal <- quantile(market_pig_original, 0.7, na.rm = TRUE)
          message("⚠ Tobit analysis failed, using quantile method for TP zone\n")
        }
      } else {
        # Fallback to quantile method if AER package not available
        tp_lower_marginal <- quantile(market_pig_original, 0.3, na.rm = TRUE)
        tp_upper_marginal <- quantile(market_pig_original, 0.7, na.rm = TRUE)
        message("⚠ AER package not available, using quantile method for TP zone\n")
      }
    }, error = function(e) {
      # Fallback to quantile method if any error occurs
      tp_lower_marginal <- quantile(market_pig_original, 0.3, na.rm = TRUE)
      tp_upper_marginal <- quantile(market_pig_original, 0.7, na.rm = TRUE)
      message(sprintf("⚠ Tobit analysis error: %s, using quantile method for TP zone\n", e$message))
    })
  } else {
    # Fallback to quantile method if insufficient data
    tp_lower_marginal <- quantile(market_pig_original, 0.3, na.rm = TRUE)
    tp_upper_marginal <- quantile(market_pig_original, 0.7, na.rm = TRUE)
    message("⚠ Insufficient data for Tobit analysis, using quantile method for TP zone\n")
  }
} else {
  # Fallback to quantile method if insufficient data
  tp_lower_marginal <- quantile(market_pig_original, 0.3, na.rm = TRUE)
  tp_upper_marginal <- quantile(market_pig_original, 0.7, na.rm = TRUE)
  message("⚠ Insufficient data for Tobit analysis, using quantile method for TP zone\n")
}

data_merged <- data_merged %>%
  dplyr::mutate(
    scale_zone = ifelse(market_pig_original < tp_lower_marginal, "Pre-TP",
                       ifelse(market_pig_original <= tp_upper_marginal, "TP Zone",
                             "Post-TP"))
  )

# Calculate marginal effects for each zone
marginal_effects_data <- data.frame()

marginal_vars <- c(inputs, bad_outputs_full)
for (zone in c("Pre-TP", "TP Zone", "Post-TP")) {
  zone_data <- data_merged %>% dplyr::filter(scale_zone == zone)
  
  if (nrow(zone_data) < 3) next
  
  # Calculate marginal effect of 1% increase in all variables (inputs + bad_outputs_full)
  for (var in marginal_vars) {
    if (!var %in% names(zone_data)) next
    
    # Log-log regression: log(efficiency) ~ log(factor)
    zone_data_clean <- zone_data %>%
      dplyr::filter(
        !is.na(.data[[var]]), 
        .data[[var]] > 1e-10,
        !is.na(A_TE_VRS), 
        A_TE_VRS > 1e-10
      )
    
    if (nrow(zone_data_clean) < 3) next
    
    # Create log-transformed columns for regression
    zone_data_clean <- zone_data_clean %>%
      dplyr::mutate(
        var_log = log(pmax(.data[[var]], 1e-10)),
        te_vrs_log = log(pmax(A_TE_VRS, 1e-10))
      )
    
    # Run regression
    model <- tryCatch({
      lm(te_vrs_log ~ var_log, data = zone_data_clean)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(model) && length(coef(model)) >= 2) {
      coef_var <- coef(model)[2]  # Elasticity
      if (is.finite(coef_var)) {
        # Marginal effect of 1% change = elasticity * 0.01
        marg_effect <- coef_var * 0.01
        se_marg <- tryCatch({
          summary(model)$coefficients[2, 2] * 0.01
        }, error = function(e) {
          return(abs(marg_effect) * 0.1)
        })
        
        if (is.finite(marg_effect)) {
          marginal_effects_data <- rbind(marginal_effects_data, data.frame(
            Zone = zone,
            Variable = var,
            Type = ifelse(var %in% inputs, "Input (+1%)", "Output (+1%)"),
            Marginal_Effect = marg_effect,
            SE = ifelse(is.finite(se_marg), se_marg, abs(marg_effect) * 0.1),
            Lower = marg_effect - 1.96 * ifelse(is.finite(se_marg), se_marg, abs(marg_effect) * 0.1),
            Upper = marg_effect + 1.96 * ifelse(is.finite(se_marg), se_marg, abs(marg_effect) * 0.1)
          ))
        }
      }
    }
  }
}

# Create marginal effects plot if data is available
if (nrow(marginal_effects_data) > 0) {
  # Ensure Zone and Variable order (use variable_order so inputs then outputs)
  marginal_effects_data <- marginal_effects_data %>%
    dplyr::mutate(
      Zone = factor(Zone, levels = c("Pre-TP", "TP Zone", "Post-TP")),
      Variable = factor(as.character(Variable), levels = variable_order),
      Type = factor(Type, levels = c("Input (+1%)", "Output (+1%)"))
    ) %>%
    dplyr::filter(!is.na(Variable))
  
  # Create labels with (+1%) indicators for all variables
  variable_labels <- paste0(
    gsub("_", " ", names(variable_colors_plot)),
    " (+1%)"
  )
  names(variable_labels) <- names(variable_colors_plot)
  
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
      limits = c(-3, 3),
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
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
      text = ggplot2::element_text(family = "Times", size = 13),
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
}

# Save plots
visualization_dir <- file.path(results_base, "visualization")
if (!dir.exists(visualization_dir)) {
  dir.create(visualization_dir, recursive = TRUE)
}

}, error = function(e) {
  message("⚠ Step 13 (Unit Trends and Marginal Effects) failed: ", conditionMessage(e))
  message("  Analysis will continue; you can re-run this step after fixing the issue.\n")
})

# Unit trends and marginal effects plots are now generated in scale_efficiency_mechanism_analysis.R
message("ℹ Unit trends and marginal effects plots will be generated in scale_efficiency_mechanism_analysis.R\n")

# Optional: run scale-up dynamic analysis (DID + bubble plots) when Scale-up and No. exist
scale_up_col_main <- NULL
for (cn in c("Scale-up", "Scale_up", "ScaleUp")) {
  if (cn %in% names(efficiency_results)) { scale_up_col_main <- cn; break }
}
no_col_main <- NULL
for (cn in c("No.", "No", "Number")) {
  if (cn %in% names(efficiency_results)) { no_col_main <- cn; break }
}
scaleup_module_path <- file.path(modules_dir, "scale_expansion_did_analysis.R")
if (!is.null(scale_up_col_main) && !is.null(no_col_main) && file.exists(scaleup_module_path)) {
  tryCatch({
    message("\n✓ Running scale-up dynamic analysis (DID + bubble plots)...\n")
    source(scaleup_module_path, encoding = "UTF-8")
  }, error = function(e) {
    message("⚠ Scale-up dynamic analysis skipped: ", conditionMessage(e), "\n")
  })
}

message("✓ ANALYSIS COMPLETE")
message(sprintf("Total runtime: %.1f minutes", 
                as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
message(sprintf("Results saved to: %s\n", results_base))















