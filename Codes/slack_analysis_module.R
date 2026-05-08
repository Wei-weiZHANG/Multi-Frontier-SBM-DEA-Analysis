# ============================================================================
# Slack Analysis Module: Identify Sources of Inefficiency
# ============================================================================
# 
# This module analyzes the sources of inefficiency by calculating slack variables
# for inputs and bad outputs. It identifies which variables contribute most
# to efficiency gaps between individual frontiers and the meta-frontier.
#
# ============================================================================

# Required packages for visualization
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!require("viridis", quietly = TRUE)) install.packages("viridis")
if (!require("scales", quietly = TRUE)) install.packages("scales")
if (!require("readr", quietly = TRUE)) install.packages("readr")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")

library(ggplot2)
library(patchwork)
library(viridis)
library(scales)
library(readr)
library(dplyr)
library(tidyr)

# SBM function with slack extraction
# Args:
#   X: Input matrix (n × m)
#   Y: Good output matrix (n × s_g)
#   bad: Bad output matrix (n × s_b) - ALL bad outputs (full set)
#   constrained_bad_indices: Indices of bad outputs to be constrained (penalized in objective)
#                            If NULL, all bad outputs are constrained (original behavior)
#   bad_weights: Weight vector for bad outputs (length = s_b)
#                If NULL, all constrained bad outputs use standard weight (1.0)
#   RTS: Returns to scale ("vrs" or "crs")
#   orientation: "output" (only output-oriented implemented)
sbm_efficiency_with_slack <- function(X, Y, bad = NULL, constrained_bad_indices = NULL, 
                                      bad_weights = NULL, RTS = "vrs", orientation = "output") {
  
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
  
  normalization_denom <- m + s_b
  
  # Initialize results
  eff_scores <- numeric(n)
  input_slack <- matrix(0, nrow = n, ncol = m)
  bad_output_slack <- if (s_b > 0) matrix(0, nrow = n, ncol = s_b) else matrix(0, nrow = n, ncol = 0)
  good_output_slack <- matrix(0, nrow = n, ncol = s_g)
  
  # Progress indicator for large datasets
  show_progress <- n > 20
  progress_interval <- if (n > 100) 10 else if (n > 50) 5 else 1
  
  for (k in seq_len(n)) {
    # Progress update
    if (show_progress && (k %% progress_interval == 0 || k == n)) {
      cat(sprintf("\r    Progress: %d/%d DMUs (%.1f%%)", k, n, 100*k/n))
      flush.console()  # Force output to appear immediately
      if (k == n) cat("\n")
    }
    x_k <- pmax(X[k, ], 1e-8)
    y_k <- pmax(Y[k, ], 1e-8)
    b_k <- if (s_b > 0) pmax(bad[k, ], 1e-8) else numeric(0)
    
    num_vars <- 1 + n + m + s_g + s_b
    
    if (orientation == "output") {
      obj_coef <- numeric(num_vars)
      obj_coef[1] <- 1
      
      normalization_denom <- m + s_b
      for (i in 1:m) {
        obj_coef[1 + n + i] <- -1 / (normalization_denom * x_k[i])
      }
      
      # Bad output penalties: ALL bad outputs are penalized (fixed variable set)
      # Weights differentiate core vs non-core constraints
      if (s_b > 0) {
        for (i in 1:s_b) {
          if (i %in% constrained_bad_indices) {
            # Use weight parameter if provided (default: 1.0 for standard weight)
            if (!is.null(bad_weights) && length(bad_weights) >= i && !is.na(bad_weights[i]) && bad_weights[i] > 0) {
              weight_i <- bad_weights[i]
            } else {
              weight_i <- 1.0  # Default standard weight
            }
            
            # Weighted penalty term: weight × standard penalty
            obj_coef[1 + n + m + s_g + i] <- -weight_i / (normalization_denom * b_k[i])
          } else {
            # Should not happen with fixed variable set, but handle gracefully
            obj_coef[1 + n + m + s_g + i] <- 0
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
        for (j in 1:s_b) {
          constr[row_idx, 2:(1+n)] <- bad[, j]
          constr[row_idx, 1 + n + m + s_g + j] <- 1
          constr[row_idx, 1] <- -b_k[j]
          dir_vec[row_idx] <- "=="
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
      
      # Solve LP with proper status handling and optimization settings
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
        list(status = 99, objval = NA, solution = rep(NA, num_vars), message = e$message)
      })
      
      # Handle LP status: 0=success, 1=infeasible, 2=unbounded, 99=error
      if (lp_result$status == 0 && !is.na(lp_result$objval)) {
        eff_scores[k] <- pmax(0, pmin(1, lp_result$objval))
        solution <- if (!is.null(lp_result$solution) && length(lp_result$solution) == num_vars) {
          lp_result$solution
        } else {
          rep(NA, num_vars)  # Invalid solution
        }
        
        # Extract slack variables from solution
        # Solution structure: [theta, lambda_1...lambda_n, s_x_1...s_x_m, s_y_1...s_y_sg, s_b_1...s_b_sb]
        
        # Input slacks (position: 1+n+1 to 1+n+m)
        if (length(solution) >= 1 + n + m) {
          input_slack[k, ] <- pmax(0, solution[(1 + n + 1):(1 + n + m)])
        }
        
        # Good output slacks (position: 1+n+m+1 to 1+n+m+s_g)
        if (length(solution) >= 1 + n + m + s_g) {
          good_output_slack[k, ] <- pmax(0, solution[(1 + n + m + 1):(1 + n + m + s_g)])
        }
        
        # Bad output slacks (position: 1+n+m+s_g+1 to 1+n+m+s_g+s_b)
        if (s_b > 0 && length(solution) >= 1 + n + m + s_g + s_b) {
          bad_output_slack[k, ] <- pmax(0, solution[(1 + n + m + s_g + 1):(1 + n + m + s_g + s_b)])
        }
      } else {
        # Handle failure cases
        eff_scores[k] <- NA_real_
        solution <- rep(NA, num_vars)
        
        # Log status for debugging (only first few failures to avoid spam)
        if (k <= 5) {
          if (lp_result$status == 1) {
            warning(sprintf("DMU %d: LP infeasible", k))
          } else if (lp_result$status == 2) {
            warning(sprintf("DMU %d: LP unbounded", k))
          } else if (lp_result$status == 99) {
            warning(sprintf("DMU %d: LP error - %s", k, 
                           ifelse(!is.null(lp_result$message), lp_result$message, "unknown")))
          }
        }
      }
    } else {
      stop("Input-oriented SBM not implemented. Use orientation='output'")
    }
  }
  
  # Return results
  result_list <- list(
    efficiency = eff_scores,
    input_slack = input_slack,
    bad_output_slack = bad_output_slack,
    good_output_slack = good_output_slack
  )
  
  return(result_list)
}

# Scaling function
scale_matrix <- function(M) {
  scaling_factors <- apply(M, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
  sweep(M, 2, scaling_factors, "/")
}

# Scale group definition (same as Monte Carlo sensitivity: market pig heads per year)
SCALE_BREAKS <- c(0, 5000, 10000, 20000, 50000, Inf)
SCALE_LABELS <- c("<5k", "5-10k", "10-20k", "20-50k", ">50k")

# Main slack analysis function
# Modified to support fixed variable set method
# Now reports slack for ALL bad outputs with weight > 0 (not just core constraints)
# Slack is also aggregated by scale group (same grouping as Monte Carlo).
analyze_slack_by_frontier <- function(efficiency_results, df_clean, inputs, good_output, 
                                       bad_outputs_list, results_base,
                                       bad_outputs_full = NULL, constrained_bad_indices_list = NULL,
                                       bad_weights_list = NULL) {
  
  message(paste(rep("=", 78), collapse = ""))
  message("SLACK ANALYSIS: Sources of Inefficiency by Frontier")
  message("Using Fixed Variable Set Method")
  message(paste(rep("=", 78), collapse = ""))
  
  # Add scale_bin from raw scale variable (Market_pig original, not log or scaled)
  # Prefer Market_pig_original when present (standardized_data); otherwise use good_output (raw)
  scale_col_raw <- NULL
  if ("Market_pig_original" %in% names(df_clean)) {
    scale_col_raw <- "Market_pig_original"
  } else {
    scale_col_raw <- if (length(good_output) == 1L) good_output else good_output[1L]
  }
  if (!"scale_bin" %in% names(df_clean) && !is.null(scale_col_raw) && scale_col_raw %in% names(df_clean)) {
    df_clean$scale_bin <- cut(
      as.numeric(df_clean[[scale_col_raw]]),
      breaks = SCALE_BREAKS,
      labels = SCALE_LABELS,
      include.lowest = TRUE,
      right = FALSE
    )
    df_clean$scale_bin <- factor(df_clean$scale_bin, levels = SCALE_LABELS, ordered = TRUE)
    message(sprintf("  Scale groups from raw variable '%s' (not log/scaled): %s", scale_col_raw, paste(SCALE_LABELS, collapse = ", ")))
  } else if ("scale_bin" %in% names(df_clean)) {
    df_clean$scale_bin <- factor(df_clean$scale_bin, levels = SCALE_LABELS, ordered = TRUE)
  }
  
  # Scale inputs and outputs
  X_inputs <- as.matrix(df_clean[, inputs])
  Y_good <- as.matrix(df_clean[, good_output, drop = FALSE])
  X_scaled <- scale_matrix(X_inputs)
  Y_scaled <- scale_matrix(Y_good)
  
  # Fixed Variable Set Method: Use full bad output set for all frontiers
  if (!is.null(bad_outputs_full) && !is.null(constrained_bad_indices_list)) {
    # Use fixed variable set method
    Z_bad_full <- as.matrix(df_clean[, bad_outputs_full, drop = FALSE])
    Z_scaled_full <- scale_matrix(Z_bad_full)
    
    message("\nFixed Variable Set Method:")
    message(sprintf("  All frontiers use the same bad output set: %s", 
                    paste(bad_outputs_full, collapse = ", ")))
    
    # Define frontiers with constrained indices
    frontiers <- list(
      M = list(name = "Mortality", constrained_indices = constrained_bad_indices_list$E),
      L = list(name = "Local Environmental", constrained_indices = constrained_bad_indices_list$L),
      G = list(name = "Global", constrained_indices = constrained_bad_indices_list$G),
      A = list(name = "Aggregated", constrained_indices = constrained_bad_indices_list$A)
    )
    
    # For slack analysis, we still need to know which bad outputs are constrained
    # to properly label the slack variables
    constrained_bad_names_list <- list(
      E = bad_outputs_full[constrained_bad_indices_list$E],
      L = bad_outputs_full[constrained_bad_indices_list$L],
      G = bad_outputs_full[constrained_bad_indices_list$G],
      A = bad_outputs_full[constrained_bad_indices_list$A]
    )
  } else {
    # Fallback to original method
    message("\nUsing original method (variable set may differ by frontier)")
    Z_bad_full <- NULL
    Z_scaled_full <- NULL
    
    frontiers <- list(
      M = list(name = "Mortality", bad = bad_outputs_list$E),
      L = list(name = "Local Environmental", bad = bad_outputs_list$L),
      G = list(name = "Global", bad = bad_outputs_list$G),
      A = list(name = "Aggregated", bad = bad_outputs_list$A)
    )
    constrained_bad_names_list <- bad_outputs_list
  }
  
  # Initialize results storage
  slack_results_list <- list()
  
  # Calculate slack for each frontier
  for (fid in names(frontiers)) {
    message(sprintf("\nAnalyzing %s-Frontier (%s)...", fid, frontiers[[fid]]$name))
    
    if (!is.null(Z_scaled_full) && !is.null(constrained_bad_indices_list)) {
      # Fixed Variable Set Method: Use full bad output set with constrained indices
      message(sprintf("  Using fixed variable set: all %d bad outputs", length(bad_outputs_full)))
      message(sprintf("  Constrained bad outputs: %s (indices: %s)",
                      paste(bad_outputs_full[frontiers[[fid]]$constrained_indices], collapse = ", "),
                      paste(frontiers[[fid]]$constrained_indices, collapse = ", ")))
      
      # Get weights for this frontier (if available)
      bad_weights <- NULL
      if (!is.null(bad_weights_list) && fid %in% names(bad_weights_list)) {
        bad_weights <- bad_weights_list[[fid]]
        # Ensure weights are in the correct order matching bad_outputs_full
        # Handle both named and unnamed vectors
        if (length(bad_weights) == length(bad_outputs_full)) {
          # If named vector, reorder to match bad_outputs_full
          if (!is.null(names(bad_weights))) {
            bad_weights <- bad_weights[bad_outputs_full]
          }
          # Ensure no NA values
          if (any(is.na(bad_weights))) {
            warning(sprintf("  ⚠ Some weights are NA for %s-frontier, replacing with 0", fid))
            bad_weights[is.na(bad_weights)] <- 0
          }
          message(sprintf("  Using weights: %s", 
                         paste(sprintf("%s=%.3f", bad_outputs_full, bad_weights), collapse = ", ")))
        } else {
          warning(sprintf("  ⚠ Weight length mismatch for %s-frontier: expected %d, got %d", 
                         fid, length(bad_outputs_full), length(bad_weights)))
          bad_weights <- NULL
        }
      }
      
      # Calculate efficiency and slack (VRS) with constrained indices and weights
      message("  Calculating VRS efficiency with slack extraction...")
      result_vrs <- sbm_efficiency_with_slack(X_scaled, Y_scaled, bad = Z_scaled_full,
                                              constrained_bad_indices = frontiers[[fid]]$constrained_indices,
                                              bad_weights = bad_weights,
                                              RTS = "vrs")
      
      # Determine which bad outputs to report: ALL with weight > 0 (not just core constraints)
      # Initialize to constrained outputs as default
      weighted_bad_names <- constrained_bad_names_list[[fid]]
      weighted_bad_indices <- constrained_bad_indices_list[[fid]]
      
      if (!is.null(bad_weights) && length(bad_weights) == length(bad_outputs_full) && 
          !any(is.na(bad_weights))) {
        # Report all bad outputs with weight > 0 (exclude NA values)
        weighted_bad_indices_temp <- which(!is.na(bad_weights) & bad_weights > 0)
        if (length(weighted_bad_indices_temp) > 0) {
          weighted_bad_names <- bad_outputs_full[weighted_bad_indices_temp]
          weighted_bad_indices <- weighted_bad_indices_temp
          message(sprintf("  Reporting slack for %d bad outputs with weight > 0: %s",
                         length(weighted_bad_names), paste(weighted_bad_names, collapse = ", ")))
        } else {
          # No weights > 0, use constrained outputs (already set above)
          message(sprintf("  No weights > 0 - reporting slack for core constraints only: %s",
                         paste(weighted_bad_names, collapse = ", ")))
        }
      } else {
        # Fallback: use constrained bad output names (already set above)
        message(sprintf("  Weights not available - reporting slack for core constraints only: %s",
                       paste(weighted_bad_names, collapse = ", ")))
      }
    } else {
      # Original method: Use frontier-specific bad output set
      Z_bad <- as.matrix(df_clean[, frontiers[[fid]]$bad, drop = FALSE])
      Z_scaled <- scale_matrix(Z_bad)
      
      # Calculate efficiency and slack (VRS)
      result_vrs <- sbm_efficiency_with_slack(X_scaled, Y_scaled, bad = Z_scaled, RTS = "vrs")
      
      constrained_bad_names <- frontiers[[fid]]$bad
      # For original method, use constrained names as weighted names
      weighted_bad_names <- constrained_bad_names
      weighted_bad_indices <- seq_along(constrained_bad_names)
    }
    
    # Store results (include scale_bin for aggregation by scale group)
    slack_df <- data.frame(
      obs_id = seq_len(nrow(df_clean)),
      frontier = fid,
      efficiency = result_vrs$efficiency,
      scale_bin = df_clean$scale_bin,
      stringsAsFactors = FALSE
    )
    
    # Add input slacks
    for (i in seq_along(inputs)) {
      slack_df[[paste0("input_slack_", inputs[i])]] <- result_vrs$input_slack[, i]
    }
    
    # Add bad output slacks
    # Note: In fixed variable set method, bad_output_slack has columns for ALL bad outputs
    # Now we report slack for ALL bad outputs with weight > 0 (not just core constraints)
    if (ncol(result_vrs$bad_output_slack) > 0) {
      if (!is.null(Z_scaled_full) && !is.null(constrained_bad_indices_list) && 
          fid %in% names(constrained_bad_indices_list)) {
        # Fixed Variable Set Method: Report slack for ALL bad outputs with weight > 0
        # Use weighted_bad_names and weighted_bad_indices determined above
        if (length(weighted_bad_names) > 0 && length(weighted_bad_indices) > 0 && 
            length(weighted_bad_names) == length(weighted_bad_indices)) {
          for (i in seq_along(weighted_bad_names)) {
            bad_idx <- weighted_bad_indices[i]
            if (!is.na(bad_idx) && bad_idx >= 1 && bad_idx <= ncol(result_vrs$bad_output_slack)) {
              slack_df[[paste0("bad_slack_", weighted_bad_names[i])]] <- result_vrs$bad_output_slack[, bad_idx]
            }
          }
        } else {
          warning(sprintf("  ⚠ Cannot extract slack for %s-frontier: weighted_bad_names/indices mismatch", fid))
        }
      } else {
        # Original method: Report all bad output slacks (sequential mapping)
        for (i in seq_along(constrained_bad_names)) {
          if (i <= ncol(result_vrs$bad_output_slack)) {
            slack_df[[paste0("bad_slack_", constrained_bad_names[i])]] <- result_vrs$bad_output_slack[, i]
          }
        }
      }
    }
    
    # Add good output slack
    slack_df[["good_output_slack"]] <- result_vrs$good_output_slack[, 1]
    
    slack_results_list[[fid]] <- slack_df
    
    # Summary statistics
    message(sprintf("  ✓ Efficiency calculated with slack extraction"))
    message(sprintf("  ✓ Mean input slack per variable: %.4f", 
                    mean(rowSums(result_vrs$input_slack, na.rm = TRUE), na.rm = TRUE)))
    if (ncol(result_vrs$bad_output_slack) > 0) {
      message(sprintf("  ✓ Mean bad output slack per variable: %.4f",
                      mean(rowSums(result_vrs$bad_output_slack, na.rm = TRUE), na.rm = TRUE)))
    }
  }
  
  # Combine all slack results
  slack_all <- dplyr::bind_rows(slack_results_list)
  
  # Save detailed slack results
  write_csv(slack_all, 
            file.path(results_base, "efficiency_calculation/slack_analysis_detailed.csv"))
  message("\n✓ Detailed slack results saved")
  
  # Aggregate mean slack by frontier and scale_bin (for heatmap by scale group)
  if ("scale_bin" %in% names(slack_all)) {
    message("\nAggregating slack by frontier and scale group...")
    slack_cols <- c(
      grep("^input_slack_", names(slack_all), value = TRUE),
      grep("^bad_slack_", names(slack_all), value = TRUE)
    )
    slack_long <- slack_all %>%
      dplyr::select(obs_id, frontier, scale_bin, dplyr::all_of(slack_cols)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(slack_cols), names_to = "slack_var", values_to = "slack_value")
    slack_long <- slack_long %>%
      dplyr::mutate(
        variable_type = ifelse(grepl("^input_slack_", slack_var), "Input", "Bad Output"),
        variable_name = gsub("^input_slack_|^bad_slack_", "", slack_var)
      )
    slack_summary_by_scale <- slack_long %>%
      dplyr::group_by(frontier, scale_bin, variable_type, variable_name) %>%
      dplyr::summarise(mean_slack = mean(slack_value, na.rm = TRUE), .groups = "drop")
    slack_summary_by_scale$frontier <- as.character(slack_summary_by_scale$frontier)
    slack_summary_by_scale$frontier <- trimws(toupper(slack_summary_by_scale$frontier))
    slack_summary_by_scale$frontier[slack_summary_by_scale$frontier == "M"] <- "M"
    slack_summary_by_scale$frontier <- factor(slack_summary_by_scale$frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)
    slack_summary_by_scale$scale_bin <- factor(slack_summary_by_scale$scale_bin, levels = SCALE_LABELS, ordered = TRUE)
    write_csv(slack_summary_by_scale,
              file.path(results_base, "efficiency_calculation/slack_analysis_summary_by_scale.csv"))
    message("✓ Slack summary by scale group saved")
  } else {
    slack_summary_by_scale <- data.frame()
  }
  
  # Aggregate analysis: Identify main sources of inefficiency
  message("\nAggregating slack analysis by frontier...")
  
  # For each frontier, calculate mean slack by variable
  slack_summary_list <- list()
  
  for (fid in names(frontiers)) {
    slack_df <- slack_results_list[[fid]]
    
    # Calculate mean slack for each input
    input_slack_cols <- grep("^input_slack_", names(slack_df), value = TRUE)
    input_slack_means <- sapply(input_slack_cols, function(col) {
      mean(slack_df[[col]], na.rm = TRUE)
    })
    
    # Calculate mean slack for each bad output
    bad_slack_cols <- grep("^bad_slack_", names(slack_df), value = TRUE)
    bad_slack_means <- if (length(bad_slack_cols) > 0) {
      sapply(bad_slack_cols, function(col) {
        mean(slack_df[[col]], na.rm = TRUE)
      })
    } else {
      numeric(0)
    }
    
    # Combine into summary
    summary_df <- data.frame(
      frontier = fid,
      variable_type = c(rep("Input", length(input_slack_means)),
                       rep("Bad Output", length(bad_slack_means))),
      variable_name = c(gsub("^input_slack_", "", names(input_slack_means)),
                       gsub("^bad_slack_", "", names(bad_slack_means))),
      mean_slack = c(input_slack_means, bad_slack_means),
      stringsAsFactors = FALSE
    )
    
    slack_summary_list[[fid]] <- summary_df
  }
  
  slack_summary <- dplyr::bind_rows(slack_summary_list)
  
  # Save summary
  write_csv(slack_summary,
            file.path(results_base, "efficiency_calculation/slack_analysis_summary.csv"))
  message("✓ Slack summary saved")
  
  # Analysis: Compare slack between inefficient farms and meta-frontier
  message("\nAnalyzing sources of inefficiency relative to meta-frontier...")
  
  # Identify farms that are inefficient in each frontier compared to meta-frontier
  meta_inefficiency_analysis <- data.frame()
  
  for (fid in names(frontiers)) {
    frontier_eff_col <- paste0(fid, "_TE_VRS")
    meta_eff_col <- "Meta_TE_VRS"
    
    if (frontier_eff_col %in% names(efficiency_results) && 
        meta_eff_col %in% names(efficiency_results)) {
      
      # Find inefficient farms (frontier efficiency < meta-frontier efficiency)
      inefficient <- efficiency_results[[frontier_eff_col]] < efficiency_results[[meta_eff_col]]
      inefficient_ids <- which(inefficient)
      
      if (length(inefficient_ids) > 0) {
        slack_df_frontier <- slack_results_list[[fid]]
        
        # Get slack for inefficient farms
        inefficient_slack <- slack_df_frontier[inefficient_ids, ]
        
        # Calculate mean slack by variable for inefficient farms
        input_slack_cols <- grep("^input_slack_", names(inefficient_slack), value = TRUE)
        bad_slack_cols <- grep("^bad_slack_", names(inefficient_slack), value = TRUE)
        
        for (col in input_slack_cols) {
          var_name <- gsub("^input_slack_", "", col)
          meta_inefficiency_analysis <- rbind(meta_inefficiency_analysis, data.frame(
            frontier = fid,
            variable_type = "Input",
            variable_name = var_name,
            mean_slack = mean(inefficient_slack[[col]], na.rm = TRUE),
            n_inefficient = length(inefficient_ids),
            efficiency_gap = mean(efficiency_results[[meta_eff_col]][inefficient_ids] - 
                                  efficiency_results[[frontier_eff_col]][inefficient_ids], 
                                  na.rm = TRUE),
            stringsAsFactors = FALSE
          ))
        }
        
        for (col in bad_slack_cols) {
          var_name <- gsub("^bad_slack_", "", col)
          meta_inefficiency_analysis <- rbind(meta_inefficiency_analysis, data.frame(
            frontier = fid,
            variable_type = "Bad Output",
            variable_name = var_name,
            mean_slack = mean(inefficient_slack[[col]], na.rm = TRUE),
            n_inefficient = length(inefficient_ids),
            efficiency_gap = mean(efficiency_results[[meta_eff_col]][inefficient_ids] - 
                                  efficiency_results[[frontier_eff_col]][inefficient_ids], 
                                  na.rm = TRUE),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Save meta-frontier inefficiency analysis
  if (nrow(meta_inefficiency_analysis) > 0) {
    write_csv(meta_inefficiency_analysis,
              file.path(results_base, "efficiency_calculation/meta_frontier_inefficiency_sources.csv"))
    message("✓ Meta-frontier inefficiency analysis saved")
    
    # Print key findings
    message("\nKey Findings:")
    for (fid in names(frontiers)) {
      fid_data <- meta_inefficiency_analysis[meta_inefficiency_analysis$frontier == fid, ]
      if (nrow(fid_data) > 0) {
        # Find top 3 variables with highest slack
        top_vars <- fid_data[order(-fid_data$mean_slack), ][seq_len(min(3, nrow(fid_data))), ]
        message(sprintf("\n%s-Frontier inefficiency sources (top contributors):", fid))
        for (i in seq_len(nrow(top_vars))) {
          message(sprintf("  %d. %s (%s): mean slack = %.4f, efficiency gap = %.3f",
                         i, top_vars$variable_name[i], top_vars$variable_type[i],
                         top_vars$mean_slack[i], top_vars$efficiency_gap[1]))
        }
      }
    }
  }
  
  message("\n✓ Slack analysis completed\n")
  
  return(list(
    detailed = slack_all,
    summary = slack_summary,
    summary_by_scale = slack_summary_by_scale,
    meta_inefficiency = meta_inefficiency_analysis
  ))
}

# ============================================================================
# Visualization Function for Slack Analysis
# ============================================================================

# Publication theme for slack visualization
theme_slack_publication <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "gray30"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      legend.position = "right",  # Unified legend position: right side
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Create comprehensive slack visualization
visualize_slack_analysis <- function(slack_results, results_base, inputs, bad_outputs_list) {
  
  message("Creating slack analysis visualizations...")
  
  slack_summary <- slack_results$summary
  slack_summary_by_scale <- slack_results$summary_by_scale
  meta_inefficiency <- slack_results$meta_inefficiency
  
  # Standardize frontier labels and ensure order: M-L-G-A
  slack_summary$frontier <- as.character(slack_summary$frontier)
  slack_summary$frontier <- trimws(toupper(slack_summary$frontier))
  slack_summary$frontier[slack_summary$frontier == "M"] <- "M"
  slack_summary <- slack_summary %>%
    dplyr::filter(frontier %in% c("M", "L", "G", "A"))
  slack_summary$frontier <- factor(slack_summary$frontier,
                                   levels = c("M", "L", "G", "A"),
                                   ordered = TRUE)
  
  if (nrow(meta_inefficiency) > 0) {
    meta_inefficiency$frontier <- as.character(meta_inefficiency$frontier)
    meta_inefficiency$frontier <- trimws(toupper(meta_inefficiency$frontier))
    meta_inefficiency$frontier[meta_inefficiency$frontier == "M"] <- "M"
    meta_inefficiency <- meta_inefficiency %>%
      dplyr::filter(frontier %in% c("M", "L", "G", "A"))
    meta_inefficiency$frontier <- factor(meta_inefficiency$frontier,
                                        levels = c("M", "L", "G", "A"),
                                        ordered = TRUE)
  }
  
  # ============================================================================
  # Panel 1: Top Contributing Variables to Inefficiency (Heatmap)
  # ============================================================================
  
  # Calculate relative contribution (normalize by maximum)
  slack_panel2 <- slack_summary %>%
    dplyr::group_by(variable_type) %>%
    dplyr::mutate(
      variable_name = ifelse(variable_name == "Waste_water", "Wastewater", variable_name),
      max_slack_in_type = max(mean_slack, na.rm = TRUE),
      relative_contribution = ifelse(max_slack_in_type > 0, 
                                     mean_slack / max_slack_in_type, 
                                     0)
    ) %>%
    dplyr::ungroup()
  
  # Calculate text color based on fill value (dark background = white text, light = black)
  # Use a threshold based on the color scale midpoint
  slack_panel2_with_colors <- slack_panel2 %>%
    dplyr::mutate(
      # Calculate relative position in color scale (0-1)
      # For viridis_c option="C" reversed, darker colors = higher values
      color_scale_pos = (mean_slack - min(mean_slack, na.rm = TRUE)) / 
                        (max(mean_slack, na.rm = TRUE) - min(mean_slack, na.rm = TRUE)),
      # Use threshold: if in darker half (top 50%), use white text; otherwise black
      text_color = ifelse(color_scale_pos > 0.5, "white", "black")
    )
  
  p_slack_heatmap <- slack_panel2_with_colors %>%
    ggplot(aes(x = frontier, y = variable_name, fill = mean_slack)) +
    geom_tile(color = "white", linewidth = 0.5, alpha = 0.9) +
    # Use dynamic text color based on background brightness
    # Map text_color in aes and use scale_color_identity to apply colors directly
    geom_text(aes(label = sprintf("%.3f", mean_slack), color = text_color),
              size = 2.8, fontface = "bold") +
    scale_color_identity(guide = "none") +  # Use the color values directly from data, no legend
    scale_fill_viridis_c(option = "C", direction = -1,
                        name = "Mean slack",
                        guide = guide_colorbar(title.position = "top", barwidth = 0.8, barheight = 8)) +
    facet_wrap(~ variable_type, scales = "free_y", ncol = 1) +
    labs(
      title = "Slack values by frontier and variable (heatmap)",
      x = "Frontier",
      y = "Variable",
      fill = "Mean slack"
    ) +
    theme_slack_publication() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.background = element_rect(fill = "gray90", color = "gray70"),
        strip.text = element_text(face = "bold"),
        legend.position = "right"  # Ensure legend on right side
      )
  
  # ============================================================================
  # Panel 2: Sources of Inefficiency Relative to Meta-Frontier
  # ============================================================================
  
  if (nrow(meta_inefficiency) > 0) {
    # Get top 3 variables per frontier
    top_vars_by_frontier <- meta_inefficiency %>%
      dplyr::group_by(frontier) %>%
      dplyr::arrange(-mean_slack) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        variable_name = ifelse(variable_name == "Waste_water", "Wastewater", variable_name),
        rank = rep(1:3, times = length(unique(frontier))),
        variable_label = paste0(variable_name, " (", variable_type, ")")
      )
    
    p_meta_inefficiency <- top_vars_by_frontier %>%
      ggplot(aes(x = reorder(variable_label, mean_slack), 
                 y = mean_slack, 
                 fill = frontier)) +
      geom_col(alpha = 0.85, width = 0.75) +
      geom_text(aes(label = sprintf("%.3f", mean_slack)),
                hjust = -0.1, size = 2.8, fontface = "bold") +
      coord_flip() +
      scale_fill_viridis_d(option = "D", breaks = c("M", "L", "G", "A")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(
        title = "Top 3 sources of inefficiency vs meta-frontier",
        subtitle = "Variables with highest slack for inefficient farms (frontier efficiency < meta-frontier efficiency)",
        x = "Variable (top 3 per frontier)",
        y = "Mean slack",
        fill = "Frontier"
      ) +
      facet_wrap(~ frontier, scales = "free_y", ncol = 2) +
      theme_slack_publication() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.background = element_rect(fill = "gray90", color = "gray70"),
        strip.text = element_text(face = "bold"),
        legend.position = "right"  # Ensure legend on right side
      )
  } else {
    p_meta_inefficiency <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No meta-frontier inefficiency data available",
               size = 5) +
      theme_void()
  }
  
  # ============================================================================
  # Panel 3: New Heatmap - Variables (x-axis) × Frontiers (y-axis)
  # ============================================================================
  
  # Prepare data for the new heatmap (transposed: variables on x, frontiers on y)
  slack_heatmap_data <- slack_summary %>%
    dplyr::mutate(
      variable_name = ifelse(variable_name == "Waste_water", "Wastewater", variable_name),
      variable_name = factor(variable_name, 
                            levels = unique(variable_name[order(variable_type, variable_name)]),
                            ordered = TRUE),
      frontier = as.character(frontier),
      frontier = trimws(toupper(frontier)),
      frontier = ifelse(frontier == "M", "M", frontier),
      frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)
    ) %>%
    dplyr::filter(!is.na(frontier)) %>%
    dplyr::arrange(variable_type, variable_name, frontier)
  
  # Use fixed black text for all data labels
  slack_heatmap_data <- slack_heatmap_data %>%
    dplyr::mutate(
      color_scale_pos = (mean_slack - min(mean_slack, na.rm = TRUE)) / 
                        (max(mean_slack, na.rm = TRUE) - min(mean_slack, na.rm = TRUE)),
      text_color = "black"
    )
  
  # Create the new heatmap: Variables (x) × Frontiers (y)
  # Reorder variable_type to have Input on left, Output on right
  slack_heatmap_data$variable_type <- factor(slack_heatmap_data$variable_type, 
                                           levels = c("Input", "Bad Output"), 
                                           ordered = TRUE)
  
  # Sequential palette aligned with Fig1 frontier style (warm-neutral to L-blue)
  heatmap_colors <- c("#FFF4E8", "#F6E7DE", "#E7EEF3", "#C9DCEA", "#8FB8CF")
  p_slack_heatmap_transposed <- slack_heatmap_data %>%
    ggplot(aes(x = variable_name, y = frontier, fill = mean_slack)) +
    geom_tile(color = "white", linewidth = 0.8, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.3f", mean_slack), color = text_color),
              size = 11.5, fontface = "bold") +
    scale_color_identity(guide = "none") +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
    scale_fill_gradientn(colors = heatmap_colors,
                       name = "Mean slack",
                       guide = guide_colorbar(title.position = "top",
                                             barwidth = 1.2,
                                             barheight = 15)) +
    facet_wrap(~ variable_type, scales = "free_x", ncol = 2) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = "Frontier",
      fill = "Mean slack"
    ) +
    theme_minimal(base_size = 30) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.9, size = 32, color = "black"),
      axis.text.y = element_text(size = 32, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 34),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
      strip.background = element_rect(fill = "gray90", color = "gray70", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = 32),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 32),
      legend.text = element_text(size = 28),
      plot.margin = margin(20, 20, 20, 20, "pt")
    )
  
  # ============================================================================
  # Panel 4: Heatmap by scale group (frontier A only; variable × scale_bin, facet by variable_type)
  # ============================================================================
  
  if (!is.null(slack_summary_by_scale) && nrow(slack_summary_by_scale) > 0) {
    heatmap_scale_colors <- c("#FFF4E8", "#F6E7DE", "#E7EEF3", "#C9DCEA", "#8FB8CF")
    # Fig4d: only frontier A slacks, grouped by scale
    # Keep full scale_bin levels so Y-axis shows all categories even when only one has data
    slack_scale_heat <- slack_summary_by_scale %>%
      dplyr::filter(frontier == "A") %>%
      dplyr::mutate(
        scale_bin = factor(scale_bin, levels = SCALE_LABELS, ordered = TRUE),
        color_scale_pos = (mean_slack - min(mean_slack, na.rm = TRUE)) /
                          (max(mean_slack, na.rm = TRUE) - min(mean_slack, na.rm = TRUE) + 1e-12),
        text_color = "black"
      )
    if (nrow(slack_scale_heat) > 0) {
      p_slack_heatmap_by_scale <- slack_scale_heat %>%
        ggplot(aes(x = variable_name, y = scale_bin, fill = mean_slack)) +
        geom_tile(color = "white", linewidth = 0.5, alpha = 0.9) +
        geom_text(aes(label = sprintf("%.3f", mean_slack), color = text_color),
                  size = 2.5, fontface = "bold") +
        scale_color_identity(guide = "none") +
        scale_fill_gradientn(colors = heatmap_scale_colors,
                           name = "Mean slack",
                           guide = guide_colorbar(title.position = "top", barwidth = 0.8, barheight = 8)) +
        scale_y_discrete(drop = FALSE) +
        facet_wrap(~ variable_type, scales = "free_x", ncol = 2) +
        labs(
          title = "Slack values by scale group (frontier A only)",
          subtitle = "Mean slack by scale category (market pig heads/a), frontier A",
          x = "Variable",
          y = "Scale category",
          fill = "Mean slack"
        ) +
        theme_slack_publication() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "gray90", color = "gray70"),
          strip.text = element_text(face = "bold"),
          legend.position = "right"
        )
    } else {
      p_slack_heatmap_by_scale <- NULL
    }
  } else {
    p_slack_heatmap_by_scale <- NULL
  }
  
  # ============================================================================
  # Save figures separately (each with legend on right side)
  # ============================================================================
  
  # Save Heatmap figure separately
  ggsave(file.path(results_base, "visualization/Fig4a_slack_heatmap.png"),
         p_slack_heatmap, width = 10, height = 8, dpi = 300)
  
  message("✓ Figure 4a saved: Fig4a_slack_heatmap.png")
  message("  Layout: Heatmap with legend on right side\n")
  
  # Save new transposed heatmap (Variables × Frontiers, 32×10 inches)
  ggsave(file.path(results_base, "visualization/Fig4c_slack_heatmap_variables_frontiers.png"),
         p_slack_heatmap_transposed, width = 32, height = 10, dpi = 300)
  
  message("✓ Figure 4c saved: Fig4c_slack_heatmap_variables_frontiers.png")
  message("  Layout: Variables (x-axis) × Frontiers (y-axis), 32×10 inches\n")
  
  # Save slack heatmap by scale group
  if (!is.null(p_slack_heatmap_by_scale)) {
    ggsave(file.path(results_base, "visualization/Fig4d_slack_heatmap_by_scale.png"),
           p_slack_heatmap_by_scale, width = 14, height = 8, dpi = 300)
    message("✓ Figure 4d saved: Fig4d_slack_heatmap_by_scale.png")
    message("  Layout: Slack by scale category and variable (frontier A only, facet: variable type)\n")
  }
  
  # Save Meta-frontier analysis figure separately
  if (nrow(meta_inefficiency) > 0) {
    ggsave(file.path(results_base, "visualization/Fig4b_meta_frontier_inefficiency.png"),
           p_meta_inefficiency, width = 10, height = 8, dpi = 300)
    
    message("✓ Figure 4b saved: Fig4b_meta_frontier_inefficiency.png")
    message("  Layout: Top 3 sources of inefficiency with legend on right side\n")
    
    return(list(heatmap = p_slack_heatmap,
                heatmap_transposed = p_slack_heatmap_transposed,
                heatmap_by_scale = p_slack_heatmap_by_scale,
                meta_inefficiency = p_meta_inefficiency))
  } else {
    message("⚠ Meta-frontier inefficiency data not available - only heatmap saved\n")
    return(list(heatmap = p_slack_heatmap,
                heatmap_transposed = p_slack_heatmap_transposed,
                heatmap_by_scale = p_slack_heatmap_by_scale,
                meta_inefficiency = NULL))
  }
}










