# ============================================================================
# PRE-EXPERIMENT MODULE: Optimal Weight Multiplier Search for Robust Frontier Differentiation
# ============================================================================
# 
# This module systematically tests different analysis/non-core weight multipliers
# to identify the optimal multiplier using correlation plateau detection:
# 1. Detects correlation plateau: region where correlation changes stabilize
# 2. Selects the multiplier at the start of the largest plateau region
# 3. This ensures robust weight configuration without arbitrary correlation thresholds
# 
# Method: Correlation Plateau Detection
# - Calculates correlation change between adjacent multipliers
# - Identifies plateau regions (correlation change < 0.001 threshold)
# - Selects the first multiplier in the largest plateau region
# - No fixed correlation threshold required (e.g., 0.95)
#
# Usage:
#   run_weight_ratio_pre_experiment(
#     df_clean = df_clean,
#     X_scaled = X_scaled,
#     Y_scaled = Y_scaled,
#     Z_scaled_full = Z_scaled_full,
#     constrained_bad_M = constrained_bad_M,
#     constrained_bad_L = constrained_bad_L,
#     constrained_bad_G = constrained_bad_G,
#     constrained_bad_A = constrained_bad_A,
#     sbm_efficiency = sbm_efficiency,
#     normalize_weights = normalize_weights,
#     results_base = results_base,
#     bad_outputs_full = bad_outputs_full,
#     skip_pre_experiment = FALSE  # Set to TRUE to skip and use default weights
#   )
#
# Returns: List with optimal multiplier and weight configurations
# ============================================================================

run_weight_ratio_pre_experiment <- function(
  df_clean,
  X_scaled,
  Y_scaled,
  Z_scaled_full,
  constrained_bad_M,
  constrained_bad_L,
  constrained_bad_G,
  constrained_bad_A,
  sbm_efficiency,
  normalize_weights,
  results_base,
  bad_outputs_full,
  skip_pre_experiment = FALSE
) {
  
  # Define theme_publication for plots (if not already defined)
  if (!exists("theme_publication", envir = .GlobalEnv)) {
    theme_publication <- function(base_size = 16) {
      ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5, color = "gray30"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "gray70", fill = NA, linewidth = 0.5),
          legend.position = "bottom",
          legend.title = ggplot2::element_text(face = "bold")
        )
    }
  } else {
    theme_publication <- get("theme_publication", envir = .GlobalEnv)
  }
  
  # Ensure required packages are loaded
  if (!require("ggplot2", quietly = TRUE)) {
    stop("ERROR: ggplot2 package is required but not loaded. Please install and load ggplot2.")
  }
  if (!require("dplyr", quietly = TRUE)) {
    stop("ERROR: dplyr package is required but not loaded. Please install and load dplyr.")
  }
  if (!require("tidyr", quietly = TRUE)) {
    stop("ERROR: tidyr package is required but not loaded. Please install and load tidyr.")
  }
  if (!require("patchwork", quietly = TRUE)) {
    stop("ERROR: patchwork package is required but not loaded. Please install and load patchwork.")
  }
  if (!require("writexl", quietly = TRUE)) {
    install.packages("writexl", repos = "https://cloud.r-project.org")
    library(writexl)
  }

  # Pre-experiment is required - no skipping allowed
  if (skip_pre_experiment) {
    stop("ERROR: skip_pre_experiment=TRUE is not allowed. Pre-experiment is required to determine optimal weights using correlation plateau method.")
  } else {
    message(paste(rep("=", 78), collapse = ""))
    message("PRE-EXPERIMENT: OPTIMAL WEIGHT MULTIPLIER SEARCH")
    message(paste(rep("=", 78), collapse = ""))

# Function to generate weights based on multiplier
# core_multiplier: weight multiplier for core constraints (relative to non-core)
# Returns normalized weight vectors for all frontiers
# REVISED: Use asymmetric base weights to create larger effective differences after normalization
generate_weights_by_multiplier <- function(core_multiplier) {
  # PROBLEM IDENTIFIED: Standardization to sum=5 compresses weight differences
  # SOLUTION: Use asymmetric base weights (non-core < 1.0) to amplify differences after normalization
  # 
  # Strategy: non-core = 0.5, core = core_multiplier
  # This creates larger effective ratios after normalization
  # Example: multiplier=10, non-core=0.5, core=10
  #   E-frontier: [0.5,0.5,0.5,10,0.5] → sum=12 → normalized: [0.208, 0.208, 0.208, 4.167, 0.208]
  #   Effective ratio: 4.167/0.208 = 20x (vs original 10x before normalization)
  
  non_core_base <- 0.5  # Reduced from 1.0 to amplify differences after normalization
  
  # Fixed bad output set follows bad_outputs_full from the main script (Water_footprint removed).
  target_sum <- length(bad_outputs_full)
  
  # E-frontier: Dead_pig is core
  weights_E_raw <- c(
    Waste_water = non_core_base,
    Manure = non_core_base,
    Dead_pig = core_multiplier,  # Core constraint
    Carbon_emission = non_core_base,
    Eutrophication_potential = non_core_base
  )
  
  # L-frontier: Local environmental pollutants are core
  weights_L_raw <- c(
    Waste_water = core_multiplier,      # Core
    Manure = core_multiplier,           # Core
    Dead_pig = non_core_base,
    Carbon_emission = non_core_base,
    Eutrophication_potential = non_core_base
  )
  
  # G-frontier: Carbon_emission + Eutrophication_potential are core
  weights_G_raw <- c(
    Waste_water = non_core_base,
    Manure = non_core_base,
    Dead_pig = non_core_base,
    Carbon_emission = core_multiplier,   # Core
    Eutrophication_potential = core_multiplier   # Core
  )
  
  # A-frontier: All are core (equal weight)
  weights_A_raw <- rep(core_multiplier, target_sum)
  names(weights_A_raw) <- bad_outputs_full
  
  # Normalize weights to sum = target_sum
  normalize <- function(w) {
    w * (target_sum / sum(w))
  }
  
  normalized_weights <- list(
    E = normalize(weights_E_raw),
    L = normalize(weights_L_raw),
    G = normalize(weights_G_raw),
    A = normalize(weights_A_raw)
  )
  
  # Calculate effective ratio for debugging (optional)
  # E-frontier effective ratio: core_weight / non_core_weight
  if (core_multiplier > 0 && non_core_base > 0) {
    eff_ratio_M <- normalized_weights$E["Dead_pig"] / normalized_weights$E["Waste_water"]
    eff_ratio_L <- normalized_weights$L["Waste_water"] / normalized_weights$L["Dead_pig"]
    eff_ratio_G <- normalized_weights$G["Carbon_emission"] / normalized_weights$G["Waste_water"]
    
    # Only print for first few multipliers to avoid spam
    if (core_multiplier <= 3.0) {
      message(sprintf("  Multiplier %.2fx: Effective ratios - M:%.2fx, L:%.2fx, G:%.2fx", 
                      core_multiplier, eff_ratio_M, eff_ratio_L, eff_ratio_G))
    }
  }
  
  return(normalized_weights)
}

    # Function to calculate efficiency correlations for given multiplier
    test_multiplier <- function(core_multiplier, n_subsample = NULL) {
      message(sprintf("Testing multiplier: %.2fx", core_multiplier))
      
      # Generate weights
      weight_list <- generate_weights_by_multiplier(core_multiplier)
      
      # Subsample if requested (for faster testing)
      if (!is.null(n_subsample) && n_subsample < nrow(df_clean)) {
        set.seed(123)
        sample_indices <- sample(seq_len(nrow(df_clean)), n_subsample)
        X_test <- X_scaled[sample_indices, ]
        Y_test <- Y_scaled[sample_indices, , drop = FALSE]
        Z_test <- Z_scaled_full[sample_indices, , drop = FALSE]
      } else {
        X_test <- X_scaled
        Y_test <- Y_scaled
        Z_test <- Z_scaled_full
      }
  
  # Calculate efficiencies for all frontiers
  eff_results <- data.frame(
    obs_id = seq_len(nrow(X_test))
  )
  
  # Define frontier configurations
  frontiers_config <- list(
    M = list(constrained = constrained_bad_M, weights = weight_list$E),
    L = list(constrained = constrained_bad_L, weights = weight_list$L),
    G = list(constrained = constrained_bad_G, weights = weight_list$G),
    A = list(constrained = constrained_bad_A, weights = weight_list$A)
  )
  
  for (fid in names(frontiers_config)) {
    config <- frontiers_config[[fid]]
    
    # Calculate VRS efficiency
    te_vrs <- sbm_efficiency(
      X = X_test,
      Y = Y_test,
      bad = Z_test,
      constrained_bad_indices = config$constrained,
      bad_weights = config$weights,
      RTS = "vrs"
    )
    
    eff_results[[paste0(fid, "_TE_VRS")]] <- te_vrs
  }
  
  # Calculate correlations
  cor_matrix <- cor(eff_results[, c("M_TE_VRS", "L_TE_VRS", "G_TE_VRS", "A_TE_VRS")], 
                    use = "pairwise.complete.obs")
  
  # Extract key correlations
  cor_ML <- cor_matrix["M_TE_VRS", "L_TE_VRS"]
  cor_MG <- cor_matrix["M_TE_VRS", "G_TE_VRS"]
  cor_MA <- cor_matrix["M_TE_VRS", "A_TE_VRS"]
  cor_LG <- cor_matrix["L_TE_VRS", "G_TE_VRS"]
  cor_LA <- cor_matrix["L_TE_VRS", "A_TE_VRS"]
  cor_GA <- cor_matrix["G_TE_VRS", "A_TE_VRS"]
  
  # Calculate efficiency statistics
  eff_stats <- data.frame(
    multiplier = core_multiplier,
    # Correlation statistics
    cor_min = min(cor_matrix[upper.tri(cor_matrix)]),
    cor_max = max(cor_matrix[upper.tri(cor_matrix)]),
    cor_mean = mean(cor_matrix[upper.tri(cor_matrix)]),
    cor_median = median(cor_matrix[upper.tri(cor_matrix)]),
    # Key pairwise correlations
    cor_ML = cor_ML,
    cor_MG = cor_MG,
    cor_MA = cor_MA,
    cor_LG = cor_LG,
    cor_LA = cor_LA,
    cor_GA = cor_GA,
    # Efficiency distribution statistics
    M_mean = mean(eff_results$M_TE_VRS, na.rm = TRUE),
    M_sd = sd(eff_results$M_TE_VRS, na.rm = TRUE),
    M_min = min(eff_results$M_TE_VRS, na.rm = TRUE),
    M_max = max(eff_results$M_TE_VRS, na.rm = TRUE),
    L_mean = mean(eff_results$L_TE_VRS, na.rm = TRUE),
    L_sd = sd(eff_results$L_TE_VRS, na.rm = TRUE),
    L_min = min(eff_results$L_TE_VRS, na.rm = TRUE),
    L_max = max(eff_results$L_TE_VRS, na.rm = TRUE),
    G_mean = mean(eff_results$G_TE_VRS, na.rm = TRUE),
    G_sd = sd(eff_results$G_TE_VRS, na.rm = TRUE),
    G_min = min(eff_results$G_TE_VRS, na.rm = TRUE),
    G_max = max(eff_results$G_TE_VRS, na.rm = TRUE),
    A_mean = mean(eff_results$A_TE_VRS, na.rm = TRUE),
    A_sd = sd(eff_results$A_TE_VRS, na.rm = TRUE),
    A_min = min(eff_results$A_TE_VRS, na.rm = TRUE),
    A_max = max(eff_results$A_TE_VRS, na.rm = TRUE),
    # Differentiation metrics
    mean_diff_ML = mean(abs(eff_results$M_TE_VRS - eff_results$L_TE_VRS), na.rm = TRUE),
    mean_diff_MG = mean(abs(eff_results$M_TE_VRS - eff_results$G_TE_VRS), na.rm = TRUE),
    mean_diff_MA = mean(abs(eff_results$M_TE_VRS - eff_results$A_TE_VRS), na.rm = TRUE),
    max_diff = max(c(
      max(abs(eff_results$M_TE_VRS - eff_results$L_TE_VRS), na.rm = TRUE),
      max(abs(eff_results$M_TE_VRS - eff_results$G_TE_VRS), na.rm = TRUE),
      max(abs(eff_results$M_TE_VRS - eff_results$A_TE_VRS), na.rm = TRUE)
    )),
    n_valid = sum(!is.na(eff_results$M_TE_VRS))
  )
  
  return(eff_stats)
}

    # Run pre-experiment with different multipliers
    message("Running systematic multiplier test...")
    message("REVISED STRATEGY: Using asymmetric base weights (non-core=0.5, core=multiplier)")
    message("  This amplifies effective weight differences after normalization")
    message("Testing multipliers from 1.0x to 50.0x (in 1.0x increments)")
    message("  Note: Step size = 1.0 for faster testing. If plateau not found, reduce step size.")

    # Step size = 1.0: 1.0x to 50.0x with 1.0x increments
    # This provides full range coverage with reasonable resolution
    # If plateau is not found, can reduce step size (e.g., 0.5) for finer resolution
    multipliers <- seq(1.0, 50.0, by = 1.0)
    results_list <- list()

    # Use subsample for faster testing (optional)
    # Set to NULL to use full dataset (more accurate but slower)
    n_subsample <- min(200, nrow(df_clean))  # Use 200 or all if less

for (i in seq_along(multipliers)) {
  multiplier <- multipliers[i]
  message(sprintf("  [%d/%d] Testing %.2fx multiplier...", 
                  i, length(multipliers), multiplier))
  
  # Add timeout mechanism (30 seconds per multiplier)
  start_time_iter <- Sys.time()
  result <- tryCatch({
    test_multiplier(multiplier, n_subsample = n_subsample)
  }, error = function(e) {
    message(sprintf("    ERROR: %s", e$message))
    return(NULL)
  })
  
  # Check if iteration took too long
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time_iter, units = "secs"))
  if (elapsed_time > 30) {
    message(sprintf("    WARNING: Multiplier %.2fx took %.1f seconds (may indicate numerical issues)", 
                    multiplier, elapsed_time))
  }
  
  if (!is.null(result)) {
    results_list[[i]] <- result
  }
  
  # Progress update (more frequent for better feedback)
  if (i %% 3 == 0 || i == length(multipliers)) {
    message(sprintf("    Progress: %d/%d complete (%.1f%%)", 
                    i, length(multipliers), 100*i/length(multipliers)))
  }
}

# Combine results
pre_exp_results <- do.call(rbind, results_list)

    # Save pre-experiment results
    readr::write_csv(pre_exp_results, 
                    file.path(results_base, "diagnostics/pre_experiment_multiplier_search.csv"))
    message(sprintf("✓ Pre-experiment results saved for %d multipliers", nrow(pre_exp_results)))

    # Plateau detection (before plots so we can annotate plateau position)
    # Stricter threshold: correlation change < 0.0005 (0.05%) to avoid overly long plateaus (e.g. 7x-50x)
    plateau_threshold <- 0.0001
    # Prefer the first (earliest) adequate plateau, not the longest one that runs to max multiplier
    min_plateau_points <- 3L
    pre_exp_results <- pre_exp_results %>% dplyr::arrange(multiplier)
    pre_exp_results <- pre_exp_results %>%
      dplyr::mutate(
        cor_change = abs(cor_mean - dplyr::lag(cor_mean)),
        cor_change_pct = cor_change / pmax(dplyr::lag(cor_mean), 0.001) * 100,
        in_plateau = !is.na(cor_change) & cor_change < plateau_threshold,
        plateau_group = cumsum(!in_plateau | is.na(cor_change))
      )
    plateau_summary <- pre_exp_results %>%
      dplyr::filter(in_plateau) %>%
      dplyr::group_by(plateau_group) %>%
      dplyr::summarise(
        start_multiplier = min(multiplier),
        end_multiplier = max(multiplier),
        n_points = n(),
        mean_cor = mean(cor_mean),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_points >= min_plateau_points) %>%
      dplyr::arrange(start_multiplier)
    best_plateau <- if (nrow(plateau_summary) > 0) plateau_summary[1, ] else NULL

    # ============================================================================
    # VISUALIZATION OF PRE-EXPERIMENT RESULTS
    # ============================================================================

    message("\nGenerating visualization of pre-experiment results...")

    # Y-axis limits from data (minimal padding to avoid excess blank space)
    y_cor_range <- range(pre_exp_results$cor_min, pre_exp_results$cor_max, na.rm = TRUE)
    y_cor_lim <- c(max(0, y_cor_range[1] - 0.02), min(1, y_cor_range[2] + 0.02))

    # Plot 1: Correlation vs Multiplier (no reference lines; annotate plateau)
    p_cor_vs_mult <- ggplot2::ggplot(pre_exp_results, ggplot2::aes(x = multiplier)) +
      ggplot2::geom_line(ggplot2::aes(y = cor_mean, color = "Mean Correlation"), linewidth = 1.2) +
      ggplot2::geom_point(ggplot2::aes(y = cor_mean), size = 2) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = cor_min, ymax = cor_max), 
                  alpha = 0.2, fill = "steelblue")
    if (!is.null(best_plateau)) {
      p_cor_vs_mult <- p_cor_vs_mult +
        ggplot2::geom_vline(xintercept = best_plateau$start_multiplier, 
                   linetype = "dashed", color = "darkgreen", linewidth = 0.9) +
        ggplot2::annotate("rect", 
                 xmin = best_plateau$start_multiplier, xmax = best_plateau$end_multiplier,
                 ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "darkgreen") +
        ggplot2::annotate("text", 
                 x = best_plateau$start_multiplier, y = y_cor_lim[2] - 0.02,
                 label = sprintf("Plateau\n(%.2fx\u2013%.2fx)", best_plateau$start_multiplier, best_plateau$end_multiplier),
                 color = "darkgreen", hjust = -0.05, size = 3.2, fontface = "bold")
    }
    p_cor_vs_mult <- p_cor_vs_mult +
      ggplot2::scale_x_continuous(breaks = seq(1, 10, 1)) +
      ggplot2::scale_y_continuous(limits = y_cor_lim, expand = ggplot2::expansion(mult = 0.03)) +
      ggplot2::scale_color_manual(values = c("Mean Correlation" = "steelblue")) +
      ggplot2::labs(
        title = "Frontier differentiation vs weight multiplier",
        subtitle = "Correlation plateau detection: stable correlation region identified",
        x = "Core/non-core weight multiplier",
        y = "Correlation coefficient",
        color = "Metric"
      ) +
      theme_publication() +
      ggplot2::theme(legend.position = "bottom")

    # Plot 2: Efficiency Range vs Multiplier
    eff_range_data <- pre_exp_results %>%
      dplyr::select(multiplier, 
             M_min, M_max, L_min, L_max, 
             G_min, G_max, A_min, A_max) %>%
      tidyr::pivot_longer(cols = -multiplier, 
                 names_to = c("frontier", "stat"), 
                 names_sep = "_") %>%
      tidyr::pivot_wider(names_from = stat, values_from = value) %>%
      dplyr::mutate(frontier = dplyr::if_else(frontier == "M", "M", frontier)) %>%
      dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A")))

    y_eff_range <- range(c(eff_range_data$min, eff_range_data$max), na.rm = TRUE)
    y_eff_lim <- c(max(0, y_eff_range[1] - 0.03), min(1, y_eff_range[2] + 0.03))

    p_eff_range <- ggplot2::ggplot(eff_range_data, ggplot2::aes(x = multiplier)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = min, ymax = max, fill = frontier), 
              alpha = 0.3) +
      ggplot2::geom_line(ggplot2::aes(y = (min + max)/2, color = frontier), 
            linewidth = 1, alpha = 0.8)
    if (!is.null(best_plateau)) {
      p_eff_range <- p_eff_range +
        ggplot2::geom_vline(xintercept = best_plateau$start_multiplier, 
                   linetype = "dashed", color = "darkgreen", linewidth = 0.9) +
        ggplot2::annotate("text", x = best_plateau$start_multiplier, y = y_eff_lim[2] - 0.02,
                 label = sprintf("Plateau (%.2fx)", best_plateau$start_multiplier),
                 color = "darkgreen", hjust = -0.05, size = 3)
    }
    p_eff_range <- p_eff_range +
      ggplot2::scale_x_continuous(breaks = seq(1, 10, 1)) +
      ggplot2::scale_y_continuous(limits = y_eff_lim, expand = ggplot2::expansion(mult = 0.03)) +
      ggplot2::scale_fill_viridis_d(option = "D", 
                       labels = c("M" = "M-Frontier", 
                                  "L" = "L-Frontier", 
                                  "G" = "G-Frontier", 
                                  "A" = "A-Frontier")) +
      ggplot2::scale_color_viridis_d(option = "D",
                        labels = c("M" = "M-Frontier", 
                                   "L" = "L-Frontier", 
                                   "G" = "G-Frontier", 
                                   "A" = "A-Frontier")) +
      ggplot2::labs(
        title = "Efficiency value range vs weight multiplier",
        subtitle = "Range = min to max efficiency values",
        x = "Core/non-core weight multiplier",
        y = "Efficiency value",
        fill = "Frontier",
        color = "Frontier"
      ) +
      theme_publication() +
      ggplot2::theme(legend.position = "bottom")

    # Plot 3: Mean Efficiency Differences vs Multiplier
    diff_data <- pre_exp_results %>%
      dplyr::select(multiplier, mean_diff_ML, mean_diff_MG, mean_diff_MA) %>%
      tidyr::pivot_longer(cols = -multiplier, 
                 names_to = "comparison", 
                 values_to = "mean_diff") %>%
      dplyr::mutate(comparison = factor(comparison,
                             levels = c("mean_diff_ML", "mean_diff_MG", "mean_diff_MA"),
                             labels = c("M vs L", "M vs G", "M vs A")))

    y_diff_range <- range(diff_data$mean_diff, na.rm = TRUE)
    y_diff_lim <- c(0, min(1, y_diff_range[2] * 1.08))

    p_diff_vs_mult <- ggplot2::ggplot(diff_data, ggplot2::aes(x = multiplier, y = mean_diff, 
                                         color = comparison, 
                                         linetype = comparison)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_point(size = 2)
    if (!is.null(best_plateau)) {
      p_diff_vs_mult <- p_diff_vs_mult +
        ggplot2::geom_vline(xintercept = best_plateau$start_multiplier, 
                   linetype = "dashed", color = "darkgreen", linewidth = 0.9) +
        ggplot2::annotate("text", x = best_plateau$start_multiplier, y = y_diff_lim[2] * 0.95,
                 label = sprintf("Plateau (%.2fx)", best_plateau$start_multiplier),
                 color = "darkgreen", hjust = -0.05, size = 3)
    }
    p_diff_vs_mult <- p_diff_vs_mult +
      ggplot2::scale_x_continuous(breaks = seq(1, 10, 1)) +
      ggplot2::scale_y_continuous(limits = y_diff_lim, expand = ggplot2::expansion(mult = 0.03)) +
      ggplot2::scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
      ggplot2::labs(
        title = "Mean absolute efficiency difference vs multiplier",
        subtitle = "Larger difference = better frontier differentiation",
        x = "Core/non-core weight multiplier",
        y = "Mean absolute difference",
        color = "Frontier comparison",
        linetype = "Frontier comparison"
      ) +
      theme_publication() +
      ggplot2::theme(legend.position = "bottom")

    # Combine plots
    p_pre_exp_summary <- (p_cor_vs_mult / p_eff_range / p_diff_vs_mult) +
      patchwork::plot_annotation(
        title = "Pre-Experiment: Optimal Weight Multiplier Search",
        subtitle = sprintf("Testing %d multipliers from 1.0x to 50.0x | n = %d observations",
                      length(multipliers), ifelse(is.null(n_subsample), nrow(df_clean), n_subsample)),
        caption = "Goal: Identify correlation plateau region and select optimal multiplier at plateau start",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5))
      )

    ggplot2::ggsave(file.path(results_base, "diagnostics/pre_experiment_multiplier_analysis.png"),
           p_pre_exp_summary, width = 12, height = 14, dpi = 300)

    message("✓ Pre-experiment visualization saved")

    # ============================================================================
    # OPTIMAL MULTIPLIER IDENTIFICATION: Correlation Plateau Detection
    # ============================================================================

    message("\nIdentifying optimal multiplier using correlation plateau detection...")
    message(sprintf("  Plateau criterion: cor_change < %.4f (stricter); select first plateau with >= %d points.",
                    plateau_threshold, min_plateau_points))

    # plateau_summary and best_plateau already computed before visualization
    if (nrow(plateau_summary) > 0) {
      # Select the first (earliest) point in the largest plateau
      best_plateau <- plateau_summary[1, ]
      best_multiplier <- best_plateau$start_multiplier
      
      message("✓ Correlation plateau detected:")
      message(sprintf("  Plateau range: %.2fx to %.2fx (%d points)", 
                      best_plateau$start_multiplier, 
                      best_plateau$end_multiplier,
                      best_plateau$n_points))
      message(sprintf("  Mean correlation in plateau: %.4f", best_plateau$mean_cor))
      message(sprintf("  Selected multiplier (plateau start): %.2fx", best_multiplier))
      
      # Get statistics for selected multiplier
      selected_stats <- pre_exp_results %>%
        dplyr::filter(multiplier == best_multiplier) %>%
        dplyr::slice(1)
      
      if (nrow(selected_stats) > 0) {
        message(sprintf("  Mean correlation: %.4f", selected_stats$cor_mean))
        message(sprintf("  Correlation range: [%.4f, %.4f]", 
                        selected_stats$cor_min, selected_stats$cor_max))
        message(sprintf("  Efficiency range: [%.3f, %.3f]", 
                        min(selected_stats$M_min, selected_stats$L_min, 
                            selected_stats$G_min, selected_stats$A_min),
                        max(selected_stats$M_max, selected_stats$L_max, 
                            selected_stats$G_max, selected_stats$A_max)))
      }
      
      # Update weight configuration with optimal multiplier
      optimal_weights <- generate_weights_by_multiplier(best_multiplier)
      
      message("\n✓ Optimal weights generated:")
      for (fid in names(optimal_weights)) {
        w <- optimal_weights[[fid]]
        message(sprintf("  %s-frontier: [%s]", 
                        fid, paste(sprintf("%.2f", w), collapse = ", ")))
      }
      
      # Save optimal configuration
      optimal_config <- list(
        multiplier = best_multiplier,
        weights = optimal_weights,
        statistics = selected_stats,
        plateau_info = best_plateau
      )
      
      saveRDS(optimal_config, 
              file.path(results_base, "diagnostics/optimal_weight_configuration.rds"))
      
      message("\n✓ Optimal configuration saved to optimal_weight_configuration.rds")
      
    } else {
      warning("⚠ No correlation plateau detected. Using first multiplier with correlation < 0.99")
      message("  This may indicate insufficient multiplier range or data issues.")
      
      # Fallback: use first multiplier with reasonable correlation
      fallback_candidates <- pre_exp_results %>%
        dplyr::filter(cor_mean < 0.99) %>%
        dplyr::arrange(multiplier) %>%
        head(1)
      
      if (nrow(fallback_candidates) > 0) {
        best_multiplier <- fallback_candidates$multiplier[1]
        message(sprintf("  Using fallback multiplier: %.2fx (correlation: %.4f)", 
                        best_multiplier, fallback_candidates$cor_mean[1]))
      } else {
        # Last resort: use median multiplier
        best_multiplier <- median(pre_exp_results$multiplier)
        message(sprintf("  Using median multiplier: %.2fx", best_multiplier))
      }
      # Ensure optimal_weights exists for Excel export (fallback path does not set it)
      if (!exists("optimal_weights") || is.null(optimal_weights)) {
        optimal_weights <- generate_weights_by_multiplier(best_multiplier)
      }
    }

    # Export final weight results to Excel (multiplier + detailed weights per frontier/variable)
    tryCatch({
      summary_df <- data.frame(
        Parameter = c("Selected_multiplier", "Plateau_start", "Plateau_end", "Plateau_n_points", "Plateau_mean_correlation"),
        Value = c(
          as.character(best_multiplier),
          if (!is.null(best_plateau)) as.character(best_plateau$start_multiplier) else NA_character_,
          if (!is.null(best_plateau)) as.character(best_plateau$end_multiplier) else NA_character_,
          if (!is.null(best_plateau)) as.character(best_plateau$n_points) else NA_character_,
          if (!is.null(best_plateau)) as.character(round(best_plateau$mean_cor, 6)) else NA_character_
        ),
        stringsAsFactors = FALSE
      )
      weights_rows <- list()
      for (fid in names(optimal_weights)) {
        w <- optimal_weights[[fid]]
        for (v in names(w)) {
          weights_rows[[length(weights_rows) + 1L]] <- data.frame(
            Frontier = fid,
            Variable = v,
            Weight = as.numeric(w[v]),
            stringsAsFactors = FALSE
          )
        }
      }
      weights_df <- do.call(rbind, weights_rows)
      excel_path <- file.path(results_base, "diagnostics", "optimal_weight_results.xlsx")
      writexl::write_xlsx(list(Summary = summary_df, Weights_by_Frontier = weights_df), excel_path)
      message(sprintf("✓ Optimal weight results saved to: %s", excel_path))
    }, error = function(e) {
      message(sprintf("⚠ Could not write optimal weight Excel: %s", e$message))
    })

    message(paste(rep("=", 78), collapse = ""))
    message("PRE-EXPERIMENT COMPLETE")
    message(paste(rep("=", 78), collapse = ""))
  }
  
  # Now update the main weight configuration with the optimal multiplier
  # (This replaces the original hard-coded weights)
  # ============================================================================
  # UPDATE MAIN WEIGHT CONFIGURATION WITH OPTIMAL MULTIPLIER
  # ============================================================================

  message("\nUpdating main weight configuration with optimal multiplier...")

  # Create weight vectors (returns these for use in main analysis)
  # Use same asymmetric base weights as in pre-experiment (non-core=0.5, core=multiplier)
  non_core_base <- 0.5  # Match the pre-experiment strategy
  
  bad_weights_M_raw <- c(
    Waste_water = non_core_base,
    Manure = non_core_base,
    Dead_pig = best_multiplier,  # Core constraint
    Carbon_emission = non_core_base,
    Eutrophication_potential = non_core_base
  )

  bad_weights_L_raw <- c(
    Waste_water = best_multiplier,      # Core
    Manure = best_multiplier,           # Core
    Dead_pig = non_core_base,
    Carbon_emission = non_core_base,
    Eutrophication_potential = non_core_base
  )

  bad_weights_G_raw <- c(
    Waste_water = non_core_base,
    Manure = non_core_base,
    Dead_pig = non_core_base,
    Carbon_emission = best_multiplier,   # Core
    Eutrophication_potential = best_multiplier   # Core
  )

  bad_weights_A_raw <- rep(best_multiplier, length(bad_outputs_full))  # All core
  names(bad_weights_A_raw) <- bad_outputs_full

  # Normalize weights (target_sum = number of bad outputs)
  target_sum <- length(bad_outputs_full)
  bad_weights_M <- normalize_weights(bad_weights_M_raw, target_sum = target_sum)
  bad_weights_L <- normalize_weights(bad_weights_L_raw, target_sum = target_sum)
  bad_weights_G <- normalize_weights(bad_weights_G_raw, target_sum = target_sum)
  bad_weights_A <- normalize_weights(bad_weights_A_raw, target_sum = target_sum)

  # Create weight list
  bad_weights_list <- list(
    E = bad_weights_M,
    L = bad_weights_L,
    G = bad_weights_G,
    A = bad_weights_A
  )

  message(sprintf("✓ Main weights updated with optimal multiplier (%.2fx)", best_multiplier))
  message("  Proceeding with main analysis...\n")

  # Return results
  return(list(
    best_multiplier = best_multiplier,
    bad_weights_M_raw = bad_weights_M_raw,
    bad_weights_L_raw = bad_weights_L_raw,
    bad_weights_G_raw = bad_weights_G_raw,
    bad_weights_A_raw = bad_weights_A_raw,
    bad_weights_M = bad_weights_M,
    bad_weights_L = bad_weights_L,
    bad_weights_G = bad_weights_G,
    bad_weights_A = bad_weights_A,
    bad_weights_list = bad_weights_list
  ))
}









