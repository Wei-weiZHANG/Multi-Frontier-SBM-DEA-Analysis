# Meta-Frontier Diagnostic Visualization Tool
# ============================================
#
# This script creates diagnostic plots to verify meta-frontier construction
#
# Run AFTER main analysis completes
#
# ============================================

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

cat("\n")
cat("==============================================================================\n")
cat("META-FRONTIER DIAGNOSTIC VISUALIZATION\n")
cat("==============================================================================\n\n")

# Display label mapping for diagnostics plots.
frontier_labels <- c(
  "M" = "M-frontier\n(mortality)",
  "L" = "L-frontier\n(local env.)",
  "G" = "G-frontier\n(global env.)",
  "A" = "A-frontier\n(all constraints)"
)
frontier_order_diag <- c("A", "M", "G", "L")

# Check if results_path was passed via environment variable (from main script)
if (nchar(Sys.getenv("CURRENT_RESULTS_DIR")) > 0) {
  results_path <- Sys.getenv("CURRENT_RESULTS_DIR")
  latest_result <- basename(results_path)  # Extract folder name for display
  cat(sprintf("Using current run results: %s\n", latest_result))
} else {
  # Fallback: find latest results folder
  results_base_dir <- "./results"
  all_results <- list.dirs(results_base_dir, full.names = FALSE, recursive = FALSE)
  # Match expected result folder format
  all_results <- all_results[grepl("^results_\\d{8}_\\d{6}$", all_results)]
  if (length(all_results) == 0) {
    stop("No results folders found. Please run the analysis first.")
  }
  latest_result <- sort(all_results, decreasing = TRUE)[1]
  results_path <- file.path(results_base_dir, latest_result)
  cat(sprintf("Using latest results folder: %s\n", latest_result))
}

# Load data
eff <- read_csv(file.path(results_path, "efficiency_calculation/all_frontiers_efficiency_with_tgr.csv"),
                show_col_types = FALSE)

cat(sprintf("Loaded results from: %s\n", latest_result))
cat(sprintf("Sample size: %d farms\n\n", nrow(eff)))

# Create diagnostics folder if it doesn't exist
diagnostics_dir <- file.path(results_path, "diagnostics")
if (!dir.exists(diagnostics_dir)) {
  dir.create(diagnostics_dir, recursive = TRUE)
  cat("✓ Created diagnostics folder\n\n")
}

# ============================================================================
# DIAGNOSTIC 1: Meta vs Individual Frontiers
# ============================================================================

cat("Creating Diagnostic 1: Meta vs Individual Frontiers...\n")

diag1_data <- eff %>%
  select(obs_id, Meta_TE_VRS, M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS) %>%
  pivot_longer(cols = c(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS),
               names_to = "frontier",
               values_to = "TE") %>%
  mutate(frontier = substr(frontier, 1, 1))

p_diag1 <- ggplot(diag1_data, aes(x = TE, y = Meta_TE_VRS)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  geom_point(alpha = 0.4, size = 2) +
  facet_wrap(~frontier, ncol = 4, labeller = labeller(frontier = c(
    "M" = frontier_labels[["M"]],
    "L" = frontier_labels[["L"]],
    "G" = frontier_labels[["G"]],
    "A" = frontier_labels[["A"]]
  ))) +
  labs(
    title = "Diagnostic 1: Meta-frontier vs group frontiers",
    subtitle = "Points should be on or above red line (meta >= group)",
    x = "Group frontier efficiency",
    y = "Meta-frontier efficiency",
    caption = "Red line: meta = group | Points above: meta > group (correct) | Points below: error"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

diag1_file <- file.path(results_path, "diagnostics/DIAGNOSTIC_meta_vs_groups.png")
ggsave(diag1_file, p_diag1, width = 12, height = 4, dpi = 300)
if (file.exists(diag1_file)) {
  cat("  ✓ Saved\n")
} else {
  warning("  ⚠ Failed to save Diagnostic 1")
}

# ============================================================================
# DIAGNOSTIC 2: TGR = 1.0 Distribution
# ============================================================================

cat("Creating Diagnostic 2: TGR=1 Distribution...\n")

# For each farm, which frontier has TGR ≈ 1.0?
tgr_matrix <- cbind(
  E = eff$M_TGR_VRS,
  L = eff$L_TGR_VRS,
  G = eff$G_TGR_VRS,
  A = eff$A_TGR_VRS
)

max_tgr_frontier <- apply(tgr_matrix, 1, function(x) {
  if (all(is.na(x))) return(NA)
  names(which.max(x))
})

tgr_1_source <- table(max_tgr_frontier)

diag2_data <- data.frame(
  Frontier = names(tgr_1_source),
  Count = as.numeric(tgr_1_source),
  Percentage = as.numeric(tgr_1_source) / sum(tgr_1_source) * 100
)

p_diag2 <- ggplot(diag2_data, aes(x = Frontier, y = Count, fill = Frontier)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            vjust = -0.3, fontface = "bold") +
  scale_fill_viridis_d(option = "D") +
  labs(
    title = "Diagnostic 2: Which Frontier Provides Meta for Each Farm?",
    subtitle = "Each farm should have TGR≈1.0 on ONE frontier (their best)",
    x = "Frontier with TGR ≈ 1.0",
    y = "Number of Farms",
    caption = "This shows which regulatory scenario enables peak performance for each farm"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Replace x-axis display with M-Frontier (Mortality) as the standard label.
p_diag2 <- p_diag2 +
  scale_x_discrete(
    breaks = frontier_order_diag[frontier_order_diag %in% diag2_data$Frontier],
    labels = frontier_labels[frontier_order_diag[frontier_order_diag %in% diag2_data$Frontier]]
  )

diag2_file <- file.path(results_path, "diagnostics/DIAGNOSTIC_tgr_1_sources.png")
ggsave(diag2_file, p_diag2, width = 10, height = 6, dpi = 300)
if (file.exists(diag2_file)) {
  cat("  ✓ Saved\n")
} else {
  warning("  ⚠ Failed to save Diagnostic 2")
}
if (file.exists(diag2_file)) {
  cat("  ✓ Saved\n")
} else {
  warning("  ⚠ Failed to save Diagnostic 2")
}

# ============================================================================
# DIAGNOSTIC 3: Meta-Frontier Source vs Efficiency Level
# ============================================================================

cat("Creating Diagnostic 3: Meta Source vs Efficiency Level...\n")

if ("Meta_Source" %in% names(eff)) {
  diag3_data <- eff %>%
    select(Meta_Source, Meta_TE_VRS, M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS) %>%
    mutate(
      Meta_Source = factor(
        Meta_Source,
        levels = c("M", "L", "G", "A"),
        labels = frontier_labels[c("M", "L", "G", "A")]
      )
    )
  
  p_diag3 <- ggplot(diag3_data, aes(x = Meta_Source, y = Meta_TE_VRS, fill = Meta_Source)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.3, alpha = 0.7) +
    scale_fill_viridis_d(option = "D") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Diagnostic 3: Meta-Frontier Efficiency by Source",
      subtitle = "Distribution of meta-frontier efficiency grouped by which frontier provides it",
      x = "Meta-Frontier Source",
      y = "Meta-Frontier Efficiency",
      caption = "Shows if certain frontiers consistently provide higher meta-frontier values"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  diag3_file <- file.path(results_path, "diagnostics/DIAGNOSTIC_meta_by_source.png")
  ggsave(diag3_file, p_diag3, width = 10, height = 6, dpi = 300)
  if (file.exists(diag3_file)) {
    cat("  ✓ Saved\n")
  } else {
    warning("  ⚠ Failed to save Diagnostic 3")
  }
}

# ============================================================================
# DIAGNOSTIC 4: TGR Heatmap by Farm
# ============================================================================

cat("Creating Diagnostic 4: TGR Heatmap...\n")

# Sample 30 farms for clearer visualization
if (nrow(eff) > 30) {
  sample_farms <- sample(seq_len(nrow(eff)), min(30, nrow(eff)))
} else {
  sample_farms <- seq_len(nrow(eff))
}

diag4_data <- eff[sample_farms, ] %>%
  select(obs_id, M_TGR_VRS, L_TGR_VRS, G_TGR_VRS, A_TGR_VRS) %>%
  pivot_longer(cols = c(M_TGR_VRS, L_TGR_VRS, G_TGR_VRS, A_TGR_VRS),
               names_to = "Frontier",
               values_to = "TGR") %>%
  mutate(
    Frontier = substr(Frontier, 1, 1),
    Farm_ID = sprintf("Farm %d", obs_id),
    Frontier_Label = frontier_labels[as.character(Frontier)]
  )

p_diag4 <- ggplot(diag4_data, aes(x = Frontier_Label, y = Farm_ID, fill = TGR)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", labels = scales::percent, limits = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f", TGR)), color = "white", size = 2.5, fontface = "bold") +
  labs(
    title = "Diagnostic 4: TGR Heatmap (Sample of Farms)",
    subtitle = "Each row should have at least ONE cell ≈ 1.00 (yellow)",
    x = "Frontier",
    y = "",
    fill = "TGR",
    caption = "Yellow (TGR≈1.0) shows which frontier provides meta for that farm"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

diag4_file <- file.path(results_path, "diagnostics/DIAGNOSTIC_tgr_heatmap.png")
ggsave(diag4_file, p_diag4, width = 8, height = 10, dpi = 300)
if (file.exists(diag4_file)) {
  cat("  ✓ Saved\n")
} else {
  warning("  ⚠ Failed to save Diagnostic 4")
}

# ============================================================================
# DIAGNOSTIC 5: Consistency Check
# ============================================================================

cat("Creating Diagnostic 5: TGR=1 Consistency Check...\n")

# Check if TGR=1 matches Meta_Source
# THEORETICAL EXPECTATION: If Meta_Source = "M", then M_TGR should = 1.0
# (because Meta_TE = max(M, L, G, A) = M_TE, so M_TGR = M_TE / Meta_TE = 1.0)
# This diagnostic verifies that Meta_Source is correctly identified

if ("Meta_Source" %in% names(eff)) {
  consistency_data <- data.frame()
  
  for (i in seq_len(nrow(eff))) {
    source <- eff$Meta_Source[i]
    if (is.na(source)) next
    
    # Get TGR on the source frontier
    tgr_on_source_raw <- eff[[paste0(source, "_TGR_VRS")]][i]
    
    # Also get the actual TE values for verification
    te_source <- eff[[paste0(source, "_TE_VRS")]][i]
    meta_te <- eff$Meta_TE_VRS[i]
    
    # Calculate what TGR SHOULD be (for verification)
    tgr_should_be <- ifelse(meta_te > 0, te_source / meta_te, NA)
    
    # CRITICAL: Ensure TGR is in [0, 1] format (not percentage)
    # TGR should be in [0, 1] range. Check format and convert if needed
    tgr_on_source <- tgr_on_source_raw
    
    if (!is.na(tgr_on_source)) {
      # Check if TGR is in percentage format (> 1.0)
      if (tgr_on_source > 1.0) {
        # Definitely in percentage format (e.g., 90.5 instead of 0.905)
        tgr_on_source <- tgr_on_source / 100
        if (i <= 3) {  # Warn for first few only
          cat(sprintf("  ⚠ Farm %d: TGR converted from percentage (%.2f%% → %.4f)\n", 
                     eff$obs_id[i], tgr_on_source_raw, tgr_on_source))
        }
      }
      
      # Additional validation: compare with calculated tgr_should_be
      if (!is.na(tgr_should_be) && abs(tgr_on_source - tgr_should_be) > 0.2) {
        # Large discrepancy - might indicate format issue
        # If raw TGR is close to should_be * 100, it's likely percentage
        if (abs(tgr_on_source_raw - tgr_should_be * 100) < abs(tgr_on_source_raw - tgr_should_be)) {
          tgr_on_source <- tgr_on_source_raw / 100
          if (i <= 3) {
            cat(sprintf("  ⚠ Farm %d: TGR format corrected via validation (%.2f → %.4f)\n",
                       eff$obs_id[i], tgr_on_source_raw, tgr_on_source))
          }
        }
      }
      
      # Final safety: bound to [0, 1] range
      tgr_on_source <- pmax(0, pmin(1, tgr_on_source))
    }
    
    # More nuanced classification (using decimal format 0-1)
    # IMPORTANT: Classification must match actual TGR values!
    status <- case_when(
      is.na(tgr_on_source) ~ "Missing",
      tgr_on_source >= 0.999 ~ "Excellent (≥99.9%)",
      tgr_on_source >= 0.99 ~ "Good (≥99%)",
      tgr_on_source >= 0.95 ~ "Acceptable (≥95%)",
      TRUE ~ "Problem (<95%)"
    )
    
    # Double-check: if TGR < 0.95, status MUST be "Problem"
    if (!is.na(tgr_on_source) && tgr_on_source < 0.95 && status != "Problem (<95%)") {
      # Force correct classification
      status <- "Problem (<95%)"
      if (i <= 3) {
        cat(sprintf("  ⚠ Farm %d: Status corrected (TGR=%.4f < 0.95 → Problem)\n",
                   eff$obs_id[i], tgr_on_source))
      }
    }
    
    consistency_data <- rbind(consistency_data, data.frame(
      obs_id = eff$obs_id[i],
      Meta_Source = source,
      TGR_on_source = tgr_on_source,
      TE_source = te_source,
      Meta_TE = meta_te,
      TGR_should_be = tgr_should_be,
      Status = status,
      stringsAsFactors = FALSE
    ))
  }
  
  # Debug: Check TGR value range
  cat("\n  DEBUG: TGR value range check:\n")
  cat(sprintf("    Min TGR: %.6f\n", min(consistency_data$TGR_on_source, na.rm = TRUE)))
  cat(sprintf("    Max TGR: %.6f\n", max(consistency_data$TGR_on_source, na.rm = TRUE)))
  cat(sprintf("    Mean TGR: %.6f\n", mean(consistency_data$TGR_on_source, na.rm = TRUE)))
  cat(sprintf("    TGR values > 1: %d\n", sum(consistency_data$TGR_on_source > 1, na.rm = TRUE)))
  cat(sprintf("    TGR values < 0.95: %d\n", sum(consistency_data$TGR_on_source < 0.95, na.rm = TRUE)))
  
  # CRITICAL FIX: Re-classify based on ACTUAL TGR values (force correct classification)
  consistency_data <- consistency_data %>%
    mutate(
      Status = case_when(
        is.na(TGR_on_source) ~ "Missing",
        TGR_on_source >= 0.999 ~ "Excellent (≥99.9%)",
        TGR_on_source >= 0.99 ~ "Good (≥99%)",
        TGR_on_source >= 0.95 ~ "Acceptable (≥95%)",
        TRUE ~ "Problem (<95%)"
      )
    )
  
  # Calculate statistics (after forced re-classification)
  n_excellent <- sum(consistency_data$Status == "Excellent (≥99.9%)", na.rm = TRUE)
  n_good <- sum(consistency_data$Status == "Good (≥99%)", na.rm = TRUE)
  n_acceptable <- sum(consistency_data$Status == "Acceptable (≥95%)", na.rm = TRUE)
  n_problem <- sum(consistency_data$Status == "Problem (<95%)", na.rm = TRUE)
  n_total <- nrow(consistency_data)
  
  # Verify classification is correct
  n_should_be_problem <- sum(consistency_data$TGR_on_source < 0.95, na.rm = TRUE)
  if (n_problem != n_should_be_problem) {
    cat(sprintf("\n  ✗ ERROR: Classification mismatch!\n"))
    cat(sprintf("    Farms with TGR < 0.95: %d\n", n_should_be_problem))
    cat(sprintf("    Farms classified as 'Problem': %d\n", n_problem))
    cat("    Forcing correction...\n")
    
    # Force correct classification
    consistency_data <- consistency_data %>%
      mutate(
        Status = ifelse(TGR_on_source < 0.95 & !is.na(TGR_on_source), 
                       "Problem (<95%)", Status)
      )
    
    # Recalculate
    n_problem <- sum(consistency_data$Status == "Problem (<95%)", na.rm = TRUE)
    n_excellent <- sum(consistency_data$Status == "Excellent (≥99.9%)", na.rm = TRUE)
    n_good <- sum(consistency_data$Status == "Good (≥99%)", na.rm = TRUE)
    n_acceptable <- sum(consistency_data$Status == "Acceptable (≥95%)", na.rm = TRUE)
  }
  
  # Verify statistics match actual data
  cat("\n  DEBUG: Status distribution:\n")
  status_table <- table(consistency_data$Status, useNA = "ifany")
  print(status_table)
  
  # Verify data before plotting
  cat("\n  DEBUG: Sample of consistency_data (first 10 rows):\n")
  print(head(consistency_data[, c("Meta_Source", "TGR_on_source", "Status")], 10))
  
  # Check for problematic farms
  problem_farms <- consistency_data %>%
    filter(TGR_on_source < 0.95) %>%
    select(obs_id, Meta_Source, TGR_on_source, Status, TE_source, Meta_TE, TGR_should_be)
  
  if (nrow(problem_farms) > 0) {
    cat(sprintf("\n  ⚠ PROBLEM FARMS (TGR < 0.95): %d farms\n", nrow(problem_farms)))
    print(problem_farms)
    cat("\n  Analysis: These farms have TGR < 95% on their meta-source frontier.\n")
    cat("  This suggests Meta_Source may be incorrectly identified for these farms.\n")
    cat("  Possible causes:\n")
    cat("    1. Numerical precision issues in which.max()\n")
    cat("    2. Multiple frontiers tied at maximum (rounding differences)\n")
    cat("    3. Data quality issues\n")
  } else {
    cat("\n  ✓ No farms with TGR < 0.95\n")
  }

  # Map frontier codes to M/L/G/A labels for display in the plot.
  consistency_data <- consistency_data %>%
    mutate(
      Meta_Source = factor(Meta_Source, levels = frontier_order_diag),
      Meta_Source_Label = frontier_labels[as.character(Meta_Source)]
    )
  
  # Create improved plot with better color scheme
  p_diag5 <- ggplot(consistency_data, aes(x = Meta_Source_Label, y = TGR_on_source, 
                                           color = Status)) +
    geom_jitter(width = 0.2, alpha = 0.7, size = 2.5) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_hline(yintercept = 0.99, linetype = "dotted", color = "orange", linewidth = 0.8) +
    geom_hline(yintercept = 0.95, linetype = "dotted", color = "yellow3", linewidth = 0.5, alpha = 0.5) +
    scale_color_manual(
      values = c(
        "Excellent (≥99.9%)" = "darkgreen",
        "Good (≥99%)" = "green3",
        "Acceptable (≥95%)" = "orange",
        "Problem (<95%)" = "red",
        "Missing" = "gray50"
      ),
      breaks = c("Excellent (≥99.9%)", "Good (≥99%)", "Acceptable (≥95%)", "Problem (<95%)", "Missing")
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0.9, 1.01), 
                       breaks = c(0.90, 0.95, 0.99, 1.00)) +
    labs(
      title = "Diagnostic 5: TGR Consistency with Meta Source",
      subtitle = sprintf("Theoretical: TGR should = 1.0 | Excellent: %d (%.0f%%) | Good: %d (%.0f%%) | Problem: %d (%.0f%%)",
                        n_excellent, n_excellent/n_total*100,
                        n_good, n_good/n_total*100,
                        n_problem, n_problem/n_total*100),
      x = "Meta-Frontier Source",
      y = "TGR on Meta-Source Frontier",
      color = "Status",
      caption = "Red dashed: TGR=1.0 (theoretical) | Orange: 99% threshold | Yellow: 95% threshold\nIf TGR < 95%, Meta_Source may be incorrectly identified"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  
  diag5_file <- file.path(results_path, "diagnostics/DIAGNOSTIC_tgr_consistency.png")
  ggsave(diag5_file, p_diag5, width = 12, height = 7, dpi = 300)
  if (file.exists(diag5_file)) {
    cat("  ✓ Saved\n")
  } else {
    warning("  ⚠ Failed to save Diagnostic 5")
  }
  
  # Detailed report
  cat("\n  TGR Consistency Report:\n")
  cat(sprintf("    Excellent (≥99.9%%): %d farms (%.1f%%)\n", n_excellent, n_excellent/n_total*100))
  cat(sprintf("    Good (≥99%%): %d farms (%.1f%%)\n", n_good, n_good/n_total*100))
  cat(sprintf("    Acceptable (≥95%%): %d farms (%.1f%%)\n", n_acceptable, n_acceptable/n_total*100))
  cat(sprintf("    Problem (<95%%): %d farms (%.1f%%)\n", n_problem, n_problem/n_total*100))
  
  if (n_problem > 0) {
    cat("\n  ⚠ WARNING: Some farms have TGR < 95% on their meta-source frontier!\n")
    cat("    This suggests Meta_Source may be incorrectly identified for these farms.\n")
    cat("    Possible causes:\n")
    cat("      1. Numerical precision issues in max() calculation\n")
    cat("      2. Multiple frontiers tied at maximum (rounding differences)\n")
    cat("      3. Data quality issues\n")
  } else if (n_excellent + n_good >= n_total * 0.95) {
    cat("\n  ✓ PASS: Most farms (≥95%%) have TGR ≥ 99% on meta-source frontier\n")
    cat("    Meta_Source identification is consistent with TGR values.\n")
  }
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("==============================================================================\n")
cat("DIAGNOSTIC PLOTS CREATED\n")
cat("==============================================================================\n\n")

cat("5 diagnostic plots saved to:\n")
cat(sprintf("%s/visualization/\n\n", results_path))

cat("Files:\n")
cat("  1. DIAGNOSTIC_meta_vs_groups.png - Meta vs each frontier (should be above line)\n")
cat("  2. DIAGNOSTIC_tgr_1_sources.png - Which frontier provides meta (bar chart)\n")
cat("  3. DIAGNOSTIC_meta_by_source.png - Meta efficiency by source\n")
cat("  4. DIAGNOSTIC_tgr_heatmap.png - TGR pattern across farms\n")
cat("  5. DIAGNOSTIC_tgr_consistency.png - TGR=1 consistency check\n\n")

cat("Open diagnostics:\n")
cat(sprintf("  shell.exec('%s/visualization')\n\n", results_path))

cat("==============================================================================\n\n")













