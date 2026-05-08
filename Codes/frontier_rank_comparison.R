# ============================================================================
# Rank Comparison Analysis: Comparing Sample Rankings Across Frontiers
# ============================================================================
#
# Purpose: Analyze how 62 samples rank differently across M, L, G, A frontiers
#
# Analysis includes:
# 1. Calculate ranks for each sample in each frontier
# 2. Rank correlation analysis (Spearman correlation)
# 3. Identify samples with largest rank changes
# 4. Create essential visualization (correlation heatmap)
#
# - Only essential figure: rank correlation heatmap
# - Only essential tables: correlation matrix and summary statistics
# - Integrated into main analysis workflow
#
# Usage:
#   - Automatically runs as part of run_pipeline.R
#   - Can also be run standalone: source("frontier_rank_comparison.R")
#
# ============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(corrplot)
library(tibble)  # For rownames_to_column

cat("\n")
cat("==============================================================================\n")
cat("        RANK COMPARISON ANALYSIS ACROSS FRONTIERS\n")
cat("==============================================================================\n\n")

# ============================================================================
# 1. Load Data
# ============================================================================

# Check if running as part of main analysis (results_base already defined)
if (exists("results_base") && exists("efficiency_results")) {
  # Running as module from main script
  results_dir <- results_base
  cat("✓ Running as integrated module (using current session data)\n")
} else {
  # Running standalone - find latest results
  results_base_dir <- "./results"
  all_results <- list.dirs(results_base_dir, full.names = FALSE, recursive = FALSE)
  all_results <- all_results[grepl("^results_\\d{8}_\\d{6}$", all_results)]
  latest_result <- sort(all_results, decreasing = TRUE)[1]
  results_dir <- file.path(results_base_dir, latest_result)
  
  cat(sprintf("✓ Loading results: %s\n", latest_result))
  
  # Load efficiency data
  eff_file <- file.path(results_dir, "efficiency_calculation/frontier_efficiency_results.csv")
  if (!file.exists(eff_file)) {
    stop("Efficiency data file not found")
  }
  efficiency_results <- read_csv(eff_file, show_col_types = FALSE)
  cat("✓ Loaded efficiency data from file\n")
}

# Get sample identifier (Farm or ID column)
id_col <- NULL
if ("Farm" %in% names(efficiency_results)) {
  id_col <- "Farm"
} else if ("ID" %in% names(efficiency_results)) {
  id_col <- "ID"
} else {
  # Create ID if not exists
  efficiency_results$Sample_ID <- paste0("Sample_", 1:nrow(efficiency_results))
  id_col <- "Sample_ID"
  cat("⚠ No Farm/ID column found, created Sample_ID\n")
}

cat(sprintf("✓ Using identifier column: %s\n", id_col))
cat(sprintf("✓ Total samples: %d\n\n", nrow(efficiency_results)))

# ============================================================================
# 2. Calculate Ranks for Each Frontier
# ============================================================================

cat("==============================================================================\n")
cat("Calculating Ranks\n")
cat("==============================================================================\n\n")

# Extract efficiency values and calculate ranks
# Higher efficiency = better rank (rank 1 = best)
rank_data <- efficiency_results %>%
  select(all_of(id_col), M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS) %>%
  filter(complete.cases(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS)) %>%
  mutate(
    Rank_M = rank(-M_TE_VRS, ties.method = "min"),  # Negative for descending order
    Rank_L = rank(-L_TE_VRS, ties.method = "min"),
    Rank_G = rank(-G_TE_VRS, ties.method = "min"),
    Rank_A = rank(-A_TE_VRS, ties.method = "min")
  )

n_samples <- nrow(rank_data)
cat(sprintf("✓ Calculated ranks for %d samples\n\n", n_samples))

# ============================================================================
# 3. Rank Correlation Analysis
# ============================================================================

cat("==============================================================================\n")
cat("Rank Correlation Analysis (Spearman Correlation)\n")
cat("==============================================================================\n\n")

# Calculate Spearman correlation matrix
rank_matrix <- rank_data %>%
  select(Rank_M, Rank_L, Rank_G, Rank_A) %>%
  as.matrix()

cor_matrix <- cor(rank_matrix, method = "spearman")

cat("Spearman Rank Correlation Matrix:\n")
cat("(Values closer to 1 indicate more similar rankings)\n\n")
print(round(cor_matrix, 3))

# Extract key correlations
cat("\nKey Pairwise Correlations:\n")
cat(sprintf("  M vs L: ρ = %.3f\n", cor_matrix["Rank_M", "Rank_L"]))
cat(sprintf("  M vs G: ρ = %.3f\n", cor_matrix["Rank_M", "Rank_G"]))
cat(sprintf("  M vs A: ρ = %.3f\n", cor_matrix["Rank_M", "Rank_A"]))
cat(sprintf("  L vs G: ρ = %.3f\n", cor_matrix["Rank_L", "Rank_G"]))
cat(sprintf("  L vs A: ρ = %.3f\n", cor_matrix["Rank_L", "Rank_A"]))
cat(sprintf("  G vs A: ρ = %.3f\n", cor_matrix["Rank_G", "Rank_A"]))

# Statistical significance tests
cat("\nStatistical Significance Tests:\n")
frontiers <- c("M", "L", "G", "A")
for (i in 1:(length(frontiers)-1)) {
  for (j in (i+1):length(frontiers)) {
    f1 <- frontiers[i]
    f2 <- frontiers[j]
    test_result <- cor.test(rank_data[[paste0("Rank_", f1)]], 
                           rank_data[[paste0("Rank_", f2)]], 
                           method = "spearman", exact = FALSE)
    sig <- ifelse(test_result$p.value < 0.001, "***",
                 ifelse(test_result$p.value < 0.01, "**",
                       ifelse(test_result$p.value < 0.05, "*", "ns")))
    cat(sprintf("  %s vs %s: ρ = %.3f, p = %.4f %s\n", 
                f1, f2, test_result$estimate, test_result$p.value, sig))
  }
}

cat("\n")

# ============================================================================
# 4. Rank Difference Analysis
# ============================================================================

cat("==============================================================================\n")
cat("Rank Difference Analysis\n")
cat("==============================================================================\n\n")

# Calculate rank differences (using A as reference)
rank_diff_data <- rank_data %>%
  mutate(
    Diff_M_A = Rank_M - Rank_A,  # Positive = worse in E than A
    Diff_L_A = Rank_L - Rank_A,
    Diff_G_A = Rank_G - Rank_A,
    Max_Rank_Change = pmax(abs(Diff_M_A), abs(Diff_L_A), abs(Diff_G_A)),
    Mean_Rank_Change = (abs(Diff_M_A) + abs(Diff_L_A) + abs(Diff_G_A)) / 3
  )

# Summary statistics
cat("Rank Differences (relative to A-Frontier (Aggregated)):\n")
cat(sprintf("  Mean |E - A|: %.2f ranks\n", mean(abs(rank_diff_data$Diff_M_A))))
cat(sprintf("  Mean |L - A|: %.2f ranks\n", mean(abs(rank_diff_data$Diff_L_A))))
cat(sprintf("  Mean |G - A|: %.2f ranks\n", mean(abs(rank_diff_data$Diff_G_A))))
cat(sprintf("  Mean max change: %.2f ranks\n", mean(rank_diff_data$Max_Rank_Change)))
cat("\n")

# Identify samples with largest rank changes
cat("Top 10 Samples with Largest Rank Changes:\n")
top_changes <- rank_diff_data %>%
  arrange(desc(Max_Rank_Change)) %>%
  head(10) %>%
  select(all_of(id_col), Rank_M, Rank_L, Rank_G, Rank_A, 
         Diff_M_A, Diff_L_A, Diff_G_A, Max_Rank_Change)

print(top_changes)
cat("\n")

# Identify samples that improved most in A vs E
cat("Top 10 Samples with Largest Improvement in A vs E:\n")
top_improvements <- rank_diff_data %>%
  arrange(desc(Diff_M_A)) %>%  # Positive = better in A
  head(10) %>%
  select(all_of(id_col), Rank_M, Rank_A, Diff_M_A, M_TE_VRS, A_TE_VRS)

print(top_improvements)
cat("\n")

# ============================================================================
# 5. Create Visualizations (Core Figure Set)
# ============================================================================

cat("==============================================================================\n")
cat("Creating Visualization\n")
cat("==============================================================================\n\n")

output_dir <- file.path(results_dir, "visualization")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Figure: Rank Correlation Heatmap (Essential visualization)
cat("Creating rank correlation heatmap...\n")
png(file.path(output_dir, "rank_correlation_heatmap.png"), 
    width = 9, height = 9, units = "in", res = 300)  # Increased width from 8 to 9 to accommodate labels

# mar: bottom, left, top, right margins
# Increase bottom and right margins to accommodate horizontal labels and ensure full display
# tl.offset: distance between labels and axis (increase for better spacing)
# tl.cex: label size (adjust if needed for readability)
corrplot(cor_matrix, method = "color", type = "upper", 
         order = "original", tl.col = "black", tl.srt = 0, 
         tl.cex = 1.2, tl.offset = 1.0,  # Adjusted offset for better fit
         addCoef.col = "black", number.cex = 0.9,
         col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(200),
         mar = c(4, 0, 2, 2))  # Added right margin (2) to prevent right-side label cutoff

dev.off()
cat("✓ Saved: rank_correlation_heatmap.png\n")

# ============================================================================
# 6. Save Results Tables (Core Result Tables)
# ============================================================================

cat("\n==============================================================================\n")
cat("Saving Results Tables\n")
cat("==============================================================================\n\n")

statistical_tests_dir <- file.path(results_dir, "statistical_tests")
if (!dir.exists(statistical_tests_dir)) {
  dir.create(statistical_tests_dir, recursive = TRUE)
}

# Save correlation matrix (essential)
cor_output_file <- file.path(statistical_tests_dir, "rank_correlation_matrix.csv")
write_csv(as.data.frame(cor_matrix) %>% rownames_to_column("Frontier"), cor_output_file)
cat(sprintf("✓ Saved: %s\n", cor_output_file))

# Save summary statistics (essential)
summary_stats <- data.frame(
  Metric = c("Mean |E - A|", "Mean |L - A|", "Mean |G - A|", 
             "Mean Max Change", "SD Max Change", "Average Correlation"),
  Value = c(
    mean(abs(rank_diff_data$Diff_M_A)),
    mean(abs(rank_diff_data$Diff_L_A)),
    mean(abs(rank_diff_data$Diff_G_A)),
    mean(rank_diff_data$Max_Rank_Change),
    sd(rank_diff_data$Max_Rank_Change),
    mean(cor_matrix[lower.tri(cor_matrix)])
  )
)

summary_file <- file.path(statistical_tests_dir, "rank_comparison_summary.csv")
write_csv(summary_stats, summary_file)
cat(sprintf("✓ Saved: %s\n", summary_file))

# ============================================================================
# 7. Summary Report
# ============================================================================

cat("\n==============================================================================\n")
cat("Summary Report\n")
cat("==============================================================================\n\n")

cat("Key Findings:\n")
cat("-------------\n\n")

# Rank consistency
mean_cor <- mean(cor_matrix[lower.tri(cor_matrix)])
cat(sprintf("1. Average rank correlation: %.3f\n", mean_cor))
if (mean_cor > 0.8) {
  cat("   → Rankings are highly consistent across frontiers\n")
} else if (mean_cor > 0.6) {
  cat("   → Rankings show moderate consistency\n")
} else {
  cat("   → Rankings show substantial variation across frontiers\n")
}

cat("\n")

# Most variable frontier
var_by_frontier <- c(
  E = sd(rank_diff_data$Diff_M_A),
  L = sd(rank_diff_data$Diff_L_A),
  G = sd(rank_diff_data$Diff_G_A)
)
most_variable <- names(var_by_frontier)[which.max(var_by_frontier)]
cat(sprintf("2. Most variable frontier (vs A): %s (SD = %.2f ranks)\n", 
            most_variable, max(var_by_frontier)))

cat("\n")

# Samples with consistent high/low ranks
consistent_high <- rank_data %>%
  filter(Rank_M <= 10 & Rank_L <= 10 & Rank_G <= 10 & Rank_A <= 10) %>%
  nrow()

consistent_low <- rank_data %>%
  filter(Rank_M > (n_samples - 10) & Rank_L > (n_samples - 10) & 
         Rank_G > (n_samples - 10) & Rank_A > (n_samples - 10)) %>%
  nrow()

cat(sprintf("3. Samples consistently in top 10: %d\n", consistent_high))
cat(sprintf("4. Samples consistently in bottom 10: %d\n", consistent_low))

cat("\n")
cat("Analysis complete!\n\n")













