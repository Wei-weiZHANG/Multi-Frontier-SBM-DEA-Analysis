# ============================================================================
# Multi-Frontier SBM Analysis
# 
# Usage: source("run_pipeline.R")

cat("\n")
cat("==============================================================================\n")
cat("    MULTI-FRONTIER SBM ANALYSIS LAUNCHER\n")
cat("==============================================================================\n\n")

start_time <- Sys.time()

# ============================================================================
# STEP 1: PRE-FLIGHT CHECKS
# ============================================================================

cat("STEP 1: Pre-flight Checks\n")
cat("---------------------------\n")

# Check R version
r_version <- as.numeric(R.version$major) + as.numeric(R.version$minor)/10
cat(sprintf("✓ R release: %s ", R.version.string))
if (r_version >= 4.0) {
  cat("(OK)\n")
} else {
  cat("(WARNING: R 4.0+ recommended)\n")
}

# Check working directory
cat(sprintf("✓ Working directory: %s\n", getwd()))

# Centralized module directory (all code modules live here)
modules_dir_candidate <- file.path(getwd(), "Multi-frontier SBM modules")
if (dir.exists(modules_dir_candidate)) {
  modules_dir <- modules_dir_candidate
} else {
  # If subdirectory doesn't exist, use current directory as modules directory
  modules_dir <- getwd()
}
cat(sprintf("✓ Modules directory: %s\n", modules_dir))
Sys.setenv(MODULES_DIR = modules_dir)

# Check data file
data_path <- "./swine_farm_data.xlsx"
if (file.exists(data_path)) {
  cat(sprintf("✓ Data file found: %s\n", basename(data_path)))
  file_info <- file.info(data_path)
  cat(sprintf("  Size: %.2f MB, Modified: %s\n", 
              file_info$size / 1024^2, file_info$mtime))
} else {
  stop(sprintf("ERROR: Data file not found at: %s", data_path))
}

# Check critical packages
critical_packages <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2", "scales", "lpSolve", "quantreg")
missing_packages <- character(0)

cat("\n✓ Checking packages:\n")
for (pkg in critical_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  ✗ %s: NOT INSTALLED\n", pkg))
    missing_packages <- c(missing_packages, pkg)
  } else {
    cat(sprintf("  ✓ %s: OK\n", pkg))
  }
}

if (length(missing_packages) > 0) {
  cat(sprintf("\n⚠ Installing missing packages: %s\n", 
              paste(missing_packages, collapse = ", ")))
  for (pkg in missing_packages) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Global plot theme for Word-friendly readability
ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 17) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 16),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 15, face = "bold"),
      legend.text = ggplot2::element_text(size = 14),
      strip.text = ggplot2::element_text(size = 15, face = "bold")
    )
)

# Check required scripts
cat("\n✓ Checking analysis scripts:\n")
required_scripts <- c(
  "multi_frontier_sbm_dea.R",
  "meta_frontier_module.R",
  "improve_scale_analysis.R",
  "scale_efficiency_mechanism_analysis.R",
  "slack_analysis_module.R",
  "frontier_rank_comparison.R",
  "temporal_effect_robustness.R",
  "scale_expansion_did_analysis.R",
  "preexperiment_weight_ratio.R"
)

scripts_status <- data.frame(
  script = required_scripts,
  exists = sapply(file.path(modules_dir, required_scripts), file.exists),
  required = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
)

for (i in seq_len(nrow(scripts_status))) {
  status_icon <- ifelse(scripts_status$exists[i], "✓", "✗")
  req_text <- ifelse(scripts_status$required[i], "(Required)", "(Optional)")
  cat(sprintf("  %s %s %s\n", status_icon, scripts_status$script[i], req_text))
}

if (any(!scripts_status$exists & scripts_status$required)) {
  stop("ERROR: Required scripts missing. Analysis cannot proceed.")
}

cat("\n")
cat("==============================================================================\n")
cat("Pre-flight checks completed successfully!\n")
cat("==============================================================================\n\n")

# ============================================================================
# STEP 2: ANALYSIS OVERVIEW
# ============================================================================

cat("STEP 2: Analysis Overview\n")
cat("--------------------------\n")
cat("This streamlined analysis will:\n")
cat("  1. Load and validate data\n")
cat("  2. Run pre-experiment to find optimal weight multipliers (if preexperiment_weight_ratio.R exists)\n")
cat("  3. Calculate efficiency for 4 frontiers (M, L, G, A) using FIXED VARIABLE SET METHOD\n")
cat("     - All frontiers use same complete variable set (5 bad outputs)\n")
cat("     - Scientific standardization of variables (scale_matrix) to avoid scale bias\n")
cat("     - Fixed normalization denominator (m + s_b = 9) for consistent scale\n")
cat("     - All bad outputs constrained with different weights (core vs non-core)\n")
cat("     - Weight standardization: sum of weights = 5 for all frontiers (consistent total penalty)\n")
cat("     - Scheme 3: Non-overlapping core (M: Dead_pig, L: Waste_water+Manure, G: Carbon_emission+Eutrophication, A: All)\n")
cat("     - Weight multipliers optimized via pre-experiment (if enabled)\n")
cat("  4. Construct meta-frontier and calculate TGR\n")
cat("  5. Generate core figures (distributions, TGR, correlations)\n")
cat("  6. Detect and validate U-shaped scale-efficiency relationships\n")
cat("  7. Create publication-ready visualizations\n")
cat("  8. Perform statistical validation (Bootstrap, Jackknife)\n")
cat("  9. Save workspace for post-analysis use\n\n")
cat("Note: By default pre-experiment runs and weights from it are used.\n")
cat("      To skip pre-experiment and use preset weights, set: Sys.setenv(SKIP_PRE_EXPERIMENT = \"TRUE\")\n\n")

cat("Expected runtime: 5-15 minutes\n")
cat("(Focused on essentials)\n\n")

# Create results directory
results_base_dir <- "./results"
results_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- file.path(results_base_dir, sprintf("results_%s", results_timestamp))

cat(sprintf("Output directory: %s\n\n", results_dir))

if (!dir.exists(results_base_dir)) {
  dir.create(results_base_dir, recursive = TRUE)
}

dir.create(results_dir, recursive = TRUE)
Sys.setenv(RESULTS_DIR = results_dir)

cat("==============================================================================\n")
cat("  Maximum Envelope Meta-Frontier\n")
cat("==============================================================================\n")
main_script <- "multi_frontier_sbm_dea.R"
cat("\n✓ Using: Maximum Envelope setting\n")
cat("  Meta-frontier: max(TE_M, TE_L, TE_G, TE_A)\n")
cat("  Method: Maximum efficiency across constrained frontiers\n\n")

# ============================================================================
# STEP 3: RUN ANALYSIS
# ============================================================================

cat("==============================================================================\n")
cat("        STARTING ANALYSIS\n")
cat("==============================================================================\n\n")

cat("STEP 3: Running Main Analysis\n")
cat("-------------------------------\n")

# Ask user whether to run pre-experiment (default: run and use its weights)
cat("\nPre-experiment option:\n")
cat("  By default the pre-experiment runs and its weights are used (recommended).\n")
cat("  The pre-experiment (weight multiplier search) can take 10-30 minutes.\n")
cat("  You may skip it only if you have already run it and want to reuse preset weights.\n\n")
skip_response <- readline("Skip pre-experiment and use preset weights? (y/n, default=n): ")

if (tolower(substr(skip_response, 1, 1)) == "y") {
  Sys.setenv(SKIP_PRE_EXPERIMENT = "TRUE")
  cat("\n✓ Pre-experiment will be skipped. Using pre-determined weights:\n")
  cat("  M-Frontier: Dead_pig=4.091, others=0.227\n")
  cat("  L-Frontier: Waste_water/Manure core, Dead_pig/Carbon non-core\n")
  cat("  G-Frontier: Carbon_emission=4.091, others=0.227\n")
  cat("  A-Frontier: All=1.000\n\n")
} else {
  Sys.setenv(SKIP_PRE_EXPERIMENT = "FALSE")
  cat("\n✓ Pre-experiment will run (may take 10-30 minutes)\n\n")
}

main_script_path <- file.path(modules_dir, main_script)
if (!file.exists(main_script_path)) {
  stop(sprintf("ERROR: Main script not found: %s", main_script_path))
}

tryCatch({
  source(main_script_path, echo = FALSE)
  analysis_success <- TRUE
  
  # Optional: Run temporal effect control analysis if Report_year column exists
  if (exists("efficiency_results")) {
    cat("\n")
    cat("Checking for Report_year column...\n")
    cat(sprintf("  Column names: %s\n", paste(names(efficiency_results), collapse = ", ")))
    
    # Check for Report_year (case-insensitive)
    year_cols <- grep("report.*year|year", names(efficiency_results), ignore.case = TRUE, value = TRUE)
    if (length(year_cols) > 0) {
      cat(sprintf("  Found year-related columns: %s\n", paste(year_cols, collapse = ", ")))
      
      # Standardize to Report_year if needed
      if (!"Report_year" %in% names(efficiency_results)) {
        efficiency_results <- efficiency_results %>%
          dplyr::rename(Report_year = !!dplyr::sym(year_cols[1]))
        cat(sprintf("  Renamed '%s' to 'Report_year'\n", year_cols[1]))
      }
    }
    
    if ("Report_year" %in% names(efficiency_results)) {
      cat("  ✓ Report_year column found!\n")
      cat(sprintf("    Range: %s\n", paste(range(efficiency_results$Report_year, na.rm = TRUE), collapse = "-")))
      
      cat("\n")
      cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
      cat("RUNNING TEMPORAL EFFECT CONTROL ANALYSIS\n")
      cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")
      
      tryCatch({
        source(file.path(modules_dir, "temporal_effect_robustness.R"), echo = FALSE)
        cat("✓ Temporal effect control analysis completed\n")
      }, error = function(e2) {
        cat("⚠ Temporal effect analysis skipped (error: ", e2$message, ")\n")
      })
    } else {
      cat("  ✗ Report_year column not found in efficiency_results\n")
      cat("  ℹ Temporal effect analysis skipped\n")
      cat("  To enable: Ensure 'Report_year' column exists in your Excel file (Sheet 3)\n")
    }
  } else {
    cat("\nℹ efficiency_results not found in workspace\n")
  }
  
  # Run rank comparison analysis
  if (exists("efficiency_results")) {
    cat("\n")
    cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
    cat("RUNNING RANK COMPARISON ANALYSIS\n")
    cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")
    
    tryCatch({
      source(file.path(modules_dir, "frontier_rank_comparison.R"), echo = FALSE)
      cat("✓ Rank comparison analysis completed\n")
    }, error = function(e3) {
      cat("⚠ Rank comparison analysis skipped (error: ", e3$message, ")\n")
    })
  }
  
  # Run scale-up dynamic analysis (DID + bubble plots; uses same data sheet and latest efficiency results)
  if (exists("efficiency_results")) {
    cat("\nChecking conditions for scale-up dynamic analysis...\n")
    scaleup_script <- file.path(modules_dir, "scale_expansion_did_analysis.R")
    cat(sprintf("  Scaleup script exists: %s\n", file.exists(scaleup_script)))
    
    # Debug: Print column names to see what's available
    cat("  Available columns in efficiency_results:\n")
    cat(sprintf("    %s\n", paste(names(efficiency_results), collapse = ", ")))
    
    scale_up_col <- grep("Scale-up|Scale_up|ScaleUp", names(efficiency_results), value = TRUE)
    no_col <- grep("^No\\.?$|^Number$", names(efficiency_results), value = TRUE)
    
    cat(sprintf("  Scale-up columns found: %s\n", ifelse(length(scale_up_col) > 0, paste(scale_up_col, collapse = ", "), "None")))
    cat(sprintf("  No./Number columns found: %s\n", ifelse(length(no_col) > 0, paste(no_col, collapse = ", "), "None")))
    
    if (file.exists(scaleup_script) && length(scale_up_col) > 0 && length(no_col) > 0) {
      cat("\n")
      cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
      cat("RUNNING SCALE-UP DYNAMIC ANALYSIS (DID + BUBBLE PLOTS)\n")
      cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")
      tryCatch({
        source(scaleup_script, encoding = "UTF-8")
        cat("✓ Scale-up dynamic analysis completed\n")
      }, error = function(e4) {
        cat("⚠ Scale-up dynamic analysis skipped (error: ", e4$message, ")\n")
      })
    } else {
      cat("  ⚠ Scale-up dynamic analysis skipped: Missing required conditions\n")
    }
  } else {
    cat("  ⚠ Scale-up dynamic analysis skipped: efficiency_results not found\n")
  }
  
  # NOTE: Tobit threshold analysis removed - no thresholds detected
  # Tobit regression analysis is now integrated into scale_efficiency_mechanism_analysis.R
  
  # Create combined figure: VRS Efficiency by Scale + TGR by Scale (2-panel horizontal layout)
  if (exists("efficiency_results")) {
    cat("\n")
    cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
    cat("CREATING COMBINED FIGURE: VRS EFFICIENCY BY SCALE + TGR BY SCALE\n")
    cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")
    
    tryCatch({
      # Load required packages
      if (!require("patchwork", quietly = TRUE)) {
        install.packages("patchwork", repos = "https://cloud.r-project.org")
        library(patchwork)
      }
      if (!require("viridis", quietly = TRUE)) {
        install.packages("viridis", repos = "https://cloud.r-project.org")
        library(viridis)
      }
      
      # Get results directory
      results_base <- Sys.getenv("RESULTS_DIR", unset = "results")
      visualization_dir <- file.path(results_base, "visualization")
      
      # Check required columns
      required_cols <- c("market_pig", "M_TE_VRS", "L_TE_VRS", "G_TE_VRS", "A_TE_VRS")
      tgr_cols <- grep("^[MLGA]_TGR_VRS$", names(efficiency_results), value = TRUE)
      if (length(tgr_cols) == 0) {
        stop("Missing TGR columns in efficiency_results. Expected columns matching pattern: [MLGA]_TGR_VRS")
      }
      missing_cols <- setdiff(required_cols, names(efficiency_results))
      if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns in efficiency_results: %s", 
                     paste(missing_cols, collapse = ", ")))
      }
      
      # Prepare data: Create unified scale bins (<5k, 5-10k, 10-20k, 20-50k, >50k)
      scale_bin_breaks <- c(0, 5000, 10000, 20000, 50000, Inf)
      scale_bin_labels <- c("<5k", "5-10k", "10-20k", "20-50k", ">50k")
      
      # Debug: Check market_pig column
      cat("DEBUG: Checking market_pig column...\n")
      cat(sprintf("  market_pig exists: %s\n", "market_pig" %in% names(efficiency_results)))
      if ("market_pig" %in% names(efficiency_results)) {
        cat(sprintf("  market_pig length: %d\n", length(efficiency_results$market_pig)))
        cat(sprintf("  market_pig NA count: %d\n", sum(is.na(efficiency_results$market_pig))))
        cat(sprintf("  market_pig min: %.2f\n", min(efficiency_results$market_pig, na.rm = TRUE)))
        cat(sprintf("  market_pig max: %.2f\n", max(efficiency_results$market_pig, na.rm = TRUE)))
      }
      
      # Create scale bins using market_pig_original (head count) so all scale categories appear
      scale_col_fig <- if ("market_pig_original" %in% names(efficiency_results)) "market_pig_original" else "market_pig"
      eff_with_scale <- efficiency_results %>%
        dplyr::mutate(
          scale_bin = cut(.data[[scale_col_fig]],
                         breaks = scale_bin_breaks,
                         labels = scale_bin_labels,
                         include.lowest = TRUE)
        ) %>%
        dplyr::filter(!is.na(scale_bin)) %>%
        dplyr::mutate(scale_bin = factor(scale_bin,
                                         levels = scale_bin_labels,
                                         ordered = TRUE))
      cat(sprintf("  eff_with_scale created (using %s), rows: %d\n", scale_col_fig, nrow(eff_with_scale)))
      cat(sprintf("  Scale bins found: %s\n", paste(unique(eff_with_scale$scale_bin), collapse = ", ")))
      
      # Panel 1: VRS Efficiency by Scale
      eff_by_scale <- eff_with_scale %>%
        tidyr::pivot_longer(
          cols = c(M_TE_VRS, L_TE_VRS, G_TE_VRS, A_TE_VRS),
          names_to = "frontier",
          values_to = "efficiency"
        ) %>%
        dplyr::mutate(frontier = stringr::str_remove(frontier, "_TE_VRS")) %>%
        dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)) %>%
        dplyr::group_by(scale_bin, frontier) %>%
        dplyr::summarise(
          n = n(),
          mean_eff = mean(efficiency, na.rm = TRUE),
          se = sd(efficiency, na.rm = TRUE) / sqrt(n),
          .groups = "drop"
        )
      
      # Calculate label positions for sample sizes (total n per scale_bin)
      eff_by_scale_labels <- eff_with_scale %>%
        dplyr::group_by(scale_bin) %>%
        dplyr::summarise(n = n(), .groups = "drop") %>%
        dplyr::mutate(label_y = 1.15)

      # Unified frontier palette (muted, coordinated tones)
      frontier_colors <- c(
        "M" = "#7EB89E",
        "L" = "#E16859",
        "G" = "#8B72A1",
        "A" = "#60BCDA"
      )
      
      p_eff_by_scale <- eff_by_scale %>%
        ggplot(aes(x = scale_bin, y = mean_eff, fill = frontier)) +
        geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.75) +
        geom_line(aes(group = frontier, color = frontier),
                  position = position_dodge(width = 0.8),
                  linewidth = 1.2, alpha = 0.9) +
        geom_errorbar(aes(ymin = mean_eff - se, ymax = mean_eff + se),
                      position = position_dodge(width = 0.8),
                      width = 0.3, linewidth = 0.4, color = "gray20", alpha = 0.8) +
        geom_point(aes(color = frontier),
                   position = position_dodge(width = 0.8),
                   size = 1.5, shape = 21, fill = "white",
                   stroke = 0.8, alpha = 1.0) +
        # Sample size labels
        geom_text(data = eff_by_scale_labels,
                  aes(x = scale_bin, y = label_y, label = sprintf("n=%d", n)),
                  inherit.aes = FALSE,
                  size = 5.6, color = "black", fontface = "bold") +
        annotate("rect", xmin = 2.5, xmax = 4.5, ymin = 0, ymax = 1.2,
                 alpha = 0.1, fill = "red") +
        annotate("text", x = 3.5, y = 1.04,
                 label = "Efficiency Trough Zone",
                 size = 5.2, color = "red", fontface = "bold.italic") +
        scale_fill_manual(values = frontier_colors,
                          breaks = c("M", "L", "G", "A"),
                          labels = c("M" = "M-Frontier",
                                     "L" = "L-Frontier",
                                     "G" = "G-Frontier",
                                     "A" = "A-Frontier")) +
        scale_color_manual(values = frontier_colors,
                           breaks = c("M", "L", "G", "A"),
                           labels = c("M" = "M-Frontier",
                                      "L" = "L-Frontier",
                                      "G" = "G-Frontier",
                                      "A" = "A-Frontier")) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           breaks = seq(0, 1, 0.2),
                           limits = c(0, 1.2),
                           expand = c(0, 0)) +
        labs(
          title = "Mean VRS efficiency by scale category",
          x = "Scale category (market pig heads/a)",
          y = "Mean VRS technical efficiency",
          fill = "Frontier",
          color = "Frontier"
        ) +
        theme_minimal(base_size = 15) +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 14, face = "bold", color = "black"),
          axis.title = element_text(size = 14),
          axis.title.y = element_text(size = 14, face = "bold", color = "black")
        )
      
      # Panel 2: TGR by Scale
      tgr_cols_available <- grep("^[MLGA]_TGR_VRS$", names(eff_with_scale), value = TRUE)
      if (length(tgr_cols_available) == 0) {
        stop("No TGR columns found matching pattern ^[MLGA]_TGR_VRS$")
      }
      
      tgr_by_scale <- eff_with_scale %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(tgr_cols_available),
          names_to = "metric",
          values_to = "tgr"
        ) %>%
        tidyr::separate(metric, into = c("frontier", "dummy1", "dummy2"), sep = "_") %>%
        dplyr::select(-dummy1, -dummy2) %>%
        dplyr::filter(!is.na(tgr)) %>%
        dplyr::mutate(frontier = factor(frontier, levels = c("M", "L", "G", "A"), ordered = TRUE)) %>%
        dplyr::mutate(tgr = pmax(0, pmin(1.05, tgr))) %>%
        dplyr::group_by(scale_bin, frontier) %>%
        dplyr::summarise(
          n = n(),
          mean_tgr = mean(tgr, na.rm = TRUE),
          se_tgr = sd(tgr, na.rm = TRUE) / sqrt(n),
          .groups = "drop"
        )
      
      # Calculate label positions for sample sizes (total n per scale_bin)
      tgr_by_scale_labels <- eff_with_scale %>%
        dplyr::group_by(scale_bin) %>%
        dplyr::summarise(n = n(), .groups = "drop") %>%
        dplyr::mutate(label_y = 1.15)
      
      # Create TGR by scale panel
      p_tgr_by_scale <- tgr_by_scale %>%
        ggplot(aes(x = scale_bin, y = mean_tgr, fill = frontier)) +
        geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.75) +
        geom_errorbar(aes(ymin = mean_tgr - se_tgr, ymax = mean_tgr + se_tgr),
                      position = position_dodge(width = 0.8),
                      width = 0.3, linewidth = 0.4, color = "gray20", alpha = 0.8) +
        geom_line(aes(group = frontier, color = frontier),
                  position = position_dodge(width = 0.8),
                  linewidth = 1.2, alpha = 0.9) +
        geom_point(aes(color = frontier),
                   position = position_dodge(width = 0.8),
                   size = 1.5, shape = 21, fill = "white",
                   stroke = 0.8, alpha = 1.0) +
        # Sample size labels
        geom_text(data = tgr_by_scale_labels,
                  aes(x = scale_bin, y = label_y, label = sprintf("n=%d", n)),
                  inherit.aes = FALSE,
                  size = 5.6, color = "black", fontface = "bold") +
        geom_hline(yintercept = 1, linetype = "dashed",
                   color = "gray40", linewidth = 0.8) +
        scale_fill_manual(values = frontier_colors,
                          breaks = c("M", "L", "G", "A"),
                          labels = c("M" = "M-Frontier",
                                     "L" = "L-Frontier",
                                     "G" = "G-Frontier",
                                     "A" = "A-Frontier")) +
        scale_color_manual(values = frontier_colors,
                           breaks = c("M", "L", "G", "A"),
                           labels = c("M" = "M-Frontier",
                                      "L" = "L-Frontier",
                                      "G" = "G-Frontier",
                                      "A" = "A-Frontier")) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           breaks = seq(0, 1, 0.2),
                           limits = c(0, 1.2),
                           expand = c(0, 0),
                           oob = scales::squish) +
        labs(
          title = "TGR by scale category (VRS)",
          x = "Scale category (market pig heads/a)",
          y = "TGR-VRS",
          fill = "Frontier",
          color = "Frontier"
        ) +
        theme_minimal(base_size = 15) +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 14, face = "bold", color = "black"),
          axis.title = element_text(size = 14),
          axis.title.y = element_text(size = 14, face = "bold", color = "black")
        )
     
      # Combine plots vertically with shared legend at bottom
      p_eff_by_scale_no_legend <- p_eff_by_scale +
        guides(color = "none", fill = "none") +
        theme(legend.position = "none")
      
      p_tgr_by_scale_with_legend <- p_tgr_by_scale +
        guides(color = "none") +
        theme(legend.position = "bottom")
      
      # Create vertical layout version
      p_combined_vertical <- (p_eff_by_scale_no_legend / p_tgr_by_scale_with_legend) +
        plot_layout(guides = "collect", heights = c(1, 1)) +
        plot_annotation(
          title = "VRS efficiency and technology gap ratios by production scale",
          subtitle = "Top: mean VRS technical efficiency | Bottom: technology gap ratio (TGR-VRS)",
          theme = theme(
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
            plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40",
                                        margin = margin(b = 10))
          )
        ) &
        theme(legend.position = "bottom",
              legend.justification = "center",
              legend.box.just = "center",
              legend.box = "horizontal",
              legend.direction = "horizontal",
              legend.margin = margin(t = 10, b = 5),
              legend.title = element_text(face = "bold", size = 14, hjust = 0.5),
              legend.text = element_text(size = 13),
              legend.title.align = 0.5) &
        guides(fill = guide_legend(override.aes = list(linetype = 0, shape = NA),
                                   nrow = 1, byrow = TRUE),
               color = "none")
      
      # Save vertical layout version
      ggsave(file.path(visualization_dir, "Fig5_VRS_Efficiency_TGR_by_Scale_vertical.png"),
             p_combined_vertical, width = 8, height = 12, dpi = 300)
      
      cat("✓ Vertical combined figure saved: Fig5_VRS_Efficiency_TGR_by_Scale_vertical.png\n")
      cat("  Size: 8×12 inches, 2-panel vertical layout\n")
      cat("  Location: ", file.path(visualization_dir, "Fig5_VRS_Efficiency_TGR_by_Scale_vertical.png"), "\n\n")
      
    }, error = function(e4) {
      cat("⚠ Combined figure creation skipped (error: ", e4$message, ")\n")
    })
  }
  
}, error = function(e) {
  cat("\n")
  cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
  cat("ERROR: Analysis failed\n")
  cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
  msg <- conditionMessage(e)
  if (is.character(msg) && nzchar(msg)) {
    cat(sprintf("Error message: %s\n\n", msg))
  } else {
    cat("Error message: (none)\n")
    cat("Condition: ", paste(format(e), collapse = " "), "\n\n")
  }
  cat("Troubleshooting:\n")
  cat("  1. Check data file format\n")
  cat("  2. Verify all required scripts are present\n")
  cat("  3. Check package installations\n\n")
  analysis_success <<- FALSE
})

# ============================================================================
# STEP 4: POST-ANALYSIS SUMMARY
# ============================================================================

if (exists("analysis_success") && analysis_success) {
  cat("\n")
  cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
  cat("        ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")
  
  runtime <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("Total runtime: %.1f minutes\n\n", as.numeric(runtime)))
  
  # ============================================================================
  # FINAL: OPEN RESULTS
  # ============================================================================
  
  cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")
  response <- readline("Open results folder now? (y/n): ")
  
  if (tolower(substr(response, 1, 1)) == "y") {
    if (.Platform$OS.type == "windows") {
      shell.exec(results_dir)
    } else {
      system(paste("open", results_dir))
    }
    cat("\n")
  }
  
  # Quick view of core figures
  response2 <- readline("Open core figures now? (y/n): ")
  
  if (tolower(substr(response2, 1, 1)) == "y") {
    if (.Platform$OS.type == "windows") {
      core_paths <- c(
        file.path(results_dir, "visualization/Fig1_efficiency_distributions_VRS_CRS.png"),
        file.path(results_dir, "visualization/Fig2_TGR_comparison_VRS_CRS.png"),
        file.path(results_dir, "visualization/Fig3_correlation_matrices_VRS_CRS.png"),
        file.path(results_dir, "visualization/Fig4_slack_analysis.png"),
        file.path(results_dir, "visualization/Fig5_VRS_Efficiency_TGR_by_Scale.png"),
        file.path(results_dir, "visualization/efficiency_by_scale_BINNED_barplot.png")
      )
      
      for (path in core_paths) {
        if (file.exists(path)) {
          shell.exec(path)
          Sys.sleep(0.5)
        }
      }
    }
  }
  
} else {
  cat("\n")
  cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
  cat("        ANALYSIS INCOMPLETE\n")
  cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")
  cat("The analysis did not complete successfully.\n")
  cat("Please review error messages above.\n\n")
}

cat("\n")
cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
cat("SESSION COMPLETE\n")
cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
cat(sprintf("Total elapsed time: %.1f minutes\n", 
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat("\n")










