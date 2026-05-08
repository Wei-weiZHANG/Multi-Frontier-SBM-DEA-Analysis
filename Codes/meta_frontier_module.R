# Meta-Frontier Construction Module
# ==================================
#
# This module implements TRUE meta-frontier construction following
# O'Donnell et al. (2008) with proper handling of heterogeneous bad outputs
# across different frontier specifications.
#
# CHALLENGE:
# ----------
# Each frontier uses different bad output sets (ORIGINAL METHOD):
#   M: {Dead_pig}
#   L: {Waste_water, Manure, Dead_pig}
#   G: {Global indicators}
#   A: {All 5 bad outputs}
#
# SOLUTION (FIXED VARIABLE SET METHOD):
# --------------------------------------
# Use UNIVERSAL bad output set (union of all) with conditional constraints:
#   - Union = {Waste_water, Manure, Dead_pig, Carbon_emission}
#   - All frontiers use the SAME complete variable set
#   - For each frontier, only relevant bad outputs are CONSTRAINED (penalized)
#   - Others are "free disposable" (appear in constraints but not penalized)
#   - This eliminates systematic bias from variable count differences
#
# ==================================

#' Construct true meta-frontier using MAXIMUM envelope approach (O'Donnell et al. 2008)
#' 
#' CORE IMPLEMENTATION: This function actually executes the meta-frontier calculation
#' by taking the maximum efficiency across all group frontiers for each DMU.
#' 
#' @param efficiency_results Dataframe with group frontier efficiencies (M_TE_VRS, L_TE_VRS, etc.)
#' @param RTS "vrs" or "crs" (default: "vrs")
#' @return efficiency_results with Meta_TE columns added
#' @details
#' O'Donnell (2008) requires: Meta_TE = max(TE_M, TE_L, TE_G, TE_A) for each DMU
#' This function implements this calculation directly, not just as a description.
construct_meta_frontier <- function(efficiency_results, RTS = "vrs") {
  
  message("\n")
  message(paste(rep("=", 78), collapse = ""))
  message("META-FRONTIER CONSTRUCTION (O'Donnell et al. 2008)")
  message(paste(rep("=", 78), collapse = ""))
  message("")
  message("Method: MAXIMUM Efficiency Envelope")
  message("  Formula: Meta_TE = max(TE_M, TE_L, TE_G, TE_A) for each DMU")
  message("  Rationale: Meta-frontier = BEST technology across all constraint scenarios\n")
  
  # Check required columns
  required_cols <- paste0(c("M", "L", "G", "A"), "_TE_", toupper(RTS))
  missing_cols <- setdiff(required_cols, names(efficiency_results))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  te_matrix <- as.matrix(efficiency_results[, required_cols])
  meta_te <- apply(te_matrix, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  
  # Store result
  meta_col <- paste0("Meta_TE_", toupper(RTS))
  efficiency_results[[meta_col]] <- meta_te
  
  # Report statistics
  message(sprintf("  Meta-frontier %s: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]",
                  toupper(RTS),
                  mean(meta_te, na.rm=TRUE), sd(meta_te, na.rm=TRUE),
                  min(meta_te, na.rm=TRUE), max(meta_te, na.rm=TRUE)))
  
  message("  ✓ Meta-frontier calculated as maximum envelope\n")
  
  return(efficiency_results)
}

#' Calculate meta-frontier efficiency using MAXIMUM envelope approach
#' 
#' @param efficiency_results Dataframe with all frontier efficiencies
#' @return efficiency_results with Meta_TE and TGR columns added
calculate_meta_frontier_maximum <- function(efficiency_results) {
  
  message("STEP 1: Calculating Meta-Frontier Efficiency (Maximum Envelope)")
  message("---------------------------------------------------------------")
  message("  Formula: Meta_TE = max(TE_M, TE_L, TE_G, TE_A) for each DMU")
  message("  Meta-frontier = BEST technology across all constraint scenarios\n")
  
  # Check if required columns exist
  required_cols <- c("M_TE_VRS", "L_TE_VRS", "G_TE_VRS", "A_TE_VRS",
                    "M_TE_CRS", "L_TE_CRS", "G_TE_CRS", "A_TE_CRS")
  missing_cols <- setdiff(required_cols, names(efficiency_results))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  te_matrix_vrs <- cbind(
    E = efficiency_results$M_TE_VRS,
    L = efficiency_results$L_TE_VRS,
    G = efficiency_results$G_TE_VRS,
    A = efficiency_results$A_TE_VRS
  )
  
  meta_te_vrs <- apply(te_matrix_vrs, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  
  # CRS Meta-Frontier (same approach)
  te_matrix_crs <- cbind(
    E = efficiency_results$M_TE_CRS,
    L = efficiency_results$L_TE_CRS,
    G = efficiency_results$G_TE_CRS,
    A = efficiency_results$A_TE_CRS
  )
  
  meta_te_crs <- apply(te_matrix_crs, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  
  relative_diff <- abs(meta_te_vrs - efficiency_results$A_TE_VRS) / 
                    pmax(abs(meta_te_vrs), abs(efficiency_results$A_TE_VRS), 1e-10)
  n_meta_eq_a_vrs <- sum(relative_diff < 1e-6, na.rm = TRUE)
  n_total <- sum(!is.na(meta_te_vrs))
  pct_from_a <- n_meta_eq_a_vrs / n_total * 100
  
  message(sprintf("  DEBUG: Meta-frontier = A-Frontier for %d/%d farms (%.1f%%)",
                  n_meta_eq_a_vrs, n_total, pct_from_a))
  message(sprintf("  DEBUG: Meta-frontier from OTHER frontiers: %.1f%%", 100 - pct_from_a))
  
  # Store in results
  efficiency_results$Meta_TE_VRS <- meta_te_vrs
  efficiency_results$Meta_TE_CRS <- meta_te_crs
  
  # Report statistics
  message(sprintf("  Meta-frontier VRS: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]",
                  mean(meta_te_vrs, na.rm=TRUE), sd(meta_te_vrs, na.rm=TRUE),
                  min(meta_te_vrs, na.rm=TRUE), max(meta_te_vrs, na.rm=TRUE)))
  message(sprintf("  Meta-frontier CRS: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]",
                  mean(meta_te_crs, na.rm=TRUE), sd(meta_te_crs, na.rm=TRUE),
                  min(meta_te_crs, na.rm=TRUE), max(meta_te_crs, na.rm=TRUE)))
  
  # Compare with individual frontiers to verify maximum
  message("\n  Verification (Meta should be ≥ all frontiers):")
  message(sprintf("    Mean Meta_TE_VRS (%.4f) vs Mean A_TE_VRS (%.4f): %s",
                  mean(meta_te_vrs, na.rm=TRUE),
                  mean(efficiency_results$A_TE_VRS, na.rm=TRUE),
                  ifelse(mean(meta_te_vrs, na.rm=TRUE) >= mean(efficiency_results$A_TE_VRS, na.rm=TRUE), 
                         "✓ Meta ≥ A", "✗ ERROR")))
  message(sprintf("    Mean Meta_TE_VRS (%.4f) vs Mean M_TE_VRS (%.4f): %s",
                  mean(meta_te_vrs, na.rm=TRUE),
                  mean(efficiency_results$M_TE_VRS, na.rm=TRUE),
                  ifelse(mean(meta_te_vrs, na.rm=TRUE) >= mean(efficiency_results$M_TE_VRS, na.rm=TRUE), 
                         "✓ Meta ≥ E", "✗ ERROR")))
  
  # Identify which frontier provides meta-frontier for each DMU
  eff_matrix <- cbind(
    E = efficiency_results$M_TE_VRS,
    L = efficiency_results$L_TE_VRS,
    G = efficiency_results$G_TE_VRS,
    A = efficiency_results$A_TE_VRS
  )
  
  # Vectorized source identification (more efficient for large samples)
  # Use which.max with proper NA handling
  frontier_source_vrs <- apply(eff_matrix, 1, function(x) {
    if (all(is.na(x))) return(NA_integer_)
    valid_x <- x[!is.na(x)]
    if (length(valid_x) == 0) return(NA_integer_)
    # Find index of maximum (handles ties by taking first)
    max_idx <- which.max(x)
    return(max_idx)
  })
  
  frontier_names <- c("M", "L", "G", "A")
  source_counts <- table(frontier_names[frontier_source_vrs], useNA = "ifany")
  
  message("\n  Meta-frontier sources (which frontier provides BEST technology for each farm):")
  for (fname in names(source_counts)) {
    if (!is.na(fname)) {
      message(sprintf("    %s-Frontier: %d farms (%.1f%%)", 
                      fname, source_counts[fname],
                      source_counts[fname]/sum(source_counts, na.rm=TRUE)*100))
    }
  }
  
  # Additional insight
  message("\n  Interpretation:")
  total_count <- sum(source_counts, na.rm=TRUE)
  if (total_count > 0) {
    if ("M" %in% names(source_counts) && source_counts["M"] / total_count > 0.5) {
      message("  → M-Frontier dominant: Most farms perform best under minimal constraints")
    }
    if ("A" %in% names(source_counts) && source_counts["A"] / total_count > 0.2) {
      message("  → Significant A-Frontier presence: Some farms excel even with all constraints")
    }
  }
  
  efficiency_results$Meta_Source <- frontier_names[frontier_source_vrs]
  
  message("\n  INTERPRETATION:")
  message("  - Meta-frontier = UPPER envelope (best possible technology)")
  message("  - For each farm, meta-frontier = MAXIMUM efficiency across scenarios")
  message("  - Represents 'best achievable with optimal technology choice'")
  message("  - Meta_TE ≥ TE_i for all frontiers i (meta-frontier dominates)\n")
  
  return(efficiency_results)
}

#' Calculate Technology Gap Ratios relative to true meta-frontier (VECTORIZED)
#' 
#' @param efficiency_results Dataframe with Meta_TE columns
#' @return efficiency_results with TGR columns added
calculate_tgr_true_meta <- function(efficiency_results) {
  
  message("STEP 2: Calculating Technology Gap Ratios (TGR)")
  message("---------------------------------------------------------------")
  message("  Formula: TGR_i = TE_i / Meta_TE (vectorized calculation)")
  message("  Where Meta_TE = max(TE_M, TE_L, TE_G, TE_A) for each farm")
  message("  Property: TGR ≤ 1 (group TE ≤ meta-frontier TE)\n")
  
  # Vectorized TGR calculation with robust zero handling
  # Use relative threshold: if meta < 1e-6, treat as numerical zero
  meta_vrs_raw <- efficiency_results$Meta_TE_VRS
  meta_crs_raw <- efficiency_results$Meta_TE_CRS
  
  meta_vrs_is_zero <- !is.na(meta_vrs_raw) & meta_vrs_raw < 1e-6
  
  if (any(meta_vrs_is_zero, na.rm = TRUE)) {
    warning(sprintf("  ⚠ %d DMUs have Meta_TE_VRS < 1e-6 (numerical zero), TGR may be unreliable",
                    sum(meta_vrs_is_zero, na.rm = TRUE)))
  }
  
  # Safe division: use 1e-8 for numerical stability, but flag true zeros
  meta_vrs_safe <- pmax(meta_vrs_raw, 1e-8)
  meta_crs_safe <- pmax(meta_crs_raw, 1e-8)
  
  for (frontier_id in c("M", "L", "G", "A")) {
    te_vrs_col <- paste0(frontier_id, "_TE_VRS")
    te_crs_col <- paste0(frontier_id, "_TE_CRS")
    
    if (!(te_vrs_col %in% names(efficiency_results))) {
      warning(sprintf("Column %s not found, skipping", te_vrs_col), call. = FALSE)
      next
    }
    
    # Vectorized TGR calculation with proper outlier detection
    te_vrs_raw <- efficiency_results[[te_vrs_col]]
    te_crs_raw <- efficiency_results[[te_crs_col]]
    
    # Calculate TGR
    tgr_vrs <- te_vrs_raw / meta_vrs_safe
    tgr_crs <- te_crs_raw / meta_crs_safe
    
    tgr_vrs_outliers <- sum(tgr_vrs > 1.01, na.rm = TRUE)
    
    if (tgr_vrs_outliers > 0) {
      warning(sprintf("  ⚠ %d TGR_VRS values > 1.01 detected (likely calculation errors)", 
                      tgr_vrs_outliers))
    }
    
    # Clip to [0, 1.01] to allow small numerical errors but flag larger ones
    tgr_vrs <- pmax(0, pmin(1.01, tgr_vrs))
    tgr_crs <- pmax(0, pmin(1.01, tgr_crs))
    
    # Store
    efficiency_results[[paste0(frontier_id, "_TGR_VRS")]] <- tgr_vrs
    efficiency_results[[paste0(frontier_id, "_TGR_CRS")]] <- tgr_crs
    
    # Statistics
    message(sprintf("%s-Frontier TGR:", frontier_id))
    message(sprintf("  VRS: mean=%.4f, median=%.4f, sd=%.4f, range=[%.4f, %.4f]",
                    mean(tgr_vrs, na.rm=TRUE), median(tgr_vrs, na.rm=TRUE),
                    sd(tgr_vrs, na.rm=TRUE),
                    min(tgr_vrs, na.rm=TRUE), max(tgr_vrs, na.rm=TRUE)))
    message(sprintf("  CRS: mean=%.4f, median=%.4f, sd=%.4f, range=[%.4f, %.4f]",
                    mean(tgr_crs, na.rm=TRUE), median(tgr_crs, na.rm=TRUE),
                    sd(tgr_crs, na.rm=TRUE),
                    min(tgr_crs, na.rm=TRUE), max(tgr_crs, na.rm=TRUE)))
    
    # Interpretation (removed misleading "often highest/lowest" language)
    mean_tgr <- mean(tgr_vrs, na.rm=TRUE)
    if (mean_tgr > 0.95) {
      message("  → Small technology gap (close to best possible)")
    } else if (mean_tgr > 0.85) {
      message("  → Moderate technology gap")
    } else {
      message("  → Large technology gap (significant improvement potential)")
    }
    message("")
  }
  
  # Summary statistics
  message("AVERAGE TGR COMPARISON:")
  message(sprintf("  M: %.4f | L: %.4f | G: %.4f | A: %.4f",
                  mean(efficiency_results$M_TGR_VRS, na.rm=TRUE),
                  mean(efficiency_results$L_TGR_VRS, na.rm=TRUE),
                  mean(efficiency_results$G_TGR_VRS, na.rm=TRUE),
                  mean(efficiency_results$A_TGR_VRS, na.rm=TRUE)))
  message("  Note: No universal ordering - each farm has TGR=1 on their best frontier\n")
  
  return(efficiency_results)
}

#' Validate meta-frontier construction
#' 
#' @param efficiency_results Dataframe with all efficiency scores
#' @param bad_weights_list Optional: List of weight vectors by frontier (for weight validation)
validate_meta_frontier <- function(efficiency_results, bad_weights_list = NULL) {
  
  message("VALIDATION: Meta-Frontier Properties")
  message("---------------------------------------------------------------")
  
  meta_vrs <- efficiency_results$Meta_TE_VRS
  e_vrs <- efficiency_results$M_TE_VRS
  l_vrs <- efficiency_results$L_TE_VRS
  g_vrs <- efficiency_results$G_TE_VRS
  a_vrs <- efficiency_results$A_TE_VRS
  
  # =========================================================================
  # CHECK 1: Meta = max(M, L, G, A) for each DMU
  # =========================================================================
  message("\n  CHECK 1: Meta-frontier equals maximum efficiency")
  
  # Calculate what meta SHOULD be (vectorized)
  te_matrix_check <- cbind(E = e_vrs, L = l_vrs, G = g_vrs, A = a_vrs)
  meta_should_be <- apply(te_matrix_check, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  
  # Check if actual matches expected (using relative tolerance for numerical precision)
  # Relative tolerance: abs(a-b) / max(|a|, |b|) < tolerance
  relative_diff <- abs(meta_vrs - meta_should_be) / 
                    pmax(abs(meta_vrs), abs(meta_should_be), 1e-10)
  matches <- sum(relative_diff < 1e-6, na.rm = TRUE)  # Relative tolerance: 1e-6
  total <- sum(!is.na(meta_vrs))
  match_pct <- matches / total * 100
  
  message(sprintf("    Meta = max(M,L,G,A): %d/%d farms (%.1f%%)",
                  matches, total, match_pct))
  
  if (match_pct > 99) {
    message("    ✓ PASS: Meta correctly calculated as maximum")
  } else {
    message("    ✗ FAIL: Meta NOT calculated as maximum")
    message(sprintf("      Only %.1f%% match expected values", match_pct))
  }
  
  # =========================================================================
  # CHECK 2: Meta ≥ all group frontiers (upper envelope property)
  # =========================================================================
  message("\n  CHECK 2: Upper envelope property (Meta ≥ all frontiers)")
  
  # Upper envelope check with relative tolerance
  # Meta should be >= all frontiers (allowing small numerical errors)
  check_m <- all((meta_vrs >= e_vrs) | (abs(meta_vrs - e_vrs) / pmax(abs(meta_vrs), abs(e_vrs), 1e-10) < 1e-6), na.rm = TRUE)
  check_l <- all((meta_vrs >= l_vrs) | (abs(meta_vrs - l_vrs) / pmax(abs(meta_vrs), abs(l_vrs), 1e-10) < 1e-6), na.rm = TRUE)
  check_g <- all((meta_vrs >= g_vrs) | (abs(meta_vrs - g_vrs) / pmax(abs(meta_vrs), abs(g_vrs), 1e-10) < 1e-6), na.rm = TRUE)
  check_a <- all((meta_vrs >= a_vrs) | (abs(meta_vrs - a_vrs) / pmax(abs(meta_vrs), abs(a_vrs), 1e-10) < 1e-6), na.rm = TRUE)
  
  message(sprintf("    Meta ≥ M-Frontier: %s", ifelse(check_m, "✓ PASS", "✗ FAIL")))
  message(sprintf("    Meta ≥ L-Frontier: %s", ifelse(check_l, "✓ PASS", "✗ FAIL")))
  message(sprintf("    Meta ≥ G-Frontier: %s", ifelse(check_g, "✓ PASS", "✗ FAIL")))
  message(sprintf("    Meta ≥ A-Frontier: %s", ifelse(check_a, "✓ PASS", "✗ FAIL")))
  
  # =========================================================================
  # CHECK 3: All TGR ∈ [0,1]
  # =========================================================================
  message("\n  CHECK 3: TGR bounds")
  
  TGR_M <- efficiency_results$M_TGR_VRS
  tgr_l <- efficiency_results$L_TGR_VRS
  tgr_g <- efficiency_results$G_TGR_VRS
  tgr_a <- efficiency_results$A_TGR_VRS
  
  check_m_bounds <- all(TGR_M >= 0 & TGR_M <= 1.001, na.rm = TRUE)
  check_l_bounds <- all(tgr_l >= 0 & tgr_l <= 1.001, na.rm = TRUE)
  check_g_bounds <- all(tgr_g >= 0 & tgr_g <= 1.001, na.rm = TRUE)
  check_a_bounds <- all(tgr_a >= 0 & tgr_a <= 1.001, na.rm = TRUE)
  
  message(sprintf("    M_TGR ∈ [0,1]: %s", ifelse(check_m_bounds, "✓ PASS", "✗ FAIL")))
  message(sprintf("    L_TGR ∈ [0,1]: %s", ifelse(check_l_bounds, "✓ PASS", "✗ FAIL")))
  message(sprintf("    G_TGR ∈ [0,1]: %s", ifelse(check_g_bounds, "✓ PASS", "✗ FAIL")))
  message(sprintf("    A_TGR ∈ [0,1]: %s", ifelse(check_a_bounds, "✓ PASS", "✗ FAIL")))
  
  # =========================================================================
  # CHECK 4: Each DMU has at least one TGR ≈ 1.0
  # =========================================================================
  message("\n  CHECK 4: Each DMU has at least one TGR ≈ 1.0")
  message("    (The frontier providing meta should have TGR=1 for that DMU)")
  
  # For each DMU, find maximum TGR across all frontiers (vectorized)
  tgr_matrix <- cbind(TGR_M, tgr_l, tgr_g, tgr_a)
  max_tgr_per_dmu <- apply(tgr_matrix, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  
  # Count how many DMUs have max TGR ≈ 1.0
  # Use adaptive threshold: 0.99 for strict check, but also report 0.995 for very strict
  # O'Donnell (2008) requires TGR=1 for DMU on meta-frontier, but numerical precision matters
  n_with_tgr_1_strict <- sum(max_tgr_per_dmu > 0.995, na.rm = TRUE)  # Very strict: > 0.995
  n_with_tgr_1 <- sum(max_tgr_per_dmu > 0.99, na.rm = TRUE)  # Standard: > 0.99
  n_total_dmus <- sum(!is.na(max_tgr_per_dmu))
  pct_with_tgr_1 <- n_with_tgr_1 / n_total_dmus * 100
  
  pct_with_tgr_1_strict <- n_with_tgr_1_strict / n_total_dmus * 100
  message(sprintf("    DMUs with max(TGR) > 0.995 (very strict): %d/%d (%.1f%%)",
                  n_with_tgr_1_strict, n_total_dmus, pct_with_tgr_1_strict))
  message(sprintf("    DMUs with max(TGR) > 0.99 (standard): %d/%d (%.1f%%)",
                  n_with_tgr_1, n_total_dmus, pct_with_tgr_1))
  
  if (pct_with_tgr_1 > 95) {
    message("    ✓ PASS: Most DMUs have TGR≈1 on their meta-source frontier")
    if (pct_with_tgr_1_strict < 80) {
      message("    ⚠ NOTE: Lower strict threshold suggests some numerical precision issues")
    }
  } else {
    message("    ✗ FAIL: Too few DMUs have TGR≈1")
    message(sprintf("      Only %.1f%% have max TGR > 0.99", pct_with_tgr_1))
  }
  
  # =========================================================================
  # CHECK 5: Meta-frontier source distribution
  # =========================================================================
  message("\n  CHECK 5: Meta-frontier source heterogeneity")
  
  if ("Meta_Source" %in% names(efficiency_results)) {
    source_dist <- table(efficiency_results$Meta_Source, useNA = "ifany")
    total_sources <- sum(source_dist, na.rm = TRUE)
    
    message("    Meta-frontier sources (which frontier provides BEST technology):")
    for (fname in names(source_dist)) {
      if (!is.na(fname)) {
        pct <- source_dist[fname] / total_sources * 100
        message(sprintf("      %s-Frontier: %d farms (%.1f%%)", fname, source_dist[fname], pct))
      }
    }
    
    # Check for unhealthy patterns
    if (length(source_dist) == 1) {
      message("    ✗ WARNING: Meta-frontier comes from ONLY ONE frontier!")
      message("      This suggests meta may not be properly constructed")
    } else if (max(source_dist) / total_sources > 0.95) {
      dominant <- names(which.max(source_dist))
      message(sprintf("    ✗ WARNING: Meta dominated by %s-Frontier (>95%%)", dominant))
      message("      Expected more heterogeneity across frontiers")
    } else {
      message("    ✓ PASS: Meta-frontier shows healthy heterogeneity")
    }
  } else {
    message("    ⚠ WARNING: Meta_Source column not found")
  }
  
  # =========================================================================
  # CHECK 6: Weight standardization validation (for fixed variable set method)
  # =========================================================================
  message("\n  CHECK 6: Weight standardization (Fixed Variable Set Method)")
  message("    Verifying that weight differences reflect constraint focus, not total penalty strength")
  
  # This check requires weight information
  # Check if weights are provided as parameter or in parent environment
  weight_check_passed <- TRUE
  weight_check_msg <- ""
  
  if (is.null(bad_weights_list)) {
    # Try to get from parent environment
    if (exists("bad_weights_list", envir = parent.frame())) {
      bad_weights_list <- get("bad_weights_list", envir = parent.frame())
    }
  }
  
  if (!is.null(bad_weights_list) && length(bad_weights_list) > 0) {
    
    # Verify weight sums are standardized (should all equal number of bad outputs)
    weight_sums <- sapply(bad_weights_list, sum)
    expected_sum <- length(bad_weights_list[[1]])  # Number of bad outputs
    
    if (all(abs(weight_sums - expected_sum) < 0.01)) {
      message(sprintf("    ✓ PASS: All frontier weight sums = %.2f (standardized)", expected_sum))
      message("      → Weight differences reflect constraint focus, not total penalty")
    } else {
      weight_check_passed <- FALSE
      weight_check_msg <- "Weight sums not standardized - differences may reflect total penalty, not constraint focus"
      message("    ✗ FAIL: Weight sums not standardized")
      for (fname in names(weight_sums)) {
        message(sprintf("      %s-Frontier: sum = %.3f (expected: %.2f)", 
                        fname, weight_sums[fname], expected_sum))
      }
    }
    
    # Verify weight ranges are within acceptable bounds [0.7, 1.6]
    all_weights <- unlist(bad_weights_list)
    if (all(all_weights >= 0.7 & all_weights <= 1.6)) {
      message(sprintf("    ✓ PASS: All weights in acceptable range [0.7, 1.6]"))
    } else {
      weight_check_passed <- FALSE
      if (weight_check_msg != "") weight_check_msg <- paste(weight_check_msg, "; ", sep = "")
      weight_check_msg <- paste(weight_check_msg, "Some weights outside [0.7, 1.6] range", sep = "")
      message(sprintf("    ✗ FAIL: Some weights outside [0.7, 1.6]: min=%.2f, max=%.2f",
                      min(all_weights), max(all_weights)))
    }
  } else {
    message("    ⚠ WARNING: Weight information not available for validation")
    message("      → Cannot verify that differences reflect constraint focus vs total penalty")
  }
  
  # =========================================================================
  # CHECK 7: TGR=1 matches meta source
  # =========================================================================
  message("\n  CHECK 7: TGR=1 consistency with meta source")
  
  if ("Meta_Source" %in% names(efficiency_results)) {
    # For each DMU, check if TGR≈1 on their meta-source frontier
    consistency_check <- 0
    
    for (i in seq_len(nrow(efficiency_results))) {
      source_frontier <- efficiency_results$Meta_Source[i]
      if (is.na(source_frontier)) next
      
      tgr_col <- paste0(source_frontier, "_TGR_VRS")
      if (tgr_col %in% names(efficiency_results)) {
        tgr_value <- efficiency_results[[tgr_col]][i]
        if (!is.na(tgr_value) && tgr_value > 0.99) {
          consistency_check <- consistency_check + 1
        }
      }
    }
    
    consistency_pct <- consistency_check / nrow(efficiency_results) * 100
    message(sprintf("    Farms where TGR≈1 on meta-source frontier: %d/%d (%.1f%%)",
                    consistency_check, nrow(efficiency_results), consistency_pct))
    
    if (consistency_pct > 95) {
      message("    ✓ PASS: Strong consistency between meta source and TGR=1")
    } else {
      message("    ✗ FAIL: Poor consistency")
      message(sprintf("      Only %.1f%% show TGR=1 on their meta-source", consistency_pct))
    }
  }
  
  # =========================================================================
  # OVERALL VALIDATION SUMMARY
  # =========================================================================
  message("\n" , paste(rep("-", 78), collapse = ""))
  
  # Check if all critical validations passed
  all_checks_pass <- (match_pct > 99) &&  # Meta = max check
                     check_m && check_l && check_g && check_a &&  # Upper envelope checks
                     check_e_bounds && check_l_bounds && check_g_bounds && check_a_bounds &&  # TGR bounds
                     (pct_with_tgr_1 > 95) &&  # Each DMU has TGR=1
                     weight_check_passed  # Weight standardization check
  
  if (all_checks_pass) {
    message("  ✓✓✓ META-FRONTIER VALIDATION: ALL CHECKS PASSED ✓✓✓")
    message("  Meta-frontier is correctly constructed as maximum envelope!")
  } else {
    message("  ✗✗✗ META-FRONTIER VALIDATION: SOME CHECKS FAILED ✗✗✗")
    message("  Review the checks above to identify issues")
    
    # Provide specific failure info
    if (match_pct <= 99) message("    → CHECK 1 failed: Meta ≠ max(M,L,G,A)")
    if (!all(check_m, check_l, check_g, check_a)) message("    → CHECK 2 failed: Upper envelope violated")
    if (!all(check_e_bounds, check_l_bounds, check_g_bounds, check_a_bounds)) message("    → CHECK 3 failed: TGR bounds violated")
    if (pct_with_tgr_1 <= 95) message("    → CHECK 4 failed: Too few DMUs with TGR=1")
    if (!weight_check_passed) message(sprintf("    → CHECK 6 failed: %s", weight_check_msg))
  }
  
  message("")
  
  invisible(all_checks_pass)
}

# ==============================================================================
# ALTERNATIVE APPROACH: Super-Union DEA
# ==============================================================================
#
# For advanced users who want to try pooling with variable bad outputs:
# 
# This function creates a "super-DEA" where:
# - All possible bad outputs are included in the model
# - For each DMU, inactive bad outputs are set to a very small value
# - This allows single DEA run with all DMUs
# 
# NOTE: This is experimental and may produce unexpected results
# ==============================================================================

#' Experimental: Meta-frontier via super-union DEA
#' 
#' WARNING: This is an experimental approach. Requires sbm_efficiency() function
#' to be available in the calling environment (from multi_frontier_sbm_dea.R).
#' 
#' @param X_scaled Scaled input matrix
#' @param Y_scaled Scaled good output matrix
#' @param df_clean Original data
#' @param bad_outputs_list List of bad outputs by frontier
#' @param all_bad_outputs Union of all bad outputs
#' @param bad_weights_list List of weight vectors by frontier (for fixed variable set method)
#' @param constrained_indices_list List of constrained bad output indices by frontier
#' @param RTS "vrs" or "crs"
#' @param sbm_efficiency_func Function to calculate SBM efficiency (must be provided)
#' @return Meta-frontier efficiency scores
meta_frontier_super_union <- function(X_scaled, Y_scaled, df_clean, 
                                      bad_outputs_list, all_bad_outputs,
                                      bad_weights_list = NULL,
                                      constrained_indices_list = NULL,
                                      RTS = "vrs",
                                      sbm_efficiency_func = NULL) {
  
  message("EXPERIMENTAL: Super-Union Meta-Frontier Construction")
  message("  Creating augmented bad output matrix with conditional activation\n")
  
  # Check if sbm_efficiency function is available
  if (is.null(sbm_efficiency_func)) {
    # Try to find it in parent environment
    if (exists("sbm_efficiency", envir = parent.frame())) {
      sbm_efficiency_func <- get("sbm_efficiency", envir = parent.frame())
    } else {
      stop("ERROR: sbm_efficiency() function not found. ",
           "This function must be defined in the calling environment ",
           "(typically from multi_frontier_sbm_dea.R)")
    }
  }
  
  n <- nrow(df_clean)
  
  # Create full bad output matrix (all bad outputs for all farms - fixed variable set)
  Z_full <- as.matrix(df_clean[, all_bad_outputs, drop = FALSE])
  
  # Scale bad outputs
  scale_matrix <- function(M) {
    scaling_factors <- apply(M, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
    sweep(M, 2, scaling_factors, "/")
  }
  Z_scaled_full <- scale_matrix(Z_full)
  
  # Strategy: Calculate efficiency for each frontier scenario
  # Then take envelope (MAXIMUM) as meta-frontier
  meta_scores <- matrix(NA, nrow = n, ncol = length(bad_outputs_list))
  frontier_names <- names(bad_outputs_list)
  
  for (i in seq_along(bad_outputs_list)) {
    frontier_name <- frontier_names[i]
    message(sprintf("  Calculating %s-Frontier efficiency...", frontier_name))
    
    # Get weights and constrained indices for this frontier
    if (!is.null(bad_weights_list) && frontier_name %in% names(bad_weights_list)) {
      weights <- bad_weights_list[[frontier_name]]
    } else {
      weights <- NULL
    }
    
    if (!is.null(constrained_indices_list) && frontier_name %in% names(constrained_indices_list)) {
      constrained_indices <- constrained_indices_list[[frontier_name]]
    } else {
      # Fallback: use all bad outputs if indices not provided
      constrained_indices <- seq_along(all_bad_outputs)
    }
    
    # Calculate efficiency under this constraint scenario
    # Use fixed variable set: all bad outputs, but with frontier-specific weights
    scores <- sbm_efficiency_func(X_scaled, Y_scaled, 
                                 bad = Z_scaled_full,
                                 constrained_bad_indices = constrained_indices,
                                 bad_weights = weights,
                                 RTS = RTS)
    meta_scores[, i] <- scores
  }
  
  # Meta-frontier = MAXIMUM across all scenarios (best possible, upper envelope)
  meta_te <- apply(meta_scores, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  
  message(sprintf("  ✓ Meta-frontier constructed from %d frontier scenarios", 
                  length(bad_outputs_list)))
  message(sprintf("  ✓ Mean meta-frontier efficiency: %.4f", mean(meta_te, na.rm=TRUE)))
  
  return(meta_te)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Scale matrix by median absolute value
scale_matrix <- function(M) {
  scaling_factors <- apply(M, 2, function(x) max(1, median(abs(x), na.rm = TRUE)))
  sweep(M, 2, scaling_factors, "/")
}

#' Generate meta-frontier analysis report
#' 
#' @param efficiency_results Dataframe with efficiency and TGR scores
#' @param output_path Path to save report
generate_meta_frontier_report <- function(efficiency_results, output_path) {
  
  report_lines <- c(
    paste(rep("=", 78), collapse = ""),
    "META-FRONTIER ANALYSIS REPORT",
    paste(rep("=", 78), collapse = ""),
    "",
    sprintf("Generated: %s", Sys.time()),
    sprintf("Sample Size: %d farms", nrow(efficiency_results)),
    "",
    "CONSTRUCTION METHOD:",
    "-------------------",
    "Maximum Envelope Approach (Correct O'Donnell et al. 2008)",
    "  Meta_TE = max(TE_M, TE_L, TE_G, TE_A)",
    "  Represents the BEST achievable performance across all technologies",
    "  Serves as UPPER envelope of all group-specific frontiers",
    "  Meta-frontier ≥ all group frontiers (by construction)",
    "",
    "META-FRONTIER EFFICIENCY:",
    "-------------------------"
  )
  
  # Meta-frontier statistics
  meta_vrs <- efficiency_results$Meta_TE_VRS
  meta_crs <- efficiency_results$Meta_TE_CRS
  
  report_lines <- c(report_lines,
    sprintf("VRS: mean=%.4f, median=%.4f, sd=%.4f",
            mean(meta_vrs, na.rm=TRUE), median(meta_vrs, na.rm=TRUE), 
            sd(meta_vrs, na.rm=TRUE)),
    sprintf("CRS: mean=%.4f, median=%.4f, sd=%.4f",
            mean(meta_crs, na.rm=TRUE), median(meta_crs, na.rm=TRUE), 
            sd(meta_crs, na.rm=TRUE)),
    "",
    "TECHNOLOGY GAP RATIOS:",
    "----------------------"
  )
  
  # TGR statistics
  for (fid in c("M", "L", "G", "A")) {
    tgr_vrs <- efficiency_results[[paste0(fid, "_TGR_VRS")]]
    tgr_crs <- efficiency_results[[paste0(fid, "_TGR_CRS")]]
    
    report_lines <- c(report_lines,
      sprintf("\n%s-Frontier:", fid),
      sprintf("  VRS TGR: mean=%.4f, median=%.4f, range=[%.4f, %.4f]",
              mean(tgr_vrs, na.rm=TRUE), median(tgr_vrs, na.rm=TRUE),
              min(tgr_vrs, na.rm=TRUE), max(tgr_vrs, na.rm=TRUE)),
      sprintf("  CRS TGR: mean=%.4f, median=%.4f, range=[%.4f, %.4f]",
              mean(tgr_crs, na.rm=TRUE), median(tgr_crs, na.rm=TRUE),
              min(tgr_crs, na.rm=TRUE), max(tgr_crs, na.rm=TRUE))
    )
  }
  
  # Meta-frontier source composition
  if ("Meta_Source" %in% names(efficiency_results)) {
    source_counts <- table(efficiency_results$Meta_Source)
    report_lines <- c(report_lines,
      "",
      "META-FRONTIER COMPOSITION:",
      "--------------------------",
      "(Which frontier provides BEST technology for each farm?)"
    )
    
    for (fname in names(source_counts)) {
      report_lines <- c(report_lines,
        sprintf("  %s-Frontier: %d farms (%.1f%%) achieve best performance here",
                fname, source_counts[fname],
                source_counts[fname]/sum(source_counts)*100)
      )
    }
    
    report_lines <- c(report_lines,
      "",
      "INTERPRETATION:",
      "  - Meta-frontier = best technology available across all constraint scenarios",
      "  - Each farm's meta comes from their highest-performing frontier",
      "  - Heterogeneous distribution shows farms excel under different regimes"
    )
  }
  
  report_lines <- c(report_lines,
    "",
    paste(rep("=", 78), collapse = ""),
    "END OF META-FRONTIER REPORT",
    paste(rep("=", 78), collapse = ""),
    ""
  )
  
  writeLines(report_lines, output_path)
  cat(paste(report_lines, collapse = "\n"))
  
  invisible(TRUE)
}

# ==============================================================================
# MODULE EXPORTS
# ==============================================================================

message("✓ Meta-frontier construction module loaded")
message("  Functions available:")
message("    - calculate_meta_frontier_maximum()  [MAIN FUNCTION]")
message("    - calculate_tgr_true_meta()")
message("    - validate_meta_frontier()")
message("    - generate_meta_frontier_report()")
message("    - meta_frontier_super_union() [experimental]")
message("")












