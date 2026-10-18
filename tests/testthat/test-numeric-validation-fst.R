# Numeric validation against an independent reference implementation
#
# WHY THIS FILE EXISTS
# ---------------------
# The audit of this package found no test comparing pgacmdr's C++
# re-implementation of the Weir & Cockerham (1984) estimators against an
# established, independently-reviewed implementation. For a scientific
# publication built around this tool, that comparison is the single most
# important piece of evidence for correctness - more important than any
# amount of code review. This file provides a first, minimal version of it.
#
# It does NOT hardcode any expected numeric value. Instead it computes the
# reference statistic at test time with the CRAN package `hierfstat`
# (Goudet, wc() function, the reference R implementation of Weir & Cockerham
# 1984) on the same toy dataset, and checks that pgacmdr's C++ output
# matches within a small numeric tolerance. This is honest validation: if
# `hierfstat` disagrees with pgacmdr, the test fails and that is exactly
# the point.
#
# HOW TO EXTEND
# -------------
# 1. Add larger / real datasets (e.g. from the package vignette) alongside
#    this toy example, ideally also cross-checked against `genepop` or
#    `adegenet::pairwise.fst()` for a second independent reference.
# 2. Do the same for FIS/FIT (this file only covers FST) and for the
#    bootstrap confidence intervals / permutation p-values, comparing
#    against `boot::boot()`-based CIs computed independently in R.
# 3. hierfstat is listed in DESCRIPTION Suggests. The test is skipped when it
#    is not installed; once your CI installs it, remove the
#    `skip_if_not_installed("hierfstat")` guard so that a missing reference
#    package fails loudly instead of silently skipping.

test_that("observed_wc84_stats_cpp() FST matches hierfstat::wc() on a toy dataset", {
  testthat::skip_if_not_installed("hierfstat")

  # ---- Toy dataset --------------------------------------------------------
  # 2 populations x 10 diploid individuals x 2 loci, one bi-allelic and one
  # tri-allelic locus, base = 100 (2-digit alleles per hierfstat convention:
  # genotype = allele1 * base + allele2). No missing data in this first pass.
  set.seed(1)
  base <- 100L
  n_per_pop <- 10L

  make_locus <- function(alleles, freqs_pop1, freqs_pop2) {
    g1 <- t(replicate(n_per_pop, sort(sample(alleles, 2, replace = TRUE, prob = freqs_pop1))))
    g2 <- t(replicate(n_per_pop, sort(sample(alleles, 2, replace = TRUE, prob = freqs_pop2))))
    rbind(g1, g2)[, 1] * base + rbind(g1, g2)[, 2]
  }

  locus1 <- make_locus(c(10L, 20L),       c(0.8, 0.2), c(0.2, 0.8))
  locus2 <- make_locus(c(10L, 20L, 30L),  c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))

  pop <- rep(1:2, each = n_per_pop)
  dat <- cbind(pop = pop, L1 = locus1, L2 = locus2)
  storage.mode(dat) <- "integer"

  # ---- Reference: hierfstat ------------------------------------------------
  hf_input <- data.frame(pop = pop, L1 = locus1, L2 = locus2)
  ref <- hierfstat::wc(hf_input, diploid = TRUE)
  ref_fst_global <- ref$FST

  # ---- pgacmdr C++ implementation --------------------------------------
  res <- pgacmdr::observed_wc84_stats_cpp(
    dat, pop_col_1based = 1L, missing_code = 0L, base = base
  )

  # Multi-locus (global) FST, Weir & Cockerham's ratio-of-sums estimator:
  # sum(a) / sum(a + b + c) across loci. pgacmdr already exposes this
  # directly as FST_overall_ratio_of_sums - use it as-is rather than
  # averaging per-locus FST values, which would NOT be the WC84 estimator.
  pgacmdr_fst_global <- res$FST_overall_ratio_of_sums

  expect_true(is.numeric(ref_fst_global))
  expect_true(!is.na(ref_fst_global))
  expect_true(is.numeric(pgacmdr_fst_global))
  expect_true(!is.na(pgacmdr_fst_global))

  # Loose tolerance on purpose for this first pass (toy sample size): the
  # goal at this stage is to catch gross implementation errors (wrong sign,
  # wrong denominator, factor-of-2 bugs), not to certify numeric precision.
  # Tighten this once the comparison is running routinely in CI.
  expect_equal(pgacmdr_fst_global, ref_fst_global, tolerance = 0.05)
})
