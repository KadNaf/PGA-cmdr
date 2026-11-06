# Bootstrap over sub-samples (populations resampled as blocks).
#
# Resampling populations with replacement duplicates some of them, which lowers
# the among-population variance of a replicate: the replicates of FST, FIT and HT
# sit below the observed value (HS is not affected). The published bias shift
# (observed - mean of the replicates; Efron & Tibshirani 1993) is applied in R by
# pga_recentre(). These tests check the bias, the re-centring, and that the
# observed values themselves are not modified.

sim_island <- function(r = 8L, n = 50L, L = 10L, K = 20L, alpha = 82, seed = 42L) {
  set.seed(seed)
  out <- matrix(0L, r * n, L + 1L)
  out[, 1] <- rep(seq_len(r), each = n)
  for (l in seq_len(L)) {
    p0 <- rgamma(K, 3); p0 <- p0 / sum(p0)
    for (i in seq_len(r)) {
      pi <- rgamma(K, alpha * p0 + 1e-3); pi <- pi / sum(pi)
      a <- matrix(sample.int(K, 2L * n, replace = TRUE, prob = pi), n, 2L)
      out[(i - 1L) * n + seq_len(n), l + 1L] <- a[, 1] * 100L + a[, 2]
    }
  }
  colnames(out) <- c("pop", paste0("L", seq_len(L)))
  out
}

test_that("BS/Ss replicates of FST and HT are lowered, and the bias shift re-centres them", {
  dat <- sim_island()
  b <- boot_popblock_wc84_parallel(dat, 1L, 0L, 100L, 2000L, 1L, 1)
  fst <- b$FST_overall_boot; fst_obs <- b$FST_overall_obs
  ht  <- b$HT_overall_boot;  ht_obs  <- b$HT_overall_obs
  hs  <- b$HS_overall_boot;  hs_obs  <- b$HS_overall_obs

  # raw replicates: a sub-sample drawn twice lowers the among-sub-sample variance
  expect_lt(mean(fst), fst_obs)
  expect_lt(mean(ht), ht_obs)
  # HS does not depend on that variance
  expect_lt(abs(mean(hs) - hs_obs), 0.25 * stats::sd(hs))

  # bias shift (observed - mean of the replicates): centred, width unchanged
  fst_s <- pga_recentre(fst, fst_obs)
  expect_equal(mean(fst_s), fst_obs, tolerance = 1e-10)
  expect_equal(stats::sd(fst_s), stats::sd(fst), tolerance = 1e-10)
  q <- stats::quantile(fst_s, c(0.025, 0.975))
  expect_true(q[[1]] <= fst_obs && fst_obs <= q[[2]])
  ht_s <- pga_recentre(ht, ht_obs)
  expect_equal(mean(ht_s), ht_obs, tolerance = 1e-10)

  # per-locus matrix, column by column
  m <- pga_recentre(b$FST_boot, b$FST_obs)
  expect_equal(unname(colMeans(m)), unname(b$FST_obs), tolerance = 1e-10)
})

test_that("pga_recentre leaves non-finite input alone", {
  expect_identical(pga_recentre(c(NA_real_, NA_real_), 0.1), c(NA_real_, NA_real_))
  expect_equal(pga_recentre(c(1, 2, 3), NA_real_), c(1, 2, 3))
})

test_that("observed values do not depend on the sub-sample correction", {
  dat <- sim_island(L = 3L)
  o <- observed_wc84_stats_cpp(dat, 1L, 0L, 100L)
  b <- boot_popblock_wc84_parallel(dat, 1L, 0L, 100L, 10L, 1L, 1)
  expect_equal(unname(b$FST_obs), unname(o$FST), tolerance = 1e-12)
  expect_equal(b$FST_overall_obs, o$FST_overall_ratio_of_sums, tolerance = 1e-12)
})
