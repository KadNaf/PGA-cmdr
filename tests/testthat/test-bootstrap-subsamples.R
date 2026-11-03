# Bootstrap over sub-samples (populations resampled as blocks).
#
# Resampling populations with replacement duplicates some of them; a duplicated
# population is counted as a distinct unit identical to another one, so the
# among-population variance of a replicate is too small by (r-1)/r. Without a
# correction the replicates of FST (and of HT) sit well below the observed value
# (on an island-model simulation, 8 populations of 50 individuals, only ~0.7 % of
# the 5000 replicates exceeded the observed FST and the 95 % interval excluded
# the estimate). These tests check that the replicates are centred on the
# observed value, and that the observed values themselves are not modified.

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

test_that("sub-sample bootstrap of FST and HT is centred on the observed value", {
  dat <- sim_island()
  b <- boot_popblock_wc84_parallel(dat, 1L, 0L, 100L, 2000L, 1L, 1)
  fst <- b$FST_overall_boot; fst_obs <- b$FST_overall_obs
  expect_lt(abs(mean(fst) / fst_obs - 1), 0.05)
  expect_gt(mean(fst > fst_obs), 0.35)
  expect_lt(mean(fst > fst_obs), 0.65)
  q <- stats::quantile(fst, c(0.025, 0.975))
  expect_true(q[[1]] <= fst_obs && fst_obs <= q[[2]])

  ht <- b$HT_overall_boot; ht_obs <- b$HT_overall_obs
  expect_lt(abs(mean(ht) - ht_obs), 0.25 * stats::sd(ht))
  q <- stats::quantile(ht, c(0.025, 0.975))
  expect_true(q[[1]] <= ht_obs && ht_obs <= q[[2]])

  # per-locus replicates are centred too
  for (l in seq_len(ncol(b$FST_boot))) {
    expect_lt(abs(mean(b$FST_boot[, l]) - b$FST_obs[l]), 0.35 * stats::sd(b$FST_boot[, l]))
  }
})

test_that("second sub-sample kernel (FST/FIT/HI/HS) is centred on the observed value", {
  dat <- sim_island()
  g <- boot_wc84_stats_popblock_cpp(dat, 1L, 0L, 100L, 2000L, 0.95, 1L, 1L)
  expect_lt(abs(mean(g$FST_boot) / mean(g$FST_obs_boot) - 1), 0.05)
  expect_lt(abs(mean(g$FIT_boot) / mean(g$FIT_obs_boot) - 1), 0.05)
})

test_that("observed values do not depend on the sub-sample correction", {
  dat <- sim_island(L = 3L)
  o <- observed_wc84_stats_cpp(dat, 1L, 0L, 100L)
  b <- boot_popblock_wc84_parallel(dat, 1L, 0L, 100L, 10L, 1L, 1)
  expect_equal(unname(b$FST_obs), unname(o$FST), tolerance = 1e-12)
  expect_equal(b$FST_overall_obs, o$FST_overall_ratio_of_sums, tolerance = 1e-12)
})
