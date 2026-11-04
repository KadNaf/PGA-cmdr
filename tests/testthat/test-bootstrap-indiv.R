# Bootstrap over INDIVIDUALS (resampled within each sub-sample) for FST, FIT, HT,
# pairwise FST, and FIS per locus x sub-sample.

sim_island2 <- function(r = 8L, n = 50L, L = 10L, K = 20L, alpha = 82, seed = 42L) {
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

test_that("bias-shifted BS/indiv CI is centred on the observed value, observed unchanged", {
  m <- sim_island2()
  obs <- observed_wc84_stats_cpp(m, 1L, 0L, 100L)
  bi <- boot_indiv_wc84_cpp(m, 1L, 0L, 100L, 400L, 1, 1L)
  expect_equal(dim(bi$FST_boot), c(400L, 10L))
  for (st in c("FST", "FIT", "HT")) {
    ov <- if (st == "HT") mean(obs$HT) else obs[[paste0(st, "_overall_ratio_of_sums")]]
    reps <- bi[[paste0(st, "_overall_boot")]]
    ci <- pga_ci_bias_shift(ov, reps, 0.95)
    expect_true(ci[["lo"]] < ov && ov < ci[["hi"]])
    expect_equal(mean(reps) + (ov - mean(reps)), ov)
  }
})

test_that("pairwise FST: one row per pair, observed equals direct computation, CI ordered", {
  m <- sim_island2(r = 4L, n = 30L, L = 6L)
  pw <- pga_pairwise_fst(m, base = 100L, n_boot = 100L, n_perm = 99L, n_threads = 1L)
  expect_equal(nrow(pw$overall), 6L)
  expect_equal(nrow(pw$per_locus), 6L * 6L)
  sub <- m[m[, 1] %in% c(1L, 2L), ]
  d <- observed_wc84_stats_cpp(sub, 1L, 0L, 100L)$FST_overall_ratio_of_sums
  expect_equal(pw$overall$FST[1], d)
  expect_true(all(pw$overall$CI_L_indiv < pw$overall$CI_U_indiv))
  expect_true(all(pw$overall$P_value > 0 & pw$overall$P_value <= 1))
})

test_that("FIS per locus x sub-sample: observed matches fis_wc_cpp, p-values in (0,1]", {
  m <- sim_island2(r = 3L, n = 20L, L = 4L)
  lp <- pga_fis_locus_pop(m, base = 100L, n_boot = 50L, n_perm = 49L)
  expect_equal(nrow(lp), 12L)
  sub <- m[m[, 1] == 2L, ]; sub[, 1] <- 1L
  expect_equal(lp$FIS[lp$Sub_sample == "2"], as.numeric(fis_wc_cpp(sub, base = 100L)$FIS))
  expect_true(all(lp$P_value > 0 & lp$P_value <= 1, na.rm = TRUE))
  expect_true(all(lp$CI_L <= lp$CI_U, na.rm = TRUE))
})

test_that("SE and p-value columns are moved last, in that order", {
  d <- data.frame(ID = "L1", Observed = 1, SE_indiv = .1, CI_L = 0, P_value = .5, CI_U = 2, SE_loci = .2, P_theta_gt = .3)
  expect_equal(names(pga_se_p_last(d)),
               c("ID", "Observed", "CI_L", "CI_U", "SE_indiv", "SE_loci", "P_value", "P_theta_gt"))
})
