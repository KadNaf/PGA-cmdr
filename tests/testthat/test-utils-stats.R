# Tests of the pure statistical helpers in R/utils_stats.R.
#
# Reference values: GENEPOP output for the default Boophilus dataset
# (inst/extdata/default_dataset.csv), as compared by the project owner:
# per-sample Ho, Hs and Fis (WC), averaged over loci (GENEPOP weights each locus by
# the number of individuals typed at it).

default_matrix <- function() {
  f <- system.file("extdata", "default_dataset.csv", package = "pgacmdr")
  if (!nzchar(f)) f <- file.path("..", "..", "inst", "extdata", "default_dataset.csv")
  testthat::skip_if_not(file.exists(f), "default dataset not available")
  d <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  loci <- c("B12", "C07", "D12", "D10", "A12", "C03")
  a <- as.matrix(d[, 9:20]); mode(a) <- "numeric"
  g <- sapply(seq_along(loci), function(j) {
    x <- a[, 2 * j - 1]; y <- a[, 2 * j]
    ifelse(x > 0 & y > 0, x * 1000L + y, 0L)
  })
  colnames(g) <- loci
  lv  <- sort(unique(d$Locality))
  mat <- cbind(Pop = match(d$Locality, lv), g)
  storage.mode(mat) <- "integer"
  attr(mat, "pop_levels") <- lv
  mat
}

genepop_ref <- data.frame(
  Population = c("Boulouparis", "Bourail", "Canala", "Gadji", "LaFoa", "Poquereux", "PortLaguerre", "Sarramea"),
  Ho  = c(0.7562, 0.6637, 0.6727, 0.6728, 0.6386, 0.6772, 0.6888, 0.6000),
  Hs  = c(0.6990, 0.7083, 0.7050, 0.6835, 0.6893, 0.7149, 0.7069, 0.7199),
  Fis = c(-0.081, 0.062, 0.045, 0.023, 0.070, 0.053, 0.026, 0.165),
  stringsAsFactors = FALSE
)

test_that("per-population Ho and Hs reproduce GENEPOP (loci weighted by n typed)", {
  mat <- default_matrix()
  w <- pga_pop_weighted_means(pga_pop_locus_stats(mat, 1000L))
  expect_equal(w$Population, genepop_ref$Population)
  expect_equal(round(w$Ho, 4), genepop_ref$Ho)
  expect_equal(round(w$Hs, 4), genepop_ref$Hs)
})

test_that("an unweighted mean over loci does NOT reproduce GENEPOP (guards the weighting)", {
  mat  <- default_matrix()
  long <- pga_pop_locus_stats(mat, 1000L)
  unw  <- tapply(long$Ho, factor(long$Population, levels = genepop_ref$Population), mean)
  expect_gt(max(abs(round(unw, 4) - genepop_ref$Ho)), 0.001)
})

test_that("overall_by_population table: Fis (WC) matches GENEPOP, detail has Ho only", {
  mat <- default_matrix()
  r <- pga_overall_by_population(mat, 1000L)
  expect_equal(round(r$overall$`Fis (WC)`, 3), genepop_ref$Fis)
  expect_named(r$detail, c("Population", "Locus", "Ho"))
  expect_equal(unique(r$detail$Locus), c("B12", "C07", "D12", "D10", "A12", "C03"))   # data order
})

test_that("the C++ kernel and nei_het_stats_cpp() give the same Hs (Nei & Chesser 1983)", {
  mat <- default_matrix()
  nei <- nei_het_stats_cpp(mat, 1L, 0L, 1000L)
  obs <- observed_wc84_stats_cpp(mat, 1L, 0L, 1000L)
  expect_equal(as.numeric(obs$HS), as.numeric(nei$Hs), tolerance = 1e-12)
  expect_equal(round(as.numeric(nei$Hs)[1], 5), 0.68164)   # locus B12, as in FSTAT/GENEPOP
})

test_that("pga_boot_over_loci(): observed values and sane intervals", {
  mat  <- default_matrix()
  comp <- wc84_locus_components_cpp(mat, 1L, 0L, 1000L)
  r <- pga_boot_over_loci(comp$A, comp$B, comp$C, comp$HS, comp$HT, n_boot = 2000L, seed = 1)
  expect_equal(r$table$Statistic, c("FST", "FIT", "FIS", "HS", "HT"))
  expect_equal(round(r$table$Observed, 4), c(0.0156, 0.0592, 0.0442, 0.7042, 0.7141))
  expect_true(all(r$table$CI_L < r$table$Observed & r$table$Observed < r$table$CI_U))
  expect_equal(dim(r$reps), c(2000L, 5L))
  # reproducible
  r2 <- pga_boot_over_loci(comp$A, comp$B, comp$C, comp$HS, comp$HT, n_boot = 2000L, seed = 1)
  expect_identical(r$reps, r2$reps)
})

test_that("pga_ci_bias_shift() centres a biased bootstrap distribution on the observed value", {
  set.seed(2)
  reps <- rnorm(5000, mean = 0.69, sd = 0.004)        # bootstrap mean 0.01 below the observed value
  plain <- pga_ci(reps)
  shifted <- pga_ci_bias_shift(0.70, reps)
  expect_lt(plain[["hi"]], 0.70)                       # plain percentile CI excludes the observed value
  expect_true(shifted[["lo"]] < 0.70 && 0.70 < shifted[["hi"]])
  expect_equal(unname(shifted[["hi"]] - shifted[["lo"]]), unname(plain[["hi"]] - plain[["lo"]]))
})

test_that("pga_hs_per_population(): observed = weighted Hs, intervals contain it", {
  mat <- default_matrix()
  h <- pga_hs_per_population(mat, 1000L, n_boot = 1000L, seed = 3)
  expect_equal(round(h$table$Observed_HS, 4), genepop_ref$Hs)
  expect_true(all(h$table$CI_L_indiv < h$table$Observed_HS & h$table$Observed_HS < h$table$CI_U_indiv))
  expect_true(all(h$table$CI_L_loci  < h$table$CI_U_loci))
  expect_equal(dim(h$reps_indiv), c(1000L, 8L))
})

test_that("pga_fis_pop_loci_boot(): observed FIS per population equals wc_fis_by_pop_wc84()", {
  mat <- default_matrix()
  f <- pga_fis_pop_loci_boot(mat, 1000L, n_boot = 1000L, seed = 4)
  ref <- wc_fis_by_pop_wc84(mat, 0L, 1000L)
  expect_equal(f$table$Observed_FIS[1:8], as.numeric(ref), tolerance = 1e-10)
  expect_equal(f$table$Population[9], "Overall")
})

test_that("table builders give the expected columns and Overall rows", {
  ft <- data.frame(ID = c("L1", "L2", "Overall"), Observed_FST = c(.01, .02, .015),
                   Boot_Mean = 0, Boot_Median = 0, P_value = c(.001, .01, .002),
                   CI_L = c(0, 0, .005), CI_U = c(.02, .04, .03))
  lb <- data.frame(Statistic = c("FST", "FIT", "FIS", "HS", "HT"), Observed = 1:5 / 10, Boot_Mean = 0,
                   SE = .01, CI_L = 0, CI_U = 1)
  s1 <- pga_fst_section(ft, lb)
  expect_equal(s1$ID, c("L1", "L2", "Overall (BS/Ss)", "Overall (BS/Loci)"))
  expect_named(s1, c("ID", "Observed_FST", "P_value", "CI_L", "CI_U"))
  so <- data.frame(Statistic = c("FIS", "FST", "FIT", "HS", "HT"), CI_L_subs = 0, CI_U_subs = 1)
  ml <- pga_multilocus_section(lb, so)
  expect_equal(ml$Statistic, c("FIS", "FST", "FIT", "HS", "HT"))                      # FIS, FST, FIT order
  expect_named(ml, c("Statistic", "Observed", "SE_loci", "CI_L_loci", "CI_U_loci", "CI_L_subs", "CI_U_subs"))
})


test_that("tie-aware counting: values equal up to rounding noise are counted as ties", {
  null <- c(0.5, 1, 1 + 1e-13, 1 - 1e-13, 2)
  expect_equal(pga_n_ge(null, 1), 4)      # 1, 1+1e-13, 1-1e-13 (tie) and 2
  expect_equal(pga_n_le(null, 1), 4)      # 0.5, and the three ties
  expect_equal(pga_n_gt(null, 1), 1)      # only 2 is strictly greater
  expect_equal(pga_n_abs_ge(c(-1, 1 + 1e-13, 0.2), -1), 2)
})

test_that("one-sided FIS p-values of a single sample match FSTAT's (ties counted in both tails)", {
  # C07 in Sarramea (24 individuals): FSTAT gives 0.4274 (deficit) and 0.7430 (excess);
  # without the tie tolerance the excess p-value used to come out around 0.59.
  mat <- default_matrix()
  p <- match("Sarramea", attr(mat, "pop_levels")); j <- match("C07", colnames(mat)[-1])
  sub <- mat[mat[, 1] == p, , drop = FALSE]; sub[, 1] <- 1L
  obs  <- fis_wc_cpp(sub, 1000L)$FIS[j]
  set.seed(1)
  perm <- batch_permute_wc_fis(dat = sub, pop_col_1based = 1L, base = 1000L, B = 6000L)[, j]
  def <- (pga_n_ge(perm, obs) + 1) / (length(perm) + 1)
  exc <- (pga_n_le(perm, obs) + 1) / (length(perm) + 1)
  expect_equal(def, 0.4274, tolerance = 0.03)
  expect_equal(exc, 0.7430, tolerance = 0.03)
  expect_gt(def + exc, 1.1)                # ties are counted on both sides
})

test_that("sub-sample bootstrap is reproducible: same result whatever the number of threads", {
  mat <- default_matrix()
  a <- boot_wc84_stats_popblock_cpp(mat, 1L, 0L, 1000L, 500L, 0.95, 1L, 1L)
  b <- boot_wc84_stats_popblock_cpp(mat, 1L, 0L, 1000L, 500L, 0.95, 1L, 2L)
  c <- boot_wc84_stats_popblock_cpp(mat, 1L, 0L, 1000L, 500L, 0.95, 1L, 2L)
  expect_identical(a$FIT_boot, b$FIT_boot)    # used to differ: the stream depended on the thread number
  expect_identical(b$FIT_boot, c$FIT_boot)
})
