# Bootstrap over individuals for the Null Alleles module (closures of the module,
# extracted from the source file; skipped when the source tree is not available).

load_na_env <- function() {
  f <- testthat::test_path("..", "..", "R", "server_null_alleles.R")
  if (!file.exists(f)) return(NULL)
  ex  <- parse(f, keep.source = FALSE)
  find_body <- function(e) {
    if (is.call(e) && identical(e[[1]], as.name("moduleServer"))) return(e[[3]])
    if (is.call(e)) for (x in as.list(e)) { r <- find_body(x); if (!is.null(r)) return(r) }
    NULL
  }
  b <- NULL; for (e in ex) { b <- find_body(e); if (!is.null(b)) break }
  want <- c("em_freena", "weir_components_allele", "compute_fst_global_full",
            "compute_fst_pairwise", "boot_indiv_na", "summarise_boot_indiv", "%||%")
  env <- new.env(parent = asNamespace("pgacmdr"))
  for (s in as.list(b[[3]])[-1]) {
    if (is.call(s) && identical(s[[1]], as.name("<-")) && is.name(s[[2]]) && as.character(s[[2]]) %in% want)
      assign(as.character(s[[2]]), eval(s[[3]], env), env)
  }
  for (n in ls(env, all.names = TRUE)) if (is.function(env[[n]])) environment(env[[n]]) <- env
  env
}

sim_raw <- function(r = 4L, n = 25L, L = 3L, K = 6L, seed = 3L) {
  set.seed(seed)
  rows <- list()
  for (l in seq_len(L)) for (i in seq_len(r)) {
    p <- rgamma(K, 2); p <- p / sum(p)
    a <- matrix(sample.int(K, 2L * n, replace = TRUE, prob = p), n, 2L) + 100L
    rows[[length(rows) + 1L]] <- data.frame(Population = as.character(i), Individual = as.character(seq_len(n)),
                                            Marker = paste0("L", l), gt = a[, 1] * 1000L + a[, 2],
                                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

test_that("BS/Indiv for null alleles: shapes, reproducibility across cores, CI contains observed", {
  skip_if(is.null(env <- load_na_env()), "source tree not available")
  raw <- sim_raw(); markers <- unique(raw$Marker); pops <- as.character(1:4)
  treats <- setNames(rep("absent", length(markers)), markers)
  em <- list()
  for (l in markers) { em[[l]] <- list()
    for (p in pops) em[[l]][[p]] <- env$em_freena(raw$gt[raw$Marker == l & raw$Population == p], 1000L) }
  g  <- env$compute_fst_global_full(em); pw <- env$compute_fst_pairwise(em)
  b1 <- env$boot_indiv_na(raw, markers, pops, treats, 1000L, 40L, 11L, n_cores = 1L)
  expect_equal(length(b1$glob_ena), 40L)
  expect_equal(dim(b1$loc_raw), c(40L, 3L))
  expect_equal(dim(b1$pair_ena), c(40L, 6L))
  if (.Platform$OS.type != "windows") {
    b2 <- env$boot_indiv_na(raw, markers, pops, treats, 1000L, 40L, 11L, n_cores = 2L)
    expect_identical(b1$glob_ena, b2$glob_ena)
  }
  s <- env$summarise_boot_indiv(b1, g, pw, 0.05)
  expect_equal(s$global$Locus[1], "GLOBAL_MULTILOCUS")
  expect_true(all(s$global$CI_lo_ENA_indiv <= s$global$CI_hi_ENA_indiv, na.rm = TRUE))
  expect_true(s$global$CI_lo_ENA_indiv[1] <= g$global_ena && g$global_ena <= s$global$CI_hi_ENA_indiv[1])
  expect_equal(nrow(s$pair), 6L)
})
