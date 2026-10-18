# ============================================================================
# utils_stats.R — pure statistical helpers (no Shiny, no database)
#
# Everything in this file works on plain R objects so that it can be unit-tested
# outside the application (see tests/testthat/test-utils_stats.R). The Shiny
# modules (server_general_stats.R) only call these functions and format the
# output.
#
# Conventions
#   mat   : integer genotype matrix as built by hf_mat_r(): column 1 = population
#           code (1..K), the other columns = loci, genotype = allele1 * base +
#           allele2, 0 (or NA) = missing. attr(mat, "pop_levels") = population
#           names (position = population code).
#   base  : packing base (1000 by default).
# ============================================================================


# Percentile confidence interval (type-7 quantiles, the R default) of the finite
# values of `x`. Returns c(lo, hi), NA when fewer than 2 finite values.
pga_ci <- function(x, conf_level = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(lo = NA_real_, hi = NA_real_))
  a <- (1 - conf_level) / 2
  q <- stats::quantile(x, probs = c(a, 1 - a), type = 7, names = FALSE)
  c(lo = q[1], hi = q[2])
}


# Bias-shifted percentile CI, for statistics whose bootstrap distribution is
# systematically displaced from the observed value.
#
# Why this exists: resampling INDIVIDUALS with replacement creates duplicated
# genotypes, so every resample is slightly less diverse than the original sample.
# For gene diversity (HS) the bootstrap mean is therefore lower than the observed
# value by about HS/(2n) (n = sample size), and a plain percentile CI ends up
# hugging, or even excluding, the observed value (e.g. upper bound = observed HS).
# The CI is shifted by (observed - bootstrap mean) so that it is centred like the
# estimator; its width (the sampling variability) is unchanged.
pga_ci_bias_shift <- function(obs, reps, conf_level = 0.95) {
  reps <- reps[is.finite(reps)]
  if (!is.finite(obs) || length(reps) < 2L) return(c(lo = NA_real_, hi = NA_real_))
  ci <- pga_ci(reps, conf_level)
  shift <- obs - mean(reps)
  c(lo = ci[["lo"]] + shift, hi = ci[["hi"]] + shift)
}


# Bootstrap over LOCI: loci are resampled with replacement and the multilocus
# Weir & Cockerham (1984) ratio-of-sums statistics are recomputed.
#   A, B, C : per-locus WC84 variance components (among populations, among
#             individuals within populations, within individuals)
#   HS, HT  : per-locus gene diversities (their multilocus value is the mean)
# Returns list(table = data.frame(Statistic, Observed, SE, CI_L, CI_U),
#              reps  = matrix n_boot x 5 with columns FST, FIT, FIS, HS, HT).
# SE is the standard deviation of the bootstrap replicates (bootstrap standard
# error).
pga_boot_over_loci <- function(A, B, C, HS, HT, n_boot = 1000L, conf_level = 0.95, seed = 1) {
  nm <- c("FST", "FIT", "FIS", "HS", "HT")
  ok <- is.finite(A) & is.finite(B) & is.finite(C) & is.finite(HS) & is.finite(HT)
  A <- A[ok]; B <- B[ok]; C <- C[ok]; HS <- HS[ok]; HT <- HT[ok]
  L <- length(A)
  empty <- list(
    table = data.frame(Statistic = nm, Observed = NA_real_, SE = NA_real_,
                       CI_L = NA_real_, CI_U = NA_real_, stringsAsFactors = FALSE),
    reps  = matrix(NA_real_, nrow = n_boot, ncol = 5L, dimnames = list(NULL, nm))
  )
  if (L < 2L) return(empty)

  abc <- sum(A + B + C); bc <- sum(B + C)
  obs <- c(FST = sum(A) / abc, FIT = sum(A + B) / abc, FIS = sum(B) / bc,
           HS = mean(HS), HT = mean(HT))

  set.seed(as.integer(seed %% .Machine$integer.max))
  n_boot <- as.integer(n_boot)
  idx <- matrix(sample.int(L, L * n_boot, replace = TRUE), nrow = n_boot, ncol = L)
  rs  <- function(v) rowSums(matrix(v[idx], nrow = n_boot, ncol = L))
  sa <- rs(A); sb <- rs(B); sc <- rs(C)
  dabc <- sa + sb + sc; dbc <- sb + sc
  reps <- cbind(
    FST = ifelse(dabc > 0, sa / dabc, NA_real_),
    FIT = ifelse(dabc > 0, (sa + sb) / dabc, NA_real_),
    FIS = ifelse(dbc  > 0, sb / dbc, NA_real_),
    HS  = rs(HS) / L,
    HT  = rs(HT) / L
  )
  tab <- do.call(rbind, lapply(nm, function(s) {
    ci <- pga_ci(reps[, s], conf_level)
    data.frame(Statistic = s, Observed = obs[[s]],
               SE = stats::sd(reps[, s], na.rm = TRUE),
               CI_L = ci[["lo"]], CI_U = ci[["hi"]], stringsAsFactors = FALSE)
  }))
  list(table = tab, reps = reps)
}


# Per population x locus: number of typed individuals, observed heterozygosity and
# unbiased gene diversity (Nei & Chesser 1983):
#     Hs = n/(n-1) * ( 1 - sum(p^2) - Ho/(2n) )
# Loci are returned in the column order of `mat` (= order of the data file) and
# populations in the order of their code.
# Returns data.frame(Population, Locus, N, Ho, Hs); Hs is NA when N < 2.
pga_pop_locus_stats <- function(mat, base = 1000L, missing_code = 0L) {
  pop    <- as.integer(mat[, 1L])
  levels <- attr(mat, "pop_levels")
  loci   <- colnames(mat)[-1L]
  pops   <- sort(unique(pop[is.finite(pop) & pop > 0L]))
  out <- vector("list", length(pops) * length(loci)); k <- 0L
  for (p in pops) {
    idx   <- which(pop == p)
    pname <- if (!is.null(levels) && p >= 1L && p <= length(levels)) as.character(levels[[p]]) else as.character(p)
    for (j in seq_along(loci)) {
      g  <- as.integer(mat[idx, j + 1L])
      a1 <- g %/% base; a2 <- g %% base
      ok <- !is.na(g) & g != missing_code & g > 0L & a1 > 0L & a2 > 0L
      a1 <- a1[ok]; a2 <- a2[ok]; n <- length(a1)
      Ho <- if (n > 0L) mean(a1 != a2) else NA_real_
      Hs <- NA_real_
      if (n >= 2L) {
        cnt <- tabulate(c(a1, a2), nbins = max(a1, a2))
        Hs  <- (n / (n - 1)) * (1 - sum((cnt / (2 * n))^2) - Ho / (2 * n))
      }
      k <- k + 1L
      out[[k]] <- data.frame(Population = pname, Locus = loci[j], N = n, Ho = Ho, Hs = Hs,
                             stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}


# Per-population Ho and Hs averaged over loci, WEIGHTED by the number of
# individuals typed at each locus (loci with fewer than 2 typed individuals are
# left out). This is the convention of GENEPOP ("statistics are averaged per
# sample over all loci with at least two individuals typed"); an unweighted mean
# over loci differs from it by up to ~0.005 when loci have different amounts of
# missing data.
# `long` = output of pga_pop_locus_stats(). Returns data.frame(Population, Ho, Hs)
# in the order of first appearance of the populations in `long`.
pga_pop_weighted_means <- function(long) {
  pops <- unique(long$Population)
  do.call(rbind, lapply(pops, function(p) {
    d <- long[long$Population == p & is.finite(long$Hs) & long$N >= 2L, , drop = FALSE]
    data.frame(Population = p,
               Ho = if (nrow(d)) sum(d$N * d$Ho) / sum(d$N) else NA_real_,
               Hs = if (nrow(d)) sum(d$N * d$Hs) / sum(d$N) else NA_real_,
               stringsAsFactors = FALSE)
  }))
}


# HS per population with its two bootstrap confidence intervals.
#   Observed HS  = weighted mean over loci of the per-locus unbiased Hs
#                  (see pga_pop_weighted_means).
#   BS INDIVIDUALS: individuals of the population are resampled with replacement
#                  (the same resampled individuals are used at every locus); the
#                  percentile CI is shifted by (observed - bootstrap mean) to
#                  remove the known downward bias of this resampling (see
#                  pga_ci_bias_shift).
#   BS LOCI      : the loci are resampled with replacement.
# Resampling individuals is done with multinomial weights (equivalent to drawing n
# individuals with replacement) so that allele counts of all replicates come from
# one matrix product: this is orders of magnitude faster than a replicate loop.
# Returns list(table, reps_indiv, reps_loci); the replicate matrices are
# n_boot x (number of populations) and let the user inspect the bootstrap
# distributions.
pga_hs_per_population <- function(mat, base = 1000L, n_boot = 1000L, conf_level = 0.95,
                                  seed = 1, missing_code = 0L) {
  n_boot <- as.integer(n_boot)
  pop    <- as.integer(mat[, 1L])
  levels <- attr(mat, "pop_levels")
  L      <- ncol(mat) - 1L
  pops   <- sort(unique(pop[is.finite(pop) & pop > 0L]))
  pnames <- vapply(pops, function(p)
    if (!is.null(levels) && p >= 1L && p <= length(levels)) as.character(levels[[p]]) else as.character(p),
    character(1))
  set.seed(as.integer(seed %% .Machine$integer.max))

  reps_i <- matrix(NA_real_, n_boot, length(pops), dimnames = list(NULL, pnames))
  reps_l <- reps_i
  rows <- vector("list", length(pops))

  for (ip in seq_along(pops)) {
    idx <- which(pop == pops[ip]); n <- length(idx)
    G   <- matrix(as.integer(mat[idx, -1L, drop = FALSE]), nrow = n)
    A1  <- G %/% as.integer(base); A2 <- G %% as.integer(base)
    typed <- !(is.na(G) | G == missing_code | G <= 0L | A1 <= 0L | A2 <= 0L)

    # multinomial resampling weights: n_boot x n (rows = replicates)
    W <- t(stats::rmultinom(n_boot, size = n, prob = rep(1 / n, n)))
    storage.mode(W) <- "double"

    num <- matrix(0, n_boot, L); den <- matrix(0, n_boot, L)   # per replicate x locus
    obs_hs <- obs_n <- numeric(L)
    for (j in seq_len(L)) {
      tj <- typed[, j]
      if (!any(tj)) { obs_n[j] <- 0; obs_hs[j] <- NA_real_; next }
      a1 <- A1[tj, j]; a2 <- A2[tj, j]
      K  <- max(a1, a2)
      X  <- matrix(0, sum(tj), K)                              # copies of allele k carried by individual i
      X[cbind(seq_along(a1), a1)] <- X[cbind(seq_along(a1), a1)] + 1
      X[cbind(seq_along(a2), a2)] <- X[cbind(seq_along(a2), a2)] + 1
      het <- as.numeric(a1 != a2)

      # observed
      nn <- length(a1)
      if (nn >= 2L) {
        cnt <- colSums(X); ho <- mean(het)
        obs_hs[j] <- (nn / (nn - 1)) * (1 - sum((cnt / (2 * nn))^2) - ho / (2 * nn))
      } else obs_hs[j] <- NA_real_
      obs_n[j] <- nn

      # replicates
      Wj  <- W[, tj, drop = FALSE]
      nb  <- rowSums(Wj)                                       # typed individuals in each replicate
      cb  <- Wj %*% X                                          # n_boot x K allele counts
      hob <- as.vector(Wj %*% het) / nb
      ok  <- nb >= 2
      hsb <- rep(NA_real_, n_boot)
      hsb[ok] <- (nb[ok] / (nb[ok] - 1)) * (1 - rowSums((cb[ok, , drop = FALSE] / (2 * nb[ok]))^2) - hob[ok] / (2 * nb[ok]))
      num[, j] <- ifelse(ok, nb * hsb, 0)
      den[, j] <- ifelse(ok, nb, 0)
    }
    use  <- is.finite(obs_hs) & obs_n >= 2
    obs  <- if (any(use)) sum(obs_n[use] * obs_hs[use]) / sum(obs_n[use]) else NA_real_
    ri   <- ifelse(rowSums(den) > 0, rowSums(num) / rowSums(den), NA_real_)

    # bootstrap over loci (population-specific): resample loci, weighted mean
    li  <- matrix(sample.int(L, L * n_boot, replace = TRUE), nrow = n_boot, ncol = L)
    wv  <- ifelse(use, obs_n, 0); hv <- ifelse(use, obs_hs, 0)
    rl  <- rowSums(matrix(wv[li] * hv[li], nrow = n_boot)) / rowSums(matrix(wv[li], nrow = n_boot))

    ci_i <- pga_ci_bias_shift(obs, ri, conf_level)
    ci_l <- pga_ci(rl, conf_level)
    reps_i[, ip] <- ri; reps_l[, ip] <- rl
    rows[[ip]] <- data.frame(
      Population  = pnames[ip],
      N_loci      = sum(use),
      Observed_HS = obs,
      SE_indiv    = stats::sd(ri, na.rm = TRUE),
      CI_L_indiv  = ci_i[["lo"]], CI_U_indiv = ci_i[["hi"]],
      CI_L_loci   = ci_l[["lo"]], CI_U_loci  = ci_l[["hi"]],
      stringsAsFactors = FALSE)
  }
  list(table = do.call(rbind, rows), reps_indiv = reps_i, reps_loci = reps_l)
}


# FIS per population (Weir & Cockerham 1984, ratio of sums over loci) with a
# bootstrap over LOCI. For each population the per-locus variance components come
# from fis_wc_cpp() applied to that population alone; the multilocus FIS of a
# population is sum(sigma_b) / sum(sigma_b + sigma_w). Loci are resampled with
# replacement (the same resampled loci for every population within a replicate).
# The "Overall" value is the mean of the population FIS values (same definition as
# in the By Population output).
# Returns list(table = data.frame(Population, Observed_FIS, CI_L_loci, CI_U_loci)
#              incl. an "Overall" row, reps = n_boot x (pops + 1)).
pga_fis_pop_loci_boot <- function(mat, base = 1000L, n_boot = 1000L, conf_level = 0.95, seed = 1) {
  n_boot <- as.integer(n_boot)
  pop    <- as.integer(mat[, 1L])
  levels <- attr(mat, "pop_levels")
  pops   <- sort(unique(pop[is.finite(pop) & pop > 0L]))
  pnames <- vapply(pops, function(p)
    if (!is.null(levels) && p >= 1L && p <= length(levels)) as.character(levels[[p]]) else as.character(p),
    character(1))
  L <- ncol(mat) - 1L
  SB <- SW <- matrix(NA_real_, length(pops), L)
  for (ip in seq_along(pops)) {
    sub <- mat[pop == pops[ip], , drop = FALSE]
    sub[, 1L] <- 1L
    if (nrow(sub) < 2L) next
    r <- fis_wc_cpp(sub, base = as.integer(base))
    SB[ip, ] <- as.numeric(r$lsigb); SW[ip, ] <- as.numeric(r$lsigw)
  }
  fis_of <- function(sb, sw) { s <- sb + sw; ifelse(is.finite(s) & s != 0, sb / s, NA_real_) }
  obs <- fis_of(rowSums(SB, na.rm = TRUE), rowSums(SW, na.rm = TRUE))

  set.seed(as.integer(seed %% .Machine$integer.max))
  idx <- matrix(sample.int(L, L * n_boot, replace = TRUE), nrow = n_boot, ncol = L)
  reps <- matrix(NA_real_, n_boot, length(pops) + 1L, dimnames = list(NULL, c(pnames, "Overall")))
  for (ip in seq_along(pops)) {
    sb <- rowSums(matrix(SB[ip, idx], nrow = n_boot)); sw <- rowSums(matrix(SW[ip, idx], nrow = n_boot))
    reps[, ip] <- fis_of(sb, sw)
  }
  reps[, length(pops) + 1L] <- rowMeans(reps[, seq_along(pops), drop = FALSE], na.rm = TRUE)

  all_obs <- c(obs, mean(obs, na.rm = TRUE))
  tab <- do.call(rbind, lapply(seq_len(ncol(reps)), function(k) {
    ci <- pga_ci(reps[, k], conf_level)
    data.frame(Population = colnames(reps)[k], Observed_FIS = all_obs[k],
               CI_L_loci = ci[["lo"]], CI_U_loci = ci[["hi"]], stringsAsFactors = FALSE)
  }))
  list(table = tab, reps = reps)
}


# Overall (multilocus) value of a statistic for each replicate of a bootstrap over
# SUB-SAMPLES, taken as the mean over loci of the per-locus replicates. Used for
# FIT and FIS, for which the C++ routines only return per-locus replicates.
#   FIT_boot, FST_boot : n_boot x L matrices of per-locus replicates.
# FIS is derived per locus and per replicate from  (1 - FIT) = (1 - FIS) (1 - FST).
pga_overall_from_locus_reps <- function(FIT_boot, FST_boot) {
  fis <- 1 - (1 - FIT_boot) / (1 - FST_boot)
  list(FIT = rowMeans(FIT_boot, na.rm = TRUE),
       FIS = rowMeans(fis,      na.rm = TRUE))
}


# Table written in overall_by_population (General Stats):
#   overall : one row per population, Ho and Hs averaged over loci WEIGHTED by the
#             number of individuals typed at each locus (GENEPOP convention, see
#             pga_pop_weighted_means), and the Weir & Cockerham (1984) multilocus
#             FIS of that population.
#   detail  : observed heterozygosity per population and locus, loci in the order
#             of the data file (the per-locus Hs is given by gene_diversity_hs_by_pop).
pga_overall_by_population <- function(mat, base = 1000L, missing_code = 0L) {
  long  <- pga_pop_locus_stats(mat, base, missing_code)
  w     <- pga_pop_weighted_means(long)
  codes <- sort(unique(as.integer(mat[, 1L])))
  codes <- codes[is.finite(codes) & codes > 0L]
  fis   <- wc_fis_by_pop_wc84(dat = mat, pop_col = 0L, base = as.integer(base))
  overall <- data.frame(Population = w$Population, Ho = w$Ho, Hs = w$Hs,
                        `Fis (WC)` = as.numeric(fis[as.character(codes)]),
                        check.names = FALSE, stringsAsFactors = FALSE)
  list(overall = overall, detail = long[, c("Population", "Locus", "Ho")])
}


# ============================================================================
# Table builders for the Subdivision / Diversities output files (pure functions)
# ============================================================================

# Subdivision, section 1: FST per locus (bootstrap over sub-samples) with the
# permutation p-value, then two Overall rows: one per type of bootstrap.
#   final_table      : ID, Observed_FST, P_value, CI_L, CI_U (+ other columns), Overall last
#   locus_boot_table : Statistic, Observed, SE, CI_L, CI_U (bootstrap over loci)
pga_fst_section <- function(final_table, locus_boot_table) {
  ft  <- final_table
  n   <- nrow(ft)
  out <- ft[, c("ID", "Observed_FST", "P_value", "CI_L", "CI_U"), drop = FALSE]
  out$ID[n] <- "Overall (BS/Ss)"
  lb  <- locus_boot_table[locus_boot_table$Statistic == "FST", , drop = FALSE]
  rbind(out, data.frame(ID = "Overall (BS/Loci)", Observed_FST = ft$Observed_FST[n],
                        P_value = ft$P_value[n],
                        CI_L = if (nrow(lb)) lb$CI_L else NA_real_,
                        CI_U = if (nrow(lb)) lb$CI_U else NA_real_,
                        stringsAsFactors = FALSE))
}

# Multilocus estimators (FIS, FST, FIT, HS, HT): observed value, bootstrap standard
# error and CI from the bootstrap over LOCI, and CI from the bootstrap over
# SUB-SAMPLES (NA when fewer than 5 sub-samples).
pga_multilocus_section <- function(locus_boot_table, subs_overall) {
  ord <- c("FIS", "FST", "FIT", "HS", "HT")
  lb  <- locus_boot_table[match(ord, locus_boot_table$Statistic), , drop = FALSE]
  sb  <- subs_overall[match(ord, subs_overall$Statistic), , drop = FALSE]
  data.frame(Statistic = ord, Observed = lb$Observed, SE_loci = lb$SE,
             CI_L_loci = lb$CI_L, CI_U_loci = lb$CI_U,
             CI_L_subs = sb$CI_L_subs, CI_U_subs = sb$CI_U_subs,
             stringsAsFactors = FALSE)
}

# A per-locus HS / HT table as written in the files: standard error instead of the
# bootstrap mean (which has no use for the reader).
pga_diversity_section <- function(tbl, obs_col) {
  out <- data.frame(ID = tbl$ID, tbl[[obs_col]], SE = tbl$Boot_SE, CI_L = tbl$CI_L, CI_U = tbl$CI_U,
                    stringsAsFactors = FALSE)
  names(out)[2] <- obs_col
  out
}

# HS per population: one row per population (bootstrap over individuals and over
# loci) and an Overall row that carries the same two types of bootstrap.
#   hs_per_pop : table of pga_hs_per_population()
#   hs_indiv_tbl, locus_boot_table : tables of run_bootstrap_fst_analysis()
pga_hs_pop_section <- function(hs_per_pop, hs_indiv_tbl, locus_boot_table) {
  ov_i <- hs_indiv_tbl[hs_indiv_tbl$ID == "Overall", , drop = FALSE]
  ov_l <- locus_boot_table[locus_boot_table$Statistic == "HS", , drop = FALSE]
  tab  <- hs_per_pop
  tab$N_loci <- NULL
  rbind(tab, data.frame(
    Population  = "Overall",
    Observed_HS = ov_i$Observed_HS,
    SE_indiv    = ov_i$Boot_SE,
    CI_L_indiv  = ov_i$CI_L, CI_U_indiv = ov_i$CI_U,
    CI_L_loci   = ov_l$CI_L, CI_U_loci  = ov_l$CI_U,
    stringsAsFactors = FALSE))
}

# Writes a titled block of bootstrap replicate values (one row per replicate).
pga_write_replicates <- function(con, title, mat) {
  if (is.null(mat) || !length(mat)) return(invisible(NULL))
  writeLines(title, con = con, useBytes = TRUE)
  utils::write.table(data.frame(Replicate = seq_len(nrow(mat)), mat, check.names = FALSE),
                     file = con, sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines("", con = con)
  invisible(NULL)
}


# ============================================================================
# Counting of permuted statistics that are "larger or equal" to the observed one
# ============================================================================
# A permutation p-value is the proportion of randomised statistics LARGER OR EQUAL to
# the observed one: ties count (FSTAT manual, sections 7 and 8.2). Ties are frequent
# for discrete statistics (e.g. the FIS of one sample under allele permutation only
# depends on the number of heterozygotes, so it takes 15-25 distinct values), but the
# same mathematical value can come out of two computations differing by ~1e-16, and a
# plain `perm >= obs` then loses a random part of the ties: the one-sided p-values of a
# locus in one sample were off by up to 0.15 and their sum was ~1 instead of ~1.17.
# All comparisons therefore use this tolerance.
PGA_TIE_TOL <- 1e-9

# number of permuted values >= observed (ties included)
pga_n_ge <- function(null, obs, tol = PGA_TIE_TOL) sum(null >= obs - tol)
# number of permuted values <= observed (ties included)
pga_n_le <- function(null, obs, tol = PGA_TIE_TOL) sum(null <= obs + tol)
# number of permuted values strictly greater than the observed one (ties excluded)
pga_n_gt <- function(null, obs, tol = PGA_TIE_TOL) sum(null > obs + tol)
# two-sided: number of permuted |values| >= |observed| (ties included)
pga_n_abs_ge <- function(null, obs, tol = PGA_TIE_TOL) sum(abs(null) >= abs(obs) - tol)
