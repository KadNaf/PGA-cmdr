# ============================================================================
# utils_stats.R - pure statistical helpers (no Shiny, no database)
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
  se  <- if ("SE" %in% names(ft)) ft$SE else rep(NA_real_, n)
  pg <- if ("P_theta_gt" %in% names(ft)) ft$P_theta_gt else rep(NA_real_, n)
  pl <- if ("P_theta_lt" %in% names(ft)) ft$P_theta_lt else rep(NA_real_, n)
  out <- data.frame(ID = ft$ID, Observed_FST = ft$Observed_FST, SE = se,
                    P_value = ft$P_value, P_theta_gt = pg, P_theta_lt = pl,
                    CI_L = ft$CI_L, CI_U = ft$CI_U,
                    stringsAsFactors = FALSE)
  out$ID[n] <- "Overall (BS/Indiv)"
  lb  <- locus_boot_table[locus_boot_table$Statistic == "FST", , drop = FALSE]
  rbind(out, data.frame(ID = "Overall (BS/Loci)", Observed_FST = ft$Observed_FST[n],
                        SE = if (nrow(lb)) lb$SE else NA_real_,
                        P_value = ft$P_value[n], P_theta_gt = pg[n], P_theta_lt = pl[n],
                        CI_L = if (nrow(lb)) lb$CI_L else NA_real_,
                        CI_U = if (nrow(lb)) lb$CI_U else NA_real_,
                        stringsAsFactors = FALSE))
}

# Multilocus estimators (FIS, FST, FIT, HS, HT): observed value, standard error and
# CI from the bootstrap over LOCI, and standard error and CI from the bootstrap
# over INDIVIDUALS (bias-shifted).
pga_multilocus_section <- function(locus_boot_table, subs_overall,
                                   stats = c("FIS", "FST", "FIT", "HS", "HT")) {
  ord <- stats
  lb  <- locus_boot_table[match(ord, locus_boot_table$Statistic), , drop = FALSE]
  sb  <- subs_overall[match(ord, subs_overall$Statistic), , drop = FALSE]
  data.frame(Statistic = ord, Observed = lb$Observed,
             SE_loci = lb$SE, CI_L_loci = lb$CI_L, CI_U_loci = lb$CI_U,
             SE_indiv = sb$SE_indiv, CI_L_indiv = sb$CI_L_subs, CI_U_indiv = sb$CI_U_subs,
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

# HS per sub-sample: one row per sub-sample (bootstrap over individuals and over
# loci) and an Overall row that carries the same two types of bootstrap.
#   hs_per_pop : table of pga_hs_per_population()
#   hs_indiv_tbl, locus_boot_table : tables of run_bootstrap_fst_analysis()
pga_hs_pop_section <- function(hs_per_pop, hs_indiv_tbl, locus_boot_table) {
  ov_i <- hs_indiv_tbl[hs_indiv_tbl$ID == "Overall", , drop = FALSE]
  ov_l <- locus_boot_table[locus_boot_table$Statistic == "HS", , drop = FALSE]
  tab  <- hs_per_pop
  tab$N_loci <- NULL
  res <- rbind(tab, data.frame(
    Population  = "Overall",
    Observed_HS = ov_i$Observed_HS,
    SE_indiv    = ov_i$Boot_SE,
    CI_L_indiv  = ov_i$CI_L, CI_U_indiv = ov_i$CI_U,
    CI_L_loci   = ov_l$CI_L, CI_U_loci  = ov_l$CI_U,
    stringsAsFactors = FALSE))
  names(res)[names(res) == "Population"] <- "Sub_sample"
  res
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


# ---------------------------------------------------------------------------
# Pairwise FST between sub-samples (Weir & Cockerham 1984), no GPS needed.
#   For every pair of sub-samples: multilocus FST (ratio of sums) and per-locus FST,
#   CI by bootstrap over INDIVIDUALS (bias-shifted, see pga_ci_bias_shift) and over
#   LOCI, and a permutation p-value (individuals reassigned between the two
#   sub-samples; p = (b + 1)/(m + 1), one-sided "greater" on FST).
# Returns list(overall = data.frame, per_locus = data.frame).
# ---------------------------------------------------------------------------
pga_pairwise_fst <- function(mat, base = 1000L, n_boot = 1000L, n_perm = 1000L,
                             conf_level = 0.95, seed = 1, missing_code = 0L,
                             n_threads = 1L) {
  mat <- as.matrix(mat); storage.mode(mat) <- "integer"
  levels <- attr(mat, "pop_levels")
  codes  <- sort(unique(mat[, 1]))
  if (is.null(levels)) levels <- as.character(codes) else levels <- as.character(levels)
  nm <- function(k) if (k <= length(levels)) levels[k] else as.character(k)
  loc <- colnames(mat)[-1L]
  if (is.null(loc)) loc <- paste0("L", seq_len(ncol(mat) - 1L))
  ov <- list(); pl <- list(); k <- 0L
  for (i in seq_len(length(codes) - 1L)) for (j in (i + 1L):length(codes)) {
    if (j > length(codes)) next
    k <- k + 1L
    sub <- mat[mat[, 1] %in% codes[c(i, j)], , drop = FALSE]
    sub[, 1] <- ifelse(sub[, 1] == codes[i], 1L, 2L)
    n1 <- sum(sub[, 1] == 1L); n2 <- sum(sub[, 1] == 2L)
    obs <- observed_wc84_stats_cpp(sub, 1L, as.integer(missing_code), as.integer(base))
    comp <- wc84_locus_components_cpp(sub, 1L, as.integer(missing_code), as.integer(base))
    bl <- pga_boot_over_loci(comp$A, comp$B, comp$C, comp$HS, comp$HT,
                             n_boot = n_boot, conf_level = conf_level, seed = seed + k)$table
    bl <- bl[bl$Statistic == "FST", ]
    bi <- boot_indiv_wc84_cpp(sub, 1L, as.integer(missing_code), as.integer(base),
                              as.integer(n_boot), as.numeric(seed + k), as.integer(n_threads))
    ci_o <- pga_ci_bias_shift(obs$FST_overall_ratio_of_sums, bi$FST_overall_boot, conf_level)
    p_o <- NA_real_; p_l <- rep(NA_real_, length(loc))
    if (n_perm > 0) {
      pr <- batch_permute_wc84_fst_auto(sub, 1L, as.integer(missing_code), as.integer(base),
                                        B = as.integer(n_perm), n_threads = as.integer(n_threads),
                                        seed = seed + k, pval_method = "greater",
                                        perm_scheme = "permute_pop_labels", debug = FALSE)
      p_o <- as.numeric(pr$p_FST_overall); p_l <- as.numeric(pr$p_G)
    }
    ov[[k]] <- data.frame(Sub_sample_1 = nm(codes[i]), Sub_sample_2 = nm(codes[j]),
                          N1 = n1, N2 = n2, FST = obs$FST_overall_ratio_of_sums,
                          SE_indiv = stats::sd(bi$FST_overall_boot, na.rm = TRUE),
                          CI_L_indiv = ci_o[["lo"]], CI_U_indiv = ci_o[["hi"]],
                          SE_loci = bl$SE, CI_L_loci = bl$CI_L, CI_U_loci = bl$CI_U,
                          P_value = p_o, stringsAsFactors = FALSE)
    ci_l <- vapply(seq_along(loc), function(q)
      pga_ci_bias_shift(obs$FST[q], bi$FST_boot[, q], conf_level), numeric(2))
    pl[[k]] <- data.frame(Sub_sample_1 = nm(codes[i]), Sub_sample_2 = nm(codes[j]),
                          Locus = loc, FST = as.numeric(obs$FST),
                          SE_indiv = apply(bi$FST_boot, 2, stats::sd, na.rm = TRUE),
                          CI_L_indiv = ci_l[1, ], CI_U_indiv = ci_l[2, ],
                          P_value = p_l, stringsAsFactors = FALSE)
  }
  list(overall = do.call(rbind, ov), per_locus = do.call(rbind, pl))
}


# ---------------------------------------------------------------------------
# FIS (WC84) per locus x sub-sample, with CI and p-values.
#   For every sub-sample: observed FIS at each locus, CI by bootstrap over
#   INDIVIDUALS (individuals resampled with replacement inside the sub-sample),
#   two-sided permutation p-value and the two one-sided ones (FSTAT: heterozygote
#   deficit = permuted FIS >= observed; excess = permuted FIS <= observed).
# Returns a long data.frame: Locus, Sub_sample, N, FIS, SE, CI_L, CI_U,
#   P_value, P_deficit, P_excess.  CI is NA when the sub-sample has < 5 individuals.
# ---------------------------------------------------------------------------
pga_fis_locus_pop <- function(mat, base = 1000L, pop_names = NULL, n_boot = 1000L,
                              n_perm = 1000L, conf_level = 0.95, seed = 1) {
  mat <- as.matrix(mat); storage.mode(mat) <- "integer"
  loc <- colnames(mat)[-1L]
  codes <- sort(unique(mat[, 1]))
  if (is.null(pop_names)) pop_names <- as.character(codes)
  out <- vector("list", length(codes))
  for (k in seq_along(codes)) {
    sub <- mat[mat[, 1] == codes[k], , drop = FALSE]; sub[, 1] <- 1L
    n <- nrow(sub)
    obs <- as.numeric(fis_wc_cpp(sub, base = as.integer(base))$FIS)
    se <- cl <- cu <- rep(NA_real_, length(loc))
    if (n >= 5L && n_boot > 0) {
      set.seed(as.integer((seed + k) %% .Machine$integer.max))
      reps <- matrix(NA_real_, n_boot, length(loc))
      for (b in seq_len(n_boot))
        reps[b, ] <- as.numeric(fis_wc_cpp(sub[sample.int(n, n, replace = TRUE), , drop = FALSE],
                                           base = as.integer(base))$FIS)
      for (j in seq_along(loc)) {
        ci <- pga_ci(reps[, j], conf_level); cl[j] <- ci[["lo"]]; cu[j] <- ci[["hi"]]
        se[j] <- stats::sd(reps[, j], na.rm = TRUE)
      }
    }
    p2 <- pd <- pe <- rep(NA_real_, length(loc))
    if (n >= 2L && n_perm > 0) {
      set.seed(as.integer((seed + k) %% .Machine$integer.max))
      perm <- batch_permute_wc_fis(dat = sub, pop_col_1based = 1L, base = as.integer(base),
                                   B = as.integer(n_perm))
      perm <- perm[, seq_along(loc), drop = FALSE]
      for (j in seq_along(loc)) {
        if (!is.finite(obs[j])) next
        pj <- perm[, j]; pj <- pj[is.finite(pj)]
        if (!length(pj)) next
        m <- length(pj)
        p2[j] <- (pga_n_abs_ge(pj, obs[j]) + 1) / (m + 1)
        pd[j] <- (pga_n_ge(pj, obs[j]) + 1) / (m + 1)
        pe[j] <- (pga_n_le(pj, obs[j]) + 1) / (m + 1)
      }
    }
    out[[k]] <- data.frame(Locus = loc, Sub_sample = pop_names[k], N = n, FIS = obs, SE = se,
                           CI_L = cl, CI_U = cu, P_value = p2, P_deficit = pd, P_excess = pe,
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}


# Moves the standard-error columns (SE, SE_*, Boot_SE) and then the p-value columns
# (P_value, P_*, p_*) to the END of a table, so that the numeric columns that are
# plotted (estimate, CI bounds) stay contiguous in a spreadsheet.
pga_se_p_last <- function(df) {
  if (!is.data.frame(df)) return(df)
  nm <- names(df)
  is_se <- grepl("^(SE($|_)|Boot_SE$)", nm)
  is_p  <- grepl("^[Pp]_", nm) | nm %in% c("P_value", "p_value")
  df[, c(nm[!is_se & !is_p], nm[is_se], nm[is_p]), drop = FALSE]
}

# write.table() with tab separation, SE and p-value columns last.
pga_write_table <- function(x, file, ...) {
  utils::write.table(pga_se_p_last(x), file = file, sep = "\t", row.names = FALSE, quote = FALSE)
}
