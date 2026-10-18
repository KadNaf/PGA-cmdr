// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ═════════════════════════════════════════════════════════════════════════
//  Bootstrap over sub-samples (= populations, resampled as whole blocks
//  with replacement) for global + per-locus Weir & Cockerham (1984) FST,
//  raw and ENA-corrected.
//
//  Unlike bootstrap-over-loci (already fully vectorised in R — one
//  sample.int() call for all replicates at once, no per-replicate R loop),
//  this bootstrap genuinely recomputes the FST components from scratch on
//  every replicate (a different set of populations enters each time), so
//  it is the one bootstrap in this module worth moving to native code.
//
//  No EM re-run is needed here or in the R precomputation step: population
//  membership, not individual genotypes, is what gets resampled, so the
//  per-locus x per-population x per-allele summary counts (ni, nA, AA,
//  AA_corr) can be extracted ONCE in R (server_null_alleles.R,
//  .precompute_boot_subs_data()) before the bootstrap starts. This
//  function only repeats the numeric aggregation step — the same
//  Weir & Cockerham analysis-of-variance formula as
//  weir_components_allele()/compute_fst_global_full() in the R code,
//  checked line-for-line against it (2026-08-12) — for each of the
//  resampled population sets.
// ═════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List boot_subsamples_fst_cpp(List perLocusData, int nboot) {
  RNGScope scope;  // synchronise with R's set.seed()

  int L = perLocusData.size();
  if (L < 1) stop("No loci provided.");

  std::vector<NumericVector> ni_raw(L), ni_ena(L);
  std::vector<NumericMatrix> nA_raw(L), AA_raw(L), nA_ena(L), AA_corr(L);
  for (int li = 0; li < L; li++) {
    List loc = perLocusData[li];
    ni_raw[li]  = as<NumericVector>(loc["ni_raw"]);
    ni_ena[li]  = as<NumericVector>(loc["ni_ena"]);
    nA_raw[li]  = as<NumericMatrix>(loc["nA_raw"]);
    AA_raw[li]  = as<NumericMatrix>(loc["AA_raw"]);
    nA_ena[li]  = as<NumericMatrix>(loc["nA_ena"]);
    AA_corr[li] = as<NumericMatrix>(loc["AA_corr"]);
  }
  int n_pop = ni_raw[0].size();

  NumericVector bootRaw(nboot), bootEna(nboot);
  NumericMatrix bootRawLoc(nboot, L), bootEnaLoc(nboot, L);

  std::vector<int> samp(n_pop);

  // Weir & Cockerham per-allele variance-component step, shared by the raw
  // and ENA branches (only the input summary vectors differ).
  auto wcAllele = [&](const NumericVector& ni_v, const NumericMatrix& nA_m,
                      const NumericMatrix& AA_m, int a, double N, double r_eff,
                      double nc, double& s2P, double& s2I, double& s2G) -> bool {
    if (N <= 0.0 || r_eff < 2.0 || nc <= 0.0 || (N - r_eff) <= 0.0) return false;
    double snA = 0.0, s2A = 0.0, sAA = 0.0;
    for (int i = 0; i < n_pop; i++) {
      int pidx = samp[i];
      double niv = ni_v[pidx];
      double nAv = nA_m(pidx, a);
      double AAv = AA_m(pidx, a);
      snA += nAv;
      if (niv > 0.0) s2A += (nAv * nAv) / (2.0 * niv);
      sAA += AAv;
    }
    double MSG = (0.5 * snA - sAA) / N;
    double MSI = (0.5 * snA + sAA - s2A) / (N - r_eff);
    double MSP = (s2A - 0.5 * snA * snA / N) / (r_eff - 1.0);
    s2P = (MSP - MSI) / (2.0 * nc);
    s2I = 0.5 * (MSI - MSG);
    s2G = MSG;
    return true;
  };

  for (int b = 0; b < nboot; b++) {
    for (int i = 0; i < n_pop; i++)
      samp[i] = (int) std::floor(unif_rand() * n_pop);  // WITH replacement (bootstrap)

    double s1_total = 0.0, s3_total = 0.0, s1c_total = 0.0, s3c_total = 0.0;

    for (int li = 0; li < L; li++) {
      int n_allele = nA_raw[li].ncol();

      // --- raw branch ---
      double N = 0.0, N2 = 0.0; int r = 0;
      for (int i = 0; i < n_pop; i++) {
        double v = ni_raw[li][samp[i]];
        N += v; N2 += v * v; if (v > 0.0) r++;
      }
      double nc = (N > 0.0 && r > 1) ? (N - N2 / N) / (r - 1) : 0.0;
      double s1l = 0.0, s3l = 0.0;
      for (int a = 0; a < n_allele; a++) {
        double s2P, s2I, s2G;
        if (wcAllele(ni_raw[li], nA_raw[li], AA_raw[li], a, N, (double) r, nc, s2P, s2I, s2G)) {
          s1l += s2P; s3l += s2P + s2I + s2G;
        }
      }
      bootRawLoc(b, li) = (s3l != 0.0) ? (s1l / s3l) : NA_REAL;
      s1_total += s1l * nc; s3_total += s3l * nc;

      // --- ENA branch ---
      double Nc = 0.0, N2c = 0.0; int rc = 0;
      for (int i = 0; i < n_pop; i++) {
        double v = ni_ena[li][samp[i]];
        Nc += v; N2c += v * v; if (v > 0.0) rc++;
      }
      double ncc = (Nc > 0.0 && rc > 1) ? (Nc - N2c / Nc) / (rc - 1) : 0.0;
      double s1lc = 0.0, s3lc = 0.0;
      for (int a = 0; a < n_allele; a++) {
        double s2P, s2I, s2G;
        if (wcAllele(ni_ena[li], nA_ena[li], AA_corr[li], a, Nc, (double) rc, ncc, s2P, s2I, s2G)) {
          s1lc += s2P; s3lc += s2P + s2I + s2G;
        }
      }
      bootEnaLoc(b, li) = (s3lc != 0.0) ? (s1lc / s3lc) : NA_REAL;
      s1c_total += s1lc * ncc; s3c_total += s3lc * ncc;
    }

    bootRaw[b] = (s3_total > 0.0) ? (s1_total / s3_total) : NA_REAL;
    bootEna[b] = (s3c_total > 0.0) ? (s1c_total / s3c_total) : NA_REAL;
  }

  return List::create(
    _["boot_raw"]     = bootRaw,
    _["boot_ena"]     = bootEna,
    _["boot_raw_loc"] = bootRawLoc,
    _["boot_ena_loc"] = bootEnaLoc
  );
}
