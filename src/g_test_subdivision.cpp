// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ═════════════════════════════════════════════════════════════════════════
//  G-based permutation test for population subdivision (FSTAT §7.1 /
//  Goudet et al. 1996 convention): per-locus G statistic (log-likelihood
//  ratio) on the population x allele contingency table, global G = sum of
//  per-locus G (additive property), significance by permuting WHOLE
//  complete multilocus genotypes among populations.
//
//  Checked against the R implementation in server_general_stats.R
//  (.g_stat_from_flat / .g_stat_vec / the run_G_test observer,
//  2026-08-13): same contingency-table G formula (2*sum(n_ij*ln(n_ij/E_ij))),
//  same permutation scheme (one global permutation of individuals per
//  replicate, re-used across all loci, subsetting to each locus's own
//  valid rows — matches FSTAT's requirement that only complete multilocus
//  genotypes are randomised).
//
//  Runs one BATCH of permutations per call (not the whole n_perm at once)
//  so the R side can call this repeatedly and update a progress bar
//  between batches, exactly like the Null Alleles module's sub-samples
//  bootstrap C++ engine.
// ═════════════════════════════════════════════════════════════════════════

static double g_stat_pop_allele(const std::vector<int>& popv, const std::vector<int>& allele0,
                                 int n_pop, int n_allele) {
  if (n_allele < 2) return NA_REAL;
  int n_ind = (int) popv.size();
  std::vector<double> cnt((size_t) n_pop * (size_t) n_allele, 0.0);
  for (int i = 0; i < n_ind; i++) {
    int p = popv[i];
    int a1 = allele0[i];
    int a2 = allele0[i + n_ind];
    cnt[(size_t) p + (size_t) n_pop * a1] += 1.0;
    cnt[(size_t) p + (size_t) n_pop * a2] += 1.0;
  }
  std::vector<double> rs(n_pop, 0.0), cs(n_allele, 0.0);
  double n = 0.0;
  for (int a = 0; a < n_allele; a++)
    for (int p = 0; p < n_pop; p++) {
      double v = cnt[(size_t) p + (size_t) n_pop * a];
      rs[p] += v; cs[a] += v; n += v;
    }
  int rpos = 0; for (double v : rs) if (v > 0.0) rpos++;
  int cpos = 0; for (double v : cs) if (v > 0.0) cpos++;
  if (n == 0.0 || rpos < 2 || cpos < 2) return NA_REAL;
  double G = 0.0; bool any = false;
  for (int a = 0; a < n_allele; a++)
    for (int p = 0; p < n_pop; p++) {
      double c = cnt[(size_t) p + (size_t) n_pop * a];
      if (c <= 0.0) continue;
      double E = rs[p] * cs[a] / n;
      if (E <= 0.0) continue;
      G += c * std::log(c / E);
      any = true;
    }
  if (!any) return NA_REAL;
  return 2.0 * G;
}

// [[Rcpp::export]]
List g_test_subdivision_batch_cpp(IntegerVector pop_idx0_full, List loci_idx0, List loci_allele0,
                                   IntegerVector loci_n_allele, int n_pop, int n_batch) {
  RNGScope scope;  // synchronise with R's set.seed()

  int n_loci = loci_idx0.size();
  int n_ind  = pop_idx0_full.size();

  std::vector<std::vector<int> > idx0(n_loci), allele0(n_loci);
  for (int j = 0; j < n_loci; j++) {
    idx0[j]    = as<std::vector<int> >(loci_idx0[j]);
    allele0[j] = as<std::vector<int> >(loci_allele0[j]);
  }
  std::vector<int> pop0(pop_idx0_full.begin(), pop_idx0_full.end());

  NumericMatrix gNullLocus(n_batch, n_loci);
  NumericVector gNullOverall(n_batch);

  std::vector<int> perm(pop0);
  std::vector<int> popj;  // scratch: per-locus subset of the permuted population vector

  for (int b = 0; b < n_batch; b++) {
    // Fisher-Yates shuffle of the FULL individual vector, once per
    // replicate — matches R's sample.int(n_ind) applied once and then
    // subset per locus (not a fresh permutation per locus).
    for (int t = n_ind - 1; t > 0; t--) {
      int k = (int) std::floor(unif_rand() * (t + 1));
      std::swap(perm[k], perm[t]);
    }
    double sumb = 0.0;
    for (int j = 0; j < n_loci; j++) {
      int nv = (int) idx0[j].size();
      popj.resize(nv);
      for (int i = 0; i < nv; i++) popj[i] = perm[idx0[j][i]];
      double g = g_stat_pop_allele(popj, allele0[j], n_pop, loci_n_allele[j]);
      gNullLocus(b, j) = g;
      if (R_finite(g)) sumb += g;
    }
    gNullOverall[b] = sumb;
  }

  return List::create(
    _["g_null_locus"]   = gNullLocus,
    _["g_null_overall"] = gNullOverall
  );
}

// [[Rcpp::export]]
List g_test_subdivision_observed_cpp(IntegerVector pop_idx0_full, List loci_idx0, List loci_allele0,
                                      IntegerVector loci_n_allele, int n_pop) {
  int n_loci = loci_idx0.size();
  std::vector<int> pop0(pop_idx0_full.begin(), pop_idx0_full.end());

  NumericVector gObs(n_loci);
  for (int j = 0; j < n_loci; j++) {
    std::vector<int> idx0j = as<std::vector<int> >(loci_idx0[j]);
    std::vector<int> allele0j = as<std::vector<int> >(loci_allele0[j]);
    std::vector<int> popj((int) idx0j.size());
    for (size_t i = 0; i < idx0j.size(); i++) popj[i] = pop0[idx0j[i]];
    gObs[j] = g_stat_pop_allele(popj, allele0j, n_pop, loci_n_allele[j]);
  }
  double gObsOverall = 0.0;
  for (int j = 0; j < n_loci; j++) if (R_finite(gObs[j])) gObsOverall += gObs[j];

  return List::create(_["g_obs"] = gObs, _["g_obs_overall"] = gObsOverall);
}
