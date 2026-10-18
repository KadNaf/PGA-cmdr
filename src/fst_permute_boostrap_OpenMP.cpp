// ============================================================================
// fst_permute_boostrap_OpenMP.cpp
// R-free OpenMP kernels + thin Rcpp wrappers.
// ============================================================================
//
// ============================================================================
// FUNCTIONS IN THIS FILE (in order of appearance) — Weir & Cockerham (1984)
// FST/FIT/FIS core, plus Nei's gene-diversity statistics. All of Subdivision,
// Genetic Diversities, and much of General Stats route through these.
//
//   nei_het_stats_cpp()             Per-locus Ho, Hs, Ht (Nei & Chesser 1983
//                                    unbiased estimator for Hs — see the
//                                    fixed formula and comment inline below).
//   observed_wc84_stats_cpp()       Per-locus WC84 FIT/FIS/FST from the
//                                    variance-component decomposition
//                                    (a = among pops, b = among indiv within
//                                    pops, c = within individuals).
//   wc84_locus_components_cpp()     Same variance components (a, b, c),
//                                    returned per-locus for downstream
//                                    multilocus ratio-of-sums aggregation.
//   batch_permute_wc84_fst_parallel()  Permutation test for FST: shuffle
//                                    individuals across populations.
//   boot_popblock_wc84_parallel()   Bootstrap over SUB-SAMPLES (populations):
//                                    resample populations with replacement,
//                                    for a CI that reflects uncertainty in
//                                    which populations were sampled.
//   boot_indiv_hs_cpp()             Bootstrap Hs over individuals within
//                                    each population.
// ============================================================================
//
// RNG / REPRODUCIBILITY NOTE (added following a code audit, see thesis/
// publication material - keep in sync with fis_permute_bootstrap_wc.cpp and
// fit_permute_bootstrap_wc.cpp, which use a DIFFERENT scheme):
//
//   This file uses std::mt19937_64, seeded per-thread as a function of the
//   `seed` parameter (see e.g. `std::mt19937_64 rng(sb)` below). This is safe
//   to call from inside an OpenMP parallel region (unlike R's own RNG, which
//   is NOT thread-safe - see the note in fis_permute_bootstrap_wc.cpp).
//
//   Reproducibility claim: for a FIXED `seed` AND a FIXED `n_threads`
//   (equivalently, a fixed static/dynamic OpenMP scheduling and chunking),
//   results are deterministic. Changing `n_threads` between two runs with
//   the same `seed` can change which permutations/bootstrap draws are
//   assigned to which thread-local generator, and therefore CAN change the
//   exact bootstrap replicates drawn (though not their statistical
//   properties). If exact bit-for-bit reproducibility across machines with
//   different core counts is claimed in the publication, it must be
//   verified explicitly against this scheme, e.g. by re-running the same
//   analysis with n_threads = 1 and n_threads > 1 and comparing full
//   bootstrap distributions (not just summary CIs).
//
//   This is in contrast to fis_permute_bootstrap_wc.cpp / 
//   fit_permute_bootstrap_wc.cpp, whose permutation loops use R's own
//   R::unif_rand() (reproducible via set.seed() from R, but explicitly NOT
//   parallelised for that reason - see the comment there). The two
//   estimator families (FST here vs FIS/FIT there) are therefore NOT
//   currently on a unified RNG scheme; harmonizing them (e.g. by moving all
//   three to a common thread-safe, seed-per-replicate scheme) is flagged in
//   the accompanying audit as a P1 item and should be done as a dedicated,
//   test-covered change rather than inside an unrelated fix.
// ============================================================================

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <random>
#include <cstring>
#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Tolerance used when a permuted statistic is compared with the observed one. A
// permutation p-value counts the randomised statistics LARGER OR EQUAL to the observed
// one (ties included, FSTAT manual); the same mathematical value can however come out
// of two computations differing by ~1e-16, and a plain `>=` then loses a random part of
// the ties and understates the p-values (see utils_stats.R, PGA_TIE_TOL).
static const double TIE_TOL = 1e-9;

// ============================================================================
// Structs
// ============================================================================
struct WC84Comp {
  double a, b, c;
  bool ok;
};

struct WC84LocusStats {
  double a, b, c;
  double Ho, Hs, Ht;
  bool ok;
};

// ============================================================================
// Genotype codec (base encoding)
// ============================================================================
inline bool decode_gt_base(int gt, int& a1, int& a2, int base) {
  if (gt == NA_INTEGER || gt <= 0) return false;
  a1 = gt / base;
  a2 = gt % base;
  return (a1 > 0 && a2 > 0);
}

inline bool decode_gt_base_with_missing(int gt, int missing_code, int& a1, int& a2, int base) {
  if (gt == missing_code) return false;
  return decode_gt_base(gt, a1, a2, base);
}

inline int encode_gt_base(int a1, int a2, int base) {
  return a1 * base + a2;
}

// ============================================================================
// Column-major matrix access (R matrices are column-major)
// ============================================================================
static inline int mat_get_colmajor(const int* ptr, int ld, int i, int j) {
  return ptr[i + (size_t)ld * (size_t)j];
}
static inline void mat_set_colmajor(int* ptr, int ld, int i, int j, int v) {
  ptr[i + (size_t)ld * (size_t)j] = v;
}

// ============================================================================
// NA helpers
// ============================================================================
static inline bool is_na_double(double v) {
  return std::isnan(v); // NA_REAL is NaN payload
}

static inline double mean_across_loci_ptr(const double* x, int L) {
  double s = 0.0; int k = 0;
  for (int i = 0; i < L; ++i) {
    const double v = x[i];
    if (is_na_double(v)) continue;
    s += v; k++;
  }
  return (k > 0) ? (s / (double)k) : NA_REAL;
}

// ============================================================================
// Strata builder (pop -> row indices) for "within_pop_alleles" permutations
// ============================================================================
static std::vector<std::vector<int>> build_pop_rows(const IntegerMatrix& mat, int pop_col) {
  const int n = mat.nrow();
  std::unordered_map<int, std::vector<int>> pop_map;
  pop_map.reserve(64);
  
  for (int i = 0; i < n; ++i) pop_map[ mat(i, pop_col) ].push_back(i);
  
  std::vector<std::vector<int>> rows;
  rows.reserve(pop_map.size());
  for (auto& kv : pop_map) rows.push_back(std::move(kv.second));
  return rows;
}

// ============================================================================
// R-FREE permutation helpers operating on raw int* + std::mt19937_64
// ============================================================================
static inline void permute_pop_labels_once_rng_buf(
    int* dat_ptr, int n, int p, int pop_col, std::mt19937_64& rng) {
  
  (void)p;
  std::vector<int> labels((size_t)n);
  for (int i = 0; i < n; ++i) labels[i] = mat_get_colmajor(dat_ptr, n, i, pop_col);
  std::shuffle(labels.begin(), labels.end(), rng);
  for (int i = 0; i < n; ++i) mat_set_colmajor(dat_ptr, n, i, pop_col, labels[i]);
}

static inline void permute_within_pops_once_rng_buf(
    int* dat_ptr, int n, int p, int pop_col,
    int missing_code, int base,
    const std::vector<std::vector<int>>& strata,
    std::mt19937_64& rng) {
  
  std::vector<int> alleles;
  std::vector<int> valid_rows;
  alleles.reserve(256);
  valid_rows.reserve(128);
  
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    for (const auto& rows : strata) {
      if (rows.size() < 2) continue;
      
      alleles.clear();
      valid_rows.clear();
      
      for (int idx : rows) {
        int a1, a2;
        const int gt = mat_get_colmajor(dat_ptr, n, idx, j);
        if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
        alleles.push_back(a1);
        alleles.push_back(a2);
        valid_rows.push_back(idx);
      }
      
      if (alleles.size() < 2 || valid_rows.empty()) continue;
      std::shuffle(alleles.begin(), alleles.end(), rng);
      
      int k = 0;
      const int aN = (int)alleles.size();
      for (int idx : valid_rows) {
        if (k + 1 >= aN) break;
        mat_set_colmajor(dat_ptr, n, idx, j, encode_gt_base(alleles[k], alleles[k + 1], base));
        k += 2;
      }
    }
  }
}

// ============================================================================
// Pointer-based WC84 locus stats (a,b,c + Ho/Hs/Ht)
// ============================================================================
static inline WC84LocusStats wc84_locus_stats_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int locus_col, int pop_col, int missing_code, int base
) {
  (void)p;
  
  std::unordered_map<int,int> pop_to_idx;
  pop_to_idx.reserve(64);
  
  std::vector<int> pop_labels;
  std::vector<int> n_i;
  std::vector< std::unordered_map<int,int> > acount;
  std::vector< std::unordered_map<int,int> > mho;
  
  auto get_pop_index = [&](int pop)->int {
    auto it = pop_to_idx.find(pop);
    if (it != pop_to_idx.end()) return it->second;
    int idx = (int)pop_labels.size();
    pop_to_idx[pop] = idx;
    pop_labels.push_back(pop);
    n_i.push_back(0);
    acount.emplace_back();
    mho.emplace_back();
    return idx;
  };
  
  std::unordered_map<int,char> all_alleles;
  all_alleles.reserve(128);
  
  int total_n = 0;
  int total_het = 0;
  std::unordered_map<int,int> pooled_acount;
  pooled_acount.reserve(256);
  
  for (int i = 0; i < n_used; ++i) {
    int a1, a2;
    const int gt = mat_get_colmajor(dat_ptr, ld, i, locus_col);
    if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
    
    const int pop = mat_get_colmajor(dat_ptr, ld, i, pop_col);
    const int pi  = get_pop_index(pop);
    
    n_i[pi] += 1;
    total_n += 1;
    
    acount[pi][a1] += 1; acount[pi][a2] += 1;
    pooled_acount[a1] += 1; pooled_acount[a2] += 1;
    all_alleles[a1] = 1; all_alleles[a2] = 1;
    
    if (a1 != a2) {
      total_het += 1;
      mho[pi][a1] += 1;
      mho[pi][a2] += 1;
    }
  }
  
  const int r = (int)pop_labels.size();
  if (r < 2) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  double n_total = 0.0, sum_n2 = 0.0;
  for (int i = 0; i < r; ++i) { n_total += n_i[i]; sum_n2 += (double)n_i[i] * (double)n_i[i]; }
  if (n_total <= 1.0) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  const double n_bar = n_total / (double)r;
  const double n_c   = (n_total - (sum_n2 / n_total)) / (double)(r - 1);
  if (n_c <= 0.0 || n_bar <= 1.0) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  double a_sum = 0.0, b_sum = 0.0, c_sum = 0.0;
  
  for (const auto& kv : all_alleles) {
    const int allele = kv.first;
    
    double pbar_num = 0.0;
    for (int i = 0; i < r; ++i) {
      if (n_i[i] <= 0) continue;
      double pi = 0.0;
      auto it = acount[i].find(allele);
      if (it != acount[i].end()) pi = (double)it->second / (2.0 * (double)n_i[i]);
      pbar_num += (double)n_i[i] * pi;
    }
    const double pbar = pbar_num / n_total;
    
    double s2_num = 0.0;
    for (int i = 0; i < r; ++i) {
      if (n_i[i] <= 0) continue;
      double pi = 0.0;
      auto it = acount[i].find(allele);
      if (it != acount[i].end()) pi = (double)it->second / (2.0 * (double)n_i[i]);
      s2_num += (double)n_i[i] * (pi - pbar) * (pi - pbar);
    }
    const double s2 = s2_num / ((double)(r - 1) * n_bar);
    
    double hbar = 0.0;
    for (int i = 0; i < r; ++i) {
      auto it = mho[i].find(allele);
      if (it != mho[i].end()) hbar += (double)it->second;
    }
    hbar /= n_total;
    
    const double term = pbar * (1.0 - pbar);
    
    const double a = (n_bar / n_c) *
      (s2 - (1.0 / (n_bar - 1.0)) *
      (term - ((double)(r - 1) / (double)r) * s2 - 0.25 * hbar));
    
    const double b = (n_bar / (n_bar - 1.0)) *
      (term - ((double)(r - 1) / (double)r) * s2 - ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * hbar);
    
    const double c = 0.5 * hbar;
    
    a_sum += a; b_sum += b; c_sum += c;
  }
  
  const double denom = a_sum + b_sum + c_sum;
  if (!(denom > 0.0)) return {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false};
  
  const double Ho = (total_n > 0) ? ((double)total_het / (double)total_n) : NA_REAL;
  
  // Nei HS: unweighted mean of per-population unbiased gene diversities.
  // Each population contributes equally regardless of sample size.
  //
  // Unbiased estimator of Nei & Chesser (1983) (eq. 7.39 in Nei 1987), the one
  // used by FSTAT and GENEPOP:
  //     Hs_pop = n/(n-1) * ( 1 - sum(p^2) - Ho_pop/(2n) )
  // where n = number of genotyped individuals in the population and Ho_pop is
  // its observed heterozygosity. The previous simplified form
  // (2n/(2n-1)) * (1 - sum(p^2)) has no Ho term and reproduces GENETIX, NOT
  // FSTAT/GENEPOP (differences up to ~0.003 for samples of ~25 individuals).
  // This is the SAME formula as in nei_het_stats_cpp(); keep both in sync.
  double Hs_num = 0.0, Hs_den = 0.0;
  for (int pi = 0; pi < r; ++pi) {
    const int ni = n_i[pi];
    if (ni <= 1) continue;  // need n > 1 for unbiased estimate
    const double n_pi   = (double)ni;
    const double denom2 = 2.0 * n_pi;
    double sum_p2 = 0.0;
    for (auto &kv2 : acount[pi]) {
      const double pk = (double)kv2.second / denom2;
      sum_p2 += pk * pk;
    }
    // Heterozygote count of this population: each heterozygous individual
    // added 1 to mho[pi] for each of its two alleles, hence the factor 1/2.
    double het_pop = 0.0;
    for (auto &kv3 : mho[pi]) het_pop += (double)kv3.second;
    het_pop *= 0.5;
    const double Ho_pi = het_pop / n_pi;
    const double Hpop  = (n_pi / (n_pi - 1.0)) * (1.0 - sum_p2 - Ho_pi / denom2);
    Hs_num += Hpop;  // equal weight per population (NOT n_i-weighted)
    Hs_den += 1.0;
  }
  const double Hs = (Hs_den > 0.0) ? (Hs_num / Hs_den) : NA_REAL;
  
  // Nei HT: based on arithmetic mean of per-population allele frequencies.
  // Each population contributes equally regardless of sample size.
  double Ht = NA_REAL;
  {
    int r_active = 0;
    for (int pi = 0; pi < r; ++pi) if (n_i[pi] > 0) ++r_active;
    
    if (r_active > 0 && total_n > 0) {
      // arithmetic mean allele frequency: p̄_a = (1/r) * Σ_i p_{i,a}
      std::unordered_map<int,double> pbar;
      pbar.reserve(128);
      for (int pi = 0; pi < r; ++pi) {
        if (n_i[pi] <= 0) continue;
        const double denom2 = 2.0 * (double)n_i[pi];
        for (auto &kv : acount[pi]) {
          pbar[kv.first] += (double)kv.second / denom2;
        }
      }
      double sum_p2 = 0.0;
      for (auto &kv : pbar) {
        const double p = kv.second / (double)r_active;
        sum_p2 += p * p;
      }
      Ht = 1.0 - sum_p2;
      // unbiased correction using total n across all individuals
      const double denom_tot = 2.0 * (double)total_n;
      if (denom_tot > 1.0) Ht *= (denom_tot / (denom_tot - 1.0));
    }
  }
  return {a_sum, b_sum, c_sum, Ho, Hs, Ht, true};
}

// ============================================================================
// Derived pointer-based helpers (FST per locus + ratio-of-sums overall)
// ============================================================================
static inline void wc84_fst_hs_ht_all_loci_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int pop_col, int missing_code, int base,
    double* fst_out, double* hs_out, double* ht_out
) {
  int out_col = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat_ptr, n_used, ld, p, j, pop_col, missing_code, base);
    if (!s.ok) {
      fst_out[out_col] = NA_REAL;
      hs_out[out_col]  = NA_REAL;
      ht_out[out_col]  = NA_REAL;
    } else {
      const double denom = s.a + s.b + s.c;
      fst_out[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
      hs_out[out_col]  = s.Hs;
      ht_out[out_col]  = s.Ht;
    }
    out_col++;
  }
}

static inline double fst_overall_ratio_of_sums_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int pop_col, int missing_code, int base
) {
  double A = 0.0, D = 0.0;
  
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat_ptr, n_used, ld, p, j, pop_col, missing_code, base);
    if (!s.ok) continue;
    
    const double denom = s.a + s.b + s.c;
    if (!(denom > 0.0)) continue;
    
    A += s.a;
    D += denom;
  }
  
  return (D > 0.0) ? (A / D) : NA_REAL;
}

// ============================================================================
// G statistic (log-likelihood ratio) for differentiation — FSTAT's test stat
// G = 2 * Σ_pop Σ_allele O_{j,a} * ln(O_{j,a} / E_{j,a})
// where E_{j,a} = 2*n_j * (total_a / 2*N)  (expected under panmixia)
// Returns sum of G across all loci (the global test statistic).
// ============================================================================
static inline double g_stat_differentiation_locus_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int locus_col, int pop_col, int missing_code, int base
) {
  std::unordered_map<int,int> pop_to_idx;
  pop_to_idx.reserve(32);
  std::vector<int> pop_labels;
  std::vector<int> n_i;
  std::vector<std::unordered_map<int,int>> acnt;
  
  auto get_pop_idx = [&](int pop)->int {
    auto it = pop_to_idx.find(pop);
    if (it != pop_to_idx.end()) return it->second;
    int idx = (int)pop_labels.size();
    pop_to_idx[pop] = idx;
    pop_labels.push_back(pop);
    n_i.push_back(0);
    acnt.emplace_back();
    return idx;
  };
  
  std::unordered_map<int,int> total_acnt;
  total_acnt.reserve(64);
  int total_n = 0;
  
  for (int i = 0; i < n_used; ++i) {
    int a1, a2;
    const int gt = mat_get_colmajor(dat_ptr, ld, i, locus_col);
    if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
    const int pop = mat_get_colmajor(dat_ptr, ld, i, pop_col);
    const int pi  = get_pop_idx(pop);
    ++n_i[pi]; ++total_n;
    ++acnt[pi][a1]; ++acnt[pi][a2];
    ++total_acnt[a1]; ++total_acnt[a2];
  }
  
  const int r = (int)pop_labels.size();
  if (r < 2 || total_n < 2) return NA_REAL;
  
  const double two_N = 2.0 * (double)total_n;
  double G = 0.0;
  
  for (int pi = 0; pi < r; ++pi) {
    if (n_i[pi] <= 0) continue;
    const double two_nj = 2.0 * (double)n_i[pi];
    for (auto &kv : acnt[pi]) {
      const int allele = kv.first;
      const double O = (double)kv.second;
      if (O <= 0.0) continue;
      auto it = total_acnt.find(allele);
      if (it == total_acnt.end()) continue;
      const double E = two_nj * ((double)it->second / two_N);
      if (E <= 0.0) continue;
      G += O * std::log(O / E);
    }
  }
  return 2.0 * G;
}

static inline double g_stat_all_loci_ptr_ld(
    const int* dat_ptr, int n_used, int ld, int p,
    int pop_col, int missing_code, int base
) {
  double G_total = 0.0;
  bool any_ok = false;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    const double g = g_stat_differentiation_locus_ptr_ld(
      dat_ptr, n_used, ld, p, j, pop_col, missing_code, base);
    if (std::isfinite(g)) { G_total += g; any_ok = true; }
  }
  return any_ok ? G_total : NA_REAL;
}


static inline void locus_names_from_colnames(const IntegerMatrix& dat, int pop_col, CharacterVector& out);

// nei_het_stats_cpp() — per-locus Ho, Hs, Ht.
//
// Inputs:
//   dat            packed-genotype matrix: column `pop_col_1based` (1-based)
//                  holds the population code, every other column is one
//                  locus with genotype = allele1*base + allele2 (both
//                  1-based; `missing_code` marks a missing genotype).
//   pop_col_1based 1-based index of the population column (default 1, i.e.
//                  the first column).
//   missing_code   genotype value treated as missing (default 0).
//   base           the packing base used to encode allele1/allele2 into one
//                  integer (must match how the R side built `dat`).
//
// Returns a List with one numeric vector per statistic, one entry per locus
// (locus order = column order of `dat`, excluding the population column):
//   Ho   observed heterozygosity, unweighted mean of each population's
//        (heterozygote count / genotyped count).
//   Hs   expected heterozygosity within populations — Nei & Chesser (1983)
//        UNBIASED estimator: for each population p with n_p diploid
//        individuals, Hs_p = n_p/(n_p-1) * (1 - sum(allele freq^2) -
//        Ho_p/(2*n_p)), then Hs = unweighted mean of Hs_p over populations.
//        The "- Ho_p/(2*n_p)" term matters: a simpler (2n/(2n-1))*(1-sum(p^2))
//        form (no Ho term) is what GENETIX uses and is NOT what FSTAT/
//        GENEPOP report — verified against FSTAT reference output; do not
//        remove this term without re-checking against both programs.
//   Ht   total expected heterozygosity, from the unweighted mean allele
//        frequency across populations (also unbiased-corrected).
// [[Rcpp::export]]
Rcpp::List nei_het_stats_cpp(const Rcpp::IntegerMatrix& dat,
                             int pop_col_1based = 1,
                             int missing_code = 0,
                             int base = 1000) {
  
  if (dat.ncol() < 2)
    Rcpp::stop("Need at least 2 columns: pop + >=1 locus.");
  
  if (base <= 1)
    Rcpp::stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  
  if (pop_col < 0 || pop_col >= dat.ncol())
    Rcpp::stop("pop_col out of range.");
  
  const int N = dat.nrow();
  const int p = dat.ncol();
  const int L = p - 1;
  
  Rcpp::NumericVector Ho(L, NA_REAL);
  Rcpp::NumericVector Hs(L, NA_REAL);
  Rcpp::NumericVector Ht(L, NA_REAL);
  
  Rcpp::IntegerVector k_eff(L);
  std::fill(k_eff.begin(), k_eff.end(), NA_INTEGER);
  
  // map pop label -> index
  std::unordered_map<int,int> pop_to_idx;
  pop_to_idx.reserve(64);
  
  std::vector<int> pop_labels;
  pop_labels.reserve(16);
  
  auto get_pop_index = [&](int pop)->int {
    auto it = pop_to_idx.find(pop);
    if (it != pop_to_idx.end()) return it->second;
    
    int idx = (int)pop_labels.size();
    pop_to_idx[pop] = idx;
    pop_labels.push_back(pop);
    return idx;
  };
  
  // pre-scan populations
  for (int i = 0; i < N; ++i) {
    int pop = dat(i, pop_col);
    get_pop_index(pop);
  }
  
  const int r = (int)pop_labels.size();
  
  if (r < 1)
    return Rcpp::List::create(_["Ho"]=Ho, _["Hs"]=Hs, _["Ht"]=Ht);
  
  int out_col = 0;
  
  for (int j = 0; j < p; ++j) {
    
    if (j == pop_col) continue;
    
    std::vector< std::unordered_map<int,int> > acount(r);
    std::vector<int> n_i(r, 0);
    std::vector<int> het_i(r, 0);
    
    std::unordered_map<int,int> total_acount;
    total_acount.reserve(256);
    
    int total_n = 0;
    int total_het = 0;
    
    for (int i = 0; i < N; ++i) {
      
      int a1, a2;
      int gt = dat(i, j);
      
      if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base))
        continue;
      
      int pop = dat(i, pop_col);
      int pi = pop_to_idx[pop];
      
      n_i[pi] += 1;
      total_n += 1;
      
      acount[pi][a1] += 1;
      acount[pi][a2] += 1;
      
      total_acount[a1] += 1;
      total_acount[a2] += 1;
      
      if (a1 != a2) {
        het_i[pi] += 1;
        total_het += 1;
      }
    }
    
    if (total_n == 0) {
      out_col++;
      continue;
    }
    
    int k_locus = 0;
    for (int pi = 0; pi < r; ++pi) {
      if (n_i[pi] > 0) k_locus++;
    }
    
    k_eff[out_col] = k_locus;
    
    // Ho = mean heterozygosity across populations
    double Ho_sum = 0.0;
    double Ho_den = 0.0;
    
    for (int pi = 0; pi < r; ++pi) {
      if (n_i[pi] <= 0) continue;
      
      Ho_sum += (double)het_i[pi] / (double)n_i[pi];
      Ho_den += 1.0;
    }
    
    Ho[out_col] = (Ho_den > 0.0) ? (Ho_sum / Ho_den) : NA_REAL;
    
    // Hs
    double Hs_num = 0.0;
    double Hs_den = 0.0;
    
    for (int pi = 0; pi < r; ++pi) {
      
      if (n_i[pi] <= 1) continue;
      
      const double n_pi = (double)n_i[pi];
      const double denom = 2.0 * n_pi;
      
      double sum_p2 = 0.0;
      
      for (auto &kv : acount[pi]) {
        
        const double pk = (double)kv.second / denom;
        sum_p2 += pk * pk;
      }
      
      const double Ho_pi = (double)het_i[pi] / n_pi;
      
      // Nei & Chesser (1983) unbiased Hs, eq. 7.39 in Nei (1987):
      //   Hs = n/(n-1) * (1 - sum(p^2) - Ho/(2n))
      // (the previous simpler (2n/(2n-1))*(1-sum(p^2)) form omitted the
      // Ho correction term, matching GENETIX's convention rather than
      // FSTAT's/GENEPOP's.)
      double Hpop = (n_pi / (n_pi - 1.0)) * (1.0 - sum_p2 - Ho_pi / denom);
      
      Hs_num += Hpop;
      Hs_den += 1.0;
    }
    
    if (Hs_den > 0.0)
      Hs[out_col] = Hs_num / Hs_den;
    
    // Ht
      {
        int k_eff = 0;
        
        for (int pi = 0; pi < r; ++pi) {
          if (n_i[pi] > 0) k_eff++;
        }

          if (k_eff <= 0) {
            
            Ht[out_col] = NA_REAL;
            
          } else {
            
            std::unordered_map<int,double> p_sum;
            p_sum.reserve(256);
            
            for (int pi = 0; pi < r; ++pi) {
              
              if (n_i[pi] <= 0) continue;
              
              const double denom_pi = 2.0 * (double)n_i[pi];
              
              for (auto &kv : acount[pi]) {
                
                const int allele = kv.first;
                const double p_aj = (double)kv.second / denom_pi;
                
                p_sum[allele] += p_aj;
              }
            }
            
            double sum_p2 = 0.0;
            
            for (auto &kv : p_sum) {
              
              const double pbar = kv.second / (double)k_eff;
              sum_p2 += pbar * pbar;
            }
            
            double Htot = 1.0 - sum_p2;
            
            const double denom_tot = 2.0 * (double)total_n;
            
            if (denom_tot > 1.0)
              Htot *= (denom_tot / (denom_tot - 1.0));
            
            Ht[out_col] = Htot;
          }
      }
      
      out_col++;
  }
  
  // attach locus names
  Rcpp::CharacterVector cn = colnames(dat);
  Rcpp::CharacterVector loc_names(L);
  
  int kk = 0;
  const bool has_colnames = (cn.size() == p);
  
  for (int j = 0; j < p; ++j) {
    
    if (j == pop_col) continue;
    
    if (has_colnames)
      loc_names[kk] = cn[j];
    else
      loc_names[kk] = std::string("L") + std::to_string(kk + 1);
    
    kk++;
  }
  
  Ho.attr("names") = loc_names;
  Hs.attr("names") = loc_names;
  Ht.attr("names") = loc_names;
  
  return Rcpp::List::create(
    _["Ho"] = Ho,
    _["Hs"] = Hs,
    _["Ht"] = Ht,
    _["k_eff"] = k_eff,
    _["locus_names"] = loc_names
  );
}

// ============================================================================
// R-FREE OpenMP permutation kernel (FST only)
// ============================================================================
static void fst_perm_kernel_parallel_rfree(
    const int* dat_ptr, int n, int p,
    int pop_col, int missing_code, int base,
    int perm_mode,
    int B, int n_threads, uint64_t seed0,
    const std::vector<std::vector<int>>& strata,
    std::vector<double>& fst_perm_buf,
    std::vector<double>& fst_overall_buf,
    std::vector<double>& g_overall_buf,  // B — global G per permutation
    std::vector<double>& g_perm_buf      // B*L — per-locus G per permutation
) {
  const int L = p - 1;

  int T = std::max(1, n_threads);
#ifdef _OPENMP
  omp_set_num_threads(T);
#else
  T = 1;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<int> pm((size_t)n * (size_t)p);
  std::vector<double> fst_tmp((size_t)L, NA_REAL);
  std::vector<double> g_tmp((size_t)L, NA_REAL);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = seed0 + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    std::mt19937_64 rng(sb);

    std::memcpy(pm.data(), dat_ptr, sizeof(int) * (size_t)n * (size_t)p);

    if (perm_mode == 1) {
      permute_pop_labels_once_rng_buf(pm.data(), n, p, pop_col, rng);
    } else {
      permute_within_pops_once_rng_buf(pm.data(), n, p, pop_col, missing_code, base, strata, rng);
    }

    int out_col = 0;
    for (int j = 0; j < p; ++j) {
      if (j == pop_col) continue;
      WC84LocusStats s = wc84_locus_stats_ptr_ld(pm.data(), n, n, p, j, pop_col, missing_code, base);
      if (!s.ok) {
        fst_tmp[out_col] = NA_REAL;
        g_tmp[out_col]   = NA_REAL;
      } else {
        const double denom = s.a + s.b + s.c;
        fst_tmp[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
        g_tmp[out_col]   = g_stat_differentiation_locus_ptr_ld(pm.data(), n, n, p, j, pop_col, missing_code, base);
      }
      out_col++;
    }

    const double fst_overall = fst_overall_ratio_of_sums_ptr_ld(pm.data(), n, n, p, pop_col, missing_code, base);
    const double g_overall = g_stat_all_loci_ptr_ld(pm.data(), n, n, p, pop_col, missing_code, base);
    g_overall_buf[(size_t)b] = g_overall;
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) {
      fst_perm_buf[row0 + (size_t)ell] = fst_tmp[(size_t)ell];
      g_perm_buf[row0 + (size_t)ell]   = g_tmp[(size_t)ell];
    }
    fst_overall_buf[(size_t)b] = fst_overall;
  }
}
}

// ============================================================================
// R-FREE OpenMP pop-block bootstrap kernel (FST + HS + HT)
// ============================================================================
static void boot_popblock_kernel_parallel_rfree(
    const int* dat_ptr, int n, int p, int pop_col,
    int missing_code, int base,
    int B, int n_threads, uint64_t seed0,
    const std::vector<std::vector<int>>& base_rows, // population blocks (indices in original)
    std::vector<double>& fst_boot_buf,              // B*L
    std::vector<double>& fst_overall_buf,           // B
    std::vector<double>& hs_boot_buf,               // B*L
    std::vector<double>& hs_overall_buf,            // B
    std::vector<double>& ht_boot_buf,               // B*L
    std::vector<double>& ht_overall_buf             // B
) {
  const int L = p - 1;
  const int P = (int)base_rows.size();
  if (P <= 0) return;
  
  int max_block = 0;
  for (const auto& rows : base_rows) max_block = std::max(max_block, (int)rows.size());
  const int ld = std::max(1, max_block * P);
  
  int T = std::max(1, n_threads);
#ifdef _OPENMP
  omp_set_num_threads(T);
#else
  T = 1;
#endif
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<int> mbuf((size_t)ld * (size_t)p);
  std::vector<double> fst_tmp((size_t)L), hs_tmp((size_t)L), ht_tmp((size_t)L);
  
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = seed0 + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    std::mt19937_64 rng(sb);
    std::uniform_int_distribution<int> U(0, P - 1);
    
    int rr = 0;
    for (int k = 0; k < P; ++k) {
      const int pick = U(rng);
      const auto& rows = base_rows[(size_t)pick];
      const int block_id = k + 1;
      
      for (int rid : rows) {
        for (int j = 0; j < p; ++j) {
          const int v = mat_get_colmajor(dat_ptr, n, rid, j);
          mat_set_colmajor(mbuf.data(), ld, rr, j, v);
        }
        mat_set_colmajor(mbuf.data(), ld, rr, pop_col, block_id);
        rr++;
        if (rr >= ld) break;
      }
      if (rr >= ld) break;
    }
    const int n_used = rr;
    
    wc84_fst_hs_ht_all_loci_ptr_ld(
      mbuf.data(), n_used, ld, p,
      pop_col, missing_code, base,
      fst_tmp.data(), hs_tmp.data(), ht_tmp.data()
    );
    
    const double fst_overall = fst_overall_ratio_of_sums_ptr_ld(mbuf.data(), n_used, ld, p, pop_col, missing_code, base);
    const double hs_overall  = mean_across_loci_ptr(hs_tmp.data(), L);
    const double ht_overall  = mean_across_loci_ptr(ht_tmp.data(), L);
    
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) {
      fst_boot_buf[row0 + (size_t)ell] = fst_tmp[(size_t)ell];
      hs_boot_buf [row0 + (size_t)ell] = hs_tmp [(size_t)ell];
      ht_boot_buf [row0 + (size_t)ell] = ht_tmp [(size_t)ell];
    }
    fst_overall_buf[(size_t)b] = fst_overall;
    hs_overall_buf [(size_t)b] = hs_overall;
    ht_overall_buf [(size_t)b] = ht_overall;
  }
}
}

// ============================================================================
// Export: observed WC84 stats (per-locus + overall) using pointer core
// ============================================================================
static inline void locus_names_from_colnames(const IntegerMatrix& dat, int pop_col, CharacterVector& out) {
  const int p = dat.ncol();
  const int L = p - 1;
  CharacterVector cn = colnames(dat);
  
  out = CharacterVector(L);
  int k = 0;
  const bool has = (cn.size() == p);
  
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    if (has) {
      out[k] = cn[j];                        // same type: Rcpp::String proxy
    } else {
      out[k] = Rcpp::String(std::string("L") + std::to_string(k + 1));
    }
    ++k;
  }
}

// [[Rcpp::export]]
List observed_wc84_stats_cpp(const IntegerMatrix& dat,
                             int pop_col_1based = 1,
                             int missing_code = 0,
                             int base = 1000) {
  
  if (dat.ncol() < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (base <= 1) stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= dat.ncol()) stop("pop_col out of range.");
  
  const int n = dat.nrow();
  const int p = dat.ncol();
  const int L = p - 1;
  
  NumericVector FST(L, NA_REAL), FIT(L, NA_REAL), FIS(L, NA_REAL);
  NumericVector Ho(L, NA_REAL), HS(L, NA_REAL), HT(L, NA_REAL);
  
  double A_sum = 0.0, B_sum = 0.0, C_sum = 0.0;
  
  int out_col = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat.begin(), n, n, p, j, pop_col, missing_code, base);
    if (s.ok) {
      const double denom = s.a + s.b + s.c;
      const double denom_fis = s.b + s.c;
      
      FST[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
      FIT[out_col] = (denom > 0.0) ? ((s.a + s.b) / denom) : NA_REAL;
      FIS[out_col] = (denom_fis > 0.0) ? (s.b / denom_fis) : NA_REAL;
      
      Ho[out_col]  = s.Ho;
      HS[out_col]  = s.Hs;
      HT[out_col]  = s.Ht;
      
      if (denom > 0.0) {
        A_sum += s.a;
        B_sum += s.b;
        C_sum += s.c;
      }
    }
    out_col++;
  }
  
  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);
  
  FST.attr("names") = loc_names;
  FIT.attr("names") = loc_names;
  FIS.attr("names") = loc_names;
  Ho.attr("names")  = loc_names;
  HS.attr("names")  = loc_names;
  HT.attr("names")  = loc_names;
  
  const double denom_overall = A_sum + B_sum + C_sum;
  const double denom_fis_overall = B_sum + C_sum;
  
  const double FST_overall = (denom_overall > 0.0) ? (A_sum / denom_overall) : NA_REAL;
  const double FIT_overall = (denom_overall > 0.0) ? ((A_sum + B_sum) / denom_overall) : NA_REAL;
  const double FIS_overall = (denom_fis_overall > 0.0) ? (B_sum / denom_fis_overall) : NA_REAL;
  
  return List::create(
    _["FST"] = FST,
    _["FIT"] = FIT,
    _["FIS"] = FIS,
    _["HI"]  = Ho,
    _["HS"]  = HS,
    _["HT"]  = HT,
    _["FST_overall_ratio_of_sums"] = FST_overall,
    _["FIT_overall_ratio_of_sums"] = FIT_overall,
    _["FIS_overall_ratio_of_sums"] = FIS_overall,
    _["locus_names"] = loc_names
  );
}

// ============================================================================
// Export: per-locus WC84 variance components for locus bootstrap
// Returns a data frame with columns A, B, C, HS, HT per locus.
// The caller resamples rows (loci) with replacement and recomputes
// ratio-of-sums to obtain the locus-bootstrap distribution.
// ============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame wc84_locus_components_cpp(
    const Rcpp::IntegerMatrix& dat,
    int pop_col_1based = 1,
    int missing_code = 0,
    int base = 1000
) {
  if (dat.ncol() < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (base <= 1) stop("base must be > 1.");

  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= dat.ncol()) stop("pop_col out of range.");

  const int n = dat.nrow();
  const int p = dat.ncol();
  const int L = p - 1;

  NumericVector A_v(L, NA_REAL), B_v(L, NA_REAL), C_v(L, NA_REAL);
  NumericVector HS_v(L, NA_REAL), HT_v(L, NA_REAL);

  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);

  int out = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    WC84LocusStats s = wc84_locus_stats_ptr_ld(
      dat.begin(), n, n, p, j, pop_col, missing_code, base);
    if (s.ok) {
      A_v[out]  = s.a;
      B_v[out]  = s.b;
      C_v[out]  = s.c;
      HS_v[out] = s.Hs;
      HT_v[out] = s.Ht;
    }
    ++out;
  }

  return DataFrame::create(
    _["locus"] = loc_names,
    _["A"]     = A_v,
    _["B"]     = B_v,
    _["C"]     = C_v,
    _["HS"]    = HS_v,
    _["HT"]    = HT_v
  );
}

// ============================================================================
// Export: parallel permutation (FST)
// ============================================================================
static inline uint64_t seed0_from_double(double seed) {
  const uint64_t s = (uint64_t)(seed < 0 ? -seed : seed);
  return s + 0x9e3779b97f4a7c15ULL;
}

// [[Rcpp::export]]
List batch_permute_wc84_fst_parallel(const IntegerMatrix& dat,
                                           int pop_col_1based,
                                           int missing_code,
                                           int base,
                                           int B,
                                           int n_threads = 1,
                                           double seed = 1.0,
                                           std::string pval_method = "two_sided_abs",
                                           std::string perm_scheme = "within_pop_alleles") {
  
  const int n = dat.nrow();
  const int p = dat.ncol();
  if (p < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (B <= 0) stop("B must be positive.");
  if (base <= 1) stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= p) stop("pop_col out of range.");
  
  const int L = p - 1;
  
  int perm_mode = 0;
  std::vector<std::vector<int>> strata;
  if (perm_scheme == "within_pop_alleles") {
    perm_mode = 0;
    strata = build_pop_rows(dat, pop_col);
  } else if (perm_scheme == "permute_pop_labels") {
    perm_mode = 1;
  } else {
    stop("Unknown perm_scheme.");
  }
  
  // observed (pointer core)
  NumericVector fst_obs(L, NA_REAL);
  NumericVector g_obs(L, NA_REAL);
  int out_col = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    WC84LocusStats s = wc84_locus_stats_ptr_ld(dat.begin(), n, n, p, j, pop_col, missing_code, base);
    if (s.ok) {
      const double denom = s.a + s.b + s.c;
      fst_obs[out_col] = (denom > 0.0) ? (s.a / denom) : NA_REAL;
      g_obs[out_col]   = g_stat_differentiation_locus_ptr_ld(dat.begin(), n, n, p, j, pop_col, missing_code, base);
    }
    out_col++;
  }
  const double fst_overall_obs = fst_overall_ratio_of_sums_ptr_ld(dat.begin(), n, n, p, pop_col, missing_code, base);
  const double g_overall_obs = g_stat_all_loci_ptr_ld(dat.begin(), n, n, p, pop_col, missing_code, base);
  int T = std::max(1, n_threads);
#ifndef _OPENMP
  T = 1;
#endif
  
  const uint64_t seed0 = seed0_from_double(seed);
  
  std::vector<double> fst_perm_buf((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> fst_overall_buf((size_t)B, NA_REAL);
  std::vector<double> g_overall_buf((size_t)B, NA_REAL);
  std::vector<double> g_perm_buf((size_t)B * (size_t)L, NA_REAL);

  fst_perm_kernel_parallel_rfree(
    dat.begin(), n, p,
    pop_col, missing_code, base,
    perm_mode,
    B, T, seed0,
    strata,
    fst_perm_buf,
    fst_overall_buf,
    g_overall_buf,
    g_perm_buf
  );
  
  // materialize R outputs
  NumericMatrix fst_perm(B, L);
  NumericVector fst_overall_perm(B, NA_REAL);
  for (int b = 0; b < B; ++b) {
    fst_overall_perm[b] = fst_overall_buf[(size_t)b];
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) fst_perm(b, ell) = fst_perm_buf[row0 + (size_t)ell];
  }
  
  // FST-based per-locus p-values (kept for reference)
  NumericVector p_fst(L, NA_REAL);
  for (int ell = 0; ell < L; ++ell) {
    const double obs = fst_obs[ell];
    if (NumericVector::is_na(obs)) continue;

    int count = 0, valid = 0;
    for (int b = 0; b < B; ++b) {
      const double v = fst_perm(b, ell);
      if (NumericVector::is_na(v)) continue;
      valid++;

      if (pval_method == "greater") {
        if (v >= obs - TIE_TOL) count++;
      } else if (pval_method == "less") {
        if (v <= obs + TIE_TOL) count++;
      } else {
        if (std::fabs(v) >= std::fabs(obs) - TIE_TOL) count++;
      }
    }
    if (valid > 0) p_fst[ell] = (1.0 + count) / ((double)valid + 1.0);
  }

  // G-based per-locus p-values (one-sided: G >= 0 by definition)
  NumericMatrix g_perm(B, L);
  for (int b = 0; b < B; ++b) {
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) g_perm(b, ell) = g_perm_buf[row0 + (size_t)ell];
  }

  NumericVector p_G(L, NA_REAL);
  for (int ell = 0; ell < L; ++ell) {
    const double obs_g = g_obs[ell];
    if (NumericVector::is_na(obs_g) || !std::isfinite(obs_g)) continue;
    int count = 0, valid = 0;
    for (int b = 0; b < B; ++b) {
      const double v = g_perm(b, ell);
      if (NumericVector::is_na(v) || !std::isfinite(v)) continue;
      ++valid;
      if (v >= obs_g - TIE_TOL) ++count;
    }
    if (valid > 0) p_G[ell] = (1.0 + count) / ((double)valid + 1.0);
  }

  // Overall p-value uses G statistic (FSTAT convention), not theta.
  double p_fst_overall = NA_REAL;
  if (std::isfinite(g_overall_obs)) {
    int count = 0, valid = 0;
    for (int b = 0; b < B; ++b) {
      const double v = g_overall_buf[(size_t)b];
      if (!std::isfinite(v)) continue;
      ++valid;
      if (v >= g_overall_obs - TIE_TOL) ++count;  // one-sided: G_perm >= G_obs (ties included)
    }
    if (valid > 0) p_fst_overall = (1.0 + count) / ((double)valid + 1.0);
  }
  
  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);
  colnames(fst_perm) = loc_names;
  colnames(g_perm)   = loc_names;
  fst_obs.attr("names") = loc_names;
  g_obs.attr("names")   = loc_names;
  p_fst.attr("names")   = loc_names;
  p_G.attr("names")     = loc_names;
  
  NumericVector g_overall_perm(B, NA_REAL);  //check place
  for (int b = 0; b < B; ++b) g_overall_perm[b] = g_overall_buf[(size_t)b]; //check place
  
  return List::create(
    _["FST_obs"]         = fst_obs,
    _["FST_overall_obs"] = fst_overall_obs,
    _["FST_perm"]        = fst_perm,
    _["FST_overall_perm"]= fst_overall_perm,
    _["G_locus_obs"]     = g_obs,
    _["G_locus_perm"]    = g_perm,
    _["G_overall_obs"]   = g_overall_obs,
    _["G_overall_perm"]  = g_overall_perm,
    _["p_FST"]           = p_fst,
    _["p_G"]             = p_G,
    _["p_FST_overall"]   = p_fst_overall,
    _["locus_names"]     = loc_names,
    _["perm_scheme"]     = perm_scheme,
    _["pval_method"]     = pval_method,
    _["n_perm"]          = B,
    _["n_threads"]       = T,
    _["seed"]            = seed
  );
}

// ============================================================================
// Export: parallel pop-block bootstrap (FST + HS + HT)
// ============================================================================

// [[Rcpp::export]]
List boot_popblock_wc84_parallel(const IntegerMatrix& mat,
                                       int pop_col_1based = 1,
                                       int missing_code = 0,
                                       int base = 1000,
                                       int B = 1000,
                                       int n_threads = 1,
                                       double seed = 1.0) {
  
  const int p = mat.ncol();
  if (p < 2) stop("Need at least 2 columns: pop + >=1 locus.");
  if (B <= 0) stop("B must be positive.");
  if (base <= 1) stop("base must be > 1.");
  
  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= p) stop("pop_col out of range.");
  
  const int n = mat.nrow();
  const int L = p - 1;
  
  auto base_rows = build_pop_rows(mat, pop_col);
  
  // observed from pointer core
  std::vector<double> fst_obs_buf((size_t)L, NA_REAL), hs_obs_buf((size_t)L, NA_REAL), ht_obs_buf((size_t)L, NA_REAL);
  wc84_fst_hs_ht_all_loci_ptr_ld(mat.begin(), n, n, p, pop_col, missing_code, base,
                                 fst_obs_buf.data(), hs_obs_buf.data(), ht_obs_buf.data());
  
  NumericVector fst_obs(L), hs_obs(L), ht_obs(L);
  for (int i = 0; i < L; ++i) { fst_obs[i] = fst_obs_buf[i]; hs_obs[i] = hs_obs_buf[i]; ht_obs[i] = ht_obs_buf[i]; }
  
  const double fst_overall_obs = fst_overall_ratio_of_sums_ptr_ld(mat.begin(), n, n, p, pop_col, missing_code, base);
  const double hs_overall_obs  = mean_across_loci_ptr(hs_obs_buf.data(), L);
  const double ht_overall_obs  = mean_across_loci_ptr(ht_obs_buf.data(), L);
  
  int T = std::max(1, n_threads);
#ifndef _OPENMP
  T = 1;
#endif
  const uint64_t seed0 = seed0_from_double(seed);
  
  // buffers
  std::vector<double> fst_boot_buf((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> hs_boot_buf ((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> ht_boot_buf ((size_t)B * (size_t)L, NA_REAL);
  
  std::vector<double> fst_overall_buf((size_t)B, NA_REAL);
  std::vector<double> hs_overall_buf ((size_t)B, NA_REAL);
  std::vector<double> ht_overall_buf ((size_t)B, NA_REAL);
  
  boot_popblock_kernel_parallel_rfree(
    mat.begin(), n, p, pop_col,
    missing_code, base,
    B, T, seed0,
    base_rows,
    fst_boot_buf, fst_overall_buf,
    hs_boot_buf,  hs_overall_buf,
    ht_boot_buf,  ht_overall_buf
  );
  
  // materialize
  NumericMatrix fst_boot(B, L), hs_boot(B, L), ht_boot(B, L);
  NumericVector fst_overall_boot(B), hs_overall_boot(B), ht_overall_boot(B);
  
  for (int b = 0; b < B; ++b) {
    fst_overall_boot[b] = fst_overall_buf[(size_t)b];
    hs_overall_boot[b]  = hs_overall_buf[(size_t)b];
    ht_overall_boot[b]  = ht_overall_buf[(size_t)b];
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) {
      fst_boot(b, ell) = fst_boot_buf[row0 + (size_t)ell];
      hs_boot (b, ell) = hs_boot_buf [row0 + (size_t)ell];
      ht_boot (b, ell) = ht_boot_buf [row0 + (size_t)ell];
    }
  }
  
  CharacterVector loc_names;
  locus_names_from_colnames(mat, pop_col, loc_names);
  colnames(fst_boot) = loc_names;
  colnames(hs_boot)  = loc_names;
  colnames(ht_boot)  = loc_names;
  
  fst_obs.attr("names") = loc_names;
  hs_obs.attr("names")  = loc_names;
  ht_obs.attr("names")  = loc_names;
  
  return List::create(
    _["FST_obs"] = fst_obs,
    _["FST_overall_obs"] = fst_overall_obs,
    _["FST_boot"] = fst_boot,
    _["FST_overall_boot"] = fst_overall_boot,
    
    _["HS_obs"] = hs_obs,
    _["HS_overall_obs"] = hs_overall_obs,
    _["HS_boot"] = hs_boot,
    _["HS_overall_boot"] = hs_overall_boot,
    
    _["HT_obs"] = ht_obs,
    _["HT_overall_obs"] = ht_overall_obs,
    _["HT_boot"] = ht_boot,
    _["HT_overall_boot"] = ht_overall_boot,
    
    _["boot_type"] = "pop_block",
    _["n_boot"] = B,
    _["n_threads"] = T,
    _["seed"] = seed,
    _["locus_names"] = loc_names
  );
}

// ============================================================================
// Individual bootstrap for HS: resample individuals within each population
// Returns HS_boot (B x L) and HS_overall_boot (B) — loci with replacement
// would give locus-level CI, handled separately in R (pga_boot_over_loci() in utils_stats.R).
// ============================================================================
// [[Rcpp::export]]
List boot_indiv_hs_cpp(
    const IntegerMatrix& dat,
    int pop_col_1based = 1,
    int missing_code   = 0,
    int base           = 1000,
    int B              = 1000,
    double seed        = 1.0,
    int n_threads      = 1
) {
  const int n = dat.nrow();
  const int p = dat.ncol();
  if (p < 2) stop("Need at least 2 columns.");
  if (B <= 0) stop("B must be positive.");
  if (base <= 1) stop("base must be > 1.");

  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= p) stop("pop_col out of range.");
  const int L = p - 1;

  std::vector<std::vector<int>> pop_rows = build_pop_rows(dat, pop_col);
  const int P = (int)pop_rows.size();
  if (P < 1) stop("No populations found.");

  const uint64_t seed0 = seed0_from_double(seed);

  std::vector<double> hs_boot_buf((size_t)B * (size_t)L, NA_REAL);
  std::vector<double> hs_overall_buf((size_t)B, NA_REAL);

  int T = std::max(1, n_threads);
#ifndef _OPENMP
  T = 1;
#endif
#ifdef _OPENMP
  omp_set_num_threads(T);
#pragma omp parallel
#endif
  {
    std::vector<int> mbuf((size_t)n * (size_t)p);
    std::vector<double> fst_tmp((size_t)L, NA_REAL);
    std::vector<double> hs_tmp((size_t)L, NA_REAL);
    std::vector<double> ht_tmp((size_t)L, NA_REAL);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int b = 0; b < B; ++b) {
      const uint64_t sb = seed0 + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
      std::mt19937_64 rng(sb);

      int rr = 0;
      for (int k = 0; k < P; ++k) {
        const auto& rows = pop_rows[(size_t)k];
        const int nk = (int)rows.size();
        if (nk == 0) continue;
        std::uniform_int_distribution<int> U_indiv(0, nk - 1);
        for (int i = 0; i < nk; ++i) {
          const int rid = rows[(size_t)U_indiv(rng)];
          for (int j = 0; j < p; ++j)
            mat_set_colmajor(mbuf.data(), n, rr, j,
                             mat_get_colmajor(dat.begin(), n, rid, j));
          rr++;
        }
      }

      wc84_fst_hs_ht_all_loci_ptr_ld(
        mbuf.data(), rr, n, p,
        pop_col, missing_code, base,
        fst_tmp.data(), hs_tmp.data(), ht_tmp.data()
      );

      const size_t row0 = (size_t)b * (size_t)L;
      for (int ell = 0; ell < L; ++ell) hs_boot_buf[row0 + (size_t)ell] = hs_tmp[(size_t)ell];
      hs_overall_buf[(size_t)b] = mean_across_loci_ptr(hs_tmp.data(), L);
    }
  }

  NumericMatrix hs_boot(B, L);
  NumericVector hs_overall_boot(B, NA_REAL);
  for (int b = 0; b < B; ++b) {
    hs_overall_boot[b] = hs_overall_buf[(size_t)b];
    const size_t row0 = (size_t)b * (size_t)L;
    for (int ell = 0; ell < L; ++ell) hs_boot(b, ell) = hs_boot_buf[row0 + (size_t)ell];
  }

  CharacterVector loc_names;
  locus_names_from_colnames(dat, pop_col, loc_names);
  colnames(hs_boot) = loc_names;

  return List::create(
    _["HS_boot"]         = hs_boot,
    _["HS_overall_boot"] = hs_overall_boot,
    _["locus_names"]     = loc_names,
    _["boot_type"]       = "individual",
    _["n_boot"]          = B
  );
}