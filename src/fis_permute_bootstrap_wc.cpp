// src/fis_permute_bootstrap_wc.cpp
//
// FUNCTIONS IN THIS FILE (in order of appearance) — within-population FIS
// (Weir & Cockerham 1984), used by the Local Panmixia module.
//
//   boot_indiv_wc_fis_by_pop()      Bootstrap over INDIVIDUALS within one
//                                    population, for that population's FIS CI
//                                    (matrix output: rows = replicates).
//   batch_permute_wc_fis_by_pop()   Permutation test, per population.
//   fis_wc_cpp()                    Multilocus FIS (WC84 variance
//                                    components), all populations combined.
//   wc_fis_by_pop_wc84()            Per-population FIS using the TRUE WC84
//                                    variance-component decomposition (this
//                                    is the validated function reused by
//                                    General Stats' overall_by_population
//                                    output and cross-checked against FSTAT).
//   batch_permute_wc_fis()          Permutation test for multilocus FIS.
//   calculate_observed_fis()        Observed per-locus FIS (WC84), no
//                                    resampling.
//   boot_indiv_wc_fis()             Bootstrap over individuals, per locus,
//                                    across all populations.
//   boot_popblock_wc_fis()          Bootstrap over SUB-SAMPLES (populations),
//                                    per locus — needs >=5 populations to be
//                                    meaningful (see the N<5 guard at the R
//                                    call site in server_general_stats.R).
//   summarize_fis_results()         Turns a bootstrap matrix (from either
//                                    boot_indiv_wc_fis or boot_popblock_wc_fis)
//                                    into mean/CI summary fields.
//   create_results_dataframe()      Assembles the final per-locus results
//                                    data.frame from observed + summary lists.
//   wc84_per_allele_fstats_cpp()    Per-ALLELE (not just per-locus) FIT/FIS/
//                                    FST from the allele-specific variance
//                                    components — used by General Stats'
//                                    "Fstats_per_allele" output.
//
// RNG / REPRODUCIBILITY NOTE (added following a code audit; keep in sync
// with fit_permute_bootstrap_wc.cpp, which follows the same scheme, and with
// fst_permute_boostrap_OpenMP.cpp, which uses a DIFFERENT scheme):
//
//   The permutation/bootstrap loops in this file draw randomness with R's
//   own R::unif_rand() (see e.g. the `R::unif_rand()` calls below). This is
//   reproducible from R via set.seed() the same way any base-R random draw
//   is, and is the reason this file explicitly does NOT parallelise those
//   loops with OpenMP (R's RNG state is process-global and not thread-safe -
//   see the "IMPORTANT: do NOT OpenMP-parallelise this loop" comment at the
//   relevant call site). Any future change that adds OpenMP around a loop
//   using R::unif_rand() would silently produce incorrect / non-reproducible
//   results and must not be made without first switching that loop to a
//   thread-local generator (as fst_permute_boostrap_OpenMP.cpp does).
//
//   For the publication: state explicitly that FIS/FIT permutation p-values
//   and bootstrap CIs are reproducible via set.seed() but are NOT computed
//   in parallel, whereas FST bootstrap/permutation (fst_permute_boostrap_
//   OpenMP.cpp) IS parallelised and has a different reproducibility
//   contract (fixed seed AND fixed thread count - see the note in that
//   file). Harmonizing the three estimator families onto one RNG scheme is
//   tracked as a P1 item in the accompanying audit and is intentionally NOT
//   done here, since it would change numerical output and must be validated
//   against reference software (hierfstat/GENEPOP) before merging - see
//   tests/testthat/test-numeric-validation-fst.R.
//
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <cstdint>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// -----------------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------------

/*
 * =============================================================================
 *  fis_wc_cpp.cpp — Weir & Cockerham (1984) FIS for hierfstat-encoded genotypes
 * =============================================================================
 *
 *  Purpose
 *  -------
 *  Compute locus-wise FIS (inbreeding coefficient within populations) following
 *  Weir & Cockerham (1984), using genotypes encoded in the integer format required
 *  by hierfstat:   gt = allele1 * base + allele2.
 *
 *  Input expectations
 *  ------------------
 *  - ndat: IntegerMatrix with:
 *      * column 0: population codes (integers, typically 1..np)
 *      * columns 1..L: loci encoded as a1*base + a2
 *  - base: integer used for decoding (must be the same base used at encoding).
 *
 *  Missing/invalid genotypes
 *  -------------------------
 *  Genotypes are treated as missing if:
 *    - gt is NA or <= 0, or
 *    - decoded alleles are not strictly positive.
 *
 *  Algorithm overview
 *  ------------------
 *  For each locus:
 *    1) Count individuals per population (n), total counts (nt), and number of
 *       populations with data (npl).
 *    2) Enumerate alleles observed across all populations.
 *    3) Compute allele frequencies per population p(a,p).
 *    4) Count heterozygotes per allele and population (mho).
 *    5) Accumulate Weir & Cockerham (1984) variance components:
 *         - within-population (lsigw)
 *         - between-population (lsigb)
 *    6) Compute FIS = lsigb / (lsigb + lsigw).
 *
 *  Output
 *  ------
 *  A named R list with:
 *    - loc.names: locus names
 *    - n        : individuals per locus × population
 *    - nt       : total individuals per locus
 *    - alploc   : number of alleles per locus
 *    - p        : list of allele-frequency matrices (allele × population)
 *    - lsigb    : between-population variance component per locus
 *    - lsigw    : within-population variance component per locus
 *    - FIS      : locus-wise FIS
 *
 *  Notes
 *  -----
 *  - A 1-locus input is internally duplicated to satisfy downstream matrix logic.
 *  - The population codes are assumed to be consecutive integers (1..np) after
 *    factor recoding on the R side.
 *
 *  References
 *  ----------
 *  Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the analysis of
 *  population structure. Evolution 38:1358–1370.
 * =============================================================================
 */


// ------------------------------------------------------------
// Decoder with a base-aware decoder
// ------------------------------------------------------------
inline bool decode_gt_base(int gt, int& a1, int& a2, int base) {
  if (gt == NA_INTEGER || gt <= 0) return false;
  a1 = gt / base;
  a2 = gt % base;
  return (a1 > 0 && a2 > 0);
}

// forward declaration — implemented after wc_fis_locus_ptr
static Rcpp::NumericVector wc_fis_by_pop_wc84_impl(
    const Rcpp::IntegerMatrix& dat, int pop_col, int base);

static std::vector< std::vector<int> >
  build_pop_rows(const IntegerMatrix &mat, const int pop_col);

// [[Rcpp::export]]
Rcpp::NumericMatrix boot_indiv_wc_fis_by_pop(
    const Rcpp::IntegerMatrix &mat,
    const int pop_col = 0,
    const int NAcode = 0,
    const int B = 1000,
    const int base = 1000,
    const int debug = 0
) {
  Rcpp::RNGScope rng;
  
  const int N = mat.nrow();
  const int P = mat.ncol();
  
  // number of populations (unique pop codes)
  Rcpp::IntegerVector pop = mat(_, pop_col);
  Rcpp::IntegerVector pops = Rcpp::sort_unique(pop);
  const int np = pops.size();
  
  Rcpp::NumericMatrix out(B, np);
  
  // column names = population codes (as strings)
  Rcpp::CharacterVector pop_names(np);
  for (int i = 0; i < np; ++i) pop_names[i] = std::to_string((int)pops[i]);
  Rcpp::colnames(out) = pop_names;
  
  // strata: list of row indices per population (based on pop_col)
  auto strata_base = build_pop_rows(mat, pop_col);
  
  // IMPORTANT: do NOT OpenMP-parallelise this loop using R::unif_rand()
  for (int b = 0; b < B; ++b) {
    
    Rcpp::IntegerMatrix bm(N, P);
    int wr = 0;  // packed row index
    
    for (size_t s = 0; s < strata_base.size(); ++s) {
      const auto &rows = strata_base[s];
      const int n = (int)rows.size();
      if (n == 0) continue;
      
      for (int i = 0; i < n; ++i) {
        const int src = rows[(int)std::floor(R::unif_rand() * n)];
        for (int j = 0; j < P; ++j) bm(wr, j) = mat(src, j);
        wr++;
      }
    }
    
    if (wr != N) {
      Rcpp::stop("bootstrap internal error: wr (%d) != N (%d)", wr, N);
    }
    
    // POPULATION-wise FIS (WC84 variance-component ratio-of-sums)
    Rcpp::NumericVector fis_pop = wc_fis_by_pop_wc84_impl(bm, pop_col, base);

    // fill output row
    for (int pp = 0; pp < np; ++pp) out(b, pp) = fis_pop[pp];
  }
  
  return out;
}

static void permute_within_pops(IntegerMatrix& m, int pop_col, int base);

// [[Rcpp::export]]
Rcpp::NumericMatrix batch_permute_wc_fis_by_pop(
    Rcpp::IntegerMatrix dat,
    int pop_col_1based,
    int base,
    int B
) {
  Rcpp::RNGScope scope;
  
  const int pop_col = pop_col_1based - 1;
  
  // number of populations (unique pop codes)
  Rcpp::IntegerVector pop = dat(_, pop_col);
  Rcpp::IntegerVector pops = Rcpp::sort_unique(pop);
  const int np = pops.size();
  
  Rcpp::NumericMatrix out(B, np);
  
  // column names = population codes (as strings)
  Rcpp::CharacterVector pop_names(np);
  for (int i = 0; i < np; ++i) pop_names[i] = std::to_string((int)pops[i]);
  Rcpp::colnames(out) = pop_names;
  
  Rcpp::IntegerMatrix pm = Rcpp::clone(dat);
  
  for (int b = 0; b < B; ++b) {
    if (b % 100 == 0 && b > 0) Rcpp::checkUserInterrupt();
    
    std::copy(dat.begin(), dat.end(), pm.begin());
    permute_within_pops(pm, pop_col, base);
    
    // POPULATION-wise FIS (WC84 variance-component ratio-of-sums)
    Rcpp::NumericVector fis_pop = wc_fis_by_pop_wc84_impl(pm, pop_col, base);

    for (int pp = 0; pp < np; ++pp) out(b, pp) = fis_pop[pp];
  }
  
  return out;
}


// [[Rcpp::export]]
Rcpp::List fis_wc_cpp(Rcpp::IntegerMatrix ndat, int base = 1000) {
  
  // ----------------------------
  // pop, ni, dat (factor-style)
  // ----------------------------
  IntegerVector pop = ndat(_, 0); 
  IntegerVector pops = sort_unique(pop);
  int ni = pop.size();
  IntegerMatrix dat = clone(ndat);
  
  // ----------------------------
  // handle 1-locus case
  // ----------------------------
  if (dat.ncol() == 2) {
    IntegerMatrix tmp(dat.nrow(), 3);
    tmp(_,0) = dat(_,0);
    tmp(_,1) = dat(_,1);
    tmp(_,2) = dat(_,1);
    dat = tmp;
  }
  
  int nl = dat.ncol() - 1;
  CharacterVector cn = colnames(dat);
  CharacterVector loc_names(nl);
  
  if (cn.size() == dat.ncol()) {
    for (int l = 0; l < nl; ++l) loc_names[l] = cn[l + 1];
  } else {
    for (int l = 0; l < nl; ++l) loc_names[l] = "L" + std::to_string(l + 1);
  }
  
  // ----------------------------
  // ind.count(dat) → n
  // ----------------------------
  int np = pops.size();
  
  NumericMatrix n(nl, np);
  n.fill(NA_REAL);
  
  for (int l = 0; l < nl; ++l) {
    for (int p = 0; p < np; ++p) {
      int count = 0;
      bool any = false;
      
      for (int i = 0; i < ni; ++i) {
        if (pop[i] != pops[p]) continue;
        int gt = dat(i, l + 1);
        int a1, a2;
        if (!decode_gt_base(gt, a1, a2, base)) continue;
        count++;
        any = true;
      }
      
      if (any) n(l, p) = count;
    }
  }
  
  // ----------------------------
  // nt, npl
  // ----------------------------
  NumericVector nt(nl);
  NumericVector npl(nl);
  
  for (int l = 0; l < nl; ++l) {
    for (int p = 0; p < np; ++p) {
      if (!NumericVector::is_na(n(l,p))) {
        nt[l] += n(l,p);
        npl[l]++;
      }
    }
  }
  
  // ----------------------------
  // nb.alleles + allele list
  // ----------------------------
  IntegerVector alploc(nl);
  std::vector< std::vector<int> > alleles(nl);
  
  for (int l = 0; l < nl; ++l) {
    std::unordered_map<int,bool> seen;
    for (int i = 0; i < ni; ++i) {
      int gt = dat(i, l + 1);
      int a1, a2;
      if (!decode_gt_base(gt, a1, a2, base)) continue;
      seen[a1] = true;
      seen[a2] = true;
    }
    for (auto &kv : seen) alleles[l].push_back(kv.first);
    alploc[l] = alleles[l].size();
  }
  
  // ----------------------------
  // pop.freq(dat) → p
  // ----------------------------
  std::vector< NumericMatrix > p_list(nl);
  
  for (int l = 0; l < nl; ++l) {
    NumericMatrix pmat(alploc[l], np);
    
    for (int a = 0; a < alploc[l]; ++a) {
      int allele = alleles[l][a];
      for (int p = 0; p < np; ++p) {
        double denom = n(l,p);
        if (NumericVector::is_na(denom) || denom == 0) {
          pmat(a,p) = NA_REAL;
          continue;
        }
        int count = 0;
        for (int i = 0; i < ni; ++i) {
          if (pop[i] != pops[p]) continue;
          int gt = dat(i, l + 1);
          int a1, a2;
          if (!decode_gt_base(gt, a1, a2, base)) continue;
          if (a1 == allele) count++;
          if (a2 == allele) count++;
        }
        pmat(a,p) = double(count) / (2.0 * denom);      }
    }
    p_list[l] = pmat;
  }
  
  // ----------------------------
  // mho: heterozygote counts per allele x pop
  // ----------------------------
  std::vector< NumericMatrix > mho(nl);
  
  for (int l = 0; l < nl; ++l) {
    NumericMatrix mho_l(alploc[l], np);
    mho_l.fill(0.0);
    
    for (int i = 0; i < ni; ++i) {
      int pop_idx = pop[i] - 1;              // pop is 1..np after match()
      int gt = dat(i, l + 1);
      int a1, a2;
      if (!decode_gt_base(gt, a1, a2, base)) continue;
      
      if (a1 != a2) {
        for (int a = 0; a < alploc[l]; ++a) {
          int allele = alleles[l][a];
          if (a1 == allele || a2 == allele) {
            mho_l(a, pop_idx) += 1.0;
          }
        }
      }
    }
    
    mho[l] = mho_l;
  }
  
  
  
  // ----------------------------
  // WC84 variance components
  // ----------------------------
  NumericVector lsigb(nl, 0.0);
  NumericVector lsigw(nl, 0.0);
  
  for (int l = 0; l < nl; ++l) {
    NumericMatrix pmat = p_list[l];
    
    for (int a = 0; a < alploc[l]; ++a) {
      
      double SSG = 0.0;
      double SSi = 0.0;
      
      for (int p = 0; p < np; ++p) {
        if (NumericVector::is_na(n(l,p))) continue;
        
        double nal = n(l,p);
        double pap = pmat(a,p);
        double mhom = (2.0 * nal * pap - mho[l](a, p)) / 2.0;
        
        SSG += nal * pap - mhom;
        SSi += nal * (pap - 2 * pap * pap) + mhom;
      }
      
      double MSG = SSG / nt[l];
      double MSI = SSi / (nt[l] - npl[l]);
      
      lsigw[l] += MSG;
      lsigb[l] += 0.5 * (MSI - MSG);
    }
  }
  
  // ----------------------------
  // FIS
  // ----------------------------
  NumericVector FIS(nl);
  for (int l = 0; l < nl; ++l)
    FIS[l] = lsigb[l] / (lsigb[l] + lsigw[l]);
  
  // ----------------------------
  // return
  // ----------------------------
  return List::create(
    _["loc.names"] = loc_names,
    _["n"] = n,
    _["nt"] = nt,
    _["alploc"] = alploc,
    _["p"] = p_list,
    _["lsigb"] = lsigb,
    _["lsigw"] = lsigw,
    _["FIS"]   = FIS
  );
}





// ------------------------------------------------------------
// Base-aware diploid decoder (must match R encoder)
// ------------------------------------------------------------
inline bool decode_gt_base_perm(int gt, int& a1, int& a2, int base) {
  if (gt == NA_INTEGER || gt <= 0) return false;
  a1 = gt / base;
  a2 = gt % base;
  return (a1 > 0 && a2 > 0);
}

inline void fy_shuffle(std::vector<int>& x) {
  for (int i = (int)x.size() - 1; i > 0; --i) {
    int j = (int)(R::unif_rand() * (i + 1));
    std::swap(x[i], x[j]);
  }
}

// ============================================================================
// OpenMP-safe PRNG and raw-pointer kernels
// ============================================================================

struct FastRng {
  uint64_t x;
  FastRng(uint64_t seed = 1) : x(seed ? seed : 1ULL) {}
  inline uint64_t next() {
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    return x * 0x2545F4914F6CDD1DULL;
  }
  inline int unif_int(int n) {
    return (int)((next() >> 11) * ((double)n * (1.0/9007199254740992.0)));
  }
};

inline void fy_shuffle_rng(std::vector<int>& x, FastRng& rng) {
  for (int i = (int)x.size() - 1; i > 0; --i) {
    int j = rng.unif_int(i + 1);
    std::swap(x[i], x[j]);
  }
}

// WC84 FIS for one locus from a raw column-major int* matrix.
// Uses variance-component decomposition (WC84 eq. 4): FIS = b/(b+c).
// Also accumulates lsigb and lsigw into the caller's totals so that
// batch_permute_wc_fis can compute the ratio-of-sums overall FIS.
static double wc_fis_locus_ptr(
    const int* dat, int n, int locus, int pop_col, int base,
    double& lsigb_total, double& lsigw_total
) {
  std::unordered_map<int,int>   n_map;
  std::unordered_map<int, std::unordered_map<int,int>> acounts;
  std::unordered_map<int, std::unordered_map<int,int>> mho_map;
  n_map.reserve(32); acounts.reserve(32); mho_map.reserve(32);

  double nt = 0.0;
  for (int i = 0; i < n; ++i) {
    int gt = dat[i + (size_t)n * locus];
    if (gt == NA_INTEGER || gt <= 0) continue;
    int a1 = gt / base, a2 = gt % base;
    if (a1 <= 0 || a2 <= 0) continue;
    int pop = dat[i + (size_t)n * pop_col];
    n_map[pop]++;
    nt++;
    acounts[pop][a1]++;
    acounts[pop][a2]++;
    if (a1 != a2) { mho_map[pop][a1]++; mho_map[pop][a2]++; }
  }

  const int npl = (int)n_map.size();
  if (npl < 1 || nt < 2.0) return NA_REAL;
  const double denom = nt - (double)npl;
  if (denom <= 0.0) return NA_REAL;

  std::unordered_map<int,bool> all_alleles;
  for (auto& kv : acounts)
    for (auto& av : kv.second) all_alleles[av.first] = true;

  double lsigb = 0.0, lsigw = 0.0;

  for (auto& allele_kv : all_alleles) {
    const int allele = allele_kv.first;
    double SSG = 0.0, SSi = 0.0;

    for (auto& pop_kv : n_map) {
      const int  pop = pop_kv.first;
      const double ni  = pop_kv.second;
      if (ni <= 0.0) continue;

      double ai = 0.0;
      auto ait = acounts.find(pop);
      if (ait != acounts.end()) {
        auto av = ait->second.find(allele);
        if (av != ait->second.end()) ai = av->second;
      }
      const double pap = ai / (2.0 * ni);

      double mho_val = 0.0;
      auto hit = mho_map.find(pop);
      if (hit != mho_map.end()) {
        auto hv = hit->second.find(allele);
        if (hv != hit->second.end()) mho_val = hv->second;
      }
      const double mhom = (2.0 * ni * pap - mho_val) / 2.0;

      SSG += ni * pap - mhom;
      SSi += ni * (pap - 2.0 * pap * pap) + mhom;
    }

    const double MSG = SSG / nt;
    const double MSI = SSi / denom;
    lsigw += MSG;
    lsigb += 0.5 * (MSI - MSG);
  }

  lsigb_total += lsigb;
  lsigw_total += lsigw;

  const double d = lsigb + lsigw;
  return (d == 0.0) ? NA_REAL : lsigb / d;
}

// WC84 FIS per population via variance-component ratio-of-sums.
// For each population, subsets the data, recodes pop code to 1, and
// accumulates lsigb/lsigw across all loci using wc_fis_locus_ptr.
static Rcpp::NumericVector wc_fis_by_pop_wc84_impl(
    const Rcpp::IntegerMatrix& dat,
    int pop_col,
    int base
) {
  const int N = dat.nrow();
  const int P = dat.ncol();

  Rcpp::IntegerVector pop_vec = dat(_, pop_col);
  Rcpp::IntegerVector pops    = Rcpp::sort_unique(pop_vec);
  const int np = pops.size();

  std::unordered_map<int,int> pop_idx_map;
  pop_idx_map.reserve((size_t)np);
  for (int pi = 0; pi < np; ++pi) pop_idx_map[pops[pi]] = pi;

  std::vector<std::vector<int>> pop_rows(np);
  for (int i = 0; i < N; ++i) {
    auto it = pop_idx_map.find(pop_vec[i]);
    if (it != pop_idx_map.end()) pop_rows[it->second].push_back(i);
  }

  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)(P - 1));
  for (int j = 0; j < P; ++j) if (j != pop_col) locus_cols.push_back(j);
  const int L = (int)locus_cols.size();

  const int* dat_ptr = &dat[0];
  Rcpp::NumericVector FIS(np, NA_REAL);

  for (int pi = 0; pi < np; ++pi) {
    const auto& rows = pop_rows[pi];
    const int n_pop  = (int)rows.size();
    if (n_pop < 2) continue;

    // Column-major sub-matrix for this population; recode pop to 1.
    std::vector<int> sub((size_t)n_pop * (size_t)P, 0);
    for (int i = 0; i < n_pop; ++i) {
      const int src = rows[i];
      for (int j = 0; j < P; ++j)
        sub[(size_t)i + (size_t)n_pop * (size_t)j] =
          dat_ptr[src + (size_t)N * (size_t)j];
      sub[(size_t)i + (size_t)n_pop * (size_t)pop_col] = 1;
    }

    // Accumulate WC84 components across all loci (ratio-of-sums FIS).
    double total_lsigb = 0.0, total_lsigw = 0.0;
    for (int ell = 0; ell < L; ++ell)
      wc_fis_locus_ptr(sub.data(), n_pop, locus_cols[ell], pop_col, base,
                       total_lsigb, total_lsigw);

    const double d = total_lsigb + total_lsigw;
    FIS[pi] = (d > 0.0) ? total_lsigb / d : NA_REAL;
  }

  Rcpp::CharacterVector pop_names(np);
  for (int p = 0; p < np; ++p) pop_names[p] = std::to_string((int)pops[p]);
  FIS.attr("names") = pop_names;
  return FIS;
}

// [[Rcpp::export]]
Rcpp::NumericVector wc_fis_by_pop_wc84(
    const Rcpp::IntegerMatrix& dat,
    int pop_col = 0,
    int base    = 1000
) {
  return wc_fis_by_pop_wc84_impl(dat, pop_col, base);
}

// Fisher-Yates within-population permutation on raw column-major int* matrix.
static void permute_within_pops_rng_ptr(
    int* pm, int n, int p, int pop_col, int base,
    const std::vector<std::vector<int>>& strata,
    FastRng& rng
) {
  std::vector<int> alleles, valid_rows;
  alleles.reserve(256); valid_rows.reserve(128);

  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    for (const auto& rows : strata) {
      if ((int)rows.size() < 2) continue;
      alleles.clear(); valid_rows.clear();
      for (int idx : rows) {
        int gt = pm[idx + (size_t)n * j];
        if (gt == NA_INTEGER || gt <= 0) continue;
        int a1 = gt / base, a2 = gt % base;
        if (a1 <= 0 || a2 <= 0) continue;
        alleles.push_back(a1); alleles.push_back(a2);
        valid_rows.push_back(idx);
      }
      if ((int)alleles.size() < 2 || valid_rows.empty()) continue;
      fy_shuffle_rng(alleles, rng);
      int k = 0;
      const int aN = (int)alleles.size();
      for (int idx : valid_rows) {
        if (k + 1 >= aN) break;
        pm[idx + (size_t)n * j] = alleles[k] * base + alleles[k+1];
        k += 2;
      }
    }
  }
}


static Rcpp::NumericVector wc_fis_all_loci(const Rcpp::IntegerMatrix& dat,
                                           int base) {
  Rcpp::List res = fis_wc_cpp(dat, base);
  return Rcpp::as<Rcpp::NumericVector>(res["FIS"]);
}


// -----------------------------------------------------------------------------
// Permutation within population × locus
// -----------------------------------------------------------------------------

static void permute_within_pops(IntegerMatrix& m,
                                int pop_col,
                                int base) {
  const int n = m.nrow();
  const int p = m.ncol();
  
  std::unordered_map<int, std::vector<int>> pop_rows;
  for (int i = 0; i < n; ++i)
    pop_rows[m(i, pop_col)].push_back(i);
  
  std::vector<int> alleles;
  std::vector<int> rows;
  
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    
    for (auto& kv : pop_rows) {
      alleles.clear();
      rows.clear();
      
      for (int i : kv.second) {
        int a1, a2;
        if (!decode_gt_base_perm(m(i, j), a1, a2, base)) continue;
        alleles.push_back(a1);
        alleles.push_back(a2);
        rows.push_back(i);
      }
      
      if (alleles.size() < 2) continue;
      
      fy_shuffle(alleles);
      
      int k = 0;
      for (int i : rows) {
        m(i, j) = alleles[k] * base + alleles[k + 1];
        k += 2;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Exact Weir & Cockerham FIS per locus
// -----------------------------------------------------------------------------

// Compute WC-FIS for one locus
// -----------------------------------------------------------------------------
// Exact port of hierfstat::wc() — per-locus WC FIS (diploid only)
// -----------------------------------------------------------------------------

static double wc_fis_locus(const IntegerMatrix& dat,
                           int locus,
                           int pop_col,
                           int base)
{
  const int N = dat.nrow();
  
  std::unordered_map<int,int> n;      // individuals per pop
  std::unordered_map<int,int> het;    // heterozygotes per pop
  std::unordered_map<int,int> alleles;
  
  for (int i = 0; i < N; ++i) {
    int gt = dat(i,locus);
    int a1, a2;
    if (!decode_gt_base_perm(gt, a1, a2, base)) continue;
    
    if (a1 <= 0 || a2 <= 0) continue;
    
    int pop = dat(i,pop_col);
    n[pop]++;
    
    if (a1 != a2) het[pop]++;
    
    alleles[a1]++;
    alleles[a2]++;
  }
  
  const int np = n.size();
  if (np < 2) return NA_REAL;
  
  double nt = 0.0, sum_n2 = 0.0;
  for (auto& kv : n) {
    nt += kv.second;
    sum_n2 += (double)kv.second * kv.second;
  }
  
  const double nc =
    (nt - sum_n2 / nt) / (np - 1.0);
  
  if (nc <= 0.0) return NA_REAL;
  
  // HT
  double sum_p2 = 0.0;
  for (auto& kv : alleles) {
    double p = kv.second / (2.0 * nt);
    sum_p2 += p * p;
  }
  const double HT = 1.0 - sum_p2;
  
  // HS
  double HS = 0.0;
  for (auto& kv : n) {
    int pop = kv.first;
    double ni = kv.second;
    double Ho = (double)het[pop] / ni;
    HS += Ho;
  }
  HS /= np;
  
  if (HT <= 0.0) return NA_REAL;
  
  return (HT - HS) / HT;
}



// -----------------------------------------------------------------------------
// Exported function
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix batch_permute_wc_fis(IntegerMatrix dat,
                                   int pop_col_1based,
                                   int base,
                                   int B,
                                   int seed = 1) {
  const int n = dat.nrow();
  const int p = dat.ncol();
  const int pop_col = pop_col_1based - 1;
  const int L = p - 1;

  CharacterVector cn = colnames(dat);
  CharacterVector loc_names(L);
  int k = 0;
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    loc_names[k] = (cn.size() == p) ? Rcpp::String(cn[j]) : Rcpp::String("L" + std::to_string(k + 1));
    k++;
  }

  auto strata = build_pop_rows(dat, pop_col);
  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)L);
  for (int j = 0; j < p; ++j) if (j != pop_col) locus_cols.push_back(j);

  const int* dat_ptr = &dat[0];
  std::vector<std::vector<double>> out_raw(B, std::vector<double>(L + 1, NA_REAL));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = (uint64_t)seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    FastRng rng(sb);

    std::vector<int> pm((size_t)n * (size_t)p);
    std::memcpy(pm.data(), dat_ptr, sizeof(int) * (size_t)n * (size_t)p);

    permute_within_pops_rng_ptr(pm.data(), n, p, pop_col, base, strata, rng);

    double total_lsigb = 0.0, total_lsigw = 0.0;
    for (int ell = 0; ell < L; ++ell)
      out_raw[b][ell] = wc_fis_locus_ptr(pm.data(), n, locus_cols[ell], pop_col, base,
                                          total_lsigb, total_lsigw);
    out_raw[b][L] = (total_lsigb + total_lsigw > 0.0)
      ? total_lsigb / (total_lsigb + total_lsigw) : NA_REAL;
  }

  CharacterVector loc_names_ext(L + 1);
  for (int i = 0; i < L; ++i) loc_names_ext[i] = loc_names[i];
  loc_names_ext[L] = "Overall";

  NumericMatrix out(B, L + 1);
  colnames(out) = loc_names_ext;
  for (int b = 0; b < B; ++b)
    for (int ell = 0; ell <= L; ++ell)
      out(b, ell) = out_raw[b][ell];

  return out;
}








// ============================================================================
//   fis_bootstrap.cpp  — FAST BOOTSTRAP FOR WEIR & COCKERHAM FIS (Option C)
//   Only includes essential core functions:
//      1) calculate_observed_fis()
//      2) boot_indiv_wc_fis()
//      3) summarize_fis_results()
//      4) create_results_dataframe()
//   Parallelised with OpenMP when available
// ============================================================================

// ============================================================================
// GENOTYPE DECODER (very fast)
// ============================================================================
inline bool decode_geno_base(const int g, int &a1, int &a2, const int NAcode, const int base) {
  if (g == NAcode || g <= 0) return false;
  a1 = g / base;
  a2 = g % base;
  return (a1 > 0 && a2 > 0);
}


// ============================================================================
// BUILD VECTORS OF INDIVIDUAL ROWS PER POPULATION
// ============================================================================
static std::vector< std::vector<int> >
  build_pop_rows(const IntegerMatrix &mat, const int pop_col) {
    
    const int N = mat.nrow();
    std::unordered_map<int, std::vector<int>> tmp;
    
    for (int i = 0; i < N; i++) {
      tmp[ mat(i, pop_col) ].push_back(i);
    }
    
    std::vector< std::vector<int> > out;
    out.reserve(tmp.size());
    for (auto &kv : tmp) out.push_back(kv.second);
    
    return out;
  }

// ============================================================================
// 1) OBSERVED FIS (no bootstrap)
// ============================================================================
// [[Rcpp::export]]
Rcpp::NumericVector calculate_observed_fis(
    const Rcpp::IntegerMatrix &mat,
    const int pop_col = 0,
    const int NAcode = 0,
    const int base = 1000
) {
  // pop_col and NAcode are unused here because fis_wc_cpp() already handles missing
  // and assumes pop is in column 0 (hierfstat style) in your current implementation.
  // If your fis_wc_cpp assumes pop is always col 0, keep it that way consistently.
  return wc_fis_all_loci(mat, base);
}


// ============================================================================
// 2) BOOTSTRAP INDIVIDUALS WITHIN POPULATIONS (CORRECT FOR FIS)
// ============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix boot_indiv_wc_fis(
    const Rcpp::IntegerMatrix &mat,
    const int pop_col = 0,
    const int NAcode = 0,
    const int B = 1000,
    const int base = 1000,
    const int debug = 0,
    const int seed = 1
) {
  const int N = mat.nrow();
  const int P = mat.ncol();
  const int L = P - 1;

  Rcpp::CharacterVector cn = Rcpp::colnames(mat);
  Rcpp::CharacterVector loc_names_ext(L + 1);
  if (cn.size() == P) {
    for (int j = 1; j < P; ++j) loc_names_ext[j - 1] = cn[j];
  } else {
    for (int j = 0; j < L; ++j) loc_names_ext[j] = "L" + std::to_string(j + 1);
  }
  loc_names_ext[L] = "Overall";

  Rcpp::NumericMatrix out(B, L + 1);
  Rcpp::colnames(out) = loc_names_ext;

  auto strata = build_pop_rows(mat, pop_col);
  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)L);
  for (int j = 0; j < P; ++j) if (j != pop_col) locus_cols.push_back(j);

  const int* dat_ptr = &mat[0];
  std::vector<std::vector<double>> out_raw(B, std::vector<double>(L + 1, NA_REAL));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = (uint64_t)seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    FastRng rng(sb);

    std::vector<int> bm((size_t)N * (size_t)P);
    int wr = 0;

    for (const auto& rows : strata) {
      const int ns = (int)rows.size();
      if (ns == 0) continue;
      for (int i = 0; i < ns; ++i) {
        const int src = rows[(size_t)rng.unif_int(ns)];
        for (int j = 0; j < P; ++j)
          bm[(size_t)wr + (size_t)N * (size_t)j] = dat_ptr[src + (size_t)N * (size_t)j];
        wr++;
      }
    }

    if (wr != N) continue;

    double total_lsigb = 0.0, total_lsigw = 0.0;
    for (int ell = 0; ell < L; ++ell)
      out_raw[b][ell] = wc_fis_locus_ptr(bm.data(), N, locus_cols[ell], pop_col, base,
                                          total_lsigb, total_lsigw);
    out_raw[b][L] = (total_lsigb + total_lsigw > 0.0)
      ? total_lsigb / (total_lsigb + total_lsigw) : NA_REAL;
  }

  for (int b = 0; b < B; ++b)
    for (int ell = 0; ell <= L; ++ell)
      out(b, ell) = out_raw[b][ell];

  return out;
}

// ============================================================================
// 2b) BOOTSTRAP POPULATIONS (POP-BLOCK) FOR FIS
// Resamples populations with replacement — same scheme as FST pop-block boot.
// Returns B x L matrix of FIS estimates, one row per bootstrap replicate.
// ============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix boot_popblock_wc_fis(
    const Rcpp::IntegerMatrix &mat,
    const int pop_col = 0,
    const int NAcode = 0,
    const int B = 1000,
    const int base = 1000,
    const int seed = 1
) {
  const int N = mat.nrow();
  const int P = mat.ncol();
  const int L = P - 1;

  Rcpp::CharacterVector cn = Rcpp::colnames(mat);
  Rcpp::CharacterVector loc_names_ext(L + 1);
  if (cn.size() == P) {
    for (int j = 1; j < P; ++j) loc_names_ext[j - 1] = cn[j];
  } else {
    for (int j = 0; j < L; ++j) loc_names_ext[j] = "L" + std::to_string(j + 1);
  }
  loc_names_ext[L] = "Overall";

  Rcpp::NumericMatrix out(B, L + 1);
  Rcpp::colnames(out) = loc_names_ext;

  auto strata = build_pop_rows(mat, pop_col);
  const int n_pops = (int)strata.size();
  if (n_pops < 2) Rcpp::stop("boot_popblock_wc_fis: need at least 2 populations.");

  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)L);
  for (int j = 0; j < P; ++j) if (j != pop_col) locus_cols.push_back(j);

  const int* dat_ptr = &mat[0];
  std::vector<std::vector<double>> out_raw(B, std::vector<double>(L + 1, NA_REAL));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = (uint64_t)seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    FastRng rng(sb);

    std::vector<int> times((size_t)n_pops, 0);
    for (int k = 0; k < n_pops; ++k) ++times[(size_t)rng.unif_int(n_pops)];

    int wr_needed = 0;
    for (int k = 0; k < n_pops; ++k)
      wr_needed += times[(size_t)k] * (int)strata[(size_t)k].size();
    if (wr_needed < 2) continue;

    std::vector<int> bm((size_t)wr_needed * (size_t)P);
    int wr = 0;
    int new_pop_code = 1;

    for (int k = 0; k < n_pops; ++k) {
      const int mult = times[(size_t)k];
      if (mult == 0) continue;
      const auto& rows = strata[(size_t)k];
      for (int rep = 0; rep < mult; ++rep) {
        for (int rid : rows) {
          for (int j = 0; j < P; ++j)
            bm[(size_t)wr + (size_t)wr_needed * (size_t)j] = dat_ptr[rid + (size_t)N * (size_t)j];
          bm[(size_t)wr + (size_t)wr_needed * (size_t)pop_col] = new_pop_code;
          wr++;
        }
        ++new_pop_code;
      }
    }

    double total_lsigb = 0.0, total_lsigw = 0.0;
    for (int ell = 0; ell < L; ++ell)
      out_raw[b][ell] = wc_fis_locus_ptr(bm.data(), wr_needed, locus_cols[ell], pop_col, base,
                                          total_lsigb, total_lsigw);
    out_raw[b][L] = (total_lsigb + total_lsigw > 0.0)
      ? total_lsigb / (total_lsigb + total_lsigw) : NA_REAL;
  }

  for (int b = 0; b < B; ++b)
    for (int ell = 0; ell <= L; ++ell)
      out(b, ell) = out_raw[b][ell];

  return out;
}


// ============================================================================
// 3) SUMMARISE BOOTSTRAP -> MEAN + CI
// ============================================================================
// [[Rcpp::export]]
Rcpp::List summarize_fis_results(
    const Rcpp::NumericMatrix &boot,
    const double conf = 0.95
) {
  const int B = boot.nrow();
  const int Ltot = boot.ncol();

  // Detect ratio-of-sums "Overall" column (last col, added by boot_indiv/popblock/permute).
  Rcpp::CharacterVector boot_cn = Rcpp::colnames(boot);
  const bool has_overall =
    (boot_cn.size() == Ltot) &&
    (std::string(Rcpp::as<std::string>(boot_cn[Ltot - 1])) == "Overall");
  const int L = has_overall ? Ltot - 1 : Ltot;

  Rcpp::NumericVector mean(L), LCI(L), UCI(L);

  const double alpha = (1.0 - conf) / 2.0;

  // ---- per-locus summaries (columns 0..L-1 only)
  for (int ell = 0; ell < L; ++ell) {
    
    std::vector<double> v;
    v.reserve(B);
    
    for (int b = 0; b < B; ++b) {
      double x = boot(b, ell);
      if (!Rcpp::NumericVector::is_na(x)) v.push_back(x);
    }
    
    if (v.empty()) {
      mean[ell] = LCI[ell] = UCI[ell] = NA_REAL;
      continue;
    }
    
    double s = 0.0;
    for (double x : v) s += x;
    mean[ell] = s / v.size();
    
    if (v.size() > 1) {
      std::sort(v.begin(), v.end());
      int lo = std::max(0, (int)std::floor(v.size() * alpha));
      int hi = std::min((int)v.size() - 1, (int)std::floor(v.size() * (1.0 - alpha)));
      LCI[ell] = v[lo];
      UCI[ell] = v[hi];
    } else {
      LCI[ell] = UCI[ell] = NA_REAL;
    }
  }
  
  // ---- overall summaries
  // Use ratio-of-sums "Overall" column when present; fall back to rowMeans.
  std::vector<double> overall;
  overall.reserve(B);

  if (has_overall) {
    for (int b = 0; b < B; ++b) {
      double x = boot(b, L);   // column index L = last column "Overall"
      if (!Rcpp::NumericVector::is_na(x)) overall.push_back(x);
    }
  } else {
    for (int b = 0; b < B; ++b) {
      double s = 0.0;
      int k = 0;
      for (int ell = 0; ell < L; ++ell) {
        double x = boot(b, ell);
        if (!Rcpp::NumericVector::is_na(x)) { s += x; k++; }
      }
      if (k > 0) overall.push_back(s / (double)k);
    }
  }
  
  double overall_mean = NA_REAL, overall_lci = NA_REAL, overall_uci = NA_REAL;
  
  if (!overall.empty()) {
    double s = 0.0;
    for (double x : overall) s += x;
    overall_mean = s / overall.size();
    
    if (overall.size() > 1) {
      std::sort(overall.begin(), overall.end());
      int lo = std::max(0, (int)std::floor(overall.size() * alpha));
      int hi = std::min((int)overall.size() - 1, (int)std::floor(overall.size() * (1.0 - alpha)));
      overall_lci = overall[lo];
      overall_uci = overall[hi];
    }
  }
  
  return Rcpp::List::create(
    Rcpp::_["mean"] = mean,
    Rcpp::_["ci_lower"] = LCI,
    Rcpp::_["ci_upper"] = UCI,
    Rcpp::_["overall_mean"] = overall_mean,
    Rcpp::_["overall_ci_lower"] = overall_lci,
    Rcpp::_["overall_ci_upper"] = overall_uci,
    Rcpp::_["n_bootstrap"] = B,
    Rcpp::_["confidence"] = conf
  );
}



// ============================================================================
// 4) CREATE A CLEAN DATAFRAME FOR R
// ============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame create_results_dataframe(
    const Rcpp::NumericVector &obs,
    const Rcpp::List &sum,
    const Rcpp::CharacterVector &locus_names
) {
  Rcpp::NumericVector LCI  = sum["ci_lower"];
  Rcpp::NumericVector UCI  = sum["ci_upper"];
  
  // Observed overall = mean of observed locus FIS (NA-safe)
  double obs_overall = NA_REAL;
  {
    double s = 0.0;
    int k = 0;
    for (int i = 0; i < obs.size(); ++i) {
      double x = obs[i];
      if (!Rcpp::NumericVector::is_na(x)) {
        s += x;
        k++;
      }
    }
    if (k > 0) obs_overall = s / (double)k;
  }
  
  // Bootstrap overall (already computed in summarize_fis_results)
  double overall_lci = sum.containsElementNamed("overall_ci_lower")
    ? Rcpp::as<double>(sum["overall_ci_lower"]) : NA_REAL;
  double overall_uci = sum.containsElementNamed("overall_ci_upper")
    ? Rcpp::as<double>(sum["overall_ci_upper"]) : NA_REAL;
  
  // Build extended vectors (+1 for Overall)
  const int L = obs.size();
  
  Rcpp::CharacterVector Locus(L + 1);
  Rcpp::NumericVector Observed_FIS(L + 1);
  Rcpp::NumericVector CI_L(L + 1);
  Rcpp::NumericVector CI_U(L + 1);
  
  for (int i = 0; i < L; ++i) {
    Locus[i] = locus_names[i];
    Observed_FIS[i] = obs[i];
    CI_L[i] = LCI[i];
    CI_U[i] = UCI[i];
  }
  
  Locus[L] = "Overall";
  Observed_FIS[L] = obs_overall;
  CI_L[L] = overall_lci;
  CI_U[L] = overall_uci;
  
  return Rcpp::DataFrame::create(
    Rcpp::_["Locus"]        = Locus,
    Rcpp::_["Observed_FIS"] = Observed_FIS,
    Rcpp::_["CI_L"]         = CI_L,
    Rcpp::_["CI_U"]         = CI_U,
    Rcpp::_["stringsAsFactors"] = false
  );
}


// ============================================================================
// Per-allele WC84 F-statistics (FIS, FST, FIT)
// ============================================================================
// For each allele k at each locus, compute WC84 variance components a_k, b_k,
// c_k and the derived F-statistics:
//   FIS_k = b_k / (b_k + c_k)
//   FST_k = a_k / (a_k + b_k + c_k)
//   FIT_k = (a_k + b_k) / (a_k + b_k + c_k)
//
// Useful to detect amplification dropout (high FIS for a specific allele) or
// directional selection (outlier FST for a specific allele).
// ============================================================================
// [[Rcpp::export]]
Rcpp::DataFrame wc84_per_allele_fstats_cpp(
    const Rcpp::IntegerMatrix& dat,
    int pop_col   = 0,
    int base      = 1000,
    int missing_code = 0
) {
  const int N     = dat.nrow();
  const int ncols = dat.ncol();

  Rcpp::CharacterVector cn = Rcpp::colnames(dat);

  std::vector<int>         loci_cols;
  std::vector<std::string> loci_names_vec;
  loci_cols.reserve(ncols - 1);
  loci_names_vec.reserve(ncols - 1);
  for (int j = 0; j < ncols; ++j) {
    if (j == pop_col) continue;
    loci_cols.push_back(j);
    if (cn.size() == ncols)
      loci_names_vec.push_back(std::string(cn[j]));
    else
      loci_names_vec.push_back("L" + std::to_string((int)loci_cols.size()));
  }
  const int L = (int)loci_cols.size();

  // Output accumulators
  std::vector<std::string> out_locus;
  std::vector<int>         out_allele;
  std::vector<int>         out_n_geno;
  std::vector<int>         out_n_pops;
  std::vector<double>      out_freq;
  std::vector<double>      out_fis;
  std::vector<double>      out_fst;
  std::vector<double>      out_fit;

  for (int li = 0; li < L; ++li) {
    const int col = loci_cols[li];
    const std::string& loc_name = loci_names_vec[li];

    // --- collect per-population data for this locus
    std::unordered_map<int,int> pop_to_idx;
    std::vector<int> n_i;
    std::vector< std::unordered_map<int,int> > acount;
    std::vector< std::unordered_map<int,int> > mho;
    std::unordered_map<int,char> all_alleles;
    std::unordered_map<int,int>  pooled_acount;
    int total_n = 0;

    for (int i = 0; i < N; ++i) {
      int gt = dat(i, col);
      int a1, a2;
      if (!decode_gt_base(gt, a1, a2, base)) continue;

      const int pop = dat(i, pop_col);
      auto it = pop_to_idx.find(pop);
      int pi;
      if (it == pop_to_idx.end()) {
        pi = (int)n_i.size();
        pop_to_idx[pop] = pi;
        n_i.push_back(0);
        acount.emplace_back();
        mho.emplace_back();
      } else {
        pi = it->second;
      }

      n_i[pi]    += 1;
      total_n    += 1;
      acount[pi][a1] += 1;  acount[pi][a2] += 1;
      pooled_acount[a1] += 1; pooled_acount[a2] += 1;
      all_alleles[a1] = 1;  all_alleles[a2] = 1;
      if (a1 != a2) { mho[pi][a1] += 1; mho[pi][a2] += 1; }
    }

    const int r = (int)n_i.size();

    // Not enough data: emit NA rows
    if (r < 2 || total_n <= 1) {
      for (const auto& kv : all_alleles) {
        out_locus.push_back(loc_name);
        out_allele.push_back(kv.first);
        out_n_geno.push_back(total_n);
        out_n_pops.push_back(r);
        out_freq.push_back(total_n > 0 ?
          (double)pooled_acount[kv.first] / (2.0 * total_n) : NA_REAL);
        out_fis.push_back(NA_REAL);
        out_fst.push_back(NA_REAL);
        out_fit.push_back(NA_REAL);
      }
      continue;
    }

    double n_total = 0.0, sum_n2 = 0.0;
    for (int i = 0; i < r; ++i) { n_total += n_i[i]; sum_n2 += (double)n_i[i]*n_i[i]; }
    const double n_bar = n_total / (double)r;
    const double n_c   = (n_total - sum_n2 / n_total) / (double)(r - 1);

    if (n_c <= 0.0 || n_bar <= 1.0) {
      for (const auto& kv : all_alleles) {
        out_locus.push_back(loc_name);
        out_allele.push_back(kv.first);
        out_n_geno.push_back(total_n);
        out_n_pops.push_back(r);
        out_freq.push_back((double)pooled_acount[kv.first] / (2.0 * n_total));
        out_fis.push_back(NA_REAL);
        out_fst.push_back(NA_REAL);
        out_fit.push_back(NA_REAL);
      }
      continue;
    }

    // --- per-allele WC84 components
    for (const auto& kv_allele : all_alleles) {
      const int allele = kv_allele.first;

      // weighted mean allele frequency across pops
      double pbar = 0.0;
      for (int i = 0; i < r; ++i) {
        if (n_i[i] <= 0) continue;
        auto it = acount[i].find(allele);
        double pi = (it != acount[i].end()) ?
          (double)it->second / (2.0 * n_i[i]) : 0.0;
        pbar += (double)n_i[i] * pi;
      }
      pbar /= n_total;

      // between-population variance of allele frequency
      double s2_num = 0.0;
      for (int i = 0; i < r; ++i) {
        if (n_i[i] <= 0) continue;
        auto it = acount[i].find(allele);
        double pi = (it != acount[i].end()) ?
          (double)it->second / (2.0 * n_i[i]) : 0.0;
        s2_num += (double)n_i[i] * (pi - pbar) * (pi - pbar);
      }
      const double s2 = s2_num / ((double)(r - 1) * n_bar);

      // mean observed heterozygosity for this allele
      double hbar = 0.0;
      for (int i = 0; i < r; ++i) {
        auto it = mho[i].find(allele);
        if (it != mho[i].end()) hbar += (double)it->second;
      }
      hbar /= n_total;

      const double term = pbar * (1.0 - pbar);

      // WC84 variance components (a = between-pop, b = between-indiv, c = within-indiv)
      const double a = (n_bar / n_c) *
        (s2 - (1.0 / (n_bar - 1.0)) *
        (term - ((double)(r - 1) / (double)r) * s2 - 0.25 * hbar));

      const double b = (n_bar / (n_bar - 1.0)) *
        (term - ((double)(r - 1) / (double)r) * s2 -
        ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * hbar);

      const double c = 0.5 * hbar;

      const double denom_abc = a + b + c;
      const double fst = (denom_abc != 0.0) ? (a / denom_abc) : NA_REAL;
      const double fit = (denom_abc != 0.0) ? ((a + b) / denom_abc) : NA_REAL;
      const double bc  = b + c;
      const double fis = (bc != 0.0) ? (b / bc) : NA_REAL;

      out_locus.push_back(loc_name);
      out_allele.push_back(allele);
      out_n_geno.push_back(total_n);
      out_n_pops.push_back(r);
      out_freq.push_back(pbar);
      out_fis.push_back(fis);
      out_fst.push_back(fst);
      out_fit.push_back(fit);
    }
  }

  const int nrows = (int)out_locus.size();
  Rcpp::CharacterVector r_locus(nrows);
  Rcpp::IntegerVector   r_allele(nrows);
  Rcpp::IntegerVector   r_n_geno(nrows);
  Rcpp::IntegerVector   r_n_pops(nrows);
  Rcpp::NumericVector   r_freq(nrows);
  Rcpp::NumericVector   r_fis(nrows);
  Rcpp::NumericVector   r_fst(nrows);
  Rcpp::NumericVector   r_fit(nrows);

  for (int i = 0; i < nrows; ++i) {
    r_locus[i]  = out_locus[i];
    r_allele[i] = out_allele[i];
    r_n_geno[i] = out_n_geno[i];
    r_n_pops[i] = out_n_pops[i];
    r_freq[i]   = out_freq[i];
    r_fis[i]    = out_fis[i];
    r_fst[i]    = out_fst[i];
    r_fit[i]    = out_fit[i];
  }

  return Rcpp::DataFrame::create(
    Rcpp::_["Locus"]       = r_locus,
    Rcpp::_["Allele"]      = r_allele,
    Rcpp::_["N_genotyped"] = r_n_geno,
    Rcpp::_["N_pops"]      = r_n_pops,
    Rcpp::_["Freq"]        = r_freq,
    Rcpp::_["FIS"]         = r_fis,
    Rcpp::_["FST"]         = r_fst,
    Rcpp::_["FIT"]         = r_fit,
    Rcpp::_["stringsAsFactors"] = false
  );
}


