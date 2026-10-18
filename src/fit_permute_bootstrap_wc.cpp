// ============================================================================
// src/fit_permute_bootstrap_wc.cpp
// Optimized version with better memory management and performance
// ============================================================================
//
// FUNCTIONS IN THIS FILE (in order of appearance) — multilocus FIT
// (Weir & Cockerham 1984), used by the Global Panmixia module.
//
//   boot_wc84_fit_popblock_raw_cpp()   (internal, not exported) low-level
//                                       population-block bootstrap kernel for
//                                       FIT (raw variance components per
//                                       replicate).
//   boot_wc84_stats_popblock_cpp()     Bootstrap over SUB-SAMPLES
//                                       (populations) for the FIT CI — needs
//                                       >=5 populations to be meaningful
//                                       (see the N<5 guard at the R call
//                                       site in server_general_stats.R).
//   batch_permute_wc84_stats()         Permutation test for FIT.
//   batch_permute_fit_global()         Batched permutation driver, all loci
//                                       combined into the multilocus FIT.
//
// RNG / REPRODUCIBILITY NOTE (added following a code audit; see the matching
// note in fis_permute_bootstrap_wc.cpp for the full explanation): this file
// draws permutation randomness via R::unif_rand() (reproducible via
// set.seed(), NOT thread-safe, hence not parallelised for those loops -
// see "Permutation RNG: R::unif_rand() (reproducible via set.seed())"
// below). This differs from fst_permute_boostrap_OpenMP.cpp, which uses a
// thread-local std::mt19937_64 and IS parallelised. Do not harmonize the two
// schemes without re-validating numerically against reference software
// (see tests/testthat/test-numeric-validation-fst.R).
// ============================================================================
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <memory>
#include <cstdint>
#include <cstring>
#include <string>

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



// ---------------------- OPTIMIZED RNG & UTILITIES --------------------------
inline bool decode_gt_base(int gt, int& a1, int& a2, int base) {
  if (gt == NA_INTEGER || gt <= 0) return false;
  a1 = gt / base;
  a2 = gt % base;
  return (a1 > 0 && a2 > 0);
}

// Helper: decoder with missing_code guard (recommended)
inline bool decode_gt_base_with_missing(int gt, int missing_code, int& a1, int& a2, int base) {
  if (gt == missing_code) return false;
  return decode_gt_base(gt, a1, a2, base);
}

// Optimized PRNG - faster than original
struct FastRng {
  uint64_t x;
  FastRng(uint64_t seed = 1) : x(seed) {}
  
  inline uint64_t next() {
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    return x * 0x2545F4914F6CDD1DULL;
  }
  
  inline double unif01() { 
    return (next() >> 11) * (1.0/9007199254740992.0);
  }
  
  inline int unif_int(int n) {
    return static_cast<int>(unif01() * n);
  }
};


// Reusable allele counting structure to avoid allocations
struct AlleleCounter {
  std::vector<int> n_dip;
  std::vector<int> n_het;
  std::vector<std::unordered_map<int, int>> acnt;
  std::vector<std::unordered_map<int,int>> hetcnt;
  std::unordered_map<int, int> acnt_total;
  int n_dip_total = 0;
  int n_het_total = 0;
  
  void reset(size_t num_groups) {
    n_dip.assign(num_groups, 0);
    n_het.assign(num_groups, 0);
    acnt.assign(num_groups, std::unordered_map<int, int>());
    hetcnt.assign(num_groups, std::unordered_map<int,int>());
    acnt_total.clear();
    n_dip_total = 0;
    n_het_total = 0;
    
  }
};





// -----------------------------------------------------------------------------
// Per-locus stats from a matrix: FST (WC84), FIT (WC84), HI (Ho), HS (weighted Hs)
// pop_col is 0-based
// -----------------------------------------------------------------------------
static Rcpp::List wc84_fst_fit_hi_hs_per_locus_from_matrix(
    const Rcpp::IntegerMatrix& mat_int,
    const int pop_col,
    const int missing_code,
    const int base
) {
  const int N = mat_int.nrow();
  const int P = mat_int.ncol();
  if (P < 2) Rcpp::stop("Need pop + >=1 locus columns.");
  if (pop_col < 0 || pop_col >= P) Rcpp::stop("pop_col out of range.");
  
  // Build pop -> idx and rows_by_pop
  std::unordered_map<int,int> pop2idx;
  pop2idx.reserve((size_t)N);
  for (int i = 0; i < N; ++i) {
    const int p = mat_int(i, pop_col);
    if (pop2idx.find(p) == pop2idx.end()) pop2idx.emplace(p, (int)pop2idx.size());
  }
  const int R = (int)pop2idx.size();
  
  std::vector<std::vector<int>> rows_by_pop(R);
  for (int i = 0; i < N; ++i) rows_by_pop[(size_t)pop2idx[ mat_int(i, pop_col) ]].push_back(i);
  
  // locus columns
  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)P - 1);
  for (int j = 0; j < P; ++j) if (j != pop_col) locus_cols.push_back(j);
  const int L = (int)locus_cols.size();
  
  Rcpp::NumericVector FST(L, NA_REAL), FIT(L, NA_REAL), HI(L, NA_REAL), HS(L, NA_REAL);
  if (R <= 1) return Rcpp::List::create(_["FST"]=FST, _["FIT"]=FIT, _["HI"]=HI, _["HS"]=HS);
  
  for (int ell = 0; ell < L; ++ell) {
    const int j = locus_cols[(size_t)ell];
    
    std::vector<int> n_dip((size_t)R, 0), n_het((size_t)R, 0);
    std::vector<std::unordered_map<int,int>> acnt((size_t)R);
    std::vector<std::unordered_map<int,int>> hetcnt((size_t)R);
    std::unordered_map<int,int> acnt_total;
    acnt_total.reserve(64);
    
    int n_dip_total = 0;
    int n_het_total = 0;
    
    for (int gidx = 0; gidx < R; ++gidx) {
      const auto& rows = rows_by_pop[(size_t)gidx];
      int nd = 0, nh = 0;
      
      for (int rid : rows) {
        const int G = mat_int(rid, j);
        int a1 = 0, a2 = 0;
        if (!decode_gt_base_with_missing(G, missing_code, a1, a2, base)) continue;
        
        ++nd;
        if (a1 != a2) {
          ++nh;
          ++hetcnt[(size_t)gidx][a1];
          ++hetcnt[(size_t)gidx][a2];
        }
        ++acnt[(size_t)gidx][a1];
        ++acnt[(size_t)gidx][a2];
        ++acnt_total[a1];
        ++acnt_total[a2];
      }
      
      n_dip[(size_t)gidx] = nd;
      n_het[(size_t)gidx] = nh;
      n_dip_total += nd;
      n_het_total += nh;
    }
    
    // active groups
    std::vector<int> active;
    active.reserve((size_t)R);
    for (int gidx = 0; gidx < R; ++gidx) if (n_dip[(size_t)gidx] > 0) active.push_back(gidx);
    const int r = (int)active.size();
    
    if (n_dip_total == 0 || r <= 1) { FST[ell]=FIT[ell]=HI[ell]=HS[ell]=NA_REAL; continue; }
    
    // HI (observed heterozygosity)
    const double HIv = (double)n_het_total / (double)n_dip_total;
    
    // HS (within-pop gene diversity, weighted)
    double HSv = 0.0;
    const double denom_glob = 2.0 * (double)n_dip_total;
    for (int gidx : active) {
      const double nd = (double)n_dip[(size_t)gidx];
      const double denom_g = 2.0 * nd;
      
      double sum_p2 = 0.0;
      for (const auto& kv : acnt[(size_t)gidx]) {
        const double p = (double)kv.second / denom_g;
        sum_p2 += p * p;
      }
      const double Hk = 1.0 - sum_p2;
      HSv += (denom_g / denom_glob) * Hk;
    }
    
    // W&C size terms
    double n_total = 0.0, sum_n2 = 0.0;
    for (int gidx : active) {
      n_total += (double)n_dip[(size_t)gidx];
      sum_n2  += (double)n_dip[(size_t)gidx] * (double)n_dip[(size_t)gidx];
    }
    const double n_bar = n_total / (double)r;
    const double n_C   = (n_total - (sum_n2 / n_total)) / ((double)r - 1.0);
    
    if (n_bar <= 1.0 || n_C <= 0.0) {
      FST[ell]=NA_REAL; FIT[ell]=NA_REAL; HI[ell]=HIv; HS[ell]=HSv;
      continue;
    }
    
    // WC84 A,Bc,C
    double A = 0.0, Bc = 0.0, C = 0.0;
    const double denom_total_gene = 2.0 * n_total;
    
    for (const auto& ap : acnt_total) {
      const int al = ap.first;
      const double p_bar = (double)ap.second / denom_total_gene;
      
      // s2
      double s_num = 0.0;
      for (int gidx : active) {
        const double nd = (double)n_dip[(size_t)gidx];
        const double denom_g = 2.0 * nd;
        double p_i = 0.0;
        auto it = acnt[(size_t)gidx].find(al);
        if (it != acnt[(size_t)gidx].end()) p_i = (double)it->second / denom_g;
        const double diff = p_i - p_bar;
        s_num += nd * diff * diff;
      }
      const double s2 = s_num / (n_bar * ((double)r - 1.0));
      
      // hbar_j
      double hbar_num = 0.0, hbar_den = 0.0;
      for (int gidx : active) {
        const double nd = (double)n_dip[(size_t)gidx];
        if (nd <= 0.0) continue;
        hbar_den += nd;
        auto hit = hetcnt[(size_t)gidx].find(al);
        if (hit != hetcnt[(size_t)gidx].end()) hbar_num += (double)hit->second;
      }
      const double hbar_j = (hbar_den > 0.0) ? (hbar_num / hbar_den) : 0.0;
      
      const double term_common =
        (p_bar * (1.0 - p_bar)) - ((((double)r - 1.0) / (double)r) * s2);
      
      const double a_j = (n_bar / n_C) *
        (s2 - ((term_common - (hbar_j / 4.0)) / (n_bar - 1.0)));
      
      const double b_j = (n_bar / (n_bar - 1.0)) *
        (term_common - (((2.0 * n_bar) - 1.0) / (4.0 * n_bar)) * hbar_j);
      
      const double c_j = hbar_j / 2.0;
      
      A  += a_j;
      Bc += b_j;
      C  += c_j;
    }
    
    const double denom_tot = A + Bc + C;
    if (!std::isfinite(denom_tot) || denom_tot == 0.0) {
      FST[ell]=NA_REAL; FIT[ell]=NA_REAL;
    } else {
      FST[ell] = A / denom_tot;
      FIT[ell] = (A + Bc) / denom_tot;
    }
    HI[ell] = HIv;
    HS[ell] = HSv;
  }
  
  return Rcpp::List::create(_["FST"]=FST, _["FIT"]=FIT, _["HI"]=HI, _["HS"]=HS);
}



// --------------------------------------------------------------------------
// (internal: not exported to R; called by boot_wc84_fit_popblock_cpp below)
Rcpp::List boot_wc84_fit_popblock_raw_cpp(const Rcpp::IntegerMatrix mat_int,
                                          const int pop_col,
                                          const int missing_code,
                                          const int base,
                                          const int B,
                                          const int seed = 1,
                                          const int n_threads = 0) {
  
  const int N = mat_int.nrow();
  const int P = mat_int.ncol();
  if (P < 2) stop("Need at least pop + 1 locus columns.");
  if (pop_col < 0 || pop_col >= P) stop("pop_col out of range.");
  if (B <= 0) stop("B must be > 0.");
  
  // Column-major access
  const int* data_ptr = &mat_int[0];
  auto col_ptr = [&](int j) -> const int* {
    return data_ptr + static_cast<size_t>(j) * N;
  };
  
  // Pop structure
  const int* popc = col_ptr(pop_col);
  std::unordered_map<int, int> pop2idx;
  pop2idx.reserve((size_t)N);
  for (int i = 0; i < N; ++i) {
    if (pop2idx.find(popc[i]) == pop2idx.end()) pop2idx.emplace(popc[i], (int)pop2idx.size());
  }
  const int R = (int)pop2idx.size();
  if (R < 1) stop("No populations found.");
  
  std::vector<std::vector<int>> rows_by_pop(R);
  for (int i = 0; i < N; ++i) rows_by_pop[ pop2idx[ popc[i] ] ].push_back(i);
  
  // Locus columns (all except pop_col)
  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)P - 1);
  for (int j = 0; j < P; ++j) if (j != pop_col) locus_cols.push_back(j);
  const int L = (int)locus_cols.size();
  
  // Storage: [locus][rep]
  std::vector<std::vector<double>> boot_FST(L, std::vector<double>(B, NA_REAL));
  std::vector<std::vector<double>> boot_FIT(L, std::vector<double>(B, NA_REAL));
  std::vector<std::vector<double>> boot_HI (L, std::vector<double>(B, NA_REAL));
  std::vector<std::vector<double>> boot_HS (L, std::vector<double>(B, NA_REAL));
  
  int threads = n_threads;
#ifdef _OPENMP
  if (threads <= 0) threads = omp_get_max_threads();
  threads = std::max(1, threads);
#endif
  
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
{
#endif
  AlleleCounter counter;
  std::vector<int> times;
  std::vector<const std::vector<int>*> groups;
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (int b = 0; b < B; ++b) {
    
    // One random stream per REPLICATE, derived from the seed and the replicate index only.
    // It used to depend on the thread number as well (seed + thread * B + b) with a dynamic
    // schedule: which thread ran which replicate changed from one run to the next, so the
    // confidence intervals were not reproducible (and depended on the number of cores).
    FastRng rng((uint64_t)seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1));
    
    // Pop-block resampling
    times.assign((size_t)R, 0);
    for (int t = 0; t < R; ++t) ++times[(size_t)rng.unif_int(R)];
    
    groups.clear();
    groups.reserve((size_t)R);
    for (int k = 0; k < R; ++k) {
      int mult = times[(size_t)k];
      if (mult <= 0 || rows_by_pop[(size_t)k].empty()) continue;
      for (int m = 0; m < mult; ++m) groups.push_back(&rows_by_pop[(size_t)k]);
    }
    
    const int num_groups = (int)groups.size();
    if (num_groups <= 1) {
      for (int ell = 0; ell < L; ++ell) {
        boot_FST[(size_t)ell][(size_t)b] = NA_REAL;
        boot_FIT[(size_t)ell][(size_t)b] = NA_REAL;
        boot_HI [(size_t)ell][(size_t)b] = NA_REAL;
        boot_HS [(size_t)ell][(size_t)b] = NA_REAL;
      }
      continue;
    }
    
    for (int ell = 0; ell < L; ++ell) {
      const int j = locus_cols[(size_t)ell];
      const int* col = col_ptr(j);
      
      counter.reset((size_t)num_groups);
      
      // Count alleles/hets
      for (int gidx = 0; gidx < num_groups; ++gidx) {
        const std::vector<int>& rows = *groups[(size_t)gidx];
        int nd = 0, nh = 0;
        
        for (int rid : rows) {
          int G = col[rid];
          int a1 = 0, a2 = 0;
          if (!decode_gt_base_with_missing(G, missing_code, a1, a2, base)) continue;
          
          ++nd;
          if (a1 != a2) {
            ++nh;
            ++counter.hetcnt[(size_t)gidx][a1];
            ++counter.hetcnt[(size_t)gidx][a2];
          }
          
          ++counter.acnt[(size_t)gidx][a1];
          ++counter.acnt[(size_t)gidx][a2];
          ++counter.acnt_total[a1];
          ++counter.acnt_total[a2];
        }
        
        counter.n_dip[(size_t)gidx] = nd;
        counter.n_het[(size_t)gidx] = nh;
        counter.n_dip_total += nd;
        counter.n_het_total += nh;
      }
      
      // Active groups
      std::vector<int> active;
      active.reserve((size_t)num_groups);
      
      for (int gidx = 0; gidx < num_groups; ++gidx) {
        if (counter.n_dip[(size_t)gidx] > 0) {
          active.push_back(gidx);
        }
      }
      
      const int r = (int)active.size();
      if (counter.n_dip_total == 0 || r <= 1) {
        boot_FST[(size_t)ell][(size_t)b] = NA_REAL;
        boot_FIT[(size_t)ell][(size_t)b] = NA_REAL;
        boot_HI [(size_t)ell][(size_t)b] = NA_REAL;
        boot_HS [(size_t)ell][(size_t)b] = NA_REAL;
        continue;
      }
      
        
        // HI (observed heterozygosity)
        const double HI = (double)counter.n_het_total / (double)counter.n_dip_total;
        
        // HS (within-group gene diversity, weighted by group size)
        double HS = 0.0;
        const double denom_glob = 2.0 * (double)counter.n_dip_total;
        for (int gidx : active) {
          const double nd = (double)counter.n_dip[(size_t)gidx];
          const double denom_g = 2.0 * nd;
          
          double sum_p2 = 0.0;
          for (const auto& kv : counter.acnt[(size_t)gidx]) {
            const double p = (double)kv.second / denom_g;
            sum_p2 += p * p;
          }
          const double Hk = 1.0 - sum_p2;
          HS += (denom_g / denom_glob) * Hk;
        }
        
        // W&C size terms
        double n_total = 0.0, sum_n2 = 0.0;
        for (int gidx : active) {
          n_total += (double)counter.n_dip[(size_t)gidx];
          sum_n2 += (double)counter.n_dip[(size_t)gidx] * (double)counter.n_dip[(size_t)gidx];
        }
        const double n_bar = n_total / (double)r;
        const double n_C = (n_total - (sum_n2 / n_total)) / ((double)r - 1.0);
        
        if (n_bar <= 1.0 || n_C <= 0.0) {
          boot_FST[(size_t)ell][(size_t)b] = NA_REAL;
          boot_FIT[(size_t)ell][(size_t)b] = NA_REAL;
          boot_HI [(size_t)ell][(size_t)b] = HI;
          boot_HS [(size_t)ell][(size_t)b] = HS;
          continue;
        }
        
        // WC84 A,Bc,C
        double A = 0.0, Bc = 0.0, C = 0.0;
        const double denom_total_gene = 2.0 * n_total;
        
        for (const auto& allele_pair : counter.acnt_total) {
          const int al = allele_pair.first;
          
          // p_bar
          const double p_bar = (double)allele_pair.second / denom_total_gene;
          
          // s2
          double s_num = 0.0;
          for (int gidx : active) {
            const double nd = (double)counter.n_dip[(size_t)gidx];
            const double denom_g = 2.0 * nd;
            double p_i = 0.0;
            
            auto it = counter.acnt[(size_t)gidx].find(al);
            if (it != counter.acnt[(size_t)gidx].end()) p_i = (double)it->second / denom_g;
            
            const double diff = p_i - p_bar;
            s_num += nd * diff * diff;
          }
          const double s2 = s_num / (n_bar * ((double)r - 1.0));
          
          // hbar_j
          double hbar_num = 0.0;
          double hbar_den = 0.0;
          
          for (int gidx : active) {
            const double nd = (double)counter.n_dip[(size_t)gidx];
            if (nd <= 0.0) continue;
            
            hbar_den += nd;
            
            auto hit = counter.hetcnt[(size_t)gidx].find(al);
            if (hit != counter.hetcnt[(size_t)gidx].end()) hbar_num += (double)hit->second;
          }
          const double hbar_j = (hbar_den > 0.0) ? (hbar_num / hbar_den) : 0.0;
          
          const double term_common =
            (p_bar * (1.0 - p_bar)) - ((((double)r - 1.0) / (double)r) * s2);
          
          const double a_j = (n_bar / n_C) *
            (s2 - ((term_common - (hbar_j / 4.0)) / (n_bar - 1.0)));
          
          const double b_j = (n_bar / (n_bar - 1.0)) *
            (term_common - (((2.0 * n_bar) - 1.0) / (4.0 * n_bar)) * hbar_j);
          
          const double c_j = hbar_j / 2.0;
          
          A  += a_j;
          Bc += b_j;
          C  += c_j;
        }
        
        const double denom_tot = A + Bc + C;
        if (!std::isfinite(denom_tot) || denom_tot == 0.0) {
          boot_FST[(size_t)ell][(size_t)b] = NA_REAL;
          boot_FIT[(size_t)ell][(size_t)b] = NA_REAL;
        } else {
          boot_FST[(size_t)ell][(size_t)b] = A / denom_tot;
          boot_FIT[(size_t)ell][(size_t)b] = (A + Bc) / denom_tot; // WC84 FIT
        }
        
        boot_HI[(size_t)ell][(size_t)b] = HI;
        boot_HS[(size_t)ell][(size_t)b] = HS;
    } // loci
  } // reps
#ifdef _OPENMP
}
#endif

auto build_BxL = [&](const std::vector<std::vector<double>>& src) {
  Rcpp::NumericMatrix out(B, L);
  for (int ell = 0; ell < L; ++ell)
    for (int b = 0; b < B; ++b)
      out(b, ell) = src[(size_t)ell][(size_t)b];
  return out;
};

Rcpp::NumericMatrix out_FST = build_BxL(boot_FST);
Rcpp::NumericMatrix out_FIT = build_BxL(boot_FIT);
Rcpp::NumericMatrix out_HI  = build_BxL(boot_HI);
Rcpp::NumericMatrix out_HS  = build_BxL(boot_HS);

// Column names
Rcpp::CharacterVector cn = Rcpp::colnames(mat_int);
if (cn.size() == P) {
  Rcpp::CharacterVector loc_names(L);
  int k = 0;
  for (int j = 0; j < P; ++j) if (j != pop_col) loc_names[k++] = cn[j];
  Rcpp::colnames(out_FST) = loc_names;
  Rcpp::colnames(out_FIT) = loc_names;
  Rcpp::colnames(out_HI)  = loc_names;
  Rcpp::colnames(out_HS)  = loc_names;
}

return Rcpp::List::create(
  _["FST"] = out_FST,
  _["FIT"] = out_FIT,
  _["HI"]  = out_HI,
  _["HS"]  = out_HS
);
}


// [[Rcpp::export]]
Rcpp::List boot_wc84_stats_popblock_cpp(
    const Rcpp::IntegerMatrix mat_int,
    const int pop_col_1based,
    const int missing_code,
    const int base,
    const int B,
    const double conf_level = 0.95,
    const int seed = 1,
    const int n_threads = 0
) {
  const int pop_col = pop_col_1based - 1;
  if (conf_level <= 0.0 || conf_level >= 1.0) Rcpp::stop("conf_level must be in (0,1).");
  
  // 1) bootstrap replicates (your existing function, pop_col is 0-based there)
  Rcpp::List boot = boot_wc84_fit_popblock_raw_cpp(mat_int, pop_col, missing_code, base, B, seed, n_threads);
  
  Rcpp::NumericMatrix bst_FST = boot["FST"];
  Rcpp::NumericMatrix bst_FIT = boot["FIT"];
  Rcpp::NumericMatrix bst_HI  = boot["HI"];
  Rcpp::NumericMatrix bst_HS  = boot["HS"];
  
  // 2) observed stats (bootstrap-side)
  Rcpp::List obs = wc84_fst_fit_hi_hs_per_locus_from_matrix(mat_int, pop_col, missing_code, base);
  Rcpp::NumericVector fst_obs = obs["FST"];
  Rcpp::NumericVector fit_obs = obs["FIT"];
  Rcpp::NumericVector hi_obs  = obs["HI"];
  Rcpp::NumericVector hs_obs  = obs["HS"];
  
  // 3) locus names
  const int P = mat_int.ncol();
  const int L = P - 1;
  Rcpp::CharacterVector cn = Rcpp::colnames(mat_int);
  Rcpp::CharacterVector loc_names(L);
  int k = 0;
  for (int j = 0; j < P; ++j) {
    if (j == pop_col) continue;
    // fallback label "L<k>" (1-based) when the matrix has no column names
    loc_names[k] = (cn.size()==P ? cn[j] : Rcpp::String("L" + std::to_string(k + 1)));
    ++k;
  }
  
  fst_obs.names() = fit_obs.names() = hi_obs.names() = hs_obs.names() = loc_names;
  
  // 4) CI helper (per locus quantiles ignoring NA/NaN)
  auto ci_from_boot = [&](const Rcpp::NumericMatrix& M) -> Rcpp::NumericMatrix {
    const double alpha = 1.0 - conf_level;
    const double lo_p = alpha / 2.0;
    const double hi_p = 1.0 - alpha / 2.0;
    
    Rcpp::NumericMatrix ci(2, L);
    Rcpp::rownames(ci) = Rcpp::CharacterVector::create("lo", "hi");
    Rcpp::colnames(ci) = loc_names;
    
    std::vector<double> vals;
    vals.reserve((size_t)B);
    
    for (int j = 0; j < L; ++j) {
      vals.clear();
      for (int b = 0; b < B; ++b) {
        const double v = M(b, j);
        if (std::isfinite(v)) vals.push_back(v);
      }
      if (vals.empty()) {
        ci(0, j) = NA_REAL; ci(1, j) = NA_REAL;
        continue;
      }
      std::sort(vals.begin(), vals.end());
      auto q = [&](double p) {
        const double idx = p * (double)(vals.size() - 1);
        const size_t i0 = (size_t)std::floor(idx);
        const size_t i1 = (size_t)std::ceil(idx);
        if (i0 == i1) return vals[i0];
        const double w = idx - (double)i0;
        return vals[i0] * (1.0 - w) + vals[i1] * w;
      };
      ci(0, j) = q(lo_p);
      ci(1, j) = q(hi_p);
    }
    return ci;
  };
  
  Rcpp::NumericMatrix CI_FIT = ci_from_boot(bst_FIT);
  Rcpp::NumericMatrix CI_HI  = ci_from_boot(bst_HI);
  Rcpp::NumericMatrix CI_HS  = ci_from_boot(bst_HS);
  
  return Rcpp::List::create(
    _["FST_boot"] = bst_FST,
    _["FIT_boot"] = bst_FIT,
    _["HI_boot"]  = bst_HI,
    _["HS_boot"]  = bst_HS,
    _["FST_obs_boot"] = fst_obs,
    _["FIT_obs_boot"] = fit_obs,
    _["HI_obs_boot"]  = hi_obs,
    _["HS_obs_boot"]  = hs_obs,
    _["CI_FIT"] = CI_FIT,
    _["CI_HI"]  = CI_HI,
    _["CI_HS"]  = CI_HS,
    _["conf_level"] = conf_level,
    _["n_boot"] = B,
    _["n_loci"] = L,
    _["locus_names"] = loc_names
  );
}
// 



// ============================================================================
// wc84_fit_boot_permute.cpp
// - Bootstrap (pop-block) WC84: FST, FIT, HI, HS (BxL matrices)
// - Permute alleles within populations: FIT_obs, FIT_perm, pvals
// - Observed FIT per locus
//
// Decoder unified on YOUR decode_gt_base() exactly.
// Bootstrap RNG: FastRng (thread-safe, OpenMP)
// Permutation RNG: R::unif_rand() (reproducible via set.seed())
// ============================================================================


// ============================================================================
// PERMUTATION SECTION (allele shuffle within populations) [from fit_permute.cpp]
// Decoder aligned on decode_gt_base + missing_code guard.
// ============================================================================



inline void fy_shuffle(std::vector<int>& x) {
  const int n = (int)x.size();
  for (int i = n - 1; i > 0; --i) {
    int j = (int)(R::unif_rand() * (i + 1));
    std::swap(x[(size_t)i], x[(size_t)j]);
  }
}

// ============================================================================
// OpenMP-safe raw-pointer helpers (no Rcpp allocation inside parallel blocks)
// ============================================================================

inline void fy_shuffle_rng(std::vector<int>& x, FastRng& rng) {
  for (int i = (int)x.size() - 1; i > 0; --i) {
    int j = rng.unif_int(i + 1);
    std::swap(x[(size_t)i], x[(size_t)j]);
  }
}

// Build strata (row indices per population) from raw column-major int*.
static std::vector<std::vector<int>> build_pop_rows_raw(
    const int* dat, int n, int pop_col
) {
  std::unordered_map<int, std::vector<int>> tmp;
  tmp.reserve(64);
  for (int i = 0; i < n; ++i) tmp[dat[i + (size_t)n * pop_col]].push_back(i);
  std::vector<std::vector<int>> out;
  out.reserve(tmp.size());
  for (auto& kv : tmp) out.push_back(std::move(kv.second));
  return out;
}

// Within-population permutation on raw column-major int*, using FastRng.
static void permute_within_pops_rng_ptr(
    int* pm, int n, int p, int pop_col, int missing_code, int base,
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
        int a1, a2;
        if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
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

// Global permutation (all individuals, ignoring population) on raw int*, using FastRng.
static void permute_globally_rng_ptr(
    int* pm, int n, int p, int pop_col, int missing_code, int base,
    FastRng& rng
) {
  std::vector<int> alleles, valid_rows;
  alleles.reserve((size_t)n * 2); valid_rows.reserve((size_t)n);
  for (int j = 0; j < p; ++j) {
    if (j == pop_col) continue;
    alleles.clear(); valid_rows.clear();
    for (int i = 0; i < n; ++i) {
      int gt = pm[i + (size_t)n * j];
      int a1, a2;
      if (!decode_gt_base_with_missing(gt, missing_code, a1, a2, base)) continue;
      alleles.push_back(a1); alleles.push_back(a2);
      valid_rows.push_back(i);
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

// Per-locus FST/FIT/HI/HS from raw column-major int*, filling output arrays.
// locus_cols[ell] = column index in dat for output locus ell.
static void wc84_fst_fit_loci_raw(
    const int* dat, int n,
    int pop_col, int missing_code, int base,
    double* fst_out, double* fit_out, double* hi_out, double* hs_out,
    const std::vector<int>& locus_cols
) {
  const int L = (int)locus_cols.size();
  for (int ell = 0; ell < L; ++ell)
    fst_out[ell] = fit_out[ell] = hi_out[ell] = hs_out[ell] = NA_REAL;

  std::unordered_map<int,int> pop2idx;
  pop2idx.reserve(64);
  for (int i = 0; i < n; ++i) {
    int pop = dat[i + (size_t)n * pop_col];
    if (pop2idx.find(pop) == pop2idx.end()) pop2idx.emplace(pop, (int)pop2idx.size());
  }
  const int R = (int)pop2idx.size();
  if (R <= 1) return;

  std::vector<std::vector<int>> rows_by_pop(R);
  for (int i = 0; i < n; ++i)
    rows_by_pop[(size_t)pop2idx[dat[i + (size_t)n * pop_col]]].push_back(i);

  std::vector<int> n_dip(R), n_het(R);
  std::vector<std::unordered_map<int,int>> acnt(R), hetcnt(R);
  std::unordered_map<int,int> acnt_total;

  for (int ell = 0; ell < L; ++ell) {
    const int j = locus_cols[ell];
    std::fill(n_dip.begin(), n_dip.end(), 0);
    std::fill(n_het.begin(), n_het.end(), 0);
    for (int g = 0; g < R; ++g) { acnt[g].clear(); hetcnt[g].clear(); }
    acnt_total.clear();
    int n_dip_total = 0;

    for (int g = 0; g < R; ++g) {
      for (int rid : rows_by_pop[g]) {
        int a1, a2;
        if (!decode_gt_base_with_missing(dat[rid + (size_t)n * j], missing_code, a1, a2, base)) continue;
        ++n_dip[g]; ++n_dip_total;
        if (a1 != a2) { ++n_het[g]; ++hetcnt[g][a1]; ++hetcnt[g][a2]; }
        ++acnt[g][a1]; ++acnt[g][a2];
        ++acnt_total[a1]; ++acnt_total[a2];
      }
    }

    std::vector<int> active;
    for (int g = 0; g < R; ++g) if (n_dip[g] > 0) active.push_back(g);
    const int r = (int)active.size();
    if (n_dip_total == 0 || r <= 1) continue;

    double n_total = 0.0, sum_n2 = 0.0;
    for (int g : active) { n_total += n_dip[g]; sum_n2 += (double)n_dip[g] * (double)n_dip[g]; }
    const double n_bar = n_total / (double)r;
    const double n_C = (n_total - sum_n2 / n_total) / ((double)r - 1.0);
    if (n_bar <= 1.0 || n_C <= 0.0) continue;

    double HI_sum = 0.0;
    for (int g : active) {
      if (n_dip[g] > 0) HI_sum += (double)n_het[g] / (double)n_dip[g];
    }
    hi_out[ell] = HI_sum / (double)r;

    double HS_sum = 0.0;
    for (int g : active) {
      const double nd = (double)n_dip[g];
      if (nd <= 1.0) continue;
      double sum_p2 = 0.0;
      for (const auto& kv : acnt[g]) {
        const double pk = (double)kv.second / (2.0 * nd);
        sum_p2 += pk * pk;
      }
      HS_sum += (1.0 - sum_p2) * (2.0 * nd / (2.0 * nd - 1.0));
    }
    hs_out[ell] = HS_sum / (double)r;

    const double denom_total_gene = 2.0 * n_total;
    double A = 0.0, Bc = 0.0, C = 0.0;
    for (const auto& ap : acnt_total) {
      const int al = ap.first;
      const double p_bar = (double)ap.second / denom_total_gene;
      double s_num = 0.0;
      for (int g : active) {
        const double nd = (double)n_dip[g];
        double p_i = 0.0;
        auto it = acnt[g].find(al);
        if (it != acnt[g].end()) p_i = (double)it->second / (2.0 * nd);
        s_num += nd * (p_i - p_bar) * (p_i - p_bar);
      }
      const double s2 = s_num / (n_bar * ((double)r - 1.0));
      double hbar_num = 0.0, hbar_den = 0.0;
      for (int g : active) {
        const double nd = (double)n_dip[g];
        if (nd <= 0.0) continue;
        hbar_den += nd;
        auto hit = hetcnt[g].find(al);
        if (hit != hetcnt[g].end()) hbar_num += (double)hit->second;
      }
      const double hbar_j = (hbar_den > 0.0) ? (hbar_num / hbar_den) : 0.0;
      const double tc = (p_bar * (1.0 - p_bar)) - ((((double)r-1.0)/(double)r) * s2);
      A  += (n_bar / n_C) * (s2 - ((tc - (hbar_j / 4.0)) / (n_bar - 1.0)));
      Bc += (n_bar / (n_bar-1.0)) * (tc - (((2.0*n_bar-1.0)/(4.0*n_bar)) * hbar_j));
      C  += hbar_j / 2.0;
    }
    const double D = A + Bc + C;
    if (std::isfinite(D) && D > 0.0) {
      fst_out[ell] = A / D;
      fit_out[ell] = (A + Bc) / D;
    }
  }
}

// Overall FIT ratio-of-sums from raw column-major int*.
static double wc84_fit_overall_raw(
    const int* dat, int n,
    int pop_col, int missing_code, int base,
    const std::vector<int>& locus_cols
) {
  std::unordered_map<int,int> pop2idx;
  pop2idx.reserve(64);
  for (int i = 0; i < n; ++i) {
    int pop = dat[i + (size_t)n * pop_col];
    if (pop2idx.find(pop) == pop2idx.end()) pop2idx.emplace(pop, (int)pop2idx.size());
  }
  const int R = (int)pop2idx.size();
  if (R <= 1) return NA_REAL;

  std::vector<std::vector<int>> rows_by_pop(R);
  for (int i = 0; i < n; ++i)
    rows_by_pop[(size_t)pop2idx[dat[i + (size_t)n * pop_col]]].push_back(i);

  double A_sum = 0.0, B_sum = 0.0, C_sum = 0.0;

  for (int lc : locus_cols) {
    std::vector<int> n_dip(R, 0), n_het(R, 0);
    std::vector<std::unordered_map<int,int>> acnt(R), hetcnt(R);
    std::unordered_map<int,int> acnt_total;
    int n_dip_total = 0;

    for (int g = 0; g < R; ++g) {
      for (int rid : rows_by_pop[g]) {
        int a1, a2;
        if (!decode_gt_base_with_missing(dat[rid + (size_t)n * lc], missing_code, a1, a2, base)) continue;
        ++n_dip[g]; ++n_dip_total;
        if (a1 != a2) { ++n_het[g]; ++hetcnt[g][a1]; ++hetcnt[g][a2]; }
        ++acnt[g][a1]; ++acnt[g][a2];
        ++acnt_total[a1]; ++acnt_total[a2];
      }
    }

    std::vector<int> active;
    for (int g = 0; g < R; ++g) if (n_dip[g] > 0) active.push_back(g);
    const int r = (int)active.size();
    if (n_dip_total == 0 || r <= 1) continue;

    double n_total = 0.0, sum_n2 = 0.0;
    for (int g : active) { n_total += n_dip[g]; sum_n2 += (double)n_dip[g] * (double)n_dip[g]; }
    const double n_bar = n_total / (double)r;
    const double n_C = (n_total - sum_n2 / n_total) / ((double)r - 1.0);
    if (n_bar <= 1.0 || n_C <= 0.0) continue;

    double A = 0.0, Bc = 0.0, C = 0.0;
    const double denom_total_gene = 2.0 * n_total;
    for (const auto& ap : acnt_total) {
      const int al = ap.first;
      const double p_bar = (double)ap.second / denom_total_gene;
      double s_num = 0.0;
      for (int g : active) {
        const double nd = (double)n_dip[g];
        double p_i = 0.0;
        auto it = acnt[g].find(al);
        if (it != acnt[g].end()) p_i = (double)it->second / (2.0 * nd);
        s_num += nd * (p_i - p_bar) * (p_i - p_bar);
      }
      const double s2 = s_num / (n_bar * ((double)r - 1.0));
      double hbar_num = 0.0, hbar_den = 0.0;
      for (int g : active) {
        const double nd = (double)n_dip[g];
        if (nd <= 0.0) continue;
        hbar_den += nd;
        auto hit = hetcnt[g].find(al);
        if (hit != hetcnt[g].end()) hbar_num += (double)hit->second;
      }
      const double hbar_j = (hbar_den > 0.0) ? (hbar_num / hbar_den) : 0.0;
      const double tc = (p_bar * (1.0 - p_bar)) - ((((double)r-1.0)/(double)r) * s2);
      A  += (n_bar / n_C) * (s2 - ((tc - (hbar_j / 4.0)) / (n_bar - 1.0)));
      Bc += (n_bar / (n_bar-1.0)) * (tc - (((2.0*n_bar-1.0)/(4.0*n_bar)) * hbar_j));
      C  += hbar_j / 2.0;
    }
    const double D = A + Bc + C;
    if (std::isfinite(D) && D > 0.0) { A_sum += A; B_sum += Bc; C_sum += C; }
  }

  const double denom = A_sum + B_sum + C_sum;
  return (denom > 0.0) ? ((A_sum + B_sum) / denom) : NA_REAL;
}

static inline void extract_alleles_base(
    const IntegerMatrix& m,
    const int locus_col,
    const std::vector<int>& rows,
    const int missing_code,
    const int base,
    std::vector<int>& alleles,
    std::vector<int>& valid_rows
) {
  alleles.clear();
  valid_rows.clear();
  
  int a1 = 0, a2 = 0;
  
  for (int rid : rows) {
    int gt = m(rid, locus_col);
    if (gt == missing_code) continue;
    if (!decode_gt_base(gt, a1, a2, base)) continue;
    
    alleles.push_back(a1);
    alleles.push_back(a2);
    valid_rows.push_back(rid);
  }
}

static void permute_within_pops_once_base(
    IntegerMatrix& m,
    const int pop_col,
    const int missing_code,
    const int base
) {
  const int N = m.nrow();
  const int P = m.ncol();
  
  std::unordered_map<int, std::vector<int>> pop_rows;
  pop_rows.reserve((size_t)N);
  
  for (int i = 0; i < N; ++i) pop_rows[m(i, pop_col)].push_back(i);
  
  std::vector<int> alleles;
  std::vector<int> valid_rows;
  alleles.reserve((size_t)N * 2);
  valid_rows.reserve((size_t)N);
  
  for (int j = 0; j < P; ++j) {
    if (j == pop_col) continue;
    
    for (auto& kv : pop_rows) {
      const std::vector<int>& rows = kv.second;
      if ((int)rows.size() < 2) continue;
      
      extract_alleles_base(m, j, rows, missing_code, base, alleles, valid_rows);
      
      const int n_alleles = (int)alleles.size();
      if (n_alleles < 2 || valid_rows.empty()) continue;
      
      fy_shuffle(alleles);
      
      int k = 0;
      for (int rid : valid_rows) {
        if (k + 1 >= n_alleles) break;
        m(rid, j) = alleles[(size_t)k] * base + alleles[(size_t)k + 1];
        k += 2;
      }
    }
  }
}

// [[Rcpp::export]]
Rcpp::List batch_permute_wc84_stats(
    const Rcpp::IntegerMatrix dat,
    const int pop_col_1based,
    const int missing_code,
    const int base,
    const int B,
    const int seed = 1
) {
  const int N = dat.nrow();
  const int P = dat.ncol();
  if (P < 2) Rcpp::stop("Need at least 2 columns: pop + >=1 locus.");
  if (B <= 0) Rcpp::stop("B must be > 0.");

  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= P) Rcpp::stop("pop_col out of range.");
  const int L = P - 1;

  // locus names
  Rcpp::CharacterVector cn = Rcpp::colnames(dat);
  Rcpp::CharacterVector loc_names(L);
  int kk = 0;
  for (int j = 0; j < P; ++j) {
    if (j == pop_col) continue;
    if (cn.size() == P) loc_names[kk] = cn[j];
    else loc_names[kk] = Rcpp::String("L" + std::to_string(kk + 1));
    ++kk;
  }

  // observed values (Rcpp call, outside parallel block)
  Rcpp::List obs = wc84_fst_fit_hi_hs_per_locus_from_matrix(dat, pop_col, missing_code, base);
  Rcpp::NumericVector fst_obs = obs["FST"];
  Rcpp::NumericVector fit_obs = obs["FIT"];
  Rcpp::NumericVector hi_obs  = obs["HI"];
  Rcpp::NumericVector hs_obs  = obs["HS"];
  if (fst_obs.size()!=L || fit_obs.size()!=L || hi_obs.size()!=L || hs_obs.size()!=L)
    Rcpp::stop("Internal error: observed vector length mismatch.");
  fst_obs.names() = fit_obs.names() = hi_obs.names() = hs_obs.names() = loc_names;

  // locus column list and strata, built once outside parallel block
  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)L);
  for (int j = 0; j < P; ++j) if (j != pop_col) locus_cols.push_back(j);

  const int* dat_ptr = &dat[0];
  auto strata = build_pop_rows_raw(dat_ptr, N, pop_col);

  // perm storage [B][L]
  std::vector<std::vector<double>> raw_fst(B, std::vector<double>(L, NA_REAL));
  std::vector<std::vector<double>> raw_fit(B, std::vector<double>(L, NA_REAL));
  std::vector<std::vector<double>> raw_hi (B, std::vector<double>(L, NA_REAL));
  std::vector<std::vector<double>> raw_hs (B, std::vector<double>(L, NA_REAL));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = (uint64_t)seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    FastRng rng(sb);

    std::vector<int> pm((size_t)N * (size_t)P);
    std::memcpy(pm.data(), dat_ptr, sizeof(int) * (size_t)N * (size_t)P);

    permute_within_pops_rng_ptr(pm.data(), N, P, pop_col, missing_code, base, strata, rng);

    wc84_fst_fit_loci_raw(pm.data(), N, pop_col, missing_code, base,
                          raw_fst[b].data(), raw_fit[b].data(),
                          raw_hi[b].data(), raw_hs[b].data(), locus_cols);
  }

  // assemble Rcpp output matrices
  Rcpp::NumericMatrix out_fst(B, L), out_fit(B, L), out_hi(B, L), out_hs(B, L);
  Rcpp::colnames(out_fst) = loc_names;
  Rcpp::colnames(out_fit) = loc_names;
  Rcpp::colnames(out_hi)  = loc_names;
  Rcpp::colnames(out_hs)  = loc_names;
  for (int b = 0; b < B; ++b)
    for (int j = 0; j < L; ++j) {
      out_fst(b, j) = raw_fst[b][j];
      out_fit(b, j) = raw_fit[b][j];
      out_hi (b, j) = raw_hi[b][j];
      out_hs (b, j) = raw_hs[b][j];
    }
  
  // p-values:
  // - FIT: your current "two_sided_abs" on |FIT|
  // - HI/HS: two-sided around permutation mean
  Rcpp::NumericVector p_fit(L, NA_REAL), p_hi(L, NA_REAL), p_hs(L, NA_REAL);
  
  for (int j = 0; j < L; ++j) {
    // ---- FIT ----
    {
      const double obs_v = fit_obs[j];
      if (std::isfinite(obs_v)) {
        const double thr = std::fabs(obs_v);
        int ge = 0, n_ok = 0;
        for (int b = 0; b < B; ++b) {
          const double v = out_fit(b, j);
          if (!std::isfinite(v)) continue;
          ++n_ok;
          if (std::fabs(v) >= thr) ++ge;
        }
        p_fit[j] = (n_ok > 0) ? ((double)(ge + 1) / (double)(n_ok + 1)) : NA_REAL;
      }
    }
    
    // ---- HI centred ----
    {
      const double obs_v = hi_obs[j];
      if (std::isfinite(obs_v)) {
        double mu = 0.0; int n_ok = 0;
        for (int b = 0; b < B; ++b) {
          const double v = out_hi(b, j);
          if (!std::isfinite(v)) continue;
          mu += v; ++n_ok;
        }
        if (n_ok > 0) {
          mu /= (double)n_ok;
          const double thr = std::fabs(obs_v - mu);
          int ge = 0;
          for (int b = 0; b < B; ++b) {
            const double v = out_hi(b, j);
            if (!std::isfinite(v)) continue;
            if (std::fabs(v - mu) >= thr) ++ge;
          }
          p_hi[j] = (double)(ge + 1) / (double)(n_ok + 1);
        }
      }
    }
    
    // ---- HS centred ----
    {
      const double obs_v = hs_obs[j];
      if (std::isfinite(obs_v)) {
        double mu = 0.0; int n_ok = 0;
        for (int b = 0; b < B; ++b) {
          const double v = out_hs(b, j);
          if (!std::isfinite(v)) continue;
          mu += v; ++n_ok;
        }
        if (n_ok > 0) {
          mu /= (double)n_ok;
          const double thr = std::fabs(obs_v - mu);
          int ge = 0;
          for (int b = 0; b < B; ++b) {
            const double v = out_hs(b, j);
            if (!std::isfinite(v)) continue;
            if (std::fabs(v - mu) >= thr) ++ge;
          }
          p_hs[j] = (double)(ge + 1) / (double)(n_ok + 1);
        }
      }
    }
  }
  
  p_fit.names() = p_hi.names() = p_hs.names() = loc_names;
  
  return Rcpp::List::create(
    _["FST_obs"]  = fst_obs,
    _["FIT_obs"]  = fit_obs,
    _["HI_obs"]   = hi_obs,
    _["HS_obs"]   = hs_obs,
    _["FST_perm"] = out_fst,
    _["FIT_perm"] = out_fit,
    _["HI_perm"]  = out_hi,
    _["HS_perm"]  = out_hs,
    _["p_FIT"]    = p_fit,
    _["p_HI"]     = p_hi,
    _["p_HS"]     = p_hs,
    _["pval_method_FIT"] = "two_sided_abs",
    _["pval_method_HI"]  = "two_sided_centered_mean",
    _["pval_method_HS"]  = "two_sided_centered_mean",
    _["n_perm"]   = B,
    _["n_loci"]   = L,
    _["locus_names"] = loc_names
  );
}



// ============================================================================
// Global allele shuffle (across ALL individuals, ignoring population).
// This is FSTAT's "Randomising alleles overall samples" — the correct null
// for FIT.  For each locus, all alleles from all valid individuals are pooled
// and reshuffled, then redistributed back to the same valid rows.
// ============================================================================
static void permute_globally_once_base(
    IntegerMatrix& m,
    const int pop_col,
    const int missing_code,
    const int base
) {
  const int N = m.nrow();
  const int P = m.ncol();
  
  // collect all valid row indices once (pop column never changes)
  std::vector<int> all_rows;
  all_rows.reserve((size_t)N);
  for (int i = 0; i < N; ++i) all_rows.push_back(i);
  
  std::vector<int> alleles;
  std::vector<int> valid_rows;
  alleles.reserve((size_t)N * 2);
  valid_rows.reserve((size_t)N);
  
  for (int j = 0; j < P; ++j) {
    if (j == pop_col) continue;
    
    // pool alleles from ALL individuals
    extract_alleles_base(m, j, all_rows, missing_code, base, alleles, valid_rows);
    
    const int n_alleles = (int)alleles.size();
    if (n_alleles < 2 || valid_rows.empty()) continue;
    
    fy_shuffle(alleles);  // global Fisher-Yates
    
    int k = 0;
    for (int rid : valid_rows) {
      if (k + 1 >= n_alleles) break;
      m(rid, j) = alleles[(size_t)k] * base + alleles[(size_t)k + 1];
      k += 2;
    }
  }
}

// Helper: WC84 overall FIT as ratio-of-sums — sum(A+B) / sum(A+B+C) across loci.
// Needed because wc84_fst_fit_hi_hs_per_locus_from_matrix() returns per-locus
// values only and the overall is NOT the mean of those ratios.
static double wc84_fit_overall_ratio_of_sums(
    const Rcpp::IntegerMatrix& mat_int,
    const int pop_col,
    const int missing_code,
    const int base
) {
  Rcpp::List st = wc84_fst_fit_hi_hs_per_locus_from_matrix(
    mat_int, pop_col, missing_code, base);
  
  Rcpp::NumericVector fit = st["FIT"];
  Rcpp::NumericVector fst = st["FST"];
  // FIT_l = (A+B)/D,  FST_l = A/D  =>  B/D = FIT_l - FST_l
  // We need sum(A+B) / sum(D).  Recover D from FST: A = FST_l * D
  // Easier: use FIT and FST per locus to recover ratio components.
  // Since FIT_l = (A+B)/D and FST_l = A/D:
  //   A+B = FIT_l * D,  D unknown per locus.
  // The ratio-of-sums cannot be recovered from FIT alone without D.
  // Instead, compute directly from the WC84 variance components.
  // We replicate the inner loop here for the overall accumulation.
  
  // (Fall back to simple approach: compute from scratch accumulating A,B,C)
  const int N = mat_int.nrow();
  const int P = mat_int.ncol();
  
  std::unordered_map<int,int> pop2idx;
  pop2idx.reserve((size_t)N);
  for (int i = 0; i < N; ++i) {
    const int p = mat_int(i, pop_col);
    if (pop2idx.find(p) == pop2idx.end()) pop2idx.emplace(p, (int)pop2idx.size());
  }
  const int R = (int)pop2idx.size();
  if (R <= 1) return NA_REAL;
  
  std::vector<std::vector<int>> rows_by_pop(R);
  for (int i = 0; i < N; ++i)
    rows_by_pop[(size_t)pop2idx[mat_int(i, pop_col)]].push_back(i);
  
  double A_sum = 0.0, B_sum = 0.0, C_sum = 0.0;
  
  for (int j = 0; j < P; ++j) {
    if (j == pop_col) continue;
    
    std::vector<int> n_dip(R, 0), n_het(R, 0);
    std::vector<std::unordered_map<int,int>> acnt(R), hetcnt(R);
    std::unordered_map<int,int> acnt_total;
    int n_dip_total = 0, n_het_total = 0;
    
    for (int gidx = 0; gidx < R; ++gidx) {
      for (int rid : rows_by_pop[(size_t)gidx]) {
        int a1 = 0, a2 = 0;
        if (!decode_gt_base_with_missing(mat_int(rid, j), missing_code, a1, a2, base)) continue;
        ++n_dip[(size_t)gidx]; ++n_dip_total;
        if (a1 != a2) { ++n_het[(size_t)gidx]; ++n_het_total;
        ++hetcnt[(size_t)gidx][a1]; ++hetcnt[(size_t)gidx][a2]; }
        ++acnt[(size_t)gidx][a1]; ++acnt[(size_t)gidx][a2];
        ++acnt_total[a1]; ++acnt_total[a2];
      }
    }
    
    std::vector<int> active;
    for (int g = 0; g < R; ++g) if (n_dip[(size_t)g] > 0) active.push_back(g);
    const int r = (int)active.size();
    if (n_dip_total == 0 || r <= 1) continue;
    
    double n_total = 0.0, sum_n2 = 0.0;
    for (int g : active) {
      n_total += (double)n_dip[(size_t)g];
      sum_n2  += (double)n_dip[(size_t)g] * (double)n_dip[(size_t)g];
    }
    const double n_bar = n_total / (double)r;
    const double n_C   = (n_total - (sum_n2 / n_total)) / ((double)r - 1.0);
    if (n_bar <= 1.0 || n_C <= 0.0) continue;
    
    const double denom_total_gene = 2.0 * n_total;
    double A = 0.0, Bc = 0.0, C = 0.0;
    
    for (const auto& ap : acnt_total) {
      const int al = ap.first;
      const double p_bar = (double)ap.second / denom_total_gene;
      double s_num = 0.0;
      for (int g : active) {
        const double nd = (double)n_dip[(size_t)g];
        double p_i = 0.0;
        auto it = acnt[(size_t)g].find(al);
        if (it != acnt[(size_t)g].end()) p_i = (double)it->second / (2.0 * nd);
        s_num += nd * (p_i - p_bar) * (p_i - p_bar);
      }
      const double s2 = s_num / (n_bar * ((double)r - 1.0));
      double hbar_num = 0.0, hbar_den = 0.0;
      for (int g : active) {
        const double nd = (double)n_dip[(size_t)g];
        if (nd <= 0.0) continue;
        hbar_den += nd;
        auto hit = hetcnt[(size_t)g].find(al);
        if (hit != hetcnt[(size_t)g].end()) hbar_num += (double)hit->second;
      }
      const double hbar_j = (hbar_den > 0.0) ? (hbar_num / hbar_den) : 0.0;
      const double tc = (p_bar * (1.0 - p_bar)) - ((((double)r-1.0)/(double)r) * s2);
      A  += (n_bar / n_C)         * (s2 - ((tc - (hbar_j / 4.0)) / (n_bar - 1.0)));
      Bc += (n_bar / (n_bar-1.0)) * (tc - (((2.0*n_bar-1.0)/(4.0*n_bar)) * hbar_j));
      C  += hbar_j / 2.0;
    }
    const double D = A + Bc + C;
    if (std::isfinite(D) && D > 0.0) { A_sum += A; B_sum += Bc; C_sum += C; }
  }
  
  const double denom = A_sum + B_sum + C_sum;
  return (denom > 0.0) ? ((A_sum + B_sum) / denom) : NA_REAL;
}

// batch_permute_fit_global:
// Permutes alleles GLOBALLY (across all individuals, ignoring population) —
// FSTAT's "Randomising alleles overall samples" null for FIT.
// Returns per-locus FIT permutation matrix + overall ratio-of-sums + p-values.
// [[Rcpp::export]]
Rcpp::List batch_permute_fit_global(
    const Rcpp::IntegerMatrix dat,
    const int pop_col_1based,
    const int missing_code,
    const int base,
    const int B,
    const int seed = 1
) {
  const int N = dat.nrow();
  const int P = dat.ncol();
  if (P < 2) Rcpp::stop("Need at least 2 columns: pop + >=1 locus.");
  if (B <= 0) Rcpp::stop("B must be > 0.");
  const int pop_col = pop_col_1based - 1;
  if (pop_col < 0 || pop_col >= P) Rcpp::stop("pop_col out of range.");
  const int L = P - 1;

  // locus names
  Rcpp::CharacterVector cn = Rcpp::colnames(dat);
  Rcpp::CharacterVector loc_names(L);
  int kk = 0;
  for (int j = 0; j < P; ++j) {
    if (j == pop_col) continue;
    loc_names[kk] = (cn.size() == P) ? cn[j] : Rcpp::String("L" + std::to_string(kk + 1));
    ++kk;
  }

  // observed (Rcpp calls outside parallel block)
  Rcpp::List obs = wc84_fst_fit_hi_hs_per_locus_from_matrix(dat, pop_col, missing_code, base);
  Rcpp::NumericVector fit_obs = obs["FIT"];
  fit_obs.names() = loc_names;

  std::vector<int> locus_cols;
  locus_cols.reserve((size_t)L);
  for (int j = 0; j < P; ++j) if (j != pop_col) locus_cols.push_back(j);

  const int* dat_ptr = &dat[0];
  const double fit_obs_overall = wc84_fit_overall_raw(dat_ptr, N, pop_col, missing_code, base, locus_cols);

  // permutation storage
  std::vector<std::vector<double>> raw_fit(B, std::vector<double>(L, NA_REAL));
  std::vector<double> raw_fit_overall(B, NA_REAL);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < B; ++b) {
    const uint64_t sb = (uint64_t)seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1);
    FastRng rng(sb);

    std::vector<int> pm((size_t)N * (size_t)P);
    std::memcpy(pm.data(), dat_ptr, sizeof(int) * (size_t)N * (size_t)P);

    permute_globally_rng_ptr(pm.data(), N, P, pop_col, missing_code, base, rng);

    std::vector<double> fst_tmp(L, NA_REAL), hs_tmp(L, NA_REAL), hi_tmp(L, NA_REAL);
    wc84_fst_fit_loci_raw(pm.data(), N, pop_col, missing_code, base,
                          fst_tmp.data(), raw_fit[b].data(),
                          hi_tmp.data(), hs_tmp.data(), locus_cols);

    raw_fit_overall[b] = wc84_fit_overall_raw(pm.data(), N, pop_col, missing_code, base, locus_cols);
  }

  // assemble Rcpp output
  Rcpp::NumericMatrix out_fit(B, L);
  Rcpp::colnames(out_fit) = loc_names;
  Rcpp::NumericVector out_fit_overall(B, NA_REAL);
  for (int b = 0; b < B; ++b) {
    for (int j = 0; j < L; ++j) out_fit(b, j) = raw_fit[b][j];
    out_fit_overall[b] = raw_fit_overall[b];
  }

  // p-values (two-sided: |perm| >= |obs|, consistent with FIS permutation test)
  Rcpp::NumericVector p_fit(L, NA_REAL);
  for (int j = 0; j < L; ++j) {
    const double obs_v = std::abs(fit_obs[j]);
    if (!std::isfinite(obs_v)) continue;
    int ge = 0, n_ok = 0;
    for (int b = 0; b < B; ++b) {
      const double v = out_fit(b, j);
      if (!std::isfinite(v)) continue;
      ++n_ok;
      if (std::abs(v) >= obs_v - TIE_TOL) ++ge;
    }
    p_fit[j] = (n_ok > 0) ? ((double)(ge + 1) / (double)(n_ok + 1)) : NA_REAL;
  }
  p_fit.names() = loc_names;

  double p_fit_overall = NA_REAL;
  if (std::isfinite(fit_obs_overall)) {
    const double obs_abs = std::abs(fit_obs_overall);
    int ge = 0, n_ok = 0;
    for (int b = 0; b < B; ++b) {
      const double v = out_fit_overall[b];
      if (!std::isfinite(v)) continue;
      ++n_ok;
      if (std::abs(v) >= obs_abs - TIE_TOL) ++ge;
    }
    if (n_ok > 0) p_fit_overall = (double)(ge + 1) / (double)(n_ok + 1);
  }

  return Rcpp::List::create(
    _["FIT_obs"]         = fit_obs,
    _["FIT_obs_overall"] = fit_obs_overall,
    _["FIT_perm"]        = out_fit,
    _["FIT_perm_overall"]= out_fit_overall,
    _["p_FIT"]           = p_fit,
    _["p_FIT_overall"]   = p_fit_overall,
    _["pval_method"]     = "two_sided_abs_global_allele_shuffle",
    _["n_perm"]          = B,
    _["locus_names"]     = loc_names
  );
}