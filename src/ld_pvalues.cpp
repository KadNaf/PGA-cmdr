// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// TIE TOLERANCE for the permutation p-values.
// The statistic is sum(x*log(x)) over the cells of a genotype x genotype table of a
// few dozen individuals: it takes very few distinct values, so a permuted table very
// often gives EXACTLY the observed statistic (a "tie"). FSTAT defines the p-value as
// the proportion of randomised statistics LARGER OR EQUAL to the observed one, ties
// included (manual, section 8.2). The same value can however come out of two
// different summation orders differing by ~1e-13, which made the comparison
// `perm >= obs` lose a random part of the ties and understate the p-values (found by
// comparing with FSTAT). Genuinely different statistics differ by far more than this.
static const double TIE_TOL = 1e-9;

// ============================================================================
// ld_pvalues.cpp — pairwise linkage disequilibrium, used by the LD module.
//
// ld_pvalues_cpp() tests genotypic disequilibrium between every pair of loci:
// for each population, builds the genotype x genotype contingency table for
// the two loci, computes the log-likelihood-ratio G statistic
// (g_stat_from_counts), then permutes the second locus's genotypes among
// individuals WITHIN that population and recomputes G, nbperms-1 times.
// The per-population p-value is (# permuted G >= observed G + 1)/nbperms.
// The combined "All" p-value sums the observed G across all populations and
// compares it to the SUM of the permuted G's from the same permutation round
// (so populations are permuted independently but summed jointly per round).
//
// KNOWN OPEN ISSUE (2026-09, not yet resolved — see the project's own
// verification report): p-values from this function have been compared
// against FSTAT's reference output for the same dataset and diverge by up
// to ~0.28 in the worst case (mean absolute difference ~0.09 across 15 pairs),
// far more than the ~0.005 sampling noise expected at nbperms=10000. FSTAT /
// GENEPOP estimate this same G-test p-value via a MARKOV CHAIN algorithm
// (Guo & Thompson 1992's "switch" moves between tables with fixed margins),
// not by independently re-shuffling genotypes each round as done here. Both
// approaches are valid estimators of the same target quantity in theory, but
// converge differently — this may explain the gap, or there may be a second,
// separate issue in how the "All" combined test is computed (sum of G
// statistics vs Genepop's own global-test method, which has not been
// confirmed). Before touching this function, get FSTAT's PER-POPULATION
// LD p-values (not just "All") for the same pairs to isolate whether the
// per-population G-test itself matches (it likely does, being a standard
// contingency-table statistic) or whether the "All" combination step is the
// actual source of the discrepancy.
// ============================================================================

// Fast thread-safe PRNG (xorshift64 Marsaglia)
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

// --- G = 2 * sum(n_ij * ln(n_ij / E_ij)) ---
static inline double g_stat_from_counts(const std::vector<int> &cnt, int nr, int nc) {
  std::vector<int> rs(nr,0), cs(nc,0);
  int n = 0;
  for (int i=0;i<nr;i++) for (int j=0;j<nc;j++) {
    int v = cnt[i*nc + j]; rs[i]+=v; cs[j]+=v; n += v;
  }
  if (n==0) return NA_REAL;
  double G=0.0;
  for (int i=0;i<nr;i++) for (int j=0;j<nc;j++) {
    int nij = cnt[i*nc + j]; if (!nij) continue;
    double E = (double)rs[i] * (double)cs[j] / (double)n;
    if (E<=0.0) continue;
    G += nij * std::log((double)nij / E);
  }
  return 2.0*G;
}

static inline bool informative_table(const std::vector<int> &cnt, int nr, int nc) {
  int nonzero=0, nzr=0, nzc=0;
  std::vector<int> rs(nr,0), cs(nc,0);
  for (int i=0;i<nr;i++) for (int j=0;j<nc;j++) {
    int v = cnt[i*nc + j];
    if (v>0) nonzero++;
    rs[i]+=v; cs[j]+=v;
  }
  for (int i=0;i<nr;i++) if (rs[i]>0) nzr++;
  for (int j=0;j<nc;j++) if (cs[j]>0) nzc++;
  return (nonzero>1 && nzr>1 && nzc>1);
}

// Accepts a packed-integer genotype matrix (gt = a*base + b, 0 / NA = missing).
// Genotypes are decoded and normalised to min(a,b)*base + max(a,b) before
// building contingency tables, so allele order in the original data does not matter.
// [[Rcpp::export]]
DataFrame ld_pvalues_cpp(const StringVector  &Population,
                         const IntegerMatrix &geno_mat,
                         const int            base,
                         const int            nbperms = 10000,
                         const int            seed = 1) {

  const int N = geno_mat.nrow();
  const int M = geno_mat.ncol();

  // map population label → index
  std::unordered_map<std::string,int> pop2id;
  std::vector<std::string> pops;
  pops.reserve(16);
  std::vector<int> pop_idx(N, -1);
  for (int i=0; i<N; i++) {
    std::string p = as<std::string>(Population[i]);
    auto it = pop2id.find(p);
    if (it == pop2id.end()) {
      int id = (int)pops.size();
      pop2id.emplace(p, id);
      pops.push_back(p);
      pop_idx[i] = id;
    } else {
      pop_idx[i] = it->second;
    }
  }
  const int P = (int)pops.size();

  const CharacterVector mnames = colnames(geno_mat);
  const int nPairs = (M*(M-1))/2;
  std::vector<std::string> pair_names(nPairs);
  {
    int r = 0;
    for (int a=0; a<M-1; a++)
      for (int b=a+1; b<M; b++)
        pair_names[r++] = std::string(mnames[a]) + " X " + std::string(mnames[b]);
  }

  std::vector< std::vector<double> > RES(nPairs, std::vector<double>(P+1, NA_REAL));

  // pre-build pair list for OpenMP parallel iteration
  std::vector<std::pair<int,int>> pair_list;
  pair_list.reserve((size_t)nPairs);
  for (int a=0; a<M-1; a++)
    for (int b=a+1; b<M; b++)
      pair_list.push_back({a, b});

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int r = 0; r < nPairs; ++r) {
    const int a = pair_list[r].first;
    const int b = pair_list[r].second;

    // per-pair deterministic thread-safe RNG
    const uint64_t sb = (uint64_t)seed + 0x9e3779b97f4a7c15ULL * (uint64_t)(r + 1);
    FastRng rng(sb);

    std::vector< std::vector<int> > row_ids(P), col_ids(P);
    std::vector< std::vector<int> > cnts(P);
    std::vector<int>    nr_v(P, 0), nc_v(P, 0);
    std::vector<double> Gob(P, NA_REAL);
    std::vector<char>   valid(P, 0);

    for (int p=0; p<P; p++) {
      std::unordered_map<int,int> map1; map1.reserve(32);
      std::unordered_map<int,int> map2; map2.reserve(32);
      std::vector<int> rids; rids.reserve(128);
      std::vector<int> cids; cids.reserve(128);

      for (int i=0; i<N; i++) {
        if (pop_idx[i] != p) continue;
        int gt1 = geno_mat(i, a);
        int gt2 = geno_mat(i, b);
        if (IntegerVector::is_na(gt1) || IntegerVector::is_na(gt2)) continue;
        if (gt1 <= 0 || gt2 <= 0) continue;

        // decode packed int → allele pair, normalise to min/max order
        int a1 = gt1 / base, a2 = gt1 % base;
        if (a1 <= 0 || a2 <= 0) continue;
        if (a1 > a2) std::swap(a1, a2);
        int key1 = a1 * base + a2;

        int b1 = gt2 / base, b2 = gt2 % base;
        if (b1 <= 0 || b2 <= 0) continue;
        if (b1 > b2) std::swap(b1, b2);
        int key2 = b1 * base + b2;

        auto it1 = map1.find(key1);
        int id1;
        if (it1 == map1.end()) { id1 = (int)map1.size(); map1.emplace(key1, id1); }
        else id1 = it1->second;

        auto it2 = map2.find(key2);
        int id2;
        if (it2 == map2.end()) { id2 = (int)map2.size(); map2.emplace(key2, id2); }
        else id2 = it2->second;

        rids.push_back(id1);
        cids.push_back(id2);
      }

      if ((int)rids.size() < 2) continue;

      int R = (int)map1.size(), C = (int)map2.size();
      std::vector<int> T(R*C, 0);
      for (size_t k=0; k<rids.size(); k++) T[rids[k]*C + cids[k]]++;

      if (!informative_table(T, R, C)) continue;

      valid[p] = 1;
      nr_v[p] = R; nc_v[p] = C;
      cnts[p].swap(T);
      row_ids[p].swap(rids);
      col_ids[p].swap(cids);
      Gob[p] = g_stat_from_counts(cnts[p], R, C);
    }

    bool any_valid = false;
    for (int p=0; p<P; p++) if (valid[p]) { any_valid = true; break; }
    if (!any_valid) continue;

    double Gall_obs = 0.0;
    for (int p=0; p<P; p++) if (valid[p] && !R_IsNA(Gob[p])) Gall_obs += Gob[p];

    std::vector<int> ge_pop(P, 0);
    int ge_all = 0;

    if (nbperms > 0) {
      for (int bperm=0; bperm<nbperms-1; bperm++) {
        double s_all = 0.0;
        for (int p=0; p<P; p++) {
          if (!valid[p]) continue;
          const int R = nr_v[p], C = nc_v[p];
          std::vector<int> &T = cnts[p];
          std::fill(T.begin(), T.end(), 0);
          std::vector<int> sh = col_ids[p];
          for (int k=(int)sh.size()-1; k>0; --k) {
            int j = rng.unif_int(k + 1);
            std::swap(sh[k], sh[j]);
          }
          for (size_t k=0; k<row_ids[p].size(); k++)
            T[ row_ids[p][k]*C + sh[k] ]++;
          double Gp = g_stat_from_counts(T, R, C);
          s_all += Gp;
          // Ties count as "larger or equal" (FSTAT): compare with a tolerance, see TIE_TOL.
          if (Gp >= Gob[p] - TIE_TOL) ge_pop[p]++;
        }
        if (s_all >= Gall_obs - TIE_TOL) ge_all++;
      }
      for (int p=0; p<P; p++)
        RES[r][p] = valid[p] ? ((double)ge_pop[p] + 1.0) / (double)nbperms : NA_REAL;
      RES[r][P] = ((double)ge_all + 1.0) / (double)nbperms;
    }
  }

  // assemble result DataFrame
  List cols_list(P+1);
  CharacterVector cn(P+1);
  for (int p=0; p<P; p++) cn[p] = pops[p];
  cn[P] = "All";

  for (int c=0; c<P+1; c++) {
    NumericVector v(nPairs, NA_REAL);
    for (int i=0; i<nPairs; i++) v[i] = RES[i][c];
    cols_list[c] = v;
  }

  DataFrame out(cols_list);
  out.attr("names")     = cn;
  out.attr("row.names") = pair_names;
  out.attr("class")     = "data.frame";

  return out;
}
