// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <random>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// ═════════════════════════════════════════════════════════════════════════
//  Mantel permutation test replicating Genepop's own mantelTest() algorithm,
//  checked directly against Genepop's C++ source (2026-08-12):
//    - src/F_est.cpp: joint permutation of population labels via a
//      Fisher-Yates shuffle: for(t=0;t<n-1;t++){k=aleam.randInt(n-t-1)+t;
//      swap(p[k],p[t]);}, test statistic = raw sum of cross-products x*y
//      (Mantel's 1967 original statistic, not a correlation coefficient),
//      one-sided p-values = plain b/m (NO +1/+1 correction, unlike
//      vegan::mantel()/ade4::mantel.rtest()).
//    - src/MersenneTwister.h: the "MTRand" wrapper is std::mt19937 seeded
//      deterministically, with bounded draws taken via a FRESH
//      std::uniform_int_distribution<unsigned long>(0, n) object on every
//      single call (never reused/cached) -- replicated exactly below.
//    - R packages are normally compiled with Rtools' MinGW-w64 GCC on
//      Windows, the same toolchain family that most likely built
//      Genepop.exe, so this C++ implementation has the best practical
//      chance of reproducing Genepop's permutation sequence bit-for-bit --
//      something pure R code cannot do, since R's own sample()/sample.int()
//      map the shared Mersenne-Twister bitstream to bounded integers with a
//      different algorithm than C++'s std::uniform_int_distribution.
//    - X and Y are passed in as plain symmetric n x n matrices. Pearson,
//      Rousset's slope, and Spearman are all handled by this SAME engine:
//      for Spearman, R rank-transforms X and Y before calling this function
//      (rank(x), rank(y)) -- ranking then testing on the raw cross-product
//      sum is exactly equivalent to a Pearson-style Mantel test on ranks,
//      which is what "Spearman Mantel" means, matching Genepop's own
//      idxsup()/idxinf() rank-then-cross-product approach for rankBool=true.
// ═════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List mantel_genepop_cpp(NumericMatrix Xmat, NumericMatrix Ymat,
                         int nperm, double seedVal) {
  int n = Xmat.nrow();
  if (Xmat.ncol() != n || Ymat.nrow() != n || Ymat.ncol() != n)
    stop("Xmat and Ymat must be square matrices of the same size.");

  // Collect valid (finite x & y) lower-triangle pairs (row i > col j)
  std::vector<int> pi_, pj_;
  std::vector<double> xv, yv;
  pi_.reserve((size_t) n * (size_t) n / 2);
  pj_.reserve((size_t) n * (size_t) n / 2);
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      double x = Xmat(i, j), y = Ymat(i, j);
      if (R_finite(x) && R_finite(y)) {
        pi_.push_back(i); pj_.push_back(j);
        xv.push_back(x);  yv.push_back(y);
      }
    }
  }
  int npairs = (int) xv.size();
  if (npairs < 3) stop("Not enough complete pairs to run the Mantel test.");

  double obsCross = 0.0;
  for (int k = 0; k < npairs; k++) obsCross += xv[k] * yv[k];

  // Genepop's own generator: std::mt19937 seeded deterministically; bounded
  // draws via a FRESH std::uniform_int_distribution built on every single
  // call (matches MTRand::randInt(n) exactly -- Genepop's own code also
  // constructs a new distribution object per call, not once and reused).
  std::mt19937 rng(static_cast<unsigned long>(seedVal));

  std::vector<int> p(n);
  for (int i = 0; i < n; i++) p[i] = i;

  NumericVector permStats(nperm);
  const double EPS = 1e-9;
  int bpos = 0, bneg = 0;

  for (int rep = 0; rep < nperm; rep++) {
    // Fisher-Yates shuffle, identical loop structure to Genepop's own:
    //   for(t=0;t<n-1;t++){ k=aleam.randInt(n-t-1)+t; swap(p[k],p[t]); }
    for (int t = 0; t < n - 1; t++) {
      unsigned long m = (unsigned long)(n - t - 1);
      std::uniform_int_distribution<unsigned long> dist(0, m);
      unsigned long k = dist(rng) + (unsigned long) t;
      std::swap(p[(int) k], p[t]);
    }
    double s = 0.0;
    for (int k = 0; k < npairs; k++) {
      int a = p[pi_[k]], b = p[pj_[k]];
      int lo = std::min(a, b), hi = std::max(a, b);
      s += xv[k] * Ymat(hi, lo);
    }
    permStats[rep] = s;
    if (s >= obsCross - EPS) bpos++;
    if (s <= obsCross + EPS) bneg++;
  }

  double pPos = (double) bpos / (double) nperm;
  double pNeg = (double) bneg / (double) nperm;

  return List::create(
    _["obs_cross"]  = obsCross,
    _["p_pos"]      = pPos,
    _["p_neg"]      = pNeg,
    _["perm_stats"] = permStats,
    _["n_pairs"]    = npairs,
    _["n_pops"]     = n
  );
}
