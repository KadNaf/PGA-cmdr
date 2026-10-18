// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// ═════════════════════════════════════════════════════════════════════════
//  Mantel permutation test — native C++ engine for the "corrected proportion"
//  p-value formula (b+1)/(m+1) (Davison & Hinkley 1997; avoids a p-value of
//  exactly 0 or 1 from a finite number of replicates).
//    - Test statistic: raw sum of cross-products x*y over all valid pairs
//      (Mantel 1967) — for a fixed X, permuting only Y, this is an exact
//      linear function of both the Pearson correlation and the OLS
//      regression slope, so ranking on it gives an identical p-value to
//      testing directly on either of those (established analytically).
//    - Permutation: joint row/column relabelling (Fisher-Yates shuffle)
//      driven by R's OWN random-number stream (Rcpp::RNGScope + unif_rand()),
//      so results are fully reproducible with the app's usual set.seed()
//      and land on the same statistical footing as the pure-R engine —
//      this C++ path exists to run the permutation loop natively for speed,
//      not to bit-match any external program's specific RNG-to-integer
//      mapping algorithm (unlike src/mantel_genepop.cpp, which targets that
//      instead, for the separate b/m formula tab).
// ═════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List mantel_plus1_cpp(NumericMatrix Xmat, NumericMatrix Ymat, int nperm) {
  RNGScope scope;  // synchronise with R's set.seed()

  int n = Xmat.nrow();
  if (Xmat.ncol() != n || Ymat.nrow() != n || Ymat.ncol() != n)
    stop("Xmat and Ymat must be square matrices of the same size.");

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

  std::vector<int> p(n);
  for (int i = 0; i < n; i++) p[i] = i;

  NumericVector permStats(nperm);
  const double EPS = 1e-9;
  int bpos = 0, bneg = 0;

  for (int rep = 0; rep < nperm; rep++) {
    // Standard reverse Fisher-Yates shuffle, unbiased random permutation,
    // driven by R's unif_rand() so the sequence is deterministic under
    // set.seed() from the calling R session.
    for (int t = n - 1; t > 0; t--) {
      int k = (int) std::floor(unif_rand() * (t + 1));
      std::swap(p[k], p[t]);
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

  double pPos = (double) (bpos + 1) / (double) (nperm + 1);
  double pNeg = (double) (bneg + 1) / (double) (nperm + 1);

  return List::create(
    _["obs_cross"]  = obsCross,
    _["p_pos"]      = pPos,
    _["p_neg"]      = pNeg,
    _["perm_stats"] = permStats,
    _["n_pairs"]    = npairs,
    _["n_pops"]     = n
  );
}
