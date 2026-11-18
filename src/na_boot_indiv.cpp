// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// ═════════════════════════════════════════════════════════════════════════
//  Null-allele pipeline in native code (Null Alleles module).
//
//  For one data set (or one bootstrap replicate) this reproduces, step by
//  step, the R functions of server_null_alleles.R:
//    em_freena()               EM estimate of the null allele frequency (rd)
//    compute_fst_global_full() global + per-locus FST (raw and ENA)
//    compute_fst_pairwise()    pairwise FST (raw and ENA)
//    compute_dc_pairwise()     pairwise DCSE (raw and INA)
//  The R versions stay the reference (observed values); this file exists so
//  that the bootstrap over individuals (EM re-run on every replicate) takes
//  seconds instead of minutes. na_stats_cpp() computes the observed values
//  and is compared with the R functions in the tests.
// ═════════════════════════════════════════════════════════════════════════

namespace {

struct Cell {                       // one locus x one sub-sample
  int efpop = 0, n_absent = 0, n_null_homo = 0;
  double rd = 0.0;
  std::vector<int> al;              // sorted alleles
  std::vector<double> gfo, pf, hii, hix;
};

void em_cell(const std::vector<int>& g, int base, bool want_null_homo, Cell& out) {
  out = Cell();
  const int n = (int)g.size();
  out.efpop = n;
  std::vector<int> a1(n, 0), a2(n, 0);
  std::vector<char> blank(n, 0);
  int nblank = 0;
  for (int i = 0; i < n; ++i) {
    const int v = g[i];
    if (v == NA_INTEGER || v <= 0) { blank[i] = 1; ++nblank; continue; }
    a1[i] = v / base; a2[i] = v % base;
    if (a1[i] == 999 && a2[i] == 999) { blank[i] = 1; ++nblank; }
  }
  bool nh = false;
  if (want_null_homo && nblank > 0) { out.n_absent = 0; out.n_null_homo = nblank; nh = true; }
  else { out.n_absent = nblank; out.n_null_homo = 0; }
  if (n - nblank == 0) return;

  std::vector<int> al;
  for (int i = 0; i < n; ++i) if (!blank[i]) {
    if (a1[i] >= 0 && a1[i] != 999) al.push_back(a1[i]);
    if (a2[i] >= 0 && a2[i] != 999) al.push_back(a2[i]);
  }
  std::sort(al.begin(), al.end()); al.erase(std::unique(al.begin(), al.end()), al.end());
  const int N = n - out.n_absent;
  if (N == 0 || al.empty()) return;
  const int nvalid = N - out.n_null_homo;
  if (nvalid == 0) return;

  const int K = (int)al.size();
  std::vector<double> gfo(K), hii(K), hix(K);
  for (int k = 0; k < K; ++k) {
    const int a = al[k];
    int c1 = 0, c2 = 0, ho = 0, he = 0;
    for (int i = 0; i < n; ++i) if (!blank[i]) {
      c1 += (a1[i] == a); c2 += (a2[i] == a);
      if (a1[i] == a && a2[i] == a) ++ho;
      else if (a1[i] == a || a2[i] == a) ++he;
    }
    gfo[k] = (double)(c1 + c2) / (2.0 * nvalid);
    hii[k] = ho; hix[k] = he;
  }
  double hotot = 0; for (int k = 0; k < K; ++k) hotot += hii[k];
  const double Nd = N, nnh = out.n_null_homo;
  double rd = nh ? std::sqrt(nnh / Nd) : std::sqrt(1.0 / (Nd + 1.0));
  std::vector<double> p(K, 0.0), np(K, 0.0);
  for (int k = 0; k < K; ++k) {
    if (gfo[k] <= 0) { p[k] = 0.0; continue; }
    double X, Y;
    if (nh) { X = nnh + hotot - hii[k] + (Nd - nnh - hotot) - hix[k]; Y = Nd; }
    else    { X = 1.0 + hotot - hii[k] + (Nd - hotot) - hix[k];       Y = Nd + 1.0; }
    p[k] = 1.0 - std::sqrt(std::max(0.0, X / Y));
  }
  for (int iter = 0; iter < 5000; ++iter) {
    double rdi = 0.0; int re = 0;
    for (int k = 0; k < K; ++k) {
      if (gfo[k] <= 0) { np[k] = 0.0; continue; }
      const double pa = p[k], denom = pa + 2.0 * rd;
      if (denom <= 0) { np[k] = 0.0; continue; }
      const double pn = (pa + rd) / denom * (hii[k] / Nd) + hix[k] / (2.0 * Nd);
      rdi += rd / denom * (hii[k] / Nd);
      np[k] = pn;
      if (std::fabs(pn - pa) > 1e-6) ++re;
    }
    const double rd_new = nh ? rdi + (2.0 * nnh) / (2.0 * Nd) : rdi;
    if (std::fabs(rd_new - rd) > 1e-6) ++re;
    p = np; rd = std::max(0.0, rd_new);
    if (re == 0) break;
  }
  out.rd = rd; out.al = al; out.gfo = gfo; out.pf = p; out.hii = hii; out.hix = hix;
}

struct W3 { double s2P, s2I, s2G; };

inline W3 weir(int r, const double* ni, const double* nA, const double* AA) {
  double Nt = 0, N2 = 0;
  for (int k = 0; k < r; ++k) { Nt += ni[k]; N2 += ni[k] * ni[k]; }
  if (Nt == 0 || r < 2) return {0, 0, 0};
  const double nc = (Nt - N2 / Nt) / (r - 1);
  if (nc <= 0 || Nt - r <= 0) return {0, 0, 0};
  double snA = 0, s2A = 0, sAA = 0;
  for (int k = 0; k < r; ++k) {
    snA += nA[k]; sAA += AA[k];
    if (ni[k] > 0) s2A += nA[k] * nA[k] / (2 * ni[k]);
  }
  const double MSG = (0.5 * snA - sAA) / Nt;
  const double MSI = (0.5 * snA + sAA - s2A) / (Nt - r);
  const double MSP = (s2A - 0.5 * snA * snA / Nt) / (r - 1);
  return {(MSP - MSI) / (2 * nc), 0.5 * (MSI - MSG), MSG};
}

inline double dcs(double csprod) {
  if (csprod > 1.0) return NA_REAL;
  return (2.0 / M_PI) * std::sqrt(2.0 * (1.0 - csprod));
}

struct Out {                         // statistics of one data set / replicate
  double g_raw, g_ena;
  std::vector<double> loc_raw, loc_ena;   // L
  std::vector<double> pr, pe, dr, di;     // pairs
};

// cells[l][p]
void stats_from_cells(const std::vector<std::vector<Cell>>& cells, Out& o) {
  const int L = (int)cells.size(), P = (int)cells[0].size();
  const int NP = P * (P - 1) / 2;
  o.loc_raw.assign(L, NA_REAL); o.loc_ena.assign(L, NA_REAL);
  std::vector<double> s1r(NP, 0), s3r(NP, 0), s1e(NP, 0), s3e(NP, 0);
  std::vector<double> dsr(NP, 0), dsi(NP, 0);
  std::vector<int> dnr(NP, L), dni(NP, L);
  double S1 = 0, S3 = 0, S1c = 0, S3c = 0;

  for (int l = 0; l < L; ++l) {
    const std::vector<Cell>& cl = cells[l];
    std::vector<int> U;
    for (int p = 0; p < P; ++p) U.insert(U.end(), cl[p].al.begin(), cl[p].al.end());
    std::sort(U.begin(), U.end()); U.erase(std::unique(U.begin(), U.end()), U.end());
    const int K = (int)U.size();
    // dense arrays [p][k]
    std::vector<std::vector<double>> gfo(P, std::vector<double>(K, 0.0)), pf = gfo, hii = gfo;
    for (int p = 0; p < P; ++p)
      for (size_t j = 0; j < cl[p].al.size(); ++j) {
        const int k = (int)(std::lower_bound(U.begin(), U.end(), cl[p].al[j]) - U.begin());
        gfo[p][k] = cl[p].gfo[j]; pf[p][k] = cl[p].pf[j]; hii[p][k] = cl[p].hii[j];
      }
    std::vector<double> nir(P), nic(P);
    for (int p = 0; p < P; ++p) {
      nir[p] = std::max(0, cl[p].efpop - cl[p].n_absent - cl[p].n_null_homo);
      nic[p] = std::max(0, cl[p].efpop - cl[p].n_absent);
    }
    // ---- global ----
    int rr = 0, rc = 0; double Nr = 0, N2r = 0, Nc = 0, N2c = 0;
    for (int p = 0; p < P; ++p) {
      rr += nir[p] > 0; rc += nic[p] > 0;
      Nr += nir[p]; N2r += nir[p] * nir[p]; Nc += nic[p]; N2c += nic[p] * nic[p];
    }
    const double ncr = (Nr > 0 && rr > 1) ? (Nr - N2r / Nr) / (rr - 1) : 0.0;
    const double ncc = (Nc > 0 && rc > 1) ? (Nc - N2c / Nc) / (rc - 1) : 0.0;
    double s1l = 0, s3l = 0, s1lc = 0, s3lc = 0;
    std::vector<double> ni(P), nA(P), AA(P), AAc(P);
    for (int k = 0; k < K; ++k) {
      for (int p = 0; p < P; ++p) { ni[p] = nir[p]; nA[p] = gfo[p][k] * 2 * nir[p]; AA[p] = hii[p][k]; }
      W3 w = weir(P, ni.data(), nA.data(), AA.data());
      s1l += w.s2P; s3l += w.s2P + w.s2I + w.s2G;
      for (int p = 0; p < P; ++p) {
        const double d = pf[p][k] + 2.0 * cl[p].rd;
        ni[p] = nic[p]; nA[p] = pf[p][k] * 2 * nic[p]; AA[p] = hii[p][k];
        AAc[p] = (AA[p] > 0 && d > 0) ? AA[p] * (pf[p][k] / d) : 0.0;
      }
      W3 wc = weir(P, ni.data(), nA.data(), AAc.data());
      s1lc += wc.s2P; s3lc += wc.s2P + wc.s2I + wc.s2G;
    }
    o.loc_raw[l] = s3l != 0 ? s1l / s3l : NA_REAL;
    o.loc_ena[l] = s3lc != 0 ? s1lc / s3lc : NA_REAL;
    S1 += s1l * ncr; S3 += s3l * ncr; S1c += s1lc * ncc; S3c += s3lc * ncc;

    // ---- pairs ----
    int pi = 0;
    for (int i = 0; i < P - 1; ++i) for (int j = i + 1; j < P; ++j, ++pi) {
      const double nri = nir[i], nrj = nir[j], nci = nic[i], ncj = nic[j];
      const double Nrp = nri + nrj, N2rp = nri * nri + nrj * nrj;
      const double Ncp = nci + ncj, N2cp = nci * nci + ncj * ncj;
      const double ncr_p = (Nrp > 0 && nri > 0 && nrj > 0) ? (Nrp - N2rp / Nrp) : 0.0;
      const double ncc_p = (Ncp > 0 && nci > 0 && ncj > 0) ? (Ncp - N2cp / Ncp) : 0.0;
      for (int k = 0; k < K; ++k) {
        if (ncr_p > 0) {
          const double n2[2] = {nri, nrj};
          const double a2[2] = {gfo[i][k] * 2 * nri, gfo[j][k] * 2 * nrj};
          const double h2[2] = {hii[i][k], hii[j][k]};
          W3 w = weir(2, n2, a2, h2);
          s1r[pi] += w.s2P * ncr_p; s3r[pi] += (w.s2P + w.s2I + w.s2G) * ncr_p;
        }
        if (ncc_p > 0) {
          const double di = pf[i][k] + 2.0 * cl[i].rd, dj = pf[j][k] + 2.0 * cl[j].rd;
          const double n2[2] = {nci, ncj};
          const double a2[2] = {pf[i][k] * 2 * nci, pf[j][k] * 2 * ncj};
          const double h2[2] = {(hii[i][k] > 0 && di > 0) ? hii[i][k] * (pf[i][k] / di) : 0.0,
                                (hii[j][k] > 0 && dj > 0) ? hii[j][k] * (pf[j][k] / dj) : 0.0};
          W3 w = weir(2, n2, a2, h2);
          s1e[pi] += w.s2P * ncc_p; s3e[pi] += (w.s2P + w.s2I + w.s2G) * ncc_p;
        }
      }
      // DCSE
      if (nri > 0 && nrj > 0 && !cl[i].gfo.empty() && !cl[j].gfo.empty()) {
        double cs = 0; for (int k = 0; k < K; ++k) if (gfo[i][k] > 0 && gfo[j][k] > 0) cs += std::sqrt(gfo[i][k] * gfo[j][k]);
        const double d = dcs(cs);
        if (!ISNAN(d)) dsr[pi] += d; else dnr[pi]--;
      } else dnr[pi]--;
      if (nci > 0 && ncj > 0) {
        double cs = 0; for (int k = 0; k < K; ++k) if (pf[i][k] > 0 && pf[j][k] > 0) cs += std::sqrt(pf[i][k] * pf[j][k]);
        if (cl[i].rd > 0 && cl[j].rd > 0) cs += std::sqrt(cl[i].rd * cl[j].rd);
        const double d = dcs(cs);
        if (!ISNAN(d)) dsi[pi] += d; else dni[pi]--;
      } else dni[pi]--;
    }
  }
  o.g_raw = S3 > 0 ? S1 / S3 : NA_REAL;
  o.g_ena = S3c > 0 ? S1c / S3c : NA_REAL;
  o.pr.assign(NP, NA_REAL); o.pe.assign(NP, NA_REAL); o.dr.assign(NP, NA_REAL); o.di.assign(NP, NA_REAL);
  for (int q = 0; q < NP; ++q) {
    if (s3r[q] > 0) o.pr[q] = s1r[q] / s3r[q];
    if (s3e[q] > 0) o.pe[q] = s1e[q] / s3e[q];
    if (dnr[q] > 0) o.dr[q] = dsr[q] / dnr[q];
    if (dni[q] > 0) o.di[q] = dsi[q] / dni[q];
  }
}

struct Data {
  int L = 0, P = 0, base = 1000;
  std::vector<std::vector<int>> g;     // g[p] column-major n_p x L
  std::vector<int> n;                  // n_p
  std::vector<char> nullhomo;          // per locus
};

void run_one(const Data& D, const std::vector<std::vector<int>>& idx, Out& o) {
  std::vector<std::vector<Cell>> cells(D.L, std::vector<Cell>(D.P));
  std::vector<int> buf;
  for (int p = 0; p < D.P; ++p) {
    const int n = D.n[p];
    for (int l = 0; l < D.L; ++l) {
      buf.resize(n);
      const int* col = D.g[p].data() + (size_t)l * n;
      for (int i = 0; i < n; ++i) buf[i] = col[idx[p][i]];
      em_cell(buf, D.base, D.nullhomo[l] != 0, cells[l][p]);
    }
  }
  stats_from_cells(cells, o);
}

Data read_data(List G, int base, IntegerVector nullhomo) {
  Data D; D.P = G.size(); D.base = base;
  for (int p = 0; p < D.P; ++p) {
    IntegerMatrix m = G[p];
    D.n.push_back(m.nrow()); D.L = m.ncol();
    D.g.emplace_back(m.begin(), m.end());
  }
  for (int l = 0; l < D.L; ++l) D.nullhomo.push_back(nullhomo[l] != 0);
  return D;
}

List pack(const Out& o) {
  return List::create(_["global_raw"] = o.g_raw, _["global_ena"] = o.g_ena,
                      _["loc_raw"] = wrap(o.loc_raw), _["loc_ena"] = wrap(o.loc_ena),
                      _["pair_raw"] = wrap(o.pr), _["pair_ena"] = wrap(o.pe),
                      _["dcse_raw"] = wrap(o.dr), _["dcse_ina"] = wrap(o.di));
}

}  // namespace

// Observed statistics (no resampling). G: list over sub-samples of integer
// matrices (individuals x loci, base-coded genotypes, <= 0 = missing).
// [[Rcpp::export]]
List na_stats_cpp(List G, int base, IntegerVector nullhomo) {
  Data D = read_data(G, base, nullhomo);
  std::vector<std::vector<int>> idx(D.P);
  for (int p = 0; p < D.P; ++p) { idx[p].resize(D.n[p]); for (int i = 0; i < D.n[p]; ++i) idx[p][i] = i; }
  Out o; run_one(D, idx, o);
  return pack(o);
}

// Bootstrap over individuals: individuals are resampled with replacement
// within each sub-sample (the same individuals at every locus); the EM and
// all statistics are recomputed on every replicate. Replicate b uses its own
// generator seeded from (seed, b), so results do not depend on the number of
// threads.
// [[Rcpp::export]]
List na_boot_indiv_cpp(List G, int base, IntegerVector nullhomo, int B, double seed, int n_threads) {
  Data D = read_data(G, base, nullhomo);
  const int NP = D.P * (D.P - 1) / 2;
  NumericVector gr(B), ge(B);
  NumericMatrix lr(B, D.L), le(B, D.L), pr(B, NP), pe(B, NP), dr(B, NP), di(B, NP);
  std::vector<Out> outs(B);
#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
#endif
#pragma omp parallel for schedule(dynamic, 4)
  for (int b = 0; b < B; ++b) {
    std::mt19937_64 rng((uint64_t)seed * 1000003ULL + (uint64_t)b * 7919ULL + 12345ULL);
    std::vector<std::vector<int>> idx(D.P);
    for (int p = 0; p < D.P; ++p) {
      idx[p].resize(D.n[p]);
      if (D.n[p] > 0) {
        std::uniform_int_distribution<int> U(0, D.n[p] - 1);
        for (int i = 0; i < D.n[p]; ++i) idx[p][i] = U(rng);
      }
    }
    run_one(D, idx, outs[b]);
  }
  for (int b = 0; b < B; ++b) {
    gr[b] = outs[b].g_raw; ge[b] = outs[b].g_ena;
    for (int l = 0; l < D.L; ++l) { lr(b, l) = outs[b].loc_raw[l]; le(b, l) = outs[b].loc_ena[l]; }
    for (int q = 0; q < NP; ++q) {
      pr(b, q) = outs[b].pr[q]; pe(b, q) = outs[b].pe[q];
      dr(b, q) = outs[b].dr[q]; di(b, q) = outs[b].di[q];
    }
  }
  return List::create(_["glob_raw"] = gr, _["glob_ena"] = ge, _["loc_raw"] = lr, _["loc_ena"] = le,
                      _["pair_raw"] = pr, _["pair_ena"] = pe, _["dcse_raw"] = dr, _["dcse_ina"] = di);
}
