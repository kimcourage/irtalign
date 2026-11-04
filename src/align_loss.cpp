// Alignment simplicity function, gradient and second derivatives.
//
// Standard form: every item has one slope a and K intercepts c_k. For group g
// with latent mean alpha_g and log variance l_g (s_g = exp(-l_g / 2)):
//   aligned slope      a1 = a0 * s_g
//   aligned intercept  c1 = c0 - alpha_g * a0 * s_g
// The simplicity function sums w_gh * wt * f(x_g - x_h) over item parameters
// and group pairs, with the component loss f(x) = (x^2 + eps)^(power / 2).

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// The default power 0.5 avoids std::pow: with r = u^(1/4),
// u^(-3/4) = r / u and u^(-7/4) = r / u^2.
static inline double clf0(double x, double eps, double p) {
  const double u = x * x + eps;
  if (p == 0.5) return std::sqrt(std::sqrt(u));
  return std::pow(u, 0.5 * p);
}
static inline double clf1(double x, double eps, double p) {
  const double u = x * x + eps;
  if (p == 0.5) return 0.5 * x * std::sqrt(std::sqrt(u)) / u;
  return p * x * std::pow(u, 0.5 * p - 1.0);
}
// value and first derivative together (one root evaluation)
static inline void clf01(double x, double eps, double p, double& f0, double& f1) {
  const double u = x * x + eps;
  const double v = (p == 0.5) ? std::sqrt(std::sqrt(u)) : std::pow(u, 0.5 * p);
  f0 = v;
  f1 = p * x * v / u;
}
static inline double clf2(double x, double eps, double p) {
  const double u = x * x + eps;
  if (p == 0.5) return 0.5 * std::sqrt(std::sqrt(u)) / (u * u) * (eps - 0.5 * x * x);
  return p * std::pow(u, 0.5 * p - 2.0) * (eps + (p - 1.0) * x * x);
}

// [[Rcpp::export]]
List align_loss_cpp(const NumericMatrix& a0, const NumericMatrix& c0,
                    const IntegerVector& citem, const NumericMatrix& W,
                    const NumericVector& wa, const NumericVector& wc,
                    const NumericVector& alpha, const NumericVector& lpsi,
                    double eps, double power, bool gradient) {
  const int G = a0.nrow(), I = a0.ncol(), P = c0.ncol();
  std::vector<double> s(G);
  for (int g = 0; g < G; ++g) s[g] = std::exp(-0.5 * lpsi[g]);
  NumericMatrix a1(G, I);
  for (int i = 0; i < I; ++i)
    for (int g = 0; g < G; ++g) a1(g, i) = a0(g, i) * s[g];

  double F = 0.0;
  NumericVector ga(G), gl(G);

  for (int i = 0; i < I; ++i) {
    if (wa[i] == 0.0) continue;
    for (int g = 0; g < G - 1; ++g) {
      if (ISNAN(a1(g, i))) continue;
      for (int h = g + 1; h < G; ++h) {
        if (ISNAN(a1(h, i))) continue;
        const double w = W(g, h) * wa[i];
        const double d = a1(g, i) - a1(h, i);
        if (gradient) {
          double v0, v1;
          clf01(d, eps, power, v0, v1);
          F += w * v0;
          const double fp = w * v1;
          gl[g] += fp * (-0.5 * a1(g, i));
          gl[h] += fp * (0.5 * a1(h, i));
        } else {
          F += w * clf0(d, eps, power);
        }
      }
    }
  }

  for (int p = 0; p < P; ++p) {
    if (wc[p] == 0.0) continue;
    const int i = citem[p];
    for (int g = 0; g < G - 1; ++g) {
      if (ISNAN(c0(g, p)) || ISNAN(a1(g, i))) continue;
      const double xg = c0(g, p) - alpha[g] * a1(g, i);
      for (int h = g + 1; h < G; ++h) {
        if (ISNAN(c0(h, p)) || ISNAN(a1(h, i))) continue;
        const double xh = c0(h, p) - alpha[h] * a1(h, i);
        const double w = W(g, h) * wc[p];
        const double d = xg - xh;
        if (gradient) {
          double v0, v1;
          clf01(d, eps, power, v0, v1);
          F += w * v0;
          const double fp = w * v1;
          ga[g] += fp * (-a1(g, i));
          ga[h] += fp * (a1(h, i));
          gl[g] += fp * (0.5 * alpha[g] * a1(g, i));
          gl[h] += fp * (-0.5 * alpha[h] * a1(h, i));
        } else {
          F += w * clf0(d, eps, power);
        }
      }
    }
  }

  if (!gradient) return List::create(Named("value") = F);
  NumericVector grad(2 * G);
  for (int g = 0; g < G; ++g) {
    grad[g] = ga[g];
    grad[G + g] = gl[g];
  }
  return List::create(Named("value") = F, Named("gradient") = grad);
}

// Contribution of each group to each column of aligned parameters x (G x Q):
// every pair term is split equally between its two groups, so the matrix sums
// to the simplicity function.
// [[Rcpp::export]]
NumericMatrix align_contrib_cpp(const NumericMatrix& x, const NumericMatrix& W,
                                const NumericVector& wt, double eps,
                                double power) {
  const int G = x.nrow(), Q = x.ncol();
  NumericMatrix out(G, Q);  // contribution of each group (half of each pair)
  for (int q = 0; q < Q; ++q) {
    if (wt[q] == 0.0) continue;
    for (int g = 0; g < G - 1; ++g) {
      if (ISNAN(x(g, q))) continue;
      for (int h = g + 1; h < G; ++h) {
        if (ISNAN(x(h, q))) continue;
        const double v = W(g, h) * wt[q] * clf0(x(g, q) - x(h, q), eps, power);
        out(g, q) += 0.5 * v;
        out(h, q) += 0.5 * v;
      }
    }
  }
  return out;
}

// Accumulate one term w * f(D), D = x_g - x_h, into H22 (m2 x m2) and
// H21 (m2 x m1). Local variables 0..3 belong to group g, 4..7 to group h:
// alpha, log variance, slope a0, intercept c0. A negative index marks an
// absent variable.
static void add_term(double w, double D, double eps, double power,
                     const int* idx, const bool* is_m2, const double* dD,
                     const double (*d2D)[8], NumericMatrix& H22,
                     NumericMatrix& H21) {
  const double f1 = w * clf1(D, eps, power);
  const double f2 = w * clf2(D, eps, power);
  for (int r = 0; r < 8; ++r) {
    if (idx[r] < 0 || !is_m2[r]) continue;
    for (int c = 0; c < 8; ++c) {
      if (idx[c] < 0) continue;
      const double v = f2 * dD[r] * dD[c] + f1 * d2D[r][c];
      if (v == 0.0) continue;
      if (is_m2[c]) H22(idx[r], idx[c]) += v;
      else H21(idx[r], idx[c]) += v;
    }
  }
}

// Hessian with respect to (alpha, log variance) only, used by the optimizer.
// [[Rcpp::export]]
NumericMatrix align_hess22_cpp(const NumericMatrix& a0, const NumericMatrix& c0,
                               const IntegerVector& citem, const NumericMatrix& W,
                               const NumericVector& wa, const NumericVector& wc,
                               const NumericVector& alpha, const NumericVector& lpsi,
                               double eps, double power) {
  const int G = a0.nrow(), I = a0.ncol(), P = c0.ncol();
  std::vector<double> s(G);
  for (int g = 0; g < G; ++g) s[g] = std::exp(-0.5 * lpsi[g]);
  NumericMatrix a1(G, I);
  for (int i = 0; i < I; ++i)
    for (int g = 0; g < G; ++g) a1(g, i) = a0(g, i) * s[g];
  NumericMatrix H(2 * G, 2 * G);

  for (int i = 0; i < I; ++i) {
    if (wa[i] == 0.0) continue;
    for (int g = 0; g < G - 1; ++g) {
      const double xg = a1(g, i);
      if (ISNAN(xg)) continue;
      for (int h = g + 1; h < G; ++h) {
        const double xh = a1(h, i);
        if (ISNAN(xh)) continue;
        const double w = W(g, h) * wa[i];
        const double D = xg - xh;
        const double f1 = w * clf1(D, eps, power), f2 = w * clf2(D, eps, power);
        const int lg = G + g, lh = G + h;
        H(lg, lg) += f2 * 0.25 * xg * xg + f1 * 0.25 * xg;
        H(lh, lh) += f2 * 0.25 * xh * xh - f1 * 0.25 * xh;
        const double cross = -f2 * 0.25 * xg * xh;
        H(lg, lh) += cross;
        H(lh, lg) += cross;
      }
    }
  }

  for (int p = 0; p < P; ++p) {
    if (wc[p] == 0.0) continue;
    const int i = citem[p];
    for (int g = 0; g < G - 1; ++g) {
      if (ISNAN(c0(g, p)) || ISNAN(a1(g, i))) continue;
      const double ag = a1(g, i), xg = c0(g, p) - alpha[g] * ag;
      for (int h = g + 1; h < G; ++h) {
        if (ISNAN(c0(h, p)) || ISNAN(a1(h, i))) continue;
        const double ah = a1(h, i), xh = c0(h, p) - alpha[h] * ah;
        const double w = W(g, h) * wc[p];
        const double D = xg - xh;
        const double f1 = w * clf1(D, eps, power), f2 = w * clf2(D, eps, power);
        const int idx[4] = {g, G + g, h, G + h};
        const double d[4] = {-ag, 0.5 * alpha[g] * ag, ah, -0.5 * alpha[h] * ah};
        for (int r = 0; r < 4; ++r)
          for (int c = 0; c < 4; ++c) H(idx[r], idx[c]) += f2 * d[r] * d[c];
        // second derivatives of D
        H(g, G + g) += f1 * 0.5 * ag;
        H(G + g, g) += f1 * 0.5 * ag;
        H(G + g, G + g) += -f1 * 0.25 * alpha[g] * ag;
        H(h, G + h) += -f1 * 0.5 * ah;
        H(G + h, h) += -f1 * 0.5 * ah;
        H(G + h, G + h) += f1 * 0.25 * alpha[h] * ah;
      }
    }
  }
  return H;
}

// [[Rcpp::export]]
List align_hess_cpp(const NumericMatrix& a0, const NumericMatrix& c0,
                    const IntegerVector& citem, const NumericMatrix& W,
                    const NumericVector& wa, const NumericVector& wc,
                    const NumericVector& alpha, const NumericVector& lpsi,
                    double eps, double power, const IntegerMatrix& ia,
                    const IntegerMatrix& ic, int M1) {
  const int G = a0.nrow(), I = a0.ncol(), P = c0.ncol();
  std::vector<double> s(G);
  for (int g = 0; g < G; ++g) s[g] = std::exp(-0.5 * lpsi[g]);
  NumericMatrix H22(2 * G, 2 * G), H21(2 * G, M1);

  int idx[8];
  bool is_m2[8] = {true, true, false, false, true, true, false, false};
  double dD[8];
  double d2D[8][8];

  // slopes
  for (int i = 0; i < I; ++i) {
    if (wa[i] == 0.0) continue;
    for (int g = 0; g < G - 1; ++g) {
      if (ISNAN(a0(g, i))) continue;
      for (int h = g + 1; h < G; ++h) {
        if (ISNAN(a0(h, i))) continue;
        const double w = W(g, h) * wa[i];
        const double xg = a0(g, i) * s[g], xh = a0(h, i) * s[h];
        const double D = xg - xh;
        for (int r = 0; r < 8; ++r) {
          dD[r] = 0.0;
          for (int c = 0; c < 8; ++c) d2D[r][c] = 0.0;
        }
        idx[0] = -1; idx[1] = G + g; idx[2] = ia(g, i); idx[3] = -1;
        idx[4] = -1; idx[5] = G + h; idx[6] = ia(h, i); idx[7] = -1;
        // group g side: dx/dl = -x/2, dx/da0 = s
        dD[1] = -0.5 * xg; dD[2] = s[g];
        d2D[1][1] = 0.25 * xg; d2D[1][2] = d2D[2][1] = -0.5 * s[g];
        // group h side enters with a minus sign
        dD[5] = 0.5 * xh; dD[6] = -s[h];
        d2D[5][5] = -0.25 * xh; d2D[5][6] = d2D[6][5] = 0.5 * s[h];
        add_term(w, D, eps, power, idx, is_m2, dD, d2D, H22, H21);
      }
    }
  }

  // intercepts
  for (int p = 0; p < P; ++p) {
    if (wc[p] == 0.0) continue;
    const int i = citem[p];
    for (int g = 0; g < G - 1; ++g) {
      if (ISNAN(c0(g, p)) || ISNAN(a0(g, i))) continue;
      for (int h = g + 1; h < G; ++h) {
        if (ISNAN(c0(h, p)) || ISNAN(a0(h, i))) continue;
        const double w = W(g, h) * wc[p];
        const double ag = a0(g, i) * s[g], ah = a0(h, i) * s[h];
        const double xg = c0(g, p) - alpha[g] * ag;
        const double xh = c0(h, p) - alpha[h] * ah;
        const double D = xg - xh;
        for (int r = 0; r < 8; ++r) {
          dD[r] = 0.0;
          for (int c = 0; c < 8; ++c) d2D[r][c] = 0.0;
        }
        idx[0] = g; idx[1] = G + g; idx[2] = ia(g, i); idx[3] = ic(g, p);
        idx[4] = h; idx[5] = G + h; idx[6] = ia(h, i); idx[7] = ic(h, p);
        // group g side
        dD[0] = -ag;
        dD[1] = 0.5 * alpha[g] * ag;
        dD[2] = -alpha[g] * s[g];
        dD[3] = 1.0;
        d2D[0][1] = d2D[1][0] = 0.5 * ag;
        d2D[0][2] = d2D[2][0] = -s[g];
        d2D[1][1] = -0.25 * alpha[g] * ag;
        d2D[1][2] = d2D[2][1] = 0.5 * alpha[g] * s[g];
        // group h side (negated)
        dD[4] = ah;
        dD[5] = -0.5 * alpha[h] * ah;
        dD[6] = alpha[h] * s[h];
        dD[7] = -1.0;
        d2D[4][5] = d2D[5][4] = -0.5 * ah;
        d2D[4][6] = d2D[6][4] = s[h];
        d2D[5][5] = 0.25 * alpha[h] * ah;
        d2D[5][6] = d2D[6][5] = -0.5 * alpha[h] * s[h];
        add_term(w, D, eps, power, idx, is_m2, dD, d2D, H22, H21);
      }
    }
  }
  return List::create(Named("H22") = H22, Named("H21") = H21);
}
