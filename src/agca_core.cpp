// Author: Alberto Quaini

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;
using arma::rowvec;
using arma::vec;

namespace {

constexpr double kEps = 1e-12;

void check_finite_matrix(const mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("%s must contain only finite values.", name);
  }
}

vec unit_vector(vec x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("%s must contain only finite values.", name);
  }
  const double nrm = arma::norm(x, 2);
  if (!std::isfinite(nrm) || nrm <= kEps) {
    Rcpp::stop("%s must have positive Euclidean norm.", name);
  }
  return x / nrm;
}

mat normalize_rows_impl(mat x) {
  check_finite_matrix(x, "x");
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const double nrm = arma::norm(x.row(i), 2);
    if (!std::isfinite(nrm) || nrm <= kEps) {
      Rcpp::stop("All rows must have positive finite Euclidean norm.");
    }
    x.row(i) /= nrm;
  }
  return x;
}

mat orient_columns(mat x) {
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    arma::uword idx = arma::index_max(arma::abs(x.col(j)));
    if (x(idx, j) < 0.0) {
      x.col(j) *= -1.0;
    }
  }
  return x;
}

arma::vec row_norms_impl(const arma::mat& x) {
  check_finite_matrix(x, "x");
  arma::vec out(x.n_rows);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    out(i) = arma::norm(x.row(i), 2);
  }
  return out;
}

arma::uvec order_radius_decreasing(const arma::vec& radius) {
  std::vector<arma::uword> ord(radius.n_elem);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](arma::uword a, arma::uword b) {
    if (radius(a) == radius(b)) {
      return a < b;
    }
    return radius(a) > radius(b);
  });
  arma::uvec out(ord.size());
  for (arma::uword i = 0; i < ord.size(); ++i) {
    out(i) = ord[i];
  }
  return out;
}

arma::mat subset_rows(const arma::mat& x, const arma::uvec& index) {
  arma::mat out(index.n_elem, x.n_cols);
  for (arma::uword i = 0; i < index.n_elem; ++i) {
    out.row(i) = x.row(index(i));
  }
  return out;
}

Rcpp::IntegerVector to_one_based(const arma::uvec& index) {
  Rcpp::IntegerVector out(index.n_elem);
  for (arma::uword i = 0; i < index.n_elem; ++i) {
    out[i] = static_cast<int>(index(i) + 1);
  }
  return out;
}

Rcpp::List tail_directions_impl(const arma::mat& x, int k, double threshold, bool use_k) {
  check_finite_matrix(x, "x");
  if (x.n_rows < 1 || x.n_cols < 2) {
    Rcpp::stop("x must have at least one row and two columns.");
  }

  arma::vec radius = row_norms_impl(x);
  arma::uvec index;
  double selected_threshold = threshold;

  if (use_k) {
    if (k < 1 || k > static_cast<int>(x.n_rows)) {
      Rcpp::stop("k must be an integer between 1 and nrow(x).");
    }
    arma::uvec ord = order_radius_decreasing(radius);
    index = ord.head(k);
    selected_threshold = arma::min(radius.elem(index));
  } else {
    if (!std::isfinite(threshold)) {
      Rcpp::stop("threshold must be a finite number.");
    }
    index = arma::find(radius > threshold);
  }

  if (index.n_elem == 0) {
    Rcpp::stop("No observations exceed the threshold.");
  }

  arma::mat directions = normalize_rows_impl(subset_rows(x, index));
  arma::vec selected_radius = radius.elem(index);
  return Rcpp::List::create(
    Rcpp::Named("directions") = directions,
    Rcpp::Named("radius") = selected_radius,
    Rcpp::Named("index") = to_one_based(index),
    Rcpp::Named("threshold") = selected_threshold,
    Rcpp::Named("all_radius") = radius
  );
}

arma::vec principal_anchor_impl(arma::mat g, bool normalize) {
  if (normalize) {
    g = normalize_rows_impl(std::move(g));
  } else {
    check_finite_matrix(g, "g");
  }
  arma::mat moment = (g.t() * g) / static_cast<double>(g.n_rows);
  arma::vec eigval;
  arma::mat eigvec;
  if (!arma::eig_sym(eigval, eigvec, moment)) {
    Rcpp::stop("Eigenvalue decomposition failed.");
  }
  arma::uword idx = eigval.index_max();
  arma::vec mu = eigvec.col(idx);
  if (arma::accu(mu) < 0.0) {
    mu *= -1.0;
  }
  return unit_vector(std::move(mu), "principal anchor");
}

arma::vec frechet_anchor_impl(arma::mat g, bool normalize, int max_iter, double tol) {
  if (normalize) {
    g = normalize_rows_impl(std::move(g));
  } else {
    check_finite_matrix(g, "g");
  }
  if (max_iter < 1) {
    Rcpp::stop("max_iter must be a positive integer.");
  }
  if (!std::isfinite(tol) || tol <= 0.0) {
    Rcpp::stop("tol must be a positive finite number.");
  }

  arma::vec mu = unit_vector(arma::mean(g, 0).t(), "initial Frechet anchor");
  const double theta_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::vec inner = g * mu;
    inner = arma::clamp(inner, -1.0, 1.0);
    arma::vec theta = arma::acos(inner);
    arma::vec scale(theta.n_elem, arma::fill::ones);
    arma::uvec nz = arma::find(theta > theta_eps);
    scale.elem(nz) = theta.elem(nz) / arma::sin(theta.elem(nz));
    arma::mat logs = g - inner * mu.t();
    logs.each_col() %= scale;
    arma::vec update = arma::mean(logs, 0).t();
    const double update_norm = arma::norm(update, 2);
    if (update_norm < tol) {
      break;
    }
    mu = std::cos(update_norm) * mu + std::sin(update_norm) * update / update_norm;
    mu = unit_vector(std::move(mu), "Frechet anchor");
  }
  if (arma::accu(mu) < 0.0) {
    mu *= -1.0;
  }
  return mu;
}

arma::vec variation_explained_from_eigenvalues(const arma::vec& eigenvalues) {
  const double total = arma::accu(eigenvalues);
  arma::vec out(eigenvalues.n_elem);
  if (!std::isfinite(total) || total <= 0.0) {
    out.fill(NA_REAL);
    return out;
  }
  out = arma::cumsum(eigenvalues) / total;
  return out;
}

arma::vec resolve_anchor_impl(const arma::mat& g, const arma::vec& fixed_mu,
                              bool fixed_anchor, const std::string& anchor_type,
                              const Rcpp::NumericVector& anchor_vector) {
  if (fixed_anchor) {
    return unit_vector(fixed_mu, "mu");
  }
  if (anchor_type == "canonical") {
    return arma::vec(g.n_cols, arma::fill::ones) / std::sqrt(static_cast<double>(g.n_cols));
  }
  if (anchor_type == "principal") {
    return principal_anchor_impl(g, false);
  }
  if (anchor_type == "frechet") {
    return frechet_anchor_impl(g, false, 100, 1e-10);
  }
  if (anchor_type == "numeric") {
    if (anchor_vector.size() != static_cast<int>(g.n_cols)) {
      Rcpp::stop("Numeric anchor must have length equal to ncol(g).");
    }
    return unit_vector(Rcpp::as<arma::vec>(anchor_vector), "anchor");
  }
  Rcpp::stop("Unsupported bootstrap anchor.");
}

arma::vec bounded_positive_exposure_impl(const arma::mat& g, const arma::mat& weights,
                                         double power, double cap) {
  arma::mat exposure = g * weights.t();
  exposure.transform([&](double z) {
    const double positive = z > 0.0 ? z : 0.0;
    double value = std::pow(positive, power);
    if (std::isfinite(cap) && value > cap) {
      value = cap;
    }
    return value;
  });
  return arma::mean(exposure, 0).t();
}

}  // namespace

// [[Rcpp::export]]
arma::vec agca_row_norms_cpp(const arma::mat& x) {
  return row_norms_impl(x);
}

// [[Rcpp::export]]
arma::mat agca_normalize_rows_cpp(arma::mat x) {
  return normalize_rows_impl(std::move(x));
}

// [[Rcpp::export]]
Rcpp::List agca_departures_cpp(arma::mat g, arma::vec mu, bool normalize = true) {
  if (normalize) {
    g = normalize_rows_impl(std::move(g));
  } else {
    check_finite_matrix(g, "g");
  }

  mu = unit_vector(std::move(mu), "mu");
  if (g.n_cols != mu.n_elem) {
    Rcpp::stop("g and mu have incompatible dimensions.");
  }

  vec a = g * mu;
  if (arma::any(a <= 0.0)) {
    Rcpp::stop("All directions must lie in the open hemisphere centered at mu.");
  }

  mat u = g - a * mu.t();
  return Rcpp::List::create(
    Rcpp::Named("g") = g,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("anchor_coordinate") = a,
    Rcpp::Named("u") = u
  );
}

// [[Rcpp::export]]
Rcpp::List agca_core_cpp(arma::mat g, arma::vec mu, bool normalize = true) {
  Rcpp::List dep = agca_departures_cpp(std::move(g), std::move(mu), normalize);
  mat g_norm = dep["g"];
  vec mu_unit = dep["mu"];
  vec anchor_coordinate = dep["anchor_coordinate"];
  mat u = dep["u"];

  const arma::uword n = g_norm.n_rows;
  const arma::uword d = g_norm.n_cols;
  if (n < 1) {
    Rcpp::stop("g must contain at least one row.");
  }
  if (d < 2) {
    Rcpp::stop("The ambient dimension must be at least 2.");
  }

  mat basis = arma::null(mu_unit.t());
  if (basis.n_cols != d - 1) {
    Rcpp::stop("Could not construct a full tangent basis at the anchor.");
  }

  mat tangent_scores = u * basis;
  mat sigma_tangent = (tangent_scores.t() * tangent_scores) / static_cast<double>(n);
  mat eigvec;
  vec eigval;
  if (!arma::eig_sym(eigval, eigvec, sigma_tangent)) {
    Rcpp::stop("Eigenvalue decomposition failed.");
  }

  arma::uvec ord = arma::sort_index(eigval, "descend");
  eigval = eigval.elem(ord);
  eigvec = eigvec.cols(ord);
  eigval.elem(arma::find(arma::abs(eigval) < 100.0 * std::numeric_limits<double>::epsilon())).zeros();
  eigval.elem(arma::find(eigval < 0.0)).zeros();

  mat loadings = orient_columns(basis * eigvec);
  mat scores = u * loadings;
  mat sigma = (u.t() * u) / static_cast<double>(n);

  return Rcpp::List::create(
    Rcpp::Named("g") = g_norm,
    Rcpp::Named("mu") = mu_unit,
    Rcpp::Named("anchor_coordinate") = anchor_coordinate,
    Rcpp::Named("u") = u,
    Rcpp::Named("basis") = basis,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("sigma_tangent") = sigma_tangent,
    Rcpp::Named("eigenvalues") = eigval,
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("scores") = scores
  );
}

// [[Rcpp::export]]
arma::mat agca_reconstruct_cpp(
    const arma::vec& anchor_coordinate,
    const arma::vec& mu,
    const arma::mat& scores,
    const arma::mat& loadings,
    int p) {
  const arma::uword n = scores.n_rows;
  const arma::uword d = loadings.n_rows;
  if (anchor_coordinate.n_elem != n) {
    Rcpp::stop("anchor_coordinate and scores have incompatible lengths.");
  }
  if (mu.n_elem != d) {
    Rcpp::stop("mu and loadings have incompatible dimensions.");
  }
  if (p < 0 || p > static_cast<int>(loadings.n_cols)) {
    Rcpp::stop("p is outside the admissible rank range.");
  }

  mat u_hat(n, d, arma::fill::zeros);
  if (p > 0) {
    u_hat = scores.cols(0, p - 1) * loadings.cols(0, p - 1).t();
  }

  mat numerator = anchor_coordinate * mu.t() + u_hat;
  return normalize_rows_impl(std::move(numerator));
}

// [[Rcpp::export]]
arma::vec agca_residual_risk_cpp(
    const arma::mat& u,
    const arma::mat& scores,
    const arma::mat& loadings,
    int max_rank) {
  if (u.n_rows != scores.n_rows || u.n_cols != loadings.n_rows) {
    Rcpp::stop("u, scores, and loadings have incompatible dimensions.");
  }
  if (max_rank < 0 || max_rank > static_cast<int>(loadings.n_cols)) {
    Rcpp::stop("max_rank is outside the admissible rank range.");
  }

  vec out(max_rank + 1);
  for (int p = 0; p <= max_rank; ++p) {
    mat residual = u;
    if (p > 0) {
      residual -= scores.cols(0, p - 1) * loadings.cols(0, p - 1).t();
    }
    out(p) = arma::accu(arma::square(residual)) / static_cast<double>(u.n_rows);
  }
  return out;
}

// [[Rcpp::export]]
arma::vec agca_geodesic_distance_cpp(arma::mat x, arma::mat y, bool normalize = true) {
  if (normalize) {
    x = normalize_rows_impl(std::move(x));
    y = normalize_rows_impl(std::move(y));
  } else {
    check_finite_matrix(x, "x");
    check_finite_matrix(y, "y");
  }
  if (x.n_cols != y.n_cols) {
    Rcpp::stop("x and y have incompatible dimensions.");
  }
  if (y.n_rows == 1 && x.n_rows > 1) {
    y = arma::repmat(y, x.n_rows, 1);
  }
  if (x.n_rows != y.n_rows) {
    Rcpp::stop("x and y must have the same number of rows, unless y has one row.");
  }

  vec out(x.n_rows);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    double inner = arma::dot(x.row(i), y.row(i));
    inner = std::max(-1.0, std::min(1.0, inner));
    out(i) = std::acos(inner);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat agca_rank_pareto_cpp(const arma::mat& x, std::string ties_method = "average") {
  check_finite_matrix(x, "x");
  if (ties_method != "average" && ties_method != "min" && ties_method != "max" &&
      ties_method != "first" && ties_method != "last" && ties_method != "random") {
    Rcpp::stop("Unsupported ties_method.");
  }

  Rcpp::RNGScope scope;
  const arma::uword n = x.n_rows;
  arma::mat out(n, x.n_cols);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    std::vector<arma::uword> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](arma::uword a, arma::uword b) {
      if (x(a, j) == x(b, j)) {
        return a < b;
      }
      return x(a, j) < x(b, j);
    });

    arma::vec ranks(n);
    arma::uword start = 0;
    while (start < n) {
      arma::uword end = start;
      while (end + 1 < n && x(ord[end + 1], j) == x(ord[start], j)) {
        ++end;
      }
      const double min_rank = static_cast<double>(start + 1);
      const double max_rank = static_cast<double>(end + 1);

      if (ties_method == "average") {
        const double avg = 0.5 * (min_rank + max_rank);
        for (arma::uword pos = start; pos <= end; ++pos) {
          ranks(ord[pos]) = avg;
        }
      } else if (ties_method == "min") {
        for (arma::uword pos = start; pos <= end; ++pos) {
          ranks(ord[pos]) = min_rank;
        }
      } else if (ties_method == "max") {
        for (arma::uword pos = start; pos <= end; ++pos) {
          ranks(ord[pos]) = max_rank;
        }
      } else {
        std::vector<arma::uword> group;
        for (arma::uword pos = start; pos <= end; ++pos) {
          group.push_back(ord[pos]);
        }
        if (ties_method == "first") {
          std::sort(group.begin(), group.end());
        } else if (ties_method == "last") {
          std::sort(group.begin(), group.end(), std::greater<arma::uword>());
        } else {
          for (arma::uword a = group.size(); a > 1; --a) {
            arma::uword b = static_cast<arma::uword>(std::floor(R::unif_rand() * a));
            std::swap(group[a - 1], group[b]);
          }
        }
        for (arma::uword offset = 0; offset < group.size(); ++offset) {
          ranks(group[offset]) = min_rank + static_cast<double>(offset);
        }
      }
      start = end + 1;
    }

    out.col(j) = (static_cast<double>(n) + 1.0) /
      (static_cast<double>(n) + 1.0 - ranks);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat agca_cdf_to_pareto_cpp(const arma::mat& u, double eps = 1e-12) {
  check_finite_matrix(u, "u");
  if (!std::isfinite(eps) || eps <= 0.0 || eps >= 0.5) {
    Rcpp::stop("eps must be a finite number in (0, 0.5).");
  }
  if (arma::any(arma::vectorise(u) < 0.0) || arma::any(arma::vectorise(u) > 1.0)) {
    Rcpp::stop("CDF values must lie in [0, 1].");
  }
  arma::mat clipped = arma::clamp(u, eps, 1.0 - eps);
  return 1.0 / (1.0 - clipped);
}

// [[Rcpp::export]]
Rcpp::List agca_tail_directions_cpp(const arma::mat& x, int k, double threshold, bool use_k) {
  return tail_directions_impl(x, k, threshold, use_k);
}

// [[Rcpp::export]]
Rcpp::List agca_decluster_runs_cpp(const arma::mat& x, int k, double threshold,
                                   bool use_k, int run) {
  if (run < 0) {
    Rcpp::stop("run must be nonnegative.");
  }
  Rcpp::List tail = tail_directions_impl(x, k, threshold, use_k);
  arma::vec radius = tail["all_radius"];
  Rcpp::IntegerVector selected = tail["index"];
  std::vector<bool> exceed(radius.n_elem, false);
  for (int idx : selected) {
    exceed[static_cast<std::size_t>(idx - 1)] = true;
  }

  std::vector<arma::uword> exceed_index;
  for (arma::uword i = 0; i < radius.n_elem; ++i) {
    if (exceed[i]) {
      exceed_index.push_back(i);
    }
  }
  if (exceed_index.empty()) {
    Rcpp::stop("No observations exceed the threshold.");
  }

  std::vector<std::vector<arma::uword>> clusters0;
  std::vector<arma::uword> current{exceed_index[0]};
  for (std::size_t i = 1; i < exceed_index.size(); ++i) {
    const int gap = static_cast<int>(exceed_index[i] - current.back() - 1);
    if (gap <= run) {
      current.push_back(exceed_index[i]);
    } else {
      clusters0.push_back(current);
      current = std::vector<arma::uword>{exceed_index[i]};
    }
  }
  clusters0.push_back(current);

  arma::uvec representatives(clusters0.size());
  Rcpp::List clusters(clusters0.size());
  for (std::size_t c = 0; c < clusters0.size(); ++c) {
    arma::uword best = clusters0[c][0];
    for (arma::uword idx : clusters0[c]) {
      if (radius(idx) > radius(best)) {
        best = idx;
      }
    }
    representatives(c) = best;
    Rcpp::IntegerVector cluster(clusters0[c].size());
    for (std::size_t h = 0; h < clusters0[c].size(); ++h) {
      cluster[h] = static_cast<int>(clusters0[c][h] + 1);
    }
    clusters[c] = cluster;
  }

  arma::mat directions = normalize_rows_impl(subset_rows(x, representatives));
  arma::vec selected_radius = radius.elem(representatives);
  return Rcpp::List::create(
    Rcpp::Named("directions") = directions,
    Rcpp::Named("radius") = selected_radius,
    Rcpp::Named("index") = to_one_based(representatives),
    Rcpp::Named("threshold") = tail["threshold"],
    Rcpp::Named("all_radius") = radius,
    Rcpp::Named("clusters") = clusters,
    Rcpp::Named("run") = run
  );
}

// [[Rcpp::export]]
arma::vec agca_principal_anchor_cpp(arma::mat g, bool normalize = true) {
  return principal_anchor_impl(std::move(g), normalize);
}

// [[Rcpp::export]]
arma::vec agca_frechet_anchor_cpp(arma::mat g, bool normalize = true,
                                  int max_iter = 100, double tol = 1e-10) {
  return frechet_anchor_impl(std::move(g), normalize, max_iter, tol);
}

// [[Rcpp::export]]
Rcpp::DataFrame agca_bootstrap_cpp(const arma::mat& g, const arma::vec& mu, int p,
                                   int B, Rcpp::IntegerVector ranks,
                                   bool fixed_anchor, std::string anchor_type,
                                   Rcpp::NumericVector anchor_vector) {
  check_finite_matrix(g, "g");
  if (B < 1) {
    Rcpp::stop("B must be positive.");
  }
  const int max_rank = static_cast<int>(g.n_cols) - 1;
  if (p < 0 || p > max_rank) {
    Rcpp::stop("p is outside the admissible rank range.");
  }
  for (int rank : ranks) {
    if (rank < 0 || rank > max_rank) {
      Rcpp::stop("ranks must contain integers between 0 and the maximum AGCA rank.");
    }
  }

  Rcpp::RNGScope scope;
  const int n = static_cast<int>(g.n_rows);
  const int nr = ranks.size();
  Rcpp::IntegerVector replicate(B * nr);
  Rcpp::IntegerVector rank_out(B * nr);
  Rcpp::NumericVector residual_risk(B * nr);
  Rcpp::NumericVector variation_explained(B * nr);

  for (int b = 0; b < B; ++b) {
    arma::mat boot_g(n, g.n_cols);
    for (int i = 0; i < n; ++i) {
      const int idx = static_cast<int>(std::floor(R::unif_rand() * n));
      boot_g.row(i) = g.row(idx);
    }

    arma::vec boot_mu = resolve_anchor_impl(
      boot_g, mu, fixed_anchor, anchor_type, anchor_vector
    );
    Rcpp::List fit = agca_core_cpp(boot_g, boot_mu, false);
    arma::mat u = fit["u"];
    arma::mat scores = fit["scores"];
    arma::mat loadings = fit["loadings"];
    arma::vec eigenvalues = fit["eigenvalues"];
    arma::vec risk = agca_residual_risk_cpp(u, scores, loadings, max_rank);
    arma::vec ave = variation_explained_from_eigenvalues(eigenvalues);

    for (int r = 0; r < nr; ++r) {
      const int out_idx = b * nr + r;
      const int rank = ranks[r];
      replicate[out_idx] = b + 1;
      rank_out[out_idx] = rank;
      residual_risk[out_idx] = risk(rank);
      variation_explained[out_idx] = rank == 0 ? 0.0 : ave(rank - 1);
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("replicate") = replicate,
    Rcpp::Named("rank") = rank_out,
    Rcpp::Named("residual_risk") = residual_risk,
    Rcpp::Named("variation_explained") = variation_explained
  );
}

// [[Rcpp::export]]
Rcpp::List agca_functional_error_cpp(const arma::mat& g,
                                     const arma::vec& anchor_coordinate,
                                     const arma::vec& mu,
                                     const arma::mat& scores,
                                     const arma::mat& loadings,
                                     const arma::mat& weights,
                                     Rcpp::IntegerVector ranks,
                                     double power,
                                     double cap) {
  check_finite_matrix(g, "g");
  check_finite_matrix(weights, "weights");
  if (!std::isfinite(power) || power <= 0.0) {
    Rcpp::stop("power must be a positive finite number.");
  }
  if (std::isnan(cap) || cap <= 0.0) {
    Rcpp::stop("cap must be positive or Inf.");
  }
  const int max_rank = static_cast<int>(loadings.n_cols);
  for (int rank : ranks) {
    if (rank < 0 || rank > max_rank) {
      Rcpp::stop("ranks must contain integers between 0 and the maximum AGCA rank.");
    }
  }

  arma::vec original = bounded_positive_exposure_impl(g, weights, power, cap);
  arma::mat reconstructed(ranks.size(), weights.n_rows);
  arma::mat relative(ranks.size(), weights.n_rows);
  for (int i = 0; i < ranks.size(); ++i) {
    arma::mat recon = agca_reconstruct_cpp(anchor_coordinate, mu, scores, loadings, ranks[i]);
    arma::vec value = bounded_positive_exposure_impl(recon, weights, power, cap);
    reconstructed.row(i) = value.t();
    relative.row(i) = (value / original - 1.0).t();
  }

  return Rcpp::List::create(
    Rcpp::Named("original") = original,
    Rcpp::Named("reconstructed") = reconstructed,
    Rcpp::Named("relative_error") = relative
  );
}
