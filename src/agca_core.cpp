#include <RcppArmadillo.h>

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

}  // namespace

// [[Rcpp::export]]
arma::vec agca_row_norms_cpp(const arma::mat& x) {
  check_finite_matrix(x, "x");
  arma::vec out(x.n_rows);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    out(i) = arma::norm(x.row(i), 2);
  }
  return out;
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
