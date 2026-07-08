#' Simulate the 10-dimensional AGCA example design
#'
#' Generates the 10-dimensional heavy-tailed design used as the package example.
#' Variables `X1`--`X8` share a low-dimensional logistic-block extremal
#' mechanism. Variables `X9` and `X10` contain independent Pareto sources, so
#' selected extremes include near-axis regimes in addition to the shared
#' low-rank angular structure.
#'
#' @param n Number of observations.
#' @param seed Optional random seed.
#' @param theta Logistic dependence parameter in `(0, 1)`.
#' @param tau Nonnegative finite-threshold noise scale.
#' @param axis9_scale,axis10_scale Positive scales for the independent Pareto
#'   sources in variables `X9` and `X10`.
#' @return A data frame with variables `X1`, ..., `X10` and a factor `regime`
#'   giving the dominant latent source for each observation.
#' @export
#' @examples
#' x <- simulate_agca_10d(n = 500, seed = 1)
#' fit <- agca(x[paste0("X", 1:10)], k = 75, p = 3)
#' agca_rank_summary(fit)
simulate_agca_10d <- function(n = 10000L, seed = NULL, theta = 0.45, tau = 0.25,
                              axis9_scale = 1, axis10_scale = 1) {
  n <- .check_scalar_integer(n, "n", lower = 1L)
  if (!is.null(seed)) {
    if (length(seed) != 1L || is.na(seed) || seed != as.integer(seed)) {
      stop("seed must be NULL or a single integer.", call. = FALSE)
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    set.seed(as.integer(seed))
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  if (!is.finite(theta) || theta <= 0 || theta >= 1) {
    stop("theta must lie in (0, 1).", call. = FALSE)
  }
  if (!is.finite(tau) || tau < 0) {
    stop("tau must be nonnegative.", call. = FALSE)
  }
  if (!is.finite(axis9_scale) || axis9_scale <= 0 ||
      !is.finite(axis10_scale) || axis10_scale <= 0) {
    stop("axis9_scale and axis10_scale must be positive.", call. = FALSE)
  }

  rays8 <- agca_shared_rays_8d()
  z_shared <- agca_rlogistic_pareto(n, d = nrow(rays8), theta = theta)
  shared_signal <- z_shared %*% rays8
  axis9_signal <- axis9_scale * agca_rpareto(n)
  axis10_signal <- axis10_scale * agca_rpareto(n)

  x_signal <- cbind(shared_signal, axis9_signal, axis10_signal)
  x <- x_signal + tau * matrix(stats::rexp(length(x_signal)), nrow(x_signal), ncol(x_signal))
  colnames(x) <- paste0("X", seq_len(10L))

  regime_score <- cbind(
    shared_low_rank = row_norms(shared_signal),
    axis_9 = axis9_signal,
    axis_10 = axis10_signal
  )
  regime_levels <- c("shared_low_rank", "axis_9", "axis_10")
  regime <- factor(
    colnames(regime_score)[max.col(regime_score, ties.method = "first")],
    levels = regime_levels
  )

  out <- as.data.frame(x)
  out$regime <- regime
  attr(out, "parameters") <- list(
    n = n,
    seed = seed,
    theta = theta,
    tau = tau,
    axis9_scale = axis9_scale,
    axis10_scale = axis10_scale
  )
  out
}

agca_rpareto <- function(n, shape = 1, scale = 1) {
  scale * stats::runif(n)^(-1 / shape)
}

agca_rpositive_stable <- function(n, alpha) {
  u <- stats::runif(n, min = 0, max = pi)
  w <- stats::rexp(n)
  sin(alpha * u) / (sin(u)^(1 / alpha)) *
    (sin((1 - alpha) * u) / w)^((1 - alpha) / alpha)
}

agca_rlogistic_frechet <- function(n, d, theta) {
  stable <- agca_rpositive_stable(n, theta)
  expo <- matrix(stats::rexp(n * d), n, d)
  sweep(1 / expo, 1L, stable, "*")^theta
}

agca_frechet_to_pareto <- function(z) {
  1 / (-expm1(-1 / z))
}

agca_rlogistic_pareto <- function(n, d, theta) {
  agca_frechet_to_pareto(agca_rlogistic_frechet(n, d, theta))
}

agca_shared_contrast_basis_8d <- function() {
  raw <- cbind(
    gradient = seq(-1, 1, length.out = 8L),
    block = c(rep(1, 4L), rep(-1, 4L))
  )
  raw <- scale(raw, center = TRUE, scale = FALSE)
  qr.Q(qr(raw), complete = FALSE)
}

agca_shared_rays_8d <- function() {
  mu8 <- canonical_anchor(8L)
  basis <- agca_shared_contrast_basis_8d()
  scores <- rbind(
    c(-0.42, -0.26),
    c(-0.30, 0.30),
    c(0.00, -0.36),
    c(0.28, 0.24),
    c(0.48, -0.02)
  )
  rays <- matrix(mu8, nrow(scores), 8L, byrow = TRUE) + scores %*% t(basis)
  if (any(rays <= 0)) {
    stop("Shared low-rank rays left the positive orthant.", call. = FALSE)
  }
  normalize_rows(rays)
}
