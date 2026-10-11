#' Generalized Maxwell model for normalized stress relaxation
#'
#' sigma(t) = 1 - sum_i A_i * (1 - exp(-t / tau_i))
#'
#' @param time Numeric vector of times (s), relative to the start of relaxation.
#' @param A Numeric vector of relaxation amplitudes (one per element).
#' @param tau Numeric vector of relaxation time constants (s), same length as A.
#' @return Numeric vector of normalized stress.
#' @export
maxwell_model <- function(time, A, tau) {
  stopifnot(length(A) == length(tau))
  decay <- vapply(seq_along(A),
                  function(i) A[i] * (1 - exp(-time / tau[i])),
                  numeric(length(time)))
  1 - rowSums(matrix(decay, nrow = length(time)))
}

#' Log-spaced subsampling of a relaxation curve
#'
#' Dense, evenly sampled curves over-weight late time points and have strongly
#' autocorrelated residuals. Keeping points evenly spaced in log(time) balances
#' fast and slow processes and reduces (but does not remove) autocorrelation.
#'
#' @param time Numeric vector of times (> 0).
#' @param n_points Approximate number of points to keep.
#' @return Integer indices into `time`.
#' @export
log_thin <- function(time, n_points = 200) {
  keep <- which(time > 0)
  if (length(keep) <= n_points) return(keep)
  targets <- exp(seq(log(min(time[keep])), log(max(time[keep])),
                     length.out = n_points))
  idx <- vapply(targets, function(x) keep[which.min(abs(time[keep] - x))],
                integer(1))
  sort(unique(idx))
}

#' Fit an n-element Maxwell model with multiple random starts
#'
#' Amplitudes are bounded to [0, 1]; time constants are fitted on the log scale
#' and bounded by `tau_bounds`. The fit with the lowest residual sum of squares
#' across all starts is returned, with time constants sorted in ascending order
#' so that tau1 always denotes the fastest process.
#'
#' @param time Numeric vector of times (s).
#' @param stress Numeric vector of normalized stress (starts near 1).
#' @param n_elements Number of Maxwell elements.
#' @param n_starts Number of random starting points.
#' @param tau_bounds Lower and upper bounds for every tau (s).
#' @param seed Optional random seed for reproducible starting points.
#' @return A list of class `maxwell_fit` with sorted parameters, goodness of
#'   fit, and diagnostics; NULL if every start failed.
#' @export
fit_maxwell <- function(time, stress, n_elements = 2, n_starts = 20,
                        tau_bounds = c(1, 1e5), seed = NULL) {
  ok <- is.finite(time) & is.finite(stress) & time > 0
  df <- data.frame(t = time[ok], y = stress[ok])
  if (nrow(df) < 2 * n_elements + 2) stop("Not enough data points for this model.")
  if (!is.null(seed)) set.seed(seed)

  n <- n_elements
  a_names <- paste0("A", seq_len(n))
  t_names <- paste0("logtau", seq_len(n))
  # tau is fitted on the log scale: better conditioned when time constants
  # differ by orders of magnitude.
  rhs <- paste0("maxwell_model(t, c(", paste(a_names, collapse = ", "),
                "), exp(c(", paste(t_names, collapse = ", "), ")))")
  form <- stats::as.formula(paste("y ~", rhs))

  # Starting taus are drawn across the observable window, in log space
  lo_tau <- log(max(tau_bounds[1], min(df$t)))
  hi_tau <- log(min(tau_bounds[2], max(df$t)))

  fits <- vector("list", n_starts)
  rss <- rep(NA_real_, n_starts)
  for (s in seq_len(n_starts)) {
    logtau0 <- sort(stats::runif(n, lo_tau, hi_tau))
    A0 <- stats::runif(n, 0.2, 1)
    A0 <- A0 / sum(A0) * stats::runif(1, 0.5, 0.95)
    start <- stats::setNames(c(A0, logtau0), c(a_names, t_names))
    fit <- tryCatch(
      minpack.lm::nlsLM(form, data = df, start = start,
        lower = c(rep(0, n), rep(log(tau_bounds[1]), n)),
        upper = c(rep(1, n), rep(log(tau_bounds[2]), n)),
        control = minpack.lm::nls.lm.control(maxiter = 1000)),
      error = function(e) NULL)
    if (!is.null(fit)) {
      fits[[s]] <- fit
      rss[s] <- sum(stats::residuals(fit)^2)
    }
  }
  if (all(is.na(rss))) return(NULL)

  best <- which.min(rss)          # report the best start, not the last one
  best_fit <- fits[[best]]
  cf <- stats::coef(best_fit)
  A <- unname(cf[a_names])
  tau <- exp(unname(cf[t_names]))
  ord <- order(tau)               # fixed labelling: tau1 < tau2 < ...
  A <- A[ord]
  tau <- tau[ord]

  res <- stats::residuals(best_fit)
  nobs <- length(res)
  k <- 2 * n + 1                  # model parameters + residual variance
  loglik <- as.numeric(stats::logLik(best_fit))

  structure(list(
    n_elements = n,
    A = A, tau = tau, total_relaxation = sum(A),
    rss = rss[best], rmse = sqrt(mean(res^2)),
    r_squared = 1 - rss[best] / sum((df$y - mean(df$y))^2),
    loglik = loglik,
    aic = -2 * loglik + 2 * k,
    bic = -2 * loglik + log(nobs) * k,
    n_obs = nobs,
    # Diagnostics
    starts_converged = sum(!is.na(rss)),
    starts_at_best = sum(rss <= rss[best] * 1.001, na.rm = TRUE),
    residual_lag1_acf = stats::acf(res, lag.max = 1, plot = FALSE)$acf[2],
    sum_A_exceeds_1 = sum(A) > 1,
    tau_at_bound = any(abs(log(tau) - log(tau_bounds[1])) < 1e-3 |
                       abs(log(tau) - log(tau_bounds[2])) < 1e-3),
    tau_beyond_window = any(tau > max(df$t)),
    data = df, fitted = unname(stats::fitted(best_fit)), nls = best_fit
  ), class = "maxwell_fit")
}

#' @export
print.maxwell_fit <- function(x, ...) {
  cat(sprintf("%d-element Maxwell fit (n = %d points)\n", x$n_elements, x$n_obs))
  for (i in seq_len(x$n_elements))
    cat(sprintf("  A%d = %.4f   tau%d = %.2f s\n", i, x$A[i], i, x$tau[i]))
  cat(sprintf("  RMSE = %.4g, R2 = %.5f, AIC = %.1f, BIC = %.1f\n",
              x$rmse, x$r_squared, x$aic, x$bic))
  cat(sprintf("  starts reaching the best solution: %d of %d converged\n",
              x$starts_at_best, x$starts_converged))
  if (abs(x$residual_lag1_acf) > 0.5)
    cat(sprintf(paste0("  NOTE: residual lag-1 autocorrelation = %.2f. AIC, BIC and p-values\n",
                       "        assume independent errors and will favour overly complex models.\n"),
                x$residual_lag1_acf))
  if (x$sum_A_exceeds_1) cat("  WARNING: sum of A > 1 (stress predicted below zero).\n")
  if (x$tau_at_bound) cat("  WARNING: a tau sits at its bound.\n")
  if (x$tau_beyond_window)
    cat("  NOTE: a tau exceeds the recording time; it is poorly constrained by the data.\n")
  invisible(x)
}

#' Compare Maxwell models with increasing numbers of elements
#'
#' Fits 1..max_elements elements to the same data and reports information
#' criteria, a nested F-test against the next-smaller model, and diagnostics.
#'
#' @inheritParams fit_maxwell
#' @param max_elements Largest model to fit.
#' @param ... Passed to [fit_maxwell()].
#' @return A data.frame with one row per model; the fits are attached as
#'   attribute "fits".
#' @export
compare_maxwell <- function(time, stress, max_elements = 3, ...) {
  fits <- lapply(seq_len(max_elements),
                 function(n) fit_maxwell(time, stress, n_elements = n, ...))
  out <- do.call(rbind, lapply(fits, function(f) {
    if (is.null(f)) return(NULL)
    data.frame(
      n_elements = f$n_elements, rss = f$rss, rmse = f$rmse,
      aic = f$aic, bic = f$bic,
      starts_at_best = f$starts_at_best,
      residual_lag1_acf = f$residual_lag1_acf,
      warnings = paste(c(if (f$sum_A_exceeds_1) "sum(A)>1",
                         if (f$tau_at_bound) "tau at bound",
                         if (f$tau_beyond_window) "tau beyond window"),
                       collapse = "; "))
  }))
  out$delta_aic <- out$aic - min(out$aic)
  out$delta_bic <- out$bic - min(out$bic)
  out$F_p_vs_previous <- NA_real_
  for (i in seq_len(nrow(out))[-1]) {
    df1 <- 2 * (out$n_elements[i] - out$n_elements[i - 1])
    df2 <- fits[[i]]$n_obs - 2 * out$n_elements[i]
    F <- ((out$rss[i - 1] - out$rss[i]) / df1) / (out$rss[i] / df2)
    out$F_p_vs_previous[i] <- stats::pf(F, df1, df2, lower.tail = FALSE)
  }
  attr(out, "fits") <- fits
  out
}

#' Fit Maxwell models to every sample in a long data frame
#'
#' Convenience wrapper for experiments with many samples (one row per time
#' point). Returns one row per sample with sorted parameters and diagnostics.
#'
#' @param data A data.frame in long format.
#' @param id,time,stress Column names (strings) for sample id, time and
#'   normalized stress.
#' @param t_max Optional: only fit time points up to t_max.
#' @param ... Passed to [fit_maxwell()], e.g. n_elements, n_starts, seed.
#' @return A data.frame with one row per sample; the fits are attached as
#'   attribute "fits".
#' @export
fit_maxwell_by_sample <- function(data, id, time = "time", stress = "stress",
                                  t_max = Inf, ...) {
  groups <- split(data, data[[id]], drop = TRUE)
  fits <- lapply(groups, function(d) {
    d <- d[d[[time]] <= t_max, ]
    tryCatch(fit_maxwell(d[[time]], d[[stress]], ...),
             error = function(e) { message(e$message); NULL })
  })
  rows <- lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    if (is.null(f)) return(NULL)
    p <- stats::setNames(c(f$A, f$tau),
                         c(paste0("A", seq_along(f$A)), paste0("tau", seq_along(f$tau))))
    data.frame(id = nm, as.list(p), total_relaxation = f$total_relaxation,
               rmse = f$rmse, r_squared = f$r_squared, aic = f$aic, bic = f$bic,
               starts_at_best = f$starts_at_best,
               residual_lag1_acf = f$residual_lag1_acf,
               flag = paste(c(if (f$sum_A_exceeds_1) "sum(A)>1",
                              if (f$tau_at_bound) "tau at bound",
                              if (f$tau_beyond_window) "tau beyond window"),
                            collapse = "; "),
               check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  names(out)[1] <- id
  attr(out, "fits") <- fits
  out
}
