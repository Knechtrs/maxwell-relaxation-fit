test_that("one-element parameters are recovered", {
  d <- simulate_relaxation(A = 0.6, tau = 50)
  f <- fit_maxwell(d$t, d$y, n_elements = 1, seed = 1)
  expect_equal(f$A, 0.6, tolerance = 0.02)
  expect_equal(f$tau, 50, tolerance = 0.03)
})

test_that("two-element parameters are recovered across noise realisations", {
  for (s in 1:5) {
    d <- simulate_relaxation(A = c(0.3, 0.5), tau = c(20, 800), seed = s)
    f <- fit_maxwell(d$t, d$y, n_elements = 2, seed = s)
    expect_equal(f$A, c(0.3, 0.5), tolerance = 0.05)
    expect_equal(f$tau, c(20, 800), tolerance = 0.08)
  }
})

test_that("taus are always returned in ascending order (no label switching)", {
  d <- simulate_relaxation(A = c(0.4, 0.4), tau = c(5, 1500))
  for (s in 1:5) {
    f <- fit_maxwell(d$t, d$y, n_elements = 2, seed = s)
    expect_lt(f$tau[1], f$tau[2])
  }
})

test_that("reported parameters belong to the lowest-RSS start", {
  d <- simulate_relaxation(A = c(0.3, 0.5), tau = c(20, 800))
  f <- fit_maxwell(d$t, d$y, n_elements = 2, n_starts = 10, seed = 3)
  rss_from_params <- sum((f$data$y - maxwell_model(f$data$t, f$A, f$tau))^2)
  expect_equal(rss_from_params, f$rss, tolerance = 1e-8)
})

test_that("model comparison picks two elements when two are present", {
  d <- simulate_relaxation(A = c(0.3, 0.5), tau = c(20, 800))
  cmp <- compare_maxwell(d$t, d$y, max_elements = 3, seed = 1)
  expect_equal(cmp$n_elements[which.min(cmp$bic)], 2)
})

test_that("with independent noise, BIC does not over-fit a one-element curve", {
  d <- simulate_relaxation(A = 0.6, tau = 50)
  cmp <- compare_maxwell(d$t, d$y, max_elements = 2, seed = 1)
  expect_equal(cmp$n_elements[which.min(cmp$bic)], 1)
})

test_that("log_thin returns sorted, unique indices", {
  t <- seq(0.1, 3000, by = 0.1)
  idx <- log_thin(t, 200)
  expect_false(is.unsorted(idx))
  expect_equal(length(idx), length(unique(idx)))
  expect_lte(length(idx), 200)
})
