test_that("fit_maxwell_by_sample returns one row per sample with sorted taus", {
  d1 <- simulate_relaxation(A = c(0.3, 0.5), tau = c(20, 800), seed = 1)
  d2 <- simulate_relaxation(A = c(0.4, 0.4), tau = c(5, 300), seed = 2)
  long <- rbind(data.frame(gel = "fast", time = d1$t, stress = d1$y),
                data.frame(gel = "slow", time = d2$t, stress = d2$y))
  res <- fit_maxwell_by_sample(long, id = "gel", n_elements = 2, n_starts = 8, seed = 1)
  expect_equal(nrow(res), 2)
  expect_true(all(res$tau1 < res$tau2))
  expect_equal(res$tau2[res$gel == "slow"], 300, tolerance = 0.08)
})
