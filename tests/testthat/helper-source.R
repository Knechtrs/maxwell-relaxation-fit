# When run from a source checkout without installing, load the functions directly
if (!exists("fit_maxwell", mode = "function")) {
  src <- testthat::test_path("..", "..", "R", "maxwell.R")
  if (file.exists(src)) source(src)
}

simulate_relaxation <- function(A, tau, t = seq(0.5, 3000, by = 0.5),
                                sd = 0.003, seed = 1) {
  set.seed(seed)
  list(t = t, y = maxwell_model(t, A, tau) + rnorm(length(t), 0, sd))
}
