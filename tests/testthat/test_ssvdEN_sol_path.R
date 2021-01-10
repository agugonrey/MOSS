context("Function ssvdEN_sol_path")

test_that("Testing if ssvdEN_sol_path detects informaticve features", {
  set.seed(42)
  sim_data <- simulate_data()
  sim_blocks <- sim_data$sim_blocks
  X <- sim_blocks$`Block 3`

  # Extracting features labels.
  lab.feat <- sim_data$labels$lab.feat[2001:3000]

  # Selecting features with ssvdEN_sol_path.
  out <- ssvdEN_sol_path(X, n.PC = 1, alpha.f = 0.5)

  # Test accuracy at detecting informative features
  expect_lt(chisq.test(out$SVD$v != 0, lab.feat)$p.val, 1e-10)
})
