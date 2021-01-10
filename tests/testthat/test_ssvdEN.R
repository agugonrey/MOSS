context("Function ssvdEN")

test_that("Testing if 
          Number of features/subjects selected by LASSO 
          equal to sparsity degrees ", {
  set.seed(42)
  sim_data <- simulate_data()
  sim_blocks <- sim_data$sim_blocks
  X <- sim_blocks$`Block 3`

  # Sampling a number of subjects and features for a fix sparsity degree.
  s.u <- sample(seq_len(nrow(X)), 1)
  s.v <- sample(seq_len(ncol(X)), 1)

  # Lasso penalties.
  expect_equal(sum(ssvdEN(X, dg.spar.features = s.v)$v != 0), s.v)
  expect_equal(sum(ssvdEN(X, dg.spar.subjects = s.u)$u != 0), s.u)
})

test_that("Testing if Smoothing by RIDGE penalty does not
          select features/subjects ", {
  set.seed(42)
  sim_data <- simulate_data()
  sim_blocks <- sim_data$sim_blocks
  X <- sim_blocks$`Block 3`

  # Sampling a number of subjects and features for a fix sparsity degree.
  s.u <- sample(seq_len(nrow(X)), 1)
  s.v <- sample(seq_len(ncol(X)), 1)

  # Lasso penalties.
  expect_equal(
    sum(ssvdEN(X, dg.spar.features = s.v, alpha.f = 0)$v != 0),
    ncol(X)
  )
  expect_equal(
    sum(ssvdEN(X, dg.spar.subjects = s.u, alpha.s = 0)$u != 0),
    nrow(X)
  )
})

test_that("Testing if Elastic Net penalty selects more subjects/features 
          than the specified sparsity degrees", {
  set.seed(42)
  sim_data <- simulate_data()
  sim_blocks <- sim_data$sim_blocks
  X <- sim_blocks$`Block 3`

  # Sampling a number of subjects and features for a fix sparsity degree.
  s.u <- sample(seq_len(nrow(X)), 1)
  s.v <- sample(seq_len(ncol(X)), 1)

  # Lasso penalties.
  expect_gte(sum(ssvdEN(X, dg.spar.features = s.v, alpha.f = 0.5)$v != 0), s.v)
  expect_gte(sum(ssvdEN(X, dg.spar.subjects = s.u, alpha.s = 0.5)$u != 0), s.u)
})
