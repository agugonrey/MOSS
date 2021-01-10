context("Function moss")

test_that("Test if moss detects clusters of subjects in simulated data", {
  # Load required packages.
  library(ggplot2)
  library(ggthemes)
  library(viridis)

  # Simulate data.
  set.seed(43)
  sim_data <- simulate_data()

  # Extracting simulated omic blocks.
  sim_blocks <- sim_data$sim_blocks

  # Extracting subjects and features labels.
  lab.sub <- sim_data$labels$lab.sub

  # Run moss to get clusters
  out <- moss(sim_blocks[-4],
    method = "pca",
    tSNE = TRUE,
    cluster = TRUE,
    plot = TRUE
  )

  expect_equal(
    chisq.test(out$clus_plot$dbscan.res$cluster,
      lab.sub,
      simulate.p.value = T
    )$p.val < 0.005,
    TRUE
  )
})

test_that("Test if moss detects informative features in simulated data", {
  # Load required packages.
  library(ggplot2)
  library(ggthemes)
  library(viridis)

  # Simulate data.
  set.seed(43)
  sim_data <- simulate_data()

  # Extracting simulated omic blocks.
  sim_blocks <- sim_data$sim_blocks

  # Extracting subjects and features labels.
  lab.feat <- sim_data$labels$lab.feat

  # Partial least squares with sparsity (PLS).
  out <- moss(sim_blocks[-4],
    K.X = 50,
    K.Y = 2,
    method = "pls",
    nu.v = seq(1, 100, by = 2),
    nu.u = seq(1, 100, by = 2),
    alpha.v = 1,
    alpha.u = 1,
    resp.block = 3,
    plot = TRUE
  )

  # Test accuracy at detecting features with signal
  table(out$sparse$u[, 1] != 0, lab.feat[1:2000])
  table(out$sparse$v[, 1] != 0, lab.feat[2001:3000])

  # Response's features.
  expect_equal(
    chisq.test(out$sparse$u[, 1] != 0,
      lab.feat[1:2000],
      simulate.p.value = T
    )$p.val < 0.005,
    TRUE
  )

  # Predictor's features.
  expect_equal(
    chisq.test(out$sparse$v[, 1] != 0,
      lab.feat[2001:3000],
      simulate.p.value = T
    )$p.val < 0.005,
    TRUE
  )
})

test_that("Test if moss uses FBMs correctly", {
  # Load required packages.
  library(ggplot2)
  library(ggthemes)
  library(viridis)
  library(bigstatsr)

  # Simulate data.
  set.seed(43)
  sim_data <- simulate_data()

  # Extracting simulated omic blocks.
  sim_blocks <- sim_data$sim_blocks

  # Extracting subjects and features labels.
  lab.feat <- sim_data$labels$lab.feat

  # Partial least squares with sparsity (PLS) w/FBM == TRUE.
  out_FBM <- moss(sim_blocks[-4],
    K.X = 10,
    method = "pls",
    nu.v = 4,
    nu.u = 5,
    alpha.v = 1,
    alpha.u = 1,
    use.fbm = TRUE,
    resp.block = 3
  )

  # Test if FBMs were used correctly.
  expect_equal(inherits(out_FBM$B, "FBM"), TRUE)

  # Partial least squares with sparsity (PLS) w/FBM == FALSE.
  out_matrix <- moss(sim_blocks[-4],
    K.X = 10,
    method = "pls",
    nu.v = 4,
    nu.u = 5,
    alpha.v = 1,
    alpha.u = 1,
    use.fbm = FALSE,
    resp.block = 3
  )

  # Test if solutions for both methods give the same structure.
  expect_equal(out_matrix$sparse[1:3],
    out_FBM$sparse[1:3],
    tolerance = 0.004
  )
})
