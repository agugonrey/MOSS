context("Reduced rank regression")

test_that("Testing if MOSS gives same results than matrix general LM", {
  set.seed(42)
  X <- matrix(stats::rnorm(1000),10,100)
  y <- X[,1:5] %*% rep(50,5) + stats::rnorm(10,sd=0.05)
  
  b.prod <- crossprod(MASS::ginv(crossprod(X)), crossprod(X,y))
  
  b.moss <- moss(data.blocks = list(y,X),scale. = F,
                 resp.block = 1,
                 method = 'rrr',
                 K.X=10,
                 K.Y = 1)$B
  
  expect_equal(b.prod , b.moss)
})



test_that("Testing if MOSS gives same results than matrix general LM (multivariate multiple regression)", {
  sim_data <- simulate_data()
  sim_blocks <- sim_data$sim_blocks
  
  X <- sim_blocks$`Block 1`
  y <- sim_blocks$`Block 3`
  
  b.prod <- crossprod(MASS::ginv(crossprod(X)), crossprod(X,y))
  
  b.moss <- moss(data.blocks = list(y,X),scale. = F,
                 resp.block = 1,
                 method = 'rrr',
                 K.X=min(dim(X)),
                 K.Y = min(dim(X), dim(y)))$B
  expect_equal(b.prod , b.moss)
})

