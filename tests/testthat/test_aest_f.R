context("Function aest_f")

test_that("Test if aest_f generates three colors for a 'numeric' label", {
  set.seed(42)
  aest_out <- aest.f(runif(1222), n.cat = 3)
  expect_equal(nrow(aest_out), 3)
})

test_that("Test if aest_f recognizes 'categorical' labels", {
  aest_out <- aest.f(c("a", "a", "b", "b", "b"))
  expect_equal(nrow(aest_out), 2)
})

test_that("Test if the output fields are 'numeric' and 'character'", {
  aest_out <- aest.f(c("a", "a", "b", "b", "b"))
  expect_equal(inherits(aest_out$pch, "numeric"), TRUE)
  expect_equal(inherits(aest_out$col, "character"), TRUE)

  set.seed(42)
  aest_out <- aest.f(runif(1222), n.cat = 3)
  expect_equal(inherits(aest_out$pch, "numeric"), TRUE)
  expect_equal(inherits(aest_out$col, "character"), TRUE)
})
