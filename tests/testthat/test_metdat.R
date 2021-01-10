context("Function metdat")

test_that("Test if metdat does its thing", {
  x <- "this is one chunk of characters & this is another"
  expect_equal(
    metdat(x, 1:2, " & ", collapse = " and "),
    "this is one chunk of characters and this is another"
  )
})
