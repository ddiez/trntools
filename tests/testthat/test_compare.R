context("Network comparison")

# create data set.
x <- matrix(0, nrow = 3, ncol = 5,
            dimnames = list(
              paste0("R", 1:3),
              paste0("G", 1:5)
            ))
x["R1", "G1"] <- 1
x["R1", "G2"] <- -1
x["R2", "G1"] <- 1

test_that("identical networks are identical", {
  expect_equal(trn_compare(x, x, "all"), 1)
  expect_equal(trn_compare(x, x, "empty"), 1)
  expect_equal(trn_compare(x, x, "positive"), 1)
  expect_equal(trn_compare(x, x, "negative"), 1)
})
