# Regression tests: outputs must match tests/testthat/reference.rds. After a
# deliberate change, regenerate it with tests/update_reference.R.

reference <- readRDS(reference_path)

test_that("seeded simulations match the reference", {
  runs <- reference_runs(reference$params)
  expect_named(runs, names(reference$runs))
  for (nm in names(runs)) expect_equal(runs[[nm]], reference$runs[[nm]], label = nm)
})

test_that("fitted parameters match the reference", {
  fits <- reference_fits()
  expect_named(fits, names(reference$fits))
  for (nm in names(fits)) expect_equal(fits[[nm]], reference$fits[[nm]], label = nm)
})

test_that("vial larvae simulations match the reference", {
  expect_equal(reference_vial(), reference$vial)
})
