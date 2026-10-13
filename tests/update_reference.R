# Regenerate the regression-test reference ------------------------------------
#
# The regression tests (tests/testthat/test-regression.R) compare seeded
# simulations, refitted parameters and vial larvae with
# tests/testthat/reference.rds. After a deliberate model or fitting change,
# rerun the analysis scripts so output/params is current, then run this from
# the project root:
#
#   Rscript tests/update_reference.R
#
# and commit the new reference with the change that caused it.

source("tests/testthat/helper-model.R")

params <- reference_params()
reference <- list(
  params = params,
  runs = reference_runs(params),
  fits = reference_fits(),
  vial = reference_vial()
)
saveRDS(reference, reference_path)
message("Wrote ", reference_path)
