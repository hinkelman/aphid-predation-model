# Run the test suite from the project root:
#
#   Rscript tests/run_tests.R

testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
