# Unit tests for small model helpers and samplers.

params <- readRDS(reference_path)$params

test_that("bean load and pupation success", {
  expect_equal(bean_load(0, 5), 1)
  expect_equal(bean_load(4, 1), 3)
  expect_equal(pupation_success(0, 10, gamma = 0.5), 1)
  expect_lt(pupation_success(10, 0, gamma = 0.5), pupation_success(5, 0, gamma = 0.5))
  expect_gt(pupation_success(10, 10, gamma = 0.5), pupation_success(10, 0, gamma = 0.5))
})

test_that("starvation is clamped to the trial range", {
  expect_equal(clamp_starve(c(0, 10, 50)), c(starve_range[1], 10, starve_range[2]))
})

test_that("as_aphid() orders species pea, bean", {
  expect_equal(levels(as_aphid(c("bean", "pea"))), c("pea", "bean"))
  expect_true(is.na(as_aphid("aphid")))
})

test_that("aphid_age_at_mass() inverts aphid_mass() before maturity", {
  pea <- dplyr::filter(params$mass, aphid == "pea")
  age <- c(0.5, 2, 4)
  expect_equal(aphid_age_at_mass(aphid_mass(age, pea), pea), age)
  expect_equal(aphid_mass(100, pea), pea$mass_adult)
})

test_that("death and birth samplers return waiting times", {
  set.seed(1)
  pea <- params$aphid$pea
  age <- c(0, 5, 20) * minutes_per_day
  expect_true(all(sample_death_time(age, pea$death_shape, pea$death_scale) > 0))
  waits <- sample_birth_time(rep(0, 200), pea$cif_age, pea$cif_value)
  expect_true(all(waits > 0))
  # no births left past the last CIF knot
  expect_equal(sample_birth_time(max(pea$cif_age) + 1, pea$cif_age, pea$cif_value), Inf)
})

test_that("sample_aft() has the AFT median exp(lp) for log-logistic", {
  set.seed(2)
  draws <- sample_aft(rep(log(10), 2e4), sigma = 0.5, dist = "loglogistic")
  expect_equal(median(draws), 10, tolerance = 0.03)
  expect_error(sample_aft(1, 1, "gamma"), "unsupported")
})

test_that("stable age distribution is a proper distribution with growth", {
  for (sp in c("pea", "bean")) {
    sad <- stable_age_distribution(params, sp)
    expect_equal(sum(sad$prop), 1)
    expect_gt(sad$lambda[1], 1)
  }
})

test_that("apply_factors() sets direct, by-species and derived parameters", {
  p <- apply_factors(params, c(plant_area = 1234, capture_bean = 0.6, packing_pea = 5,
                               density_lag = 1.6, starvation_scale = 2, hatch_day = 4))
  expect_equal(p$plant_area, 1234)
  expect_equal(p$capture[["bean"]], 0.6)
  expect_equal(p$capture[["pea"]], params$capture[["pea"]])
  expect_equal(p$colony_packing[["pea"]], 5)
  expect_identical(p$density_lag, 2L)
  expect_equal(p$starvation_median, params$starvation_median * 2)
  expect_null(p$hatch_day)
  expect_error(apply_factors(params, c(plant_aera = 1)))
})

test_that("design values map design bounds back to factor ranges", {
  ranges <- free_parameter_ranges()
  bounds <- design_bounds(ranges)
  expect_equal(unname(design_values(setNames(bounds$binf, ranges$factor), ranges)), ranges$low)
  expect_equal(unname(design_values(setNames(bounds$bsup, ranges$factor), ranges)), ranges$high)
})

test_that("tidy_census() reshapes the census matrix", {
  census <- matrix(1:8, nrow = 2)
  out <- tidy_census(census)
  expect_named(out, c("day", "plant", "species", "n"))
  expect_equal(nrow(out), 8)
  expect_equal(out$n[out$day == 1 & out$plant == 2 & out$species == "bean"], 8)
})
