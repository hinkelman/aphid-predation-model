# Invariants of simulate_two_plants() that hold for any parameters.

params <- readRDS(reference_path)$params

test_that("the same seed gives the same run", {
  sc <- two_plant_scenario(c(pea = 10), c(bean = 10), run_days = 15)
  expect_identical(simulate_two_plants(params, sc, seed = 11), simulate_two_plants(params, sc, seed = 11))
})

test_that("census has one row per day, plant and species, with non-negative counts", {
  run <- simulate_two_plants(params, two_plant_scenario(c(pea = 10), c(bean = 10), run_days = 20), seed = 12)
  expect_equal(nrow(run$census), 21 * 4)
  expect_setequal(unique(run$census$day), 0:20)
  expect_setequal(unique(run$census$species), c("pea", "bean"))
  expect_true(all(run$census$n >= 0))
  expect_true(all(run$emigrants >= 0))
})

test_that("day-0 census holds the founders", {
  run <- simulate_two_plants(params, two_plant_scenario(c(pea = 7), c(bean = 4), run_days = 2), seed = 13)
  day0 <- dplyr::filter(run$census, day == 0)
  expect_equal(day0$n[day0$plant == 1 & day0$species == "pea"], 7)
  expect_equal(day0$n[day0$plant == 2 & day0$species == "bean"], 4)
  expect_equal(sum(day0$n), 11)
})

test_that("meals and stages are consistent with the predator's history", {
  run <- simulate_two_plants(params, two_plant_scenario(c(pea = 10), c(pea = 10)), seed = 14)
  meals <- run$meals
  expect_gt(nrow(meals), 0)
  expect_false(is.unsorted(meals$time))
  expect_true(all(meals$handling > 0))
  expect_true(all(meals$plant %in% 1:2))
  expect_true(all(meals$species %in% c("pea", "bean")))
  expect_false(is.unsorted(run$predator$stages$time))
  expect_true(all(diff(match(run$predator$stages$stage, larval_stages)) == 1))
  expect_true(run$predator$fate %in% c("alive", "died", "starved", "failed_pupation", "pupated"))
})

test_that("without a predator there are no meals and no fate", {
  run <- simulate_two_plants(params, two_plant_scenario(c(pea = 10), c(pea = 0), predator = FALSE,
                                                        run_days = 10), seed = 15)
  expect_equal(nrow(run$meals), 0)
  expect_true(is.na(run$predator$fate))
})

test_that("vial mode keeps prey numbers constant", {
  vp <- vial_params(params)
  run <- simulate_two_plants(vp, vial_scenario("M", vp, run_days = 10), seed = 16)
  totals <- run$census |>
    dplyr::filter(plant == 1) |>
    dplyr::summarise(n = sum(n), .by = c(day, species))
  expect_equal(unique(totals$n[totals$species == "pea"]), 15)
  expect_equal(unique(totals$n[totals$species == "bean"]), 15)
  expect_gt(nrow(run$meals), 0)
})
