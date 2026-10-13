# Data readers: sizes and the documented cleaning fixes.

local_project_root()

test_that("clip-cage life tables", {
  cc <- read_clip_cage()
  expect_equal(nrow(cc), 1402)
  expect_equal(levels(cc$aphid), c("pea", "bean"))
  expect_false(any(cc$aphid == "pea" & cc$id == 23L)) # sterile aphid dropped
  expect_true(any(read_clip_cage(drop_sterile = FALSE)$id == 23L))
  expect_equal(min(cc$date), as.Date("2009-12-15")) # 2010 typo fixed
  lifespans <- clip_cage_lifespans(cc)
  expect_equal(nrow(lifespans), dplyr::n_distinct(cc$aphid, cc$id))
})

test_that("behavior trials", {
  trials <- read_behavior_trials()
  expect_equal(nrow(trials), 125)
  expect_gte(min(trials$post_handle), 0.05) # floored for log-time models
  expect_true(all(trials$activity >= 0 & trials$activity <= 1, na.rm = TRUE))
  expect_equal(nrow(read_handle_age()), 132)
})

test_that("diet experiments", {
  intervals <- read_diet_experiment()
  expect_equal(nrow(intervals), 1442)
  expect_equal(levels(intervals$diet), c("P", "M", "B"))
  # single-species diets have no kills of the other species
  expect_true(all(intervals$bean[intervals$diet == "P"] == 0, na.rm = TRUE))
  expect_true(all(intervals$pea[intervals$diet == "B"] == 0, na.rm = TRUE))
  # non-feeding stages eat nothing
  expect_true(all(intervals$pea[intervals$stage %in% c("Pupa", "Adult")] == 0))
  expect_equal(nrow(read_diet_timing()), 148)
  expect_equal(nrow(read_diet_summary()), 89)
})
