# Test setup ------------------------------------------------------------------
#
# Loads the model and defines the reference runs and fits that the
# regression tests compare against tests/testthat/reference.rds. The same
# functions regenerate that file (tests/update_reference.R).

find_project_root <- function(dir = getwd()) {
  while (!file.exists(file.path(dir, "R", "load_model.R"))) {
    parent <- dirname(dir)
    if (parent == dir) stop("Could not find the project root (R/load_model.R)")
    dir <- parent
  }
  normalizePath(dir)
}

project_root <- find_project_root()
withr::with_dir(project_root, {
  source("R/load_model.R")
  source("R/vial_sim.R") # used by analysis/03 only, so not in load_model.R
})

# Data readers and default_params() use paths relative to the project root.
local_project_root <- function(env = parent.frame()) {
  withr::local_dir(project_root, .local_envir = env)
}

reference_path <- file.path(project_root, "tests", "testthat", "reference.rds")

# Seeded simulations covering vial mode, the two-plant scenarios, adult
# founders, no predator, and changed free parameters.
reference_runs <- function(params) {
  RNGkind("default", "default", "default")
  vp <- vial_params(params)
  changed <- apply_factors(params, c(capture_pea = 0.3, packing_bean = 12,
                                     starvation_scale = 2, density_lag = 2.4))
  list(
    vial_pea = simulate_two_plants(vp, vial_scenario("P", vp), seed = 1),
    vial_mixed = simulate_two_plants(vp, vial_scenario("M", vp), seed = 2),
    pea_pea = simulate_two_plants(params, two_plant_scenario(c(pea = 10), c(pea = 10)), seed = 3),
    pea_bean = simulate_two_plants(params, two_plant_scenario(c(pea = 10), c(bean = 10)), seed = 4),
    mixed_none = simulate_two_plants(params, two_plant_scenario(c(pea = 10, bean = 10), c(pea = 0)), seed = 5),
    adult_founders = simulate_two_plants(
      params, two_plant_scenario(c(pea = 5, bean = 5), c(bean = 10), founders = "adults"), seed = 6
    ),
    no_predator = simulate_two_plants(
      params, two_plant_scenario(c(pea = 10), c(bean = 10), predator = FALSE, run_days = 50), seed = 7
    ),
    changed_factors = simulate_two_plants(changed, two_plant_scenario(c(pea = 10), c(bean = 10)), seed = 8)
  )
}

# Key fitted values, refitted from data/.
reference_fits <- function() {
  withr::local_dir(project_root)
  demography <- fit_aphid_demography()
  behavior <- fit_predator_behavior()
  movement <- fit_movement()
  development <- fit_predator_development()
  list(
    aphid = aphid_params_minutes(demography),
    frailty = demography$frailty,
    handling = behavior$handling[c("coef", "sigma")],
    post_handling = behavior$post_handling[c("coef", "sigma")],
    age_exponent = behavior$age_exponent,
    rejection = coef(behavior$rejection),
    movement = list(speed = movement$speed, activity = coef(movement$activity),
                    leaf_area = movement$leaf_area),
    mass_anchors = aphid_mass_anchors(),
    thresholds = development$thresholds[c("threshold", "sigma", "v_bean_alone", "v_bean_mixed")],
    mortality = development$mortality[c("log_hazard", "l1_bean_exponent")],
    max_intake = development$max_intake,
    critical_food = development$critical$critical_food,
    pupation_delay = coef(development$critical$delay)
  )
}

# Seeded vial larvae with development refitted from data/ and a fixed
# pupation exponent.
reference_vial <- function() {
  withr::local_dir(project_root)
  development <- fit_predator_development()
  RNGkind("default", "default", "default")
  set.seed(9)
  list(
    mixed = simulate_vial_larvae(300, "M", development, pupation_gamma = 0.1),
    starved_day2 = simulate_vial_larvae(300, "P", development, switch = list(day = 2, to = "S"),
                                        pupation_gamma = 0.1)
  )
}

# Parameters for the reference runs: default_params() without the
# mortality model object (not used by the simulation, and large).
reference_params <- function() {
  withr::local_dir(project_root)
  params <- default_params()
  params$mortality$fit <- NULL
  params
}
