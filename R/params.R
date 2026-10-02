# Simulation parameters -----------------------------------------------------
#
# Assembles fitted parameters (from the analysis scripts' outputs) and free
# parameters (no data; defaults are provisional and should be explored by
# sensitivity analysis) into one list used by simulate_two_plants().
# Times in minutes unless noted.

#' Default parameters.
#'
#' Requires output/params/*.rds from analysis/01-03 (run those first).
default_params <- function(param_dir = "output/params") {
  read_param <- function(f) {
    path <- file.path(param_dir, f)
    if (!file.exists(path)) stop("Missing ", path, "; run the analysis scripts in analysis/ first.")
    readRDS(path)
  }
  aphid <- read_param("aphid_demography.rds")
  beh <- read_param("predator_behavior.rds")
  dev <- read_param("predator_development.rds")

  hb <- beh$behavior$handling
  pb <- beh$behavior$post_handling
  d <- dev$development
  # Mean handling time (min) of a size-matched pea aphid by a recently fed
  # larva (2 h, the lower end of the trials): log-logistic mean =
  # median * (pi sigma) / sin(pi sigma)
  mean_pea_handling <- exp(hb$coef[["(Intercept)"]] + 2 * hb$coef[["starve"]]) *
    (pi * hb$sigma) / sin(pi * hb$sigma)

  intake <- d$max_intake |>
    dplyr::filter(diet == "P") |>
    dplyr::arrange(stage) |>
    dplyr::pull(kills_per_day)

  list(
    # Aphids ---------------------------------------------------------------
    aphid = aphid$params, # by species: death_shape/scale, cif_age/value, frailty_shape
    mass = beh$mass_anchors, # provisional mass at age
    ref_age = beh$ref_age, # age (d) of size-matched prey in 2008 trials

    # Predator behavior (fitted) -------------------------------------------
    handling = list(
      coef = hb$coef, sigma = hb$sigma, # log-logistic AFT ~ aphid * starve
      age_exponent = beh$behavior$age_exponent
    ),
    post_handling = list(coef = pb$coef, sigma = pb$sigma), # Weibull AFT
    rejection = coef(beh$behavior$rejection), # logit ~ aphid_age (days)
    starve_range = c(2, 24), # hours; behavior predictions clamped

    # Predator development (fitted) ----------------------------------------
    threshold = d$thresholds$threshold, # units per instar
    threshold_sd = dev$threshold_sd,
    v_bean_alone = d$thresholds$v_bean_alone,
    v_bean_mixed = d$thresholds$v_bean_mixed,
    pea_memory = 24 * 60, # bean counts as "mixed" if pea eaten within this time
    max_intake = setNames(intake, larval_stages), # kills/day, ad lib
    mean_pea_handling = mean_pea_handling,
    mortality = d$mortality,
    pupation_gamma = dev$pupation_gamma,
    critical_food = list(
      location = unname(-coef(d$critical$pupation)[1] / coef(d$critical$pupation)[2]),
      scale = unname(1 / coef(d$critical$pupation)[2])
    ),
    pupation_delay = coef(d$critical$delay), # days ~ L4 food at starvation onset
    prepupa_days = dev$prepupa_days,
    starvation_onset = 24 * 60, # no meal for this long counts as stopping feeding

    # Free parameters (no data) ----------------------------------------------
    # Gut: digestion rate (per day) and capacity set so the ad lib kill rate
    # (max_intake) is reached: capacity = 1 + max_intake / digestion_rate.
    digestion_rate = 4,
    # Search: area searched per minute (cm^2) over plant area (cm^2) gives the
    # per-aphid encounter rate. ~1.1 mm/s speed x ~3 mm detection width.
    search_rate = 2,
    plant_area = 400,
    capture = c(pea = 0.3, bean = 0.8), # pea aphids drop off; bean rarely do
    # Bean handling with experience: multiplier on bean handling time =
    # m + (1 - m) * exp(-bean_meals / learn_meals). m calibrated to the vial
    # kill rates (analysis/04), learn_meals free ("slow learners").
    bean_handling_experienced = NA_real_,
    learn_meals = 20,
    rejection_time = 1, # minutes lost on a rejected aphid
    giving_up_time = 120, # mean minutes on a plant without a meal before leaving
    travel_time = 60, # minutes between plants
    # Starvation: lognormal time from last meal to death, median by instar
    # (days). Only the starved L4s of the diet-timing experiment constrain
    # this (most died, time not recorded).
    starvation_median = c(L1 = 1.5, L2 = 2, L3 = 2.5, L4 = 4),
    starvation_sdlog = 0.25
  )
}
