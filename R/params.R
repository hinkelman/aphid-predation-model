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
  mov <- beh$behavior$movement
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

    # Search and movement (2008 tracking + physical assumptions) -------------
    # Area searched per minute = speed x activity x detection width, scaled
    # by instar body size: speed and width both scale with body length, so
    # search area ~ (length / L4 length)^2.
    # Encounter rate per aphid per minute = search area / plant_area.
    speed = mov$speed, # mm/min while moving, L4 (median, 2008 tracking)
    activity = coef(mov$activity), # logit(activity) ~ aphid * starve (2008)
    # Detection width (mm), L4: body width ~2.1 mm (taxonomic description)
    # plus ~0.7 mm detection distance each side (C. septempunctata L4).
    detection_width = 3.5,
    # Body length by instar (mm): L1 1.9 and L4 7.0 (taxonomic description),
    # L2-L3 interpolated geometrically.
    larva_length = c(L1 = 1.9, L2 = 2.95, L3 = 4.57, L4 = 7.0),
    # Searchable plant surface (cm^2): a large fava plant (~6 weeks, ~50 cm,
    # 10-12 compound leaves, ~1200 cm^2 leaf area per side) - both leaf
    # surfaces plus stem.
    plant_area = 2500,
    # Leaving a plant: post-handling (leaf-leaving) times were measured on a
    # single leaf of trial_leaf_area cm^2; on a plant the larva moves on to
    # other leaves, so times scale by plant_area / trial_leaf_area x
    # leave_scaling (assumption; 1 = proportional to area). With no meal on
    # the current plant, the pea post-handling model at current hunger applies.
    trial_leaf_area = mov$leaf_area,
    # Aphid aggregation (assumptions): colony area (cm^2) = colony_min_area +
    # aphid density (mg adult-mass equivalents) / colony_density (mg/cm^2),
    # capped at plant_area. 10 cm^2 ~ a growing tip or one leaflet underside;
    # 4 mg/cm^2 ~ 1 pea or ~4.4 bean aphids per cm^2.
    colony_min_area = 10,
    colony_density = 4,
    leave_scaling = 1,
    # Plants are separate: walk down one plant, across soil, up the other.
    # Travel time = path / (speed x activity).
    plant_path = 1200, # mm (~50 cm down, ~20 cm across, ~50 cm up)
    capture = c(pea = 0.3, bean = 0.8), # pea aphids drop off; bean rarely do
    # Prey size vs larval size (assumptions): capture is multiplied by
    # 1 / (1 + (prey length / (prey_size_ratio x larval length))^size_steepness);
    # aphid length (mm) = aphid_length_coef x mass(mg)^(1/3) (adult pea ~3.75
    # mm, adult bean ~2.3 mm). prey_size_ratio = 2 is set so that an L1 vs
    # 0.9 mg prey gets ~0.9: diet-experiment L1s killed ~3.7 size-matched
    # (0.9 mg) aphids per day. An L1 vs an adult pea aphid (never offered in
    # any experiment) then gets ~0.5; an L4 vs any aphid ~1.
    prey_size_ratio = 2,
    size_steepness = 4,
    aphid_length_coef = 2.4,

    # Free parameters (no data) ----------------------------------------------
    # Aphid carrying capacity per plant, in mg of adult-mass equivalents
    # (each aphid counts its species' adult mass; both species combined).
    # Births are thinned by max(0, 1 - density / capacity); Inf turns density
    # dependence off. Scaled to the large plant: ~2 mg per cm^2 of searchable
    # surface; 5000 mg ~ 1300 pea or ~5500 bean aphids.
    aphid_capacity = 5000,
    # Gut: digestion rate (per day); gut capacity is set so the ad lib kill
    # rate (max_intake) is reached (see simulate_two_plants()).
    digestion_rate = 4,
    # Bean handling with experience: multiplier on bean handling time =
    # m + (1 - m) * exp(-bean_meals / learn_meals). m calibrated to the vial
    # kill rates (analysis/04), learn_meals free ("slow learners").
    bean_handling_experienced = if (file.exists(file.path(param_dir, "simulation.rds"))) {
      readRDS(file.path(param_dir, "simulation.rds"))$bean_handling_experienced
    } else NA_real_,
    learn_meals = 20,
    rejection_time = 1, # minutes lost on a failed attack or rejected aphid
    # Starvation: lognormal time from last meal to death, median by instar
    # (days). Only the starved L4s of the diet-timing experiment constrain
    # this (most died, time not recorded).
    starvation_median = c(L1 = 1.5, L2 = 2, L3 = 2.5, L4 = 4),
    starvation_sdlog = 0.25
  )
}
