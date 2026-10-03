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
    # Searchable plant surface (cm^2): a medium fava plant (~3-4 weeks, ~30
    # cm, 5-6 compound leaves, ~380 cm^2 leaf area per side) - both leaf
    # surfaces plus stem. Chosen so one larva can deplete a young colony and
    # give up on the plant within its larval life (a 1200-2500 cm^2 plant
    # takes days to give up on).
    plant_area = 800,
    # Leaving a plant: post-handling (leaf-leaving) times were measured on a
    # single leaf of trial_leaf_area cm^2; on a plant the larva moves on to
    # other leaves, so times scale by plant_area / trial_leaf_area x
    # leave_scaling (assumption; 1 = proportional to area). With no meal on
    # the current plant, the pea post-handling model at current hunger applies.
    trial_leaf_area = mov$leaf_area,
    # Aphid aggregation: colony area (cm^2) = colony_min_area + sum over
    # species of (count x adult mass) / colony_packing (mg/cm^2), capped at
    # plant_area. 10 cm^2 ~ a growing tip or one leaflet underside (assumed).
    # Packing: no published per-area colony densities. Body footprints (pea
    # ~3.5 x 1.6 mm, bean ~2 x 1.2 mm) cap touching adults at ~68 and ~37
    # mg/cm^2. Bean aphids form dense colonies on stems and shoot tips; pea
    # aphids form loose colonies, aggregating within ~1-3 body lengths
    # (Social aggregation in pea aphids, PMC3869777) and dispersing from the
    # colony as late nymphs. Defaults: pea 8 mg/cm^2 (~2 aphids/cm^2), bean
    # 20 mg/cm^2 (~22 aphids/cm^2).
    colony_min_area = 10,
    colony_packing = c(pea = 8, bean = 20),
    leave_scaling = 1,
    # Plants are separate: walk down one plant, across soil, up the other.
    # Travel time = path / (speed x activity).
    plant_path = 900, # mm (~30 cm down, ~30 cm across, ~30 cm up)
    # Capture per encounter. Pea aphids escape by walking off or dropping:
    # adult H. convergens consumed 1 of 72 pea aphids encountered in alfalfa
    # field arenas (Nelson & Rosenheim 2006); ladybirds trigger dropping >3x
    # as often as damsel bugs (Losey & Denno 1998). No larval or bean aphid
    # data; bean aphids rarely drop.
    capture = c(pea = 0.1, bean = 0.8),
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
    # Aphid density dependence, after the ALMaSS aphid model (Thomsen, Duan &
    # Topping 2024, Food Ecol. Syst. Model. J. 5:e123747):
    # * delayed crowding mortality: extra hazard ln(1 + k x) per day, x =
    #   aphids per gram of plant density_lag days earlier; k from their Table
    #   5 (their density-independent survivorship Sa = 0.9 is not used: it
    #   stands for field predators, and clip-cage lifespans set baseline
    #   mortality here).
    # * winged emigration: % winged offspring = slope x + gs_coef GS +
    #   intercept (Carter 1982, adapted by ALMaSS), GS = plant growth stage
    #   (Zadoks-type 1-10 scale; ~3.5 for a vegetative plant).
    crowding_k = c(pea = 0.012, bean = 0.02), # per (aphids / g)
    # Lag 1 day, not ALMaSS's 4: on a single plant, with growth ~0.38/day,
    # a 4-day lag (and 2 days) gives boom-bust cycles; ALMaSS works on 10 x 10
    # m field cells with immigration and background mortality.
    density_lag = 1L, # days
    # Crowding between species: a species' crowding density = own count +
    # competition_alpha x other species' count. ALMaSS models species
    # separately (alpha = 0); with alpha = 1 pea aphids exclude bean aphids on
    # a shared plant (bean has the higher k). Unknown; explored in the
    # sensitivity analysis.
    competition_alpha = 0.5,
    alate = list(slope = 2.603, growth_stage_coef = 0.847, intercept = -27.189,
                 growth_stage = 3.5),
    # Fresh green biomass of the plant (g): ~370 cm^2 leaf area per side /
    # SLA 25.7 mm^2/mg ~ 1.4 g dry leaf, plus stem, at ~90% water.
    plant_biomass = 20,
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
