# Two-plant discrete-event simulation ---------------------------------------
#
# One coccinellid larva foraging on two plants that hold pea and bean
# aphids. Continuous time in minutes. Aphid events (births, deaths) are kept
# in vectors and the next one is found with which.min(); the predator is a
# state machine with its own schedule. See docs/model-design.md.
#
# Plants are 1 and 2. Species are coded PEA = 1, BEAN = 2 (index into
# species_names).

species_names <- c("pea", "bean")
PEA <- 1L
BEAN <- 2L

# Aphid event types
BIRTH <- 1L
DEATH <- 2L # life-table (Weibull) death
CROWDING_DEATH <- 3L
EMIGRATION <- 4L # a winged adult leaves the plant

# Linear predictor of the 2008 behavior models, all fitted ~ aphid * starve:
# coefficients `b`, prey species code(s) and hunger `starve_h` (hours).
aphid_starve_lp <- function(b, species, starve_h) {
  is_bean <- species == BEAN
  b[["(Intercept)"]] + is_bean * b[["aphidbean"]] + starve_h * b[["starve"]] +
    is_bean * starve_h * b[["aphidbean:starve"]]
}

#' Run one simulation.
#'
#' @param params from default_params()
#' @param scenario list with
#'   * `aphids`: data frame with columns plant, species ("pea"/"bean"),
#'     n, age (days) - the initial aphids
#'   * `predator`: list(start_day, plant, stage = "L1"), or NULL for no predator
#'   * `run_days`: total length of the run
#'   * `vial`: optional list(species_share = c(pea, bean), n = aphids
#'     available) - vial mode: no aphid demography, killed aphids are
#'     replaced immediately by size-matched prey, and the larva never
#'     leaves (used for calibration against the vial experiments)
#' @return list(census, emigrants, meals, predator)
simulate_two_plants <- function(params, scenario, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params
  vial_mode <- !is.null(scenario$vial)
  run_end <- scenario$run_days * minutes_per_day
  t <- 0

  # ---- Per-species parameters (indexed by species code) ---------------------
  aphid_param <- function(field) vapply(p$aphid[species_names], \(a) a[[field]], numeric(1))
  mass_param <- function(field) unname(setNames(p$mass[[field]], as.character(p$mass$aphid))[species_names])

  death_shape <- aphid_param("death_shape")
  death_scale <- aphid_param("death_scale")
  frailty_shape <- aphid_param("frailty_shape")
  cif_age <- lapply(p$aphid[species_names], `[[`, "cif_age")
  cif_value <- lapply(p$aphid[species_names], `[[`, "cif_value")
  cif_max <- vapply(cif_value, max, numeric(1))
  cif_last_age <- vapply(cif_age, max, numeric(1))
  # The birth CIF has knots at whole days (0, 1440, 2880, ... minutes), so
  # evaluating it is O(1) arithmetic instead of approx().
  stopifnot(all(vapply(cif_age, \(a) isTRUE(all.equal(diff(a), rep(minutes_per_day, length(a) - 1))), logical(1))))

  m_adult <- mass_param("mass_adult")
  m_neonate <- mass_param("mass_neonate")
  growth_rate <- mass_param("growth_rate")
  age_mature <- mass_param("age_mature") * minutes_per_day

  mass_of <- function(species, age_days) {
    min(m_neonate[species] * exp(growth_rate[species] * age_days), m_adult[species])
  }

  # ---- Aphid store ---------------------------------------------------------
  # Parallel vectors indexed by slot. Slots of dead aphids are reused (free
  # stack), so vectors scale with the living population, not with every
  # aphid ever born. Unused slots have a_next = Inf, so which.min(a_next)
  # over the whole vector finds the next aphid event without copying.
  #
  # Each aphid has a "base" next event (BIRTH, DEATH or, for a winged adult,
  # EMIGRATION) from its own life history, and a crowding-death time a_dd
  # (redrawn daily); a_next / a_type hold whichever comes first.
  cap <- 1024L
  a_sp <- integer(cap); a_plant <- integer(cap); a_birth <- numeric(cap)
  a_fec <- numeric(cap); a_alive <- logical(cap); a_alate <- logical(cap)
  a_base_next <- rep(Inf, cap); a_base_type <- integer(cap); a_dd <- rep(Inf, cap)
  a_next <- rep(Inf, cap); a_type <- integer(cap)
  n_used <- 0L # high-water mark of slots ever used
  free <- integer(0); n_free <- 0L
  N <- matrix(0L, 2, 2) # living aphids, plant x species

  grow <- function() {
    a_sp <<- c(a_sp, integer(cap)); a_plant <<- c(a_plant, integer(cap))
    a_birth <<- c(a_birth, numeric(cap)); a_fec <<- c(a_fec, numeric(cap))
    a_alive <<- c(a_alive, logical(cap)); a_alate <<- c(a_alate, logical(cap))
    a_base_next <<- c(a_base_next, rep(Inf, cap)); a_base_type <<- c(a_base_type, integer(cap))
    a_dd <<- c(a_dd, rep(Inf, cap))
    a_next <<- c(a_next, rep(Inf, cap)); a_type <<- c(a_type, integer(cap))
    cap <<- cap * 2L
  }

  add_aphid <- function(species, plant, birth_time, alate = FALSE) {
    if (n_free > 0L) {
      i <- free[n_free]
      n_free <<- n_free - 1L
    } else {
      if (n_used == cap) grow()
      n_used <<- n_used + 1L
      i <- n_used
    }
    a_sp[i] <<- species; a_plant[i] <<- plant; a_birth[i] <<- birth_time
    a_fec[i] <<- sample_fecundity(1, frailty_shape[[species]])
    a_alive[i] <<- TRUE
    a_alate[i] <<- alate
    a_dd[i] <<- if (vial_mode) Inf else t + rexp_or_inf(1, dd_rate[plant, species])
    N[plant, species] <<- N[plant, species] + 1L
    schedule_aphid(i)
    i
  }

  remove_aphid <- function(i) {
    a_alive[i] <<- FALSE
    a_next[i] <<- Inf; a_base_next[i] <<- Inf; a_dd[i] <<- Inf
    N[a_plant[i], a_sp[i]] <<- N[a_plant[i], a_sp[i]] - 1L
    n_free <<- n_free + 1L
    free[n_free] <<- i
  }

  # In vial mode prey were replaced daily with size-matched aphids, so prey
  # age is held at the size-matched reference age.
  aphid_age_days <- function(i) {
    if (vial_mode) return(p$ref_age[[species_names[a_sp[i]]]])
    (t - a_birth[i]) / minutes_per_day
  }

  # ---- Aphid life events ---------------------------------------------------
  cif_eval <- function(sp, age) {
    if (age >= cif_last_age[sp]) return(cif_max[sp])
    x <- age / minutes_per_day
    k <- floor(x)
    v <- cif_value[[sp]]
    v[k + 1L] + (v[k + 2L] - v[k + 1L]) * (x - k)
  }

  # Time to next birth by inversion of the CIF: sample_birth_time() with
  # O(1) CIF lookups instead of approx().
  birth_wait <- function(sp, age, fec) {
    target <- cif_eval(sp, age) + rexp(1) / fec
    if (target >= cif_max[sp]) return(Inf)
    v <- cif_value[[sp]]
    k <- findInterval(target, v, checkSorted = FALSE, checkNA = FALSE) # largest knot <= target
    (k - 1 + (target - v[k]) / (v[k + 1L] - v[k])) * minutes_per_day - age
  }

  # Next event of aphid i: the earlier of its crowding death and its
  # base event.
  set_next <- function(i) {
    if (a_dd[i] < a_base_next[i]) {
      a_next[i] <<- a_dd[i]; a_type[i] <<- CROWDING_DEATH
    } else {
      a_next[i] <<- a_base_next[i]; a_type[i] <<- a_base_type[i]
    }
  }

  # Draw aphid i's base event: death (conditional Weibull) vs its next birth
  # or, if winged, its departure at maturity - winged aphids leave without
  # reproducing on the plant.
  schedule_aphid <- function(i) {
    if (vial_mode) { a_next[i] <<- Inf; return(invisible()) }
    sp <- a_sp[i]
    age <- t - a_birth[i]
    death <- sample_death_time(age, death_shape[[sp]], death_scale[[sp]])
    other <- if (a_alate[i]) max(age_mature[[sp]] - age, 0) else birth_wait(sp, age, a_fec[i])
    if (other < death) {
      a_base_next[i] <<- t + other; a_base_type[i] <<- if (a_alate[i]) EMIGRATION else BIRTH
    } else {
      a_base_next[i] <<- t + death; a_base_type[i] <<- DEATH
    }
    set_next(i)
  }

  # ---- Aphid density dependence --------------------------------------------
  # After the ALMaSS aphid model (Thomsen et al. 2024):
  # * Delayed crowding mortality: an extra death hazard ln(1 + k x) per day,
  #   x = crowding density per gram of plant `density_lag` days earlier
  #   (daily census), k by species. Crowding density for a species = own
  #   count + competition_alpha x other species' count (alpha = 0: only
  #   intraspecific crowding, as in ALMaSS; 1: all aphids count equally).
  #   Rates change only at day boundaries, when every aphid's crowding death
  #   time is redrawn (exact: the hazard is constant within a day).
  # * Winged emigration: a newborn becomes winged with probability
  #   (alate_slope x + alate_gs GS + alate_intercept) / 100 at the current
  #   crowding density (Carter 1982); winged aphids leave at maturity without
  #   reproducing.
  dd_k <- unname(p$crowding_k[species_names])
  dd_rate <- matrix(0, 2, 2) # plant x species, per minute
  alpha <- p$competition_alpha
  rexp_or_inf <- function(n, rate) { x <- rexp(n, pmax(rate, 1e-300)); x[rate <= 0] <- Inf; x }

  # crowding density (aphids per g) felt by each species, from counts c(pea, bean)
  crowding <- function(counts) (counts + alpha * rev(counts)) / p$plant_biomass

  update_dd_rates <- function(lagged) { # lagged: plant x species counts
    for (pl in 1:2) dd_rate[pl, ] <<- log1p(dd_k * crowding(lagged[pl, ])) / minutes_per_day
  }

  # Vectorized set_next() over all living aphids after redrawing a_dd.
  redraw_crowding_deaths <- function() {
    idx <- which(a_alive)
    if (!length(idx)) return(invisible())
    rate <- dd_rate[cbind(a_plant[idx], a_sp[idx])]
    a_dd[idx] <<- t + rexp_or_inf(length(idx), rate)
    dd_first <- a_dd[idx] < a_base_next[idx]
    a_next[idx] <<- ifelse(dd_first, a_dd[idx], a_base_next[idx])
    a_type[idx] <<- ifelse(dd_first, CROWDING_DEATH, a_base_type[idx])
  }

  newborn_is_alate <- function(plant, species) {
    x <- crowding(N[plant, ])[species]
    pr_alate <- (p$alate$slope * x + p$alate$growth_stage_coef * p$alate$growth_stage + p$alate$intercept) / 100
    runif(1) < pr_alate
  }

  # ---- Initial aphids ------------------------------------------------------
  # age NA = draw founder ages from the stable age distribution
  # (two_plant_scenario(founders = "stable"))
  for (r in seq_len(nrow(scenario$aphids))) {
    row <- scenario$aphids[r, ]
    sp <- match(row$species, species_names)
    ages <- if (is.na(row$age)) sample_founder_ages(row$n, p, row$species) else rep(row$age, row$n)
    for (a in ages) add_aphid(sp, row$plant, -a * minutes_per_day)
  }

  # ---- Predator state --------------------------------------------------------
  has_pred <- !is.null(scenario$predator)
  # next time of each predator event (Inf = not scheduled)
  ev <- c(start = Inf, encounter = Inf, reject_end = Inf, handle_end = Inf, gut_ready = Inf,
          leave = Inf, arrive = Inf, death = Inf, starve_death = Inf, starve_check = Inf,
          pupation = Inf)
  pr <- list(state = "waiting", stage = 1L, plant = NA_integer_, food = 0, gut = 0, t_gut = 0,
             bean_cum = 0, pea_cum = 0, last_meal = NA_real_, last_pea = -Inf, n_bean = 0L,
             hazard_budget = rexp(1), t_hazard = 0, hazard = 0,
             thresholds = p$threshold * exp(rnorm(4, 0, p$threshold_sd)),
             crit_food = rlogis(1, p$critical_food$location, p$critical_food$scale),
             prey_species = NA_integer_, prey_age = NA_real_, attack_time = NA_real_,
             starve_at_attack = NA_real_, activity = NA_real_, in_colony = FALSE,
             n_encounters = 0L, n_failed = 0L, leave_remaining = NA_real_,
             last_encounter = NA_real_, fate = NA_character_, fate_time = NA_real_)
  if (has_pred) ev["start"] <- scenario$predator$start_day * minutes_per_day

  meals <- vector("list", 2000L); n_meals <- 0L
  moves <- list()
  stage_log <- list()

  # ---- Predator: gut ---------------------------------------------------------
  digestion <- p$digestion_rate / minutes_per_day # per minute
  # Gut capacity by instar. The larva attacks only when the gut has room for
  # one more size-matched prey (gut <= C - 1); the gut empties exponentially
  # at rate k throughout, including during handling. With ad lib prey the
  # gut cycles between C - 1 (attack) and a +1 jump at the end of handling
  # (duration h); a steady cycle of length T = 1 / (kills per day) requires
  #   C - 1 = exp(-k (T - h)) / (1 - exp(-k T)).
  capacities <- vapply(seq_along(p$max_intake), function(s) {
    cycle <- minutes_per_day / p$max_intake[[s]]
    h <- min(p$mean_pea_handling, cycle)
    1 + exp(-digestion * (cycle - h)) / (1 - exp(-digestion * cycle))
  }, numeric(1))
  gut_capacity <- function() capacities[pr$stage]
  gut_now <- function() pr$gut * exp(-digestion * (t - pr$t_gut))

  # ---- Predator: mortality and starvation ----------------------------------
  # Background death: the larva carries an Exp(1) budget of cumulative
  # hazard, spent at the current hazard rate; the death event is rescheduled
  # whenever the rate changes.
  hazard_now <- function() {
    st <- if (pr$state == "prepupa") 4L else pr$stage
    h <- exp(p$mortality$log_hazard[c("L1", "L2-3", "L2-3", "L4")[st]])
    if (st == 1L) h <- h * bean_load(pr$bean_cum, pr$pea_cum)^p$mortality$l1_bean_exponent
    unname(h) / minutes_per_day
  }
  spend_hazard <- function() {
    pr$hazard_budget <<- pr$hazard_budget - pr$hazard * (t - pr$t_hazard)
    pr$t_hazard <<- t
  }
  refresh_death <- function() {
    pr$hazard <<- hazard_now()
    ev["death"] <<- t + max(pr$hazard_budget, 0) / pr$hazard
  }

  schedule_starvation <- function() {
    # Lognormal around the median: a well-fed larva redraws this after every
    # meal, so the distribution must put ~no mass at short times.
    days <- p$starvation_median[[pr$stage]] * exp(rnorm(1, 0, p$starvation_sdlog))
    ev["starve_death"] <<- pr$last_meal + days * minutes_per_day
    ev["starve_check"] <<- pr$last_meal + p$starvation_onset
  }

  # Hours since the last meal, clamped to the range of the behavior trials.
  starve_hours <- function() {
    h <- (t - pr$last_meal) / 60
    min(max(h, p$starve_range[1]), p$starve_range[2])
  }

  # ---- Predator: search and capture ------------------------------------------
  # Proportion of time moving after a meal of `species` at hunger `starve_h`
  # (2008 tracking): ~0.4 after pea; lower after bean, more so when hungry.
  activity_after <- function(species, starve_h) plogis(aphid_starve_lp(p$activity, species, starve_h))

  # Area searched per minute (mm^2): speed x activity x detection width,
  # scaled by (body length / L4 length)^2 for earlier instars.
  size_scale <- (p$larva_length / p$larva_length[["L4"]])^2
  search_area <- function() p$speed * pr$activity * p$detection_width * size_scale[[pr$stage]]
  plant_area_mm2 <- p$plant_area * 100

  # Aphid aggregation. Aphids live in colonies whose area grows with their
  # numbers: colony area = colony_min_area + sum over species of (count x
  # adult mass) / colony_packing[species], capped at the plant area (bean
  # aphids pack densely, pea aphids loosely). A larva that has found the
  # colony (hatched beside it, or has eaten on this plant - area-restricted
  # search after a meal) searches only the colony area; on arriving at a
  # plant it searches the whole plant until it eats.
  search_space_mm2 <- function() {
    if (!pr$in_colony) return(plant_area_mm2)
    colony <- p$colony_min_area +
      N[pr$plant, PEA] * m_adult[PEA] / p$colony_packing[["pea"]] +
      N[pr$plant, BEAN] * m_adult[BEAN] / p$colony_packing[["bean"]]
    min(p$plant_area, colony) * 100
  }

  # Encounters: the larva meets aphids on its plant at random, at rate
  # (area searched / search space) x number of aphids. Each encounter is an
  # attack that succeeds with probability capture[species] x size factor.
  update_encounter <- function() {
    if (pr$state != "search") { ev["encounter"] <<- Inf; return(invisible()) }
    rate <- search_area() / search_space_mm2() * (N[pr$plant, PEA] + N[pr$plant, BEAN])
    ev["encounter"] <<- if (rate > 0) t + rexp(1, rate) else Inf
  }

  # Prey size: small larvae cannot subdue large aphids. Size factor =
  # 1 / (1 + (prey length / (prey_size_ratio x larval length))^size_steepness)
  # with aphid length (mm) = aphid_length_coef x mass^(1/3).
  capture_prob <- function(species, age_days) {
    prey_len <- p$aphid_length_coef * mass_of(species, age_days)^(1 / 3)
    size <- 1 / (1 + (prey_len / (p$prey_size_ratio * p$larva_length[[pr$stage]]))^p$size_steepness)
    p$capture[[species]] * size
  }

  # Handling time (min): 2008 log-logistic model for size-matched prey,
  # scaled by (age / reference age)^b (2005 trials), and for bean aphids
  # shortened with experience: multiplier m + (1 - m) exp(-bean meals /
  # learn_meals).
  sample_handling <- function(species, age_days, starve_h) {
    lp <- aphid_starve_lp(p$handling$coef, species, starve_h) +
      p$handling$age_exponent * log(age_days / p$ref_age[[species_names[species]]])
    mult <- if (species == BEAN) {
      m <- p$bean_handling_experienced
      m + (1 - m) * exp(-pr$n_bean / p$learn_meals)
    } else 1
    sample_aft(lp, p$handling$sigma, "loglogistic") * mult
  }

  # ---- Predator: leaving a plant -------------------------------------------
  # Post-handling times were measured for leaving one leaf (trial_leaf_area).
  # Walking off a leaf on a plant leads to another leaf of the same plant, so
  # leaving the plant means giving up on ~n_leaves = plant_area /
  # trial_leaf_area x leave_scaling leaves in turn: the leave time is the sum
  # of n_leaves leaf-level draws. (Multiplying one draw by n_leaves would
  # keep the fitted early-departure spike - Weibull shape < 1 - and larvae
  # would walk off dense colonies minutes after a meal.) With no meal on this
  # plant, the pea model at current hunger applies.
  n_leaves <- max(1L, round(p$plant_area / p$trial_leaf_area * p$leave_scaling))

  # Leaf-level post-handling (leave) times, Weibull AFT; vectorized.
  sample_leave <- function(species, starve_h) {
    lp <- aphid_starve_lp(p$post_handling$coef, species, starve_h)
    sample_aft(lp, p$post_handling$sigma, "weibull")
  }

  schedule_leave <- function(species, starve_h) {
    if (vial_mode) return(invisible())
    pr$leave_remaining <<- NA_real_
    ev["leave"] <<- t + sum(sample_leave(rep(species, n_leaves), starve_h))
  }

  # ---- Predator: state transitions -------------------------------------------
  # After a meal, rejection or arrival: wait for gut room, or search.
  go_hungry_or_satiated <- function() {
    g <- gut_now()
    cap_g <- gut_capacity()
    if (g > cap_g - 1) {
      pr$state <<- "satiated"
      ev["gut_ready"] <<- t + log(g / (cap_g - 1)) / digestion
      ev["encounter"] <<- Inf
      # The leave clock measures unsuccessful searching, so it pauses while
      # the larva is satiated (2008 leave times come from larvae searching a
      # leaf with no prey left).
      if (is.finite(ev[["leave"]])) {
        pr$leave_remaining <<- ev[["leave"]] - t
        ev["leave"] <<- Inf
      }
    } else {
      pr$state <<- "search"
      ev["gut_ready"] <<- Inf
      update_encounter()
    }
  }

  finish <- function(fate) {
    spend_hazard()
    pr$fate <<- fate; pr$fate_time <<- t; pr$state <<- "done"
    ev[] <<- Inf
  }

  start_prepupa <- function(duration) {
    pr$state <<- "prepupa"
    ev[c("encounter", "gut_ready", "leave", "starve_death", "starve_check")] <<- Inf
    ev["pupation"] <<- t + duration
  }

  log_stage <- function() stage_log[[length(stage_log) + 1L]] <<- list(time = t, stage = pr$stage)

  # In vial mode a killed aphid is replaced at once by a size-matched one.
  replenish <- function(species) {
    add_aphid(species, 1L, t - p$ref_age[[species_names[species]]] * minutes_per_day)
  }

  # Value of a meal in development units: bean aphids are worth less, and
  # less still without a recent pea meal.
  food_value <- function(species) {
    if (species == PEA) 1
    else if (t - pr$last_pea <= p$pea_memory) p$v_bean_mixed
    else p$v_bean_alone
  }

  # Molting and pupation once accumulated food passes the instar threshold.
  # Returns FALSE if the larva died attempting to pupate.
  develop <- function() {
    while (pr$state != "prepupa" && pr$food >= pr$thresholds[pr$stage]) {
      if (pr$stage < 4L) {
        pr$food <<- 0
        pr$stage <<- pr$stage + 1L
        log_stage()
      } else {
        if (runif(1) > pupation_success(pr$bean_cum, pr$pea_cum, p$pupation_gamma)) {
          finish("failed_pupation")
          return(FALSE)
        }
        start_prepupa(p$prepupa_days * minutes_per_day)
      }
    }
    TRUE
  }

  # ---- Predator event handlers ---------------------------------------------
  on_start <- function() {
    pr$plant <<- scenario$predator$plant
    pr$stage <<- match(scenario$predator$stage %||% "L1", larval_stages)
    pr$last_meal <<- t
    pr$t_gut <<- t; pr$t_hazard <<- t
    pr$activity <<- activity_after(PEA, p$starve_range[1])
    pr$in_colony <<- TRUE # eggs are laid next to aphid colonies
    log_stage()
    refresh_death()
    schedule_starvation()
    ev["start"] <<- Inf
    # No leaving before the first meal: hatchlings emerge beside a colony and
    # the leave rule (fitted to fed L4s) starts with the first meal.
    go_hungry_or_satiated()
  }

  on_encounter <- function() {
    pr$n_encounters <<- pr$n_encounters + 1L
    pr$last_encounter <<- t
    pr$attack_time <<- t
    pool <- which(a_alive & a_plant == pr$plant)
    i <- if (length(pool) == 1L) pool else sample(pool, 1L)
    sp <- a_sp[i]
    age <- aphid_age_days(i)
    # failed attack (aphid escapes) or rejection after capture: aphid survives
    if (runif(1) > capture_prob(sp, age) ||
        runif(1) < plogis(p$rejection[1] + p$rejection[2] * age)) {
      pr$n_failed <<- pr$n_failed + 1L
      pr$state <<- "rejecting"
      ev["encounter"] <<- Inf
      ev["reject_end"] <<- t + p$rejection_time
      return(invisible())
    }
    pr$prey_species <<- sp
    pr$prey_age <<- age
    pr$starve_at_attack <<- starve_hours()
    remove_aphid(i)
    if (vial_mode) replenish(sp)
    pr$state <<- "handling"
    ev["encounter"] <<- Inf
    ev["leave"] <<- Inf # leave clock restarts after the meal
    ev[c("starve_death", "starve_check")] <<- Inf # eating; restarted after the meal
    ev["handle_end"] <<- t + sample_handling(sp, max(age, 1 / 24), pr$starve_at_attack)
  }

  on_reject_end <- function() {
    ev["reject_end"] <<- Inf
    pr$state <<- "search"
    go_hungry_or_satiated()
  }

  on_handle_end <- function() {
    ev["handle_end"] <<- Inf
    sp <- pr$prey_species
    m_units <- mass_of(sp, pr$prey_age) / unit_mass
    v <- food_value(sp)
    spend_hazard()
    pr$gut <<- gut_now() + m_units; pr$t_gut <<- t
    pr$food <<- pr$food + v * m_units
    if (sp == PEA) {
      pr$pea_cum <<- pr$pea_cum + m_units; pr$last_pea <<- t
    } else {
      pr$bean_cum <<- pr$bean_cum + m_units; pr$n_bean <<- pr$n_bean + 1L
    }
    pr$last_meal <<- t
    pr$in_colony <<- TRUE
    n_meals <<- n_meals + 1L
    if (n_meals > length(meals)) meals <<- c(meals, vector("list", length(meals)))
    meals[[n_meals]] <<- list(time = t, plant = pr$plant, species = sp, aphid_age = pr$prey_age,
                              handling = t - pr$attack_time, units = v * m_units, stage = pr$stage)

    if (!develop()) return(invisible())
    refresh_death()
    if (pr$state == "prepupa") return(invisible())
    schedule_starvation()
    pr$activity <<- activity_after(sp, pr$starve_at_attack)
    schedule_leave(sp, pr$starve_at_attack)
    go_hungry_or_satiated()
  }

  on_gut_ready <- function() {
    ev["gut_ready"] <<- Inf
    pr$state <<- "search"
    if (!is.na(pr$leave_remaining)) {
      ev["leave"] <<- t + pr$leave_remaining
      pr$leave_remaining <<- NA_real_
    }
    update_encounter()
  }

  on_leave <- function() {
    if (pr$state %in% c("handling", "rejecting")) { ev["leave"] <<- Inf; return(invisible()) }
    moves[[length(moves) + 1L]] <<- list(
      time = t, from = pr$plant, stage = pr$stage, activity = pr$activity,
      in_colony = pr$in_colony, since_meal = t - pr$last_meal,
      since_encounter = t - pr$last_encounter,
      aphids_here = N[pr$plant, PEA] + N[pr$plant, BEAN]
    )
    pr$state <<- "travel"
    ev[c("leave", "encounter", "gut_ready")] <<- Inf
    ev["arrive"] <<- t + p$plant_path / (p$speed * pr$activity)
  }

  on_arrive <- function() {
    ev["arrive"] <<- Inf
    pr$plant <<- 3L - pr$plant
    pr$in_colony <<- FALSE
    schedule_leave(PEA, starve_hours())
    go_hungry_or_satiated()
  }

  # A day without food: an L4 that has eaten enough commits to pupation.
  on_starve_check <- function() {
    ev["starve_check"] <<- Inf
    if (pr$stage == 4L && pr$food >= pr$crit_food) {
      delay <- (p$pupation_delay[[1]] + p$pupation_delay[[2]] * pr$food) * minutes_per_day
      start_prepupa(max(delay - p$starvation_onset, 0))
    }
  }

  run_predator_event <- function(e) {
    switch(e,
      start = on_start(),
      encounter = on_encounter(),
      reject_end = on_reject_end(),
      handle_end = on_handle_end(),
      gut_ready = on_gut_ready(),
      leave = on_leave(),
      arrive = on_arrive(),
      death = finish("died"),
      starve_death = finish("starved"),
      starve_check = on_starve_check(),
      pupation = finish("pupated")
    )
  }

  run_aphid_event <- function(i) {
    sp <- a_sp[i]
    plant <- a_plant[i]
    if (a_type[i] == BIRTH) {
      add_aphid(sp, plant, t, alate = newborn_is_alate(plant, sp))
      schedule_aphid(i)
    } else {
      if (a_type[i] == EMIGRATION) emigrants[plant, sp] <<- emigrants[plant, sp] + 1L
      remove_aphid(i)
    }
    if (pr$state == "search" && pr$plant == plant) update_encounter()
  }

  # ---- Daily census ----------------------------------------------------------
  # Columns are as.vector(N): plant1_pea, plant2_pea, plant1_bean, plant2_bean.
  n_days <- floor(scenario$run_days)
  census <- matrix(NA_integer_, n_days + 1L, 4L)
  census[1, ] <- as.vector(N)
  census_row <- 1L
  next_census <- minutes_per_day
  emigrants <- matrix(0L, 2, 2, dimnames = list(plant = 1:2, species = species_names))

  # Record the census, then update crowding mortality from the lagged census.
  take_census <- function() {
    census_row <<- census_row + 1L
    census[census_row, ] <<- as.vector(N)
    next_census <<- next_census + minutes_per_day
    if (!vial_mode) {
      lagged <- census[max(1L, census_row - p$density_lag), ]
      update_dd_rates(matrix(lagged, 2, 2))
      redraw_crowding_deaths()
    }
  }

  # ---- Main loop --------------------------------------------------------------
  repeat {
    ia <- which.min(a_next)
    ta <- a_next[ia]
    tp <- min(ev)
    if (next_census <= min(ta, tp) && next_census <= run_end) {
      t <- next_census
      take_census()
      if (census_row > n_days) break
    } else if (min(ta, tp) >= run_end) {
      break
    } else if (ta <= tp) {
      t <- ta
      run_aphid_event(ia)
    } else {
      t <- tp
      run_predator_event(names(ev)[which.min(ev)])
    }
  }

  # ---- Output -------------------------------------------------------------------
  list(
    census = tidy_census(census[seq_len(census_row), , drop = FALSE]),
    emigrants = emigrants,
    meals = if (n_meals > 0L) {
      dplyr::bind_rows(meals[seq_len(n_meals)]) |>
        dplyr::mutate(species = species_names[species])
    } else {
      tibble::tibble(time = numeric(), plant = integer(), species = character(),
                     aphid_age = numeric(), handling = numeric(), units = numeric(),
                     stage = integer())
    },
    predator = list(
      fate = if (!has_pred) NA_character_ else if (is.na(pr$fate)) "alive" else pr$fate,
      fate_day = pr$fate_time / minutes_per_day,
      n_encounters = pr$n_encounters, n_failed = pr$n_failed, final_plant = pr$plant,
      moves = dplyr::bind_rows(moves),
      stages = if (length(stage_log)) {
        dplyr::bind_rows(stage_log) |>
          dplyr::mutate(day = time / minutes_per_day, stage = larval_stages[stage])
      } else {
        tibble::tibble(time = numeric(), stage = character(), day = numeric())
      }
    )
  )
}

#' Daily census matrix (columns plant1_pea, plant2_pea, plant1_bean,
#' plant2_bean) to a long tibble of day, plant, species, n.
tidy_census <- function(census) {
  colnames(census) <- c("plant1_pea", "plant2_pea", "plant1_bean", "plant2_bean")
  tibble::as_tibble(census) |>
    dplyr::mutate(day = dplyr::row_number() - 1L) |>
    tidyr::pivot_longer(-day, names_to = c("plant", "species"), names_sep = "_", values_to = "n") |>
    dplyr::mutate(plant = as.integer(sub("plant", "", plant)))
}
