# Two-plant discrete-event simulation ---------------------------------------
#
# One coccinellid larva foraging on two plants that hold pea and bean
# aphids. Continuous time in minutes. Aphid events (births, deaths) are kept
# in vectors and the next one is found with which.min(); the predator is a
# state machine with its own schedule. See docs/model-design.md.
#
# Species codes: 1 = pea, 2 = bean. Plants: 1, 2.

species_names <- c("pea", "bean")
minutes_per_day <- 1440

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
#' @return list(census, meals, predator)
simulate_two_plants <- function(params, scenario, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params
  vial <- scenario$vial
  run_end <- scenario$run_days * minutes_per_day

  # ---- Aphid store ---------------------------------------------------------
  # Parallel vectors indexed by slot. Slots of dead aphids are reused (free
  # stack), so vectors scale with the living population, not with every
  # aphid ever born. Unused slots have a_next = Inf, so which.min(a_next)
  # over the whole vector finds the next aphid event without copying.
  cap <- 1024L
  a_sp <- integer(cap); a_plant <- integer(cap); a_birth <- numeric(cap)
  a_fec <- numeric(cap); a_next <- rep(Inf, cap); a_type <- integer(cap) # 1 birth, 2 death
  a_alive <- logical(cap)
  n_used <- 0L # high-water mark of slots ever used
  free <- integer(0); n_free <- 0L
  N <- matrix(0L, 2, 2) # plant x species
  t <- 0

  grow <- function() {
    new_cap <- cap * 2L
    a_sp <<- c(a_sp, integer(cap)); a_plant <<- c(a_plant, integer(cap))
    a_birth <<- c(a_birth, numeric(cap)); a_fec <<- c(a_fec, numeric(cap))
    a_next <<- c(a_next, rep(Inf, cap)); a_type <<- c(a_type, integer(cap))
    a_alive <<- c(a_alive, logical(cap))
    cap <<- new_cap
  }

  # Fast per-species lookups for event sampling. The birth CIF has knots at
  # whole days (0, 1440, 2880, ... minutes), so evaluating it is O(1)
  # arithmetic instead of approx().
  death_shape <- vapply(p$aphid[species_names], `[[`, numeric(1), "death_shape")
  death_scale <- vapply(p$aphid[species_names], `[[`, numeric(1), "death_scale")
  frailty_shape <- vapply(p$aphid[species_names], `[[`, numeric(1), "frailty_shape")
  cif_value <- lapply(p$aphid[species_names], `[[`, "cif_value")
  cif_max <- vapply(cif_value, max, numeric(1))
  cif_last_age <- vapply(p$aphid[species_names], \(a) max(a$cif_age), numeric(1))
  stopifnot(all(vapply(p$aphid[species_names], \(a) isTRUE(all.equal(diff(a$cif_age), rep(minutes_per_day, length(a$cif_age) - 1))), logical(1))))

  cif_eval <- function(sp, age) {
    if (age >= cif_last_age[sp]) return(cif_max[sp])
    x <- age / minutes_per_day
    k <- floor(x)
    v <- cif_value[[sp]]
    v[k + 1L] + (v[k + 2L] - v[k + 1L]) * (x - k)
  }

  # Time to next birth by inversion of the CIF (see sample_birth_time()).
  birth_wait <- function(sp, age, fec) {
    target <- cif_eval(sp, age) + rexp(1) / fec
    if (target >= cif_max[sp]) return(Inf)
    v <- cif_value[[sp]]
    k <- findInterval(target, v, checkSorted = FALSE, checkNA = FALSE) # largest knot <= target
    (k - 1 + (target - v[k]) / (v[k + 1L] - v[k])) * minutes_per_day - age
  }

  schedule_aphid <- function(i) {
    if (!is.null(vial)) { a_next[i] <<- Inf; return(invisible()) }
    sp <- a_sp[i]
    age <- t - a_birth[i]
    sh <- death_shape[[sp]]; sc <- death_scale[[sp]]
    death <- sc * ((age / sc)^sh + rexp(1))^(1 / sh) - age
    birth <- birth_wait(sp, age, a_fec[i])
    if (birth < death) { a_next[i] <<- t + birth; a_type[i] <<- 1L } else { a_next[i] <<- t + death; a_type[i] <<- 2L }
  }

  add_aphid <- function(species, plant, birth_time) {
    if (n_free > 0L) {
      i <- free[n_free]
      n_free <<- n_free - 1L
    } else {
      if (n_used == cap) grow()
      n_used <<- n_used + 1L
      i <- n_used
    }
    a_sp[i] <<- species; a_plant[i] <<- plant; a_birth[i] <<- birth_time
    a_fec[i] <<- rgamma(1, shape = frailty_shape[[species]], rate = frailty_shape[[species]])
    a_alive[i] <<- TRUE
    if (is.finite(capacity)) update_dbar(plant)
    N[plant, species] <<- N[plant, species] + 1L
    schedule_aphid(i)
    i
  }

  remove_aphid <- function(i) {
    if (is.finite(capacity)) update_dbar(a_plant[i])
    a_alive[i] <<- FALSE
    a_next[i] <<- Inf
    N[a_plant[i], a_sp[i]] <<- N[a_plant[i], a_sp[i]] - 1L
    n_free <<- n_free + 1L
    free[n_free] <<- i
  }

  # In vial mode prey were replaced daily with size-matched aphids, so prey
  # age is held at the size-matched reference age.
  aphid_age_days <- function(i) {
    if (!is.null(vial)) return(p$ref_age[[species_names[a_sp[i]]]])
    (t - a_birth[i]) / minutes_per_day
  }

  mass_of <- function(species, age_days) {
    m <- p$mass[p$mass$aphid == species_names[species], ]
    min(m$mass_neonate * exp(m$growth_rate * age_days), m$mass_adult)
  }

  # Density dependence: births are thinned with probability
  # g = max(0, 1 - Dbar / K). D is aphid density on the plant in adult-mass
  # equivalents (mg): each aphid counts its species' adult mass, both species
  # combined (counting newborns at full weight avoids a growth lag). Dbar is an
  # exponentially weighted average of D over the past ~density_lag days, so
  # crowding acts with a delay and losses (e.g. to the predator) are not
  # replaced instantly. D is constant between events, so Dbar is updated
  # exactly: Dbar <- D + (Dbar - D) exp(-dt / lag). Exact thinning because
  # g <= 1 and the unthinned birth process is each aphid's own CIF.
  m_adult <- unname(setNames(p$mass$mass_adult, as.character(p$mass$aphid))[species_names])
  capacity <- p$aphid_capacity
  lag <- p$density_lag * minutes_per_day
  dbar <- c(0, 0); t_dbar <- c(0, 0)
  density_now <- function(plant) N[plant, 1L] * m_adult[1L] + N[plant, 2L] * m_adult[2L]
  update_dbar <- function(plant) {
    if (lag <= 0) { dbar[plant] <<- density_now(plant); return(invisible()) }
    dbar[plant] <<- density_now(plant) + (dbar[plant] - density_now(plant)) * exp(-(t - t_dbar[plant]) / lag)
    t_dbar[plant] <<- t
  }
  birth_succeeds <- function(plant) {
    if (!is.finite(capacity)) return(TRUE)
    update_dbar(plant)
    runif(1) < 1 - dbar[plant] / capacity # FALSE whenever Dbar >= capacity
  }

  # initial aphids
  for (r in seq_len(nrow(scenario$aphids))) {
    row <- scenario$aphids[r, ]
    sp <- match(row$species, species_names)
    for (k in seq_len(row$n)) add_aphid(sp, row$plant, -row$age * minutes_per_day)
  }

  # ---- Predator -------------------------------------------------------------
  has_pred <- !is.null(scenario$predator)
  ev <- c(start = Inf, encounter = Inf, reject_end = Inf, handle_end = Inf, gut_ready = Inf,
          leave = Inf, arrive = Inf, death = Inf, starve_death = Inf, starve_check = Inf,
          pupation = Inf)
  pr <- list(state = "waiting", stage = 1L, plant = NA_integer_, food = 0, gut = 0, t_gut = 0,
             bean_cum = 0, pea_cum = 0, last_meal = NA_real_, last_pea = -Inf, n_bean = 0L,
             hazard_budget = rexp(1), t_hazard = 0, hazard = 0,
             thresholds = p$threshold * exp(rnorm(4, 0, p$threshold_sd)),
             crit_food = rlogis(1, p$critical_food$location, p$critical_food$scale),
             prey = NA_integer_, prey_species = NA_integer_, prey_age = NA_real_,
             starve_at_attack = NA_real_, activity = NA_real_, in_colony = FALSE,
             n_encounters = 0L, n_failed = 0L, leave_remaining = NA_real_,
             leave_species = 1L, leave_starve = 2, last_encounter = NA_real_,
             fate = NA_character_, fate_time = NA_real_)
  if (has_pred) ev["start"] <- scenario$predator$start_day * minutes_per_day

  meals <- vector("list", 2000L); n_meals <- 0L
  moves <- list()
  stage_log <- list()

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

  hazard_now <- function() {
    st <- if (pr$state == "prepupa") 4L else pr$stage
    h <- exp(p$mortality$log_hazard[c("L1", "L2-3", "L2-3", "L4")[st]])
    if (st == 1L) h <- h * (1 + pr$bean_cum / (1 + pr$pea_cum))^p$mortality$l1_bean_exponent
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

  starve_hours <- function() {
    h <- (t - pr$last_meal) / 60
    min(max(h, p$starve_range[1]), p$starve_range[2])
  }

  # Proportion of time moving after a meal of `species` at hunger `starve_h`
  # (2008 tracking): ~0.4 after pea; lower after bean, more so when hungry.
  activity_after <- function(species, starve_h) {
    b <- p$activity
    is_bean <- species == 2L
    plogis(b[["(Intercept)"]] + is_bean * b[["aphidbean"]] + starve_h * b[["starve"]] +
             is_bean * starve_h * b[["aphidbean:starve"]])
  }

  # Area searched per minute (mm^2): speed x activity x detection width,
  # scaled by (body length / L4 length)^2 for earlier instars.
  size_scale <- (p$larva_length / p$larva_length[["L4"]])^2
  search_area <- function() p$speed * pr$activity * p$detection_width * size_scale[[pr$stage]]
  plant_area_mm2 <- p$plant_area * 100

  # Aphid aggregation. Aphids live in colonies whose area grows with aphid
  # density (adult-mass equivalents, both species): colony area =
  # colony_min_area + density / colony_density, capped at the plant area. A
  # larva that has found the colony (hatched beside it, or has eaten on this
  # plant - area-restricted search after a meal) searches only the colony
  # area; on arriving at a plant it searches the whole plant until it eats.
  search_space_mm2 <- function() {
    if (!pr$in_colony) return(plant_area_mm2)
    density <- N[pr$plant, 1L] * m_adult[1L] + N[pr$plant, 2L] * m_adult[2L]
    min(p$plant_area, p$colony_min_area + density / p$colony_density) * 100
  }

  # Encounters: the larva meets aphids on its plant at random, at rate
  # (area searched / search space) x number of aphids. Each encounter is an
  # attack that succeeds with probability capture[species] x size factor.
  update_encounter <- function() {
    if (pr$state != "search") { ev["encounter"] <<- Inf; return(invisible()) }
    rate <- search_area() / search_space_mm2() * (N[pr$plant, 1L] + N[pr$plant, 2L])
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

  # Leaving a plant. Post-handling times were measured for leaving one leaf
  # (trial_leaf_area). Walking off a leaf on a plant leads to another leaf of
  # the same plant, so leaving the plant means giving up on ~n_leaves =
  # plant_area / trial_leaf_area x leave_scaling leaves in turn: the leave time
  # is the sum of n_leaves leaf-level draws. (Multiplying one draw by n_leaves
  # would keep the fitted early-departure spike - Weibull shape < 1 - and
  # larvae would walk off dense colonies minutes after a meal.) With no meal
  # on this plant, the pea model at current hunger applies.
  n_leaves <- max(1L, round(p$plant_area / p$trial_leaf_area * p$leave_scaling))
  schedule_leave <- function(species, starve_h) {
    if (!is.null(vial)) return(invisible())
    pr$leave_remaining <<- NA_real_
    ev["leave"] <<- t + sum(sample_leave(rep(species, n_leaves), starve_h))
  }

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

  sample_handling <- function(species, age_days, starve_h) {
    b <- p$handling$coef
    is_bean <- species == 2L
    lp <- b[["(Intercept)"]] + is_bean * b[["aphidbean"]] + starve_h * b[["starve"]] +
      is_bean * starve_h * b[["aphidbean:starve"]] +
      p$handling$age_exponent * log(age_days / p$ref_age[[species_names[species]]])
    mult <- if (is_bean) {
      m <- p$bean_handling_experienced
      m + (1 - m) * exp(-pr$n_bean / p$learn_meals)
    } else 1
    exp(lp + p$handling$sigma * qlogis(runif(1))) * mult
  }

  sample_leave <- function(species, starve_h) {
    b <- p$post_handling$coef
    is_bean <- species == 2L
    lp <- b[["(Intercept)"]] + is_bean * b[["aphidbean"]] + starve_h * b[["starve"]] +
      is_bean * starve_h * b[["aphidbean:starve"]]
    exp(lp + p$post_handling$sigma * log(-log(runif(length(species)))))
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

  # ---- Predator event handlers ---------------------------------------------
  on_start <- function() {
    pr$plant <<- scenario$predator$plant
    pr$stage <<- match(scenario$predator$stage %||% "L1", larval_stages)
    pr$last_meal <<- t
    pr$t_gut <<- t; pr$t_hazard <<- t
    pr$activity <<- activity_after(1L, p$starve_range[1])
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
    if (!is.null(vial)) replenish(sp)
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
    v <- if (sp == 1L) 1 else if (t - pr$last_pea <= p$pea_memory) p$v_bean_mixed else p$v_bean_alone
    spend_hazard()
    pr$gut <<- gut_now() + m_units; pr$t_gut <<- t
    pr$food <<- pr$food + v * m_units
    if (sp == 1L) { pr$pea_cum <<- pr$pea_cum + m_units; pr$last_pea <<- t }
    else { pr$bean_cum <<- pr$bean_cum + m_units; pr$n_bean <<- pr$n_bean + 1L }
    pr$last_meal <<- t
    pr$in_colony <<- TRUE
    n_meals <<- n_meals + 1L
    if (n_meals > length(meals)) meals <<- c(meals, vector("list", length(meals)))
    meals[[n_meals]] <<- list(time = t, plant = pr$plant, species = sp, aphid_age = pr$prey_age,
                              handling = t - ev_attack_time, units = v * m_units, stage = pr$stage)

    # development
    while (pr$state != "prepupa" && pr$food >= pr$thresholds[pr$stage]) {
      if (pr$stage < 4L) {
        pr$food <<- 0
        pr$stage <<- pr$stage + 1L
        log_stage()
      } else {
        if (runif(1) > (1 + pr$bean_cum / (1 + pr$pea_cum))^(-p$pupation_gamma)) {
          finish("failed_pupation")
          return(invisible())
        }
        start_prepupa(p$prepupa_days * minutes_per_day)
      }
    }
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
      aphids_here = N[pr$plant, 1L] + N[pr$plant, 2L]
    )
    pr$state <<- "travel"
    ev[c("leave", "encounter", "gut_ready")] <<- Inf
    ev["arrive"] <<- t + p$plant_path / (p$speed * pr$activity)
  }

  on_arrive <- function() {
    ev["arrive"] <<- Inf
    pr$plant <<- 3L - pr$plant
    pr$in_colony <<- FALSE
    schedule_leave(1L, starve_hours())
    go_hungry_or_satiated()
  }

  on_starve_check <- function() {
    ev["starve_check"] <<- Inf
    if (pr$stage == 4L && pr$food >= pr$crit_food) {
      delay <- (p$pupation_delay[[1]] + p$pupation_delay[[2]] * pr$food) * minutes_per_day
      start_prepupa(max(delay - p$starvation_onset, 0))
    }
  }

  replenish <- function(species) {
    add_aphid(species, 1L, t - p$ref_age[[species_names[species]]] * minutes_per_day)
  }

  # ---- Census -----------------------------------------------------------------
  n_days <- floor(scenario$run_days)
  census <- matrix(NA_integer_, n_days + 1L, 4L)
  census[1, ] <- as.vector(N)
  next_census <- minutes_per_day
  census_row <- 1L

  # ---- Main loop --------------------------------------------------------------
  ev_attack_time <- NA_real_
  repeat {
    ia <- which.min(a_next)
    ta <- a_next[ia]
    tp <- min(ev)
    t_next <- min(ta, tp, next_census, run_end)

    if (next_census <= min(ta, tp) && next_census <= run_end) {
      t <- next_census
      census_row <- census_row + 1L
      census[census_row, ] <- as.vector(N)
      next_census <- next_census + minutes_per_day
      if (census_row > n_days) break
      next
    }
    if (min(ta, tp) >= run_end) break

    if (ta <= tp) {
      t <- ta
      plant <- a_plant[ia]
      if (a_type[ia] == 1L) {
        if (birth_succeeds(plant)) add_aphid(a_sp[ia], plant, t)
        schedule_aphid(ia)
      } else {
        remove_aphid(ia)
      }
      if (pr$state == "search" && pr$plant == plant) update_encounter()
    } else {
      t <- tp
      e <- names(ev)[which.min(ev)]
      if (e == "encounter") ev_attack_time <- t
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
  }

  census_df <- tibble::as_tibble(census[seq_len(census_row), , drop = FALSE], .name_repair = "minimal")
  names(census_df) <- c("plant1_pea", "plant2_pea", "plant1_bean", "plant2_bean")
  census_df <- census_df |>
    dplyr::mutate(day = dplyr::row_number() - 1L) |>
    tidyr::pivot_longer(-day, names_to = c("plant", "species"), names_sep = "_", values_to = "n") |>
    dplyr::mutate(plant = as.integer(sub("plant", "", plant)))

  list(
    census = census_df,
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
