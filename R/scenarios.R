# Scenario helpers ----------------------------------------------------------

#' Vial scenario reproducing the diet experiment: one larva from hatching,
#' ad lib size-matched prey replaced as eaten.
vial_scenario <- function(diet, params, n_prey = 30, run_days = 40) {
  share <- switch(diet, P = c(pea = 1, bean = 0), B = c(pea = 0, bean = 1), M = c(pea = 0.5, bean = 0.5))
  aphids <- tibble::tibble(
    plant = 1L,
    species = c("pea", "bean"),
    n = round(n_prey * share),
    age = c(params$ref_age[["pea"]], params$ref_age[["bean"]])
  ) |>
    dplyr::filter(n > 0)
  list(aphids = aphids, predator = list(start_day = 0, plant = 1L, stage = "L1"),
       run_days = run_days, vial = list())
}

#' Vial geometry: a small arena where prey are found quickly and cannot drop
#' off a plant to escape.
vial_params <- function(params) {
  params$plant_area <- 50 # vial surface; prey are searched over the whole vial
  params$colony_min_area <- 50
  params$capture[] <- 1
  params
}

#' Two-plant scenario.
#'
#' @param plant1,plant2 named vectors of founding aphids per species, e.g.
#'   c(pea = 10, bean = 0). With stable-age founders, 10 aphids (~1 adult)
#'   give ~30-45 aphids by day 3, like 2 adult founders. Small founding colonies leave the aphids room to
#'   grow before reaching capacity.
#' @param predator TRUE to add one L1 larva on `predator_plant` at
#'   `predator_day`. Default day 3: the colony (2 founders) has ~20 aphids, a
#'   hatchling can survive, and one larva can still depress or deplete it
#'   (from day 4-5 the colony outgrows a single larva).
#' @param founders "stable" (default): each colony's founders have mixed ages
#'   drawn from the species' stable age distribution, so colonies don't start
#'   as one synchronized cohort. Ages are drawn inside simulate_two_plants()
#'   with the run's seed (common random numbers across parameter sets).
#'   "adults": all founders are adults of `adult_age` days.
two_plant_scenario <- function(plant1, plant2, predator = TRUE, predator_plant = 1L,
                               predator_day = 3, run_days = 35, adult_age = 8,
                               founders = c("stable", "adults")) {
  founders <- match.arg(founders)
  aphids <- dplyr::bind_rows(
    tibble::tibble(plant = 1L, species = names(plant1), n = unname(plant1)),
    tibble::tibble(plant = 2L, species = names(plant2), n = unname(plant2))
  ) |>
    dplyr::filter(n > 0) |>
    dplyr::mutate(age = if (founders == "adults") adult_age else NA_real_)
  list(
    aphids = aphids,
    predator = if (predator) list(start_day = predator_day, plant = predator_plant, stage = "L1"),
    run_days = run_days
  )
}

#' Summarize one vial run: fate and per-instar duration and kills.
summarize_vial_run <- function(run, params) {
  st <- run$predator$stages
  end_of_last <- if (identical(run$predator$fate, "pupated")) {
    (run$predator$fate_day - params$prepupa_days) * minutes_per_day
  } else NA_real_
  ends <- c(st$time[-1], end_of_last)
  kills <- run$meals |> dplyr::count(stage)
  tibble::tibble(
    fate = run$predator$fate,
    stage = st$stage,
    days = (ends - st$time) / minutes_per_day,
    kills = kills$n[match(seq_along(st$stage), kills$stage)]
  )
}

#' Stable age distribution of an aphid species (days), from the fitted
#' life table: survival l(x) from the Weibull lifespan, daily fecundity m(x)
#' from the birth CIF; lambda solves Euler-Lotka sum lambda^-x l(x) m(x) = 1
#' and c(x) ~ lambda^-x l(x). Returns a tibble of age (days) and proportion.
stable_age_distribution <- function(params, species) {
  a <- params$aphid[[species]]
  age <- seq(0, max(a$cif_age) / minutes_per_day) # whole days
  surv <- pweibull(age * minutes_per_day, a$death_shape, a$death_scale, lower.tail = FALSE)
  fec <- c(diff(a$cif_value), 0) # offspring in (x, x + 1]
  lotka <- \(lambda) sum(lambda^-(age + 1) * surv * fec) - 1
  lambda <- uniroot(lotka, c(1.0001, 5))$root
  c_x <- lambda^-age * surv
  tibble::tibble(age = age, prop = c_x / sum(c_x), lambda = lambda)
}

#' Draw founder ages (days) from the stable age distribution.
sample_founder_ages <- function(n, params, species) {
  sad <- stable_age_distribution(params, species)
  sample(sad$age, n, replace = TRUE, prob = sad$prop) + runif(n)
}
