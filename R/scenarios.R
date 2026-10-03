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
#' @param plant1,plant2 named vectors of initial adults per species, e.g.
#'   c(pea = 2, bean = 0). Small founding colonies leave the aphids room to
#'   grow on the large default plant before reaching capacity.
#' @param predator TRUE to add one L1 larva on `predator_plant` at
#'   `predator_day`. Default day 3: the colony (2 founders) has ~20 aphids, a
#'   hatchling can survive, and one larva can still depress or deplete it
#'   (from day 4-5 the colony outgrows a single larva).
#' @param adult_age age (days) of the initial adults
two_plant_scenario <- function(plant1, plant2, predator = TRUE, predator_plant = 1L,
                               predator_day = 3, run_days = 35, adult_age = 8) {
  aphids <- dplyr::bind_rows(
    tibble::tibble(plant = 1L, species = names(plant1), n = unname(plant1)),
    tibble::tibble(plant = 2L, species = names(plant2), n = unname(plant2))
  ) |>
    dplyr::filter(n > 0) |>
    dplyr::mutate(age = adult_age)
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
