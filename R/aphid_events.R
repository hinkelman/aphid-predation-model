# Aphid event-time samplers -----------------------------------------------
#
# All times and ages are in minutes. Each function is vectorized over `age`
# and returns the time until the event (not the absolute event time).

#' Time until death given survival to `age` (conditional Weibull).
#'
#' If T ~ Weibull(shape, scale) then, conditional on T > age,
#' T = scale * ((age / scale)^shape + E)^(1 / shape) with E ~ Exp(1).
sample_death_time <- function(age, shape, scale) {
  e <- rexp(length(age))
  scale * ((age / scale)^shape + e)^(1 / shape) - age
}

#' Evaluate a piecewise-linear cumulative intensity function at `age`.
#' Flat beyond the last knot.
cif_at <- function(age, cif_age, cif_value) {
  approx(cif_age, cif_value, xout = age, rule = 2)$y
}

#' Draw individual fecundity multipliers: Gamma with mean 1, variance 1/shape.
sample_fecundity <- function(n, frailty_shape) {
  rgamma(n, shape = frailty_shape, rate = frailty_shape)
}

#' Time until next birth for an aphid of age `age` (inversion of the CIF).
#'
#' The aphid's birth intensity is `fecundity` times the population CIF's
#' intensity, so with E ~ Exp(1) the next birth is at the age where the CIF
#' has risen by E / fecundity. If that would exceed the maximum of the CIF,
#' the aphid has no more births and `Inf` is returned.
sample_birth_time <- function(age, cif_age, cif_value, fecundity = 1) {
  target <- cif_at(age, cif_age, cif_value) + rexp(length(age)) / fecundity
  out <- rep(Inf, length(age))
  ok <- target < max(cif_value)

  # Largest knot with value <= target; because `target` is strictly below the
  # next knot's value, that segment has a positive slope even where the CIF
  # has flat stretches (e.g., before maturity).
  i <- findInterval(target[ok], cif_value)
  frac <- (target[ok] - cif_value[i]) / (cif_value[i + 1] - cif_value[i])
  out[ok] <- cif_age[i] + frac * (cif_age[i + 1] - cif_age[i]) - age[ok]
  out
}

#' Next event for aphids of the given ages: whichever of death or birth
#' comes first. Returns a tibble with `event` ("death"/"birth") and `time`.
sample_aphid_event <- function(age, fecundity, params) {
  death <- sample_death_time(age, params$death_shape, params$death_scale)
  birth <- sample_birth_time(age, params$cif_age, params$cif_value, fecundity)
  tibble::tibble(
    event = dplyr::if_else(birth < death, "birth", "death"),
    time = pmin(birth, death)
  )
}

#' Simulate the lifetimes of a cohort of newborn aphids (no predators, and
#' offspring are counted but not followed). Used to check that the event
#' samplers reproduce the clip-cage life tables. Each aphid draws its own
#' fecundity multiplier at birth (not inherited from its mother).
#'
#' Returns one row per event: aphid id, event type and age at event (minutes).
simulate_aphid_cohort <- function(n, params) {
  age <- rep(0, n)
  fecundity <- sample_fecundity(n, params$frailty_shape)
  alive <- rep(TRUE, n)
  events <- list()

  while (any(alive)) {
    idx <- which(alive)
    nxt <- sample_aphid_event(age[idx], fecundity[idx], params)
    age[idx] <- age[idx] + nxt$time
    events[[length(events) + 1]] <- tibble::tibble(
      id = idx, event = nxt$event, age = age[idx]
    )
    alive[idx[nxt$event == "death"]] <- FALSE
  }

  dplyr::bind_rows(events) |> dplyr::arrange(id, age)
}
