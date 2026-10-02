# Aphid demography: fitting -----------------------------------------------
#
# Fits are done in days (the resolution of the clip-cage data) and converted
# to minutes, the simulation clock unit, by `aphid_params_minutes()`.

minutes_per_day <- 1440

#' Weibull lifespan fit for one aphid species.
#'
#' Deaths were checked daily, so a death recorded on day d happened in
#' (d - 1, d]; censored aphids were last seen alive on day d. `interval = TRUE`
#' uses that interval censoring; `interval = FALSE` treats the recorded day as
#' the exact death time, as in the 2010 analysis.
#'
#' Returns shape and scale in R's `dweibull()` parameterization (days).
fit_lifespan_weibull <- function(lifespans, interval = TRUE) {
  surv <- if (interval) {
    with(lifespans, survival::Surv(
      time = lifespan - died, # last day known alive (death) or censoring day
      time2 = dplyr::if_else(died == 1, lifespan, NA_real_),
      type = "interval2"
    ))
  } else {
    with(lifespans, survival::Surv(lifespan, died))
  }

  fit <- survival::survreg(surv ~ 1, dist = "weibull")
  tibble::tibble(
    shape = 1 / fit$scale,
    scale = exp(unname(coef(fit))),
    n = nrow(lifespans),
    n_died = sum(lifespans$died)
  )
}

#' Empirical cumulative intensity function (CIF) for births.
#'
#' Following Leemis (2004): the mean number of offspring per aphid alive on
#' each day of age is summed to give the CIF at day boundaries. Between knots
#' the CIF is linear (i.e., constant birth intensity within a day). Aphid-days
#' with a missing offspring count are dropped from the at-risk set.
#'
#' Returns knots at ages 0, 1, ..., max day (days), with the number at risk and
#' an approximate 95% CI (Poisson, as in the original analysis).
fit_birth_cif <- function(clip_cage) {
  daily <- clip_cage |>
    dplyr::filter(!is.na(offspring)) |>
    dplyr::summarise(
      n_at_risk = dplyr::n(),
      births = sum(offspring),
      .by = day
    ) |>
    dplyr::arrange(day) |>
    dplyr::mutate(
      intensity = births / n_at_risk,
      cumint = cumsum(intensity)
    )

  dplyr::bind_rows(
    tibble::tibble(day = 0L, n_at_risk = daily$n_at_risk[1], births = 0L,
                   intensity = 0, cumint = 0),
    daily
  ) |>
    dplyr::mutate(
      half_width = qnorm(0.975) * sqrt(cumint / n_at_risk),
      lower = pmax(cumint - half_width, 0),
      upper = cumint + half_width
    ) |>
    dplyr::select(-half_width)
}

#' Individual variation (frailty) in fecundity.
#'
#' Each aphid has a fecundity multiplier m ~ Gamma(shape = k, rate = k), so
#' mean 1 and variance 1/k, that scales its birth intensity. Given m, offspring
#' over an observed life of L days are Poisson(m * CIF(L)), so lifetime
#' offspring are negative binomial with offset log(CIF(L)) and theta = k.
#' Censored aphids contribute through the CIF at their censoring age; those
#' censored before any reproduction (CIF = 0) carry no information and are
#' dropped.
#'
#' Returns k (`frailty_shape`) and the SD of the multiplier.
fit_fecundity_frailty <- function(lifespans, cif) {
  dat <- lifespans |>
    dplyr::mutate(
      expected = approx(cif$day, cif$cumint, xout = lifespan, rule = 2)$y
    ) |>
    dplyr::filter(expected > 0)
  fit <- MASS::glm.nb(total_offspring ~ 1 + offset(log(expected)), data = dat)

  tibble::tibble(
    frailty_shape = fit$theta,
    frailty_sd = 1 / sqrt(fit$theta),
    # intercept should be ~0: the CIF already gives the mean
    log_mean_ratio = unname(coef(fit))
  )
}

#' Fit lifespan and birth models for both species.
fit_aphid_demography <- function(clip_cage = read_clip_cage()) {
  lifespans <- clip_cage_lifespans(clip_cage)
  birth_cif <- clip_cage |>
    tidyr::nest(.by = aphid) |>
    dplyr::mutate(cif = purrr::map(data, fit_birth_cif)) |>
    tidyr::unnest(cif) |>
    dplyr::select(-data)

  list(
    lifespan = lifespans |>
      tidyr::nest(.by = aphid) |>
      dplyr::mutate(fit = purrr::map(data, fit_lifespan_weibull)) |>
      tidyr::unnest(fit) |>
      dplyr::select(-data),
    birth_cif = birth_cif,
    frailty = lifespans |>
      tidyr::nest(.by = aphid) |>
      dplyr::mutate(fit = purrr::map2(
        data, aphid,
        \(d, sp) fit_fecundity_frailty(d, dplyr::filter(birth_cif, aphid == sp))
      )) |>
      tidyr::unnest(fit) |>
      dplyr::select(-data)
  )
}

#' Convert fitted demography to simulation parameters in minutes.
#'
#' Returns a named list (by species) of
#' `list(death_shape, death_scale, cif_age, cif_value, frailty_shape)`.
aphid_params_minutes <- function(demography) {
  species <- levels(demography$lifespan$aphid)
  purrr::set_names(species) |>
    purrr::map(\(sp) {
      life <- dplyr::filter(demography$lifespan, aphid == sp)
      cif <- dplyr::filter(demography$birth_cif, aphid == sp)
      frailty <- dplyr::filter(demography$frailty, aphid == sp)
      list(
        death_shape = life$shape,
        death_scale = life$scale * minutes_per_day,
        cif_age = cif$day * minutes_per_day,
        cif_value = cif$cumint,
        frailty_shape = frailty$frailty_shape
      )
    })
}
