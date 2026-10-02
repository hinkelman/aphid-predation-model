# Predator foraging behavior: fitting and sampling --------------------------
#
# Durations are fitted as parametric accelerated failure time (AFT) models
# with survival::survreg so they can be sampled in the simulation:
#   log(T) = X %*% beta + sigma * eps
# with eps ~ logistic (log-logistic T) or minimum extreme value (Weibull T).
# Distributions were chosen by AIC (see analysis/02_predator_behavior.R).
# Times are in minutes; starvation (time since last meal) is in hours.

# Starvation range in the 2008 behavior trials; predictions are clamped to
# this range rather than extrapolated.
starve_range <- c(2, 24)

fit_aft <- function(formula, data, dist) {
  fit <- survival::survreg(formula, data = data, dist = dist)
  list(
    coef = coef(fit),
    sigma = fit$scale,
    dist = dist,
    terms = delete.response(terms(fit)),
    xlevels = fit$xlevels,
    aic = AIC(fit),
    fit = fit
  )
}

#' Linear predictor for new data from a fitted AFT model.
aft_lp <- function(model, newdata) {
  X <- model.matrix(model$terms, model.frame(model$terms, newdata, xlev = model$xlevels))
  drop(X %*% model$coef)
}

#' Draw durations from a fitted AFT model given its linear predictor.
sample_aft <- function(lp, sigma, dist) {
  u <- runif(length(lp))
  eps <- switch(dist,
    weibull = log(-log(u)),
    loglogistic = qlogis(u),
    lognormal = qnorm(u),
    stop("unsupported distribution: ", dist)
  )
  exp(lp + sigma * eps)
}

#' Survival function of a fitted AFT model at times `t` for linear predictor `lp`.
aft_survival <- function(t, lp, sigma, dist) {
  z <- (log(t) - lp) / sigma
  switch(dist,
    weibull = exp(-exp(z)),
    loglogistic = plogis(-z),
    lognormal = pnorm(-z)
  )
}

#' Fit all predator-behavior components.
#'
#' * handling: 2008 trials, log-logistic ~ aphid * starve (size-matched prey)
#' * post_handling: 2008 trials, Weibull ~ aphid * starve (time to leave leaf)
#' * age_exponent: 2005 trials, handling time scales as aphid_age^b (no
#'   species effect once age is included)
#' * rejection: 2005 trials, logistic ~ aphid_age
#' * partial: 2008 trials, proportion of meals only partially consumed
fit_predator_behavior <- function(trials = read_behavior_trials(),
                                  handle_age = read_handle_age()) {
  trials <- dplyr::mutate(trials, post_handle = pmax(post_handle, 0.05))
  handled <- dplyr::filter(handle_age, rejected == 0, !is.na(handle))

  age_model <- fit_aft(
    survival::Surv(handle, handle_end) ~ log(aphid_age), handled, "loglogistic"
  )
  rejection <- glm(rejected ~ aphid_age, family = binomial, data = handle_age)

  list(
    handling = fit_aft(
      survival::Surv(handle, handle_end) ~ aphid * starve, trials, "loglogistic"
    ),
    post_handling = fit_aft(
      survival::Surv(post_handle, left) ~ aphid * starve, trials, "weibull"
    ),
    age_exponent = unname(age_model$coef["log(aphid_age)"]),
    age_model = age_model,
    rejection = rejection,
    partial = trials |>
      dplyr::summarise(p_partial = mean(partial, na.rm = TRUE), .by = aphid)
  )
}

#' Sample handling times (minutes).
#'
#' The 2008 fit applies to size-matched prey: an adult bean aphid, or a pea
#' aphid of the same mass. For other aphids, handling is scaled by
#' (aphid_age / ref_age)^b from the 2005 age trials, where `ref_age` is the
#' age of the 2008 prey of that species (see `reference_prey_age()`).
sample_handling_time <- function(aphid, starve, aphid_age, behavior, ref_age) {
  m <- behavior$handling
  nd <- tibble::tibble(
    aphid = factor(aphid, levels = c("pea", "bean")),
    starve = pmin(pmax(starve, starve_range[1]), starve_range[2])
  )
  age_shift <- behavior$age_exponent * log(aphid_age / ref_age[as.character(aphid)])
  sample_aft(aft_lp(m, nd) + age_shift, m$sigma, m$dist)
}

#' Sample post-handling times (minutes until leaving the plant if no other
#' prey is encountered). `starve` is hunger before the meal just eaten.
sample_post_handling_time <- function(aphid, starve, behavior) {
  m <- behavior$post_handling
  nd <- tibble::tibble(
    aphid = factor(aphid, levels = c("pea", "bean")),
    starve = pmin(pmax(starve, starve_range[1]), starve_range[2])
  )
  sample_aft(aft_lp(m, nd), m$sigma, m$dist)
}

#' Probability that a captured aphid of a given age is rejected (not eaten).
p_reject <- function(aphid_age, behavior) {
  predict(behavior$rejection, tibble::tibble(aphid_age = aphid_age), type = "response")
}

#' Ages (days) of the prey used in the 2008 trials: adult bean aphids (taken
#' as the age at maturity) and pea aphids of the same mass.
reference_prey_age <- function(mass_anchors) {
  bean <- dplyr::filter(mass_anchors, aphid == "bean")
  pea <- dplyr::filter(mass_anchors, aphid == "pea")
  c(
    pea = aphid_age_at_mass(bean$mass_adult, pea),
    bean = bean$age_mature
  )
}
