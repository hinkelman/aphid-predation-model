# Predator development and survival ----------------------------------------
#
# Food is measured in "units": one size-matched aphid (~0.9 mg) fully
# eaten, i.e. a pea aphid of that size. In the simulation an aphid of mass m
# gives m / 0.9 units if it is a pea aphid, and v * m / 0.9 if it is a bean
# aphid, where v is the value per bean kill calibrated below.
#
# Fitted from the larval diet experiment (dissertation Ch. 1) using each
# larva's own observed daily kills, so development is fitted separately
# from intake and foraging behavior.

larval_stages <- c("L1", "L2", "L3", "L4")
unit_mass <- 0.9 # mg, size-matched prey in all predator experiments

#' Food eaten and outcome for each larva in each instar.
stage_intake <- function(diet_intervals) {
  diet_intervals |>
    dplyr::filter(stage %in% larval_stages) |>
    dplyr::summarise(
      pea = sum(pea),
      bean = sum(bean),
      days = dplyr::n(),
      molted = any(next_stage != stage & died == 0),
      .by = c(id, diet, stage)
    ) |>
    dplyr::mutate(stage = droplevels(stage))
}

#' Food thresholds for each instar and the value of a bean kill.
#'
#' Model for completed instars: log(pea + v_diet * bean) = log(threshold_stage)
#' + e, e ~ N(0, sigma). v is estimated separately for the bean-only diet
#' (`v_bean_alone`) and the mixed diet (`v_bean_mixed`): mixed-diet larvae
#' killed many bean aphids but got little development from them.
fit_development_thresholds <- function(diet_intervals) {
  st <- stage_intake(diet_intervals) |>
    dplyr::filter(molted, !is.na(pea), !is.na(bean))

  units <- function(v) with(st, pea + dplyr::case_when(diet == "B" ~ v[1], diet == "M" ~ v[2], TRUE ~ 0) * bean)
  rss <- function(v) sum(resid(lm(log(units(v)) ~ 0 + stage, data = st))^2)
  opt <- optim(c(0.5, 0.5), rss, method = "L-BFGS-B", lower = 0.01, upper = 2)

  fit <- lm(log(units(opt$par)) ~ 0 + stage, data = st)
  list(
    threshold = setNames(exp(coef(fit)), larval_stages),
    sigma = sigma(fit),
    v_bean_alone = opt$par[1],
    v_bean_mixed = opt$par[2],
    data = dplyr::mutate(st, units = units(opt$par), resid = resid(fit))
  )
}

#' Daily mortality hazard of larvae and pupae (complementary log-log model
#' on check intervals, so exp(linear predictor) is a hazard per day).
#'
#' Baseline hazard by stage group, plus in L1 a bean effect that increases
#' with bean eaten so far and is buffered by pea eaten so far:
#'   hazard = exp(alpha_stage) * (1 + bean / (1 + pea))^beta   (L1 only)
#' In later instars no diet effect on the daily hazard was detectable; bean
#' diets still raise mortality there by prolonging development.
fit_larval_mortality <- function(diet_intervals) {
  d <- diet_intervals |>
    dplyr::filter(stage %in% c(larval_stages, "Pupa"), censored == 0) |>
    dplyr::arrange(id, day) |>
    dplyr::mutate(
      bean_cum = cumsum(dplyr::coalesce(dplyr::lag(bean, default = 0), 0)),
      pea_cum = cumsum(dplyr::coalesce(dplyr::lag(pea, default = 0), 0)),
      .by = id
    ) |>
    dplyr::mutate(
      stage_group = mortality_stage_group(stage),
      l1 = as.numeric(stage == "L1"),
      bean_load = log1p(bean_cum / (1 + pea_cum))
    )

  fit <- glm(died ~ 0 + stage_group + l1:bean_load, family = binomial("cloglog"), data = d)
  b <- coef(fit)
  list(
    log_hazard = setNames(b[grepl("^stage_group", names(b))],
                          sub("stage_group", "", names(b)[grepl("^stage_group", names(b))])),
    l1_bean_exponent = unname(b["l1:bean_load"]),
    fit = fit
  )
}

mortality_stage_group <- function(stage) {
  factor(dplyr::case_when(
    stage == "L1" ~ "L1",
    stage %in% c("L2", "L3") ~ "L2-3",
    stage == "L4" ~ "L4",
    TRUE ~ "Pupa"
  ), levels = c("L1", "L2-3", "L4", "Pupa"))
}

#' Daily mortality hazard for a larva in `stage` that has eaten `bean_cum`
#' and `pea_cum` (units, since hatching) - vectorized.
mortality_hazard <- function(stage, bean_cum, pea_cum, mortality) {
  base <- exp(mortality$log_hazard[as.character(mortality_stage_group(stage))])
  bean_term <- dplyr::if_else(
    stage == "L1",
    (1 + bean_cum / (1 + pea_cum))^mortality$l1_bean_exponent,
    1
  )
  unname(base * bean_term)
}

#' Ad libitum kill rate (aphids per day) in each instar, by diet. These set
#' the gut-limited maximum intake in the simulation.
#'
#' In L4, feeding tapers off before pupation (pea diet: ~39 and ~41 kills on
#' L4 days 1-2, ~18 on day 3, ~5 on day 4), so the L4 rate uses only the
#' first two L4 intervals; the taper is represented by a non-feeding
#' pre-pupa period in the simulation.
fit_max_intake <- function(diet_intervals) {
  diet_intervals |>
    dplyr::filter(stage %in% larval_stages, !is.na(pea), !is.na(bean), died == 0) |>
    dplyr::arrange(id, day) |>
    dplyr::mutate(stage_day = dplyr::row_number(), .by = c(id, stage)) |>
    dplyr::filter(!(stage == "L4" & (stage_day > 2 | next_stage == "Pupa"))) |>
    dplyr::mutate(kills = pea + bean) |>
    dplyr::summarise(kills_per_day = mean(kills), n = dplyr::n(), .by = c(diet, stage)) |>
    dplyr::arrange(stage, diet)
}

#' Critical L4 food for pupation without further feeding.
#'
#' Diet-timing larvae were starved from day `timing` of L4; the food eaten in
#' L4 before starvation is taken as the mean cumulative L4 kills of pea-diet
#' larvae by that day (diet experiment). Logistic fit of pupation on food.
#' Also returns the delay from the start of starvation to pupation for
#' larvae that pupated, as a linear function of L4 food.
fit_critical_food <- function(diet_intervals, diet_timing) {
  l4_food_by_day <- diet_intervals |>
    dplyr::filter(diet == "P", stage == "L4") |>
    dplyr::arrange(id, day) |>
    dplyr::mutate(timing = dplyr::row_number(), food = cumsum(pea) - pea, .by = id) |>
    dplyr::summarise(food = mean(food), .by = timing)

  starved <- diet_timing |>
    dplyr::filter(trt == "S") |>
    dplyr::left_join(l4_food_by_day, by = "timing")

  pupation <- glm(pupated ~ food, family = binomial, data = starved)
  delay <- lm(trt_to_pupa ~ food, data = dplyr::filter(starved, pupated == 1))

  list(
    l4_food_by_day = l4_food_by_day,
    pupation = pupation,
    critical_food = unname(-coef(pupation)[1] / coef(pupation)[2]), # 50% point
    delay = delay,
    data = starved
  )
}

#' Probability that a larva reaching the L4 food threshold pupates
#' successfully, given its lifetime bean and pea kills.
#'
#' Bean-fed 4th instars often ate well past the threshold and still failed
#' to pupate (one was noted "attempting to pupate" when it died), so bean
#' toxicity is applied at the pupation attempt, with the same pea-buffered
#' load as the L1 hazard: P(success) = (1 + bean / (1 + pea))^-gamma.
#' gamma is calibrated by simulation (analysis/03_predator_development.R).
pupation_success <- function(bean_cum, pea_cum, gamma) {
  (1 + bean_cum / (1 + pea_cum))^(-gamma)
}

#' Fit all development components.
fit_predator_development <- function(diet_intervals = read_diet_experiment(),
                                     diet_timing = read_diet_timing()) {
  list(
    thresholds = fit_development_thresholds(diet_intervals),
    mortality = fit_larval_mortality(diet_intervals),
    max_intake = fit_max_intake(diet_intervals),
    critical = fit_critical_food(diet_intervals, diet_timing)
  )
}
