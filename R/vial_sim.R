# Vial simulation of larval development ------------------------------------
#
# Validation tool: simulates larvae reared individually in vials with ad
# libitum aphids, as in the diet and diet-timing experiments, using the
# fitted development components. Intake is gut-limited: kills arrive as a
# Poisson process at the observed ad lib rate for the stage. Hourly steps,
# vectorized over larvae.
#
# Stage codes: 1-4 = L1-L4, 5 = pre-pupa (non-feeding), 6 = pupated, 0 = dead.

#' @param n number of larvae
#' @param diet "P", "B" or "M" (rearing diet from hatching)
#' @param development output of fit_predator_development()
#' @param switch optional list(day = L4 day of switch (1-4), to = "B" or "S")
#'   implementing the diet-timing treatments
#' @param threshold_sd SD of individual log-threshold multipliers
#' @param prepupa_days duration of the non-feeding pre-pupa
#' @param max_days stop after this many days
#' @param pupation_gamma bean-load exponent for pupation failure
#'   (see pupation_success())
simulate_vial_larvae <- function(n, diet, development, switch = NULL,
                                 threshold_sd = 0.2, prepupa_days = 1.5,
                                 pupation_gamma = 0, max_days = 60) {
  th <- development$thresholds
  mort <- development$mortality
  crit <- development$critical
  rates <- development$max_intake |>
    dplyr::filter(diet == !!diet) |>
    dplyr::arrange(stage) |>
    dplyr::pull(kills_per_day)
  bean_share <- c(P = 0, B = 1, M = 0.55)[[diet]]
  v_bean <- if (diet == "M") th$v_bean_mixed else th$v_bean_alone

  dt <- 1 / 24
  stage <- rep(1L, n)
  food <- rep(0, n) # units in current instar
  bean_cum <- pea_cum <- rep(0, n) # kills since hatching
  mult <- exp(rnorm(n * 4, 0, threshold_sd)) |> matrix(ncol = 4)
  threshold <- sweep(mult, 2, th$threshold, `*`)
  crit_food <- rlogis(n, -coef(crit$pupation)[1] / coef(crit$pupation)[2],
                      1 / coef(crit$pupation)[2])
  stage_entry <- matrix(NA_real_, n, 6) # day each stage was entered
  stage_entry[, 1] <- 0
  l4_entry <- rep(NA_real_, n)
  switched <- rep(FALSE, n)
  starving <- rep(FALSE, n)
  prepupa_left <- rep(NA_real_, n)
  pupation_day <- rep(NA_real_, n)
  switch_day <- rep(NA_real_, n)

  t <- 0
  while (t < max_days && any(stage %in% 1:5)) {
    alive_larva <- stage %in% 1:4

    # diet switch in L4
    if (!is.null(switch)) {
      do_switch <- stage == 4 & !switched & (t - l4_entry) >= (switch$day - 1)
      switched[do_switch] <- TRUE
      switch_day[do_switch] <- t
      if (switch$to == "S") {
        starving[do_switch] <- TRUE
        # committed larvae pupate after a food-dependent delay; others die
        # (time to death is not modeled here, only the fate)
        will_pupate <- do_switch & food >= crit_food
        delay <- pmax(predict(crit$delay, data.frame(food = food[will_pupate])), prepupa_days)
        stage[will_pupate] <- 5L
        prepupa_left[will_pupate] <- delay
        stage[do_switch & !will_pupate] <- 0L
      }
    }

    # intake
    feeding <- stage %in% 1:4 & !starving
    rate <- rates[pmin(stage, 4)] * dt
    kills <- ifelse(feeding, rpois(n, rate), 0L)
    is_bean_diet <- if (!is.null(switch) && switch$to == "B") switched else rep(FALSE, n)
    p_bean <- ifelse(is_bean_diet, 1, bean_share)
    bean <- rbinom(n, kills, p_bean)
    pea <- kills - bean
    v <- ifelse(is_bean_diet, th$v_bean_alone, v_bean)
    food <- food + pea + v * bean
    bean_cum <- bean_cum + bean
    pea_cum <- pea_cum + pea

    # mortality (larvae and pre-pupae)
    at_risk <- stage %in% 1:5
    h <- mortality_hazard(
      factor(c("L1", "L2", "L3", "L4", "L4")[pmax(stage, 1)][at_risk], levels = c("L1", "L2", "L3", "L4", "Pupa")),
      bean_cum[at_risk], pea_cum[at_risk], mort
    )
    die <- rep(FALSE, n)
    die[at_risk] <- runif(sum(at_risk)) < 1 - exp(-h * dt)
    stage[die] <- 0L

    # molting
    idx <- which(stage %in% 1:4)
    molt <- idx[food[idx] >= threshold[cbind(idx, stage[idx])]]
    to_prepupa <- molt[stage[molt] == 4]
    to_next <- setdiff(molt, to_prepupa)
    stage[to_next] <- stage[to_next] + 1L
    food[to_next] <- 0
    stage_entry[cbind(to_next, stage[to_next])] <- t + dt
    l4_entry[to_next[stage[to_next] == 4]] <- t + dt
    fails <- to_prepupa[runif(length(to_prepupa)) >
      pupation_success(bean_cum[to_prepupa], pea_cum[to_prepupa], pupation_gamma)]
    to_prepupa <- setdiff(to_prepupa, fails)
    stage[fails] <- 0L
    stage[to_prepupa] <- 5L
    prepupa_left[to_prepupa] <- prepupa_days

    # pre-pupa -> pupa
    pp <- which(stage == 5)
    prepupa_left[pp] <- prepupa_left[pp] - dt
    done <- pp[prepupa_left[pp] <= 0]
    stage[done] <- 6L
    pupation_day[done] <- t + dt

    t <- t + dt
  }

  tibble::tibble(
    larva = seq_len(n),
    diet = diet,
    pupated = stage == 6L,
    L1 = stage_entry[, 2] - stage_entry[, 1],
    L2 = stage_entry[, 3] - stage_entry[, 2],
    L3 = stage_entry[, 4] - stage_entry[, 3],
    L4 = pupation_day - stage_entry[, 4],
    pupation_day = pupation_day,
    trt_to_pupa = pupation_day - switch_day,
    reached_l4 = !is.na(stage_entry[, 4]),
    switched = switched
  )
}
