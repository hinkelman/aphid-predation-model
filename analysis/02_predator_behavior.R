# Predator foraging behavior ----------------------------------------------
#
# Fits handling time, post-handling (leaf departure) time, rejection and
# partial consumption from the 2008 behavior trials and 2005 handling-by-age
# trials, and builds the provisional aphid mass-at-age model used to place
# the size-matched trial prey on each species' age axis.
#
# Outputs:
#   output/params/predator_behavior.rds
#   output/figures/behavior_*.png

library(dplyr)
library(tidyr)
library(ggplot2)

source("R/data.R")
source("R/aphid_mass.R")
source("R/predator_behavior.R")
source("R/plot_theme.R")

dir.create("output/params", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
set.seed(2008)

trials <- read_behavior_trials()
handle_age <- read_handle_age()

# Distribution choice (AIC) -----------------------------------------------

aic_table <- crossing(
  response = c("handle", "post_handle"),
  dist = c("exponential", "weibull", "lognormal", "loglogistic")
) |>
  mutate(aic = purrr::map2_dbl(response, dist, \(r, d) {
    dat <- mutate(trials, post_handle = pmax(post_handle, 0.05))
    f <- if (r == "handle") {
      survival::Surv(handle, handle_end) ~ aphid * starve
    } else {
      survival::Surv(post_handle, left) ~ aphid * starve
    }
    AIC(survival::survreg(f, data = dat, dist = d))
  })) |>
  mutate(delta_aic = aic - min(aic), .by = response) |>
  arrange(response, aic)
aic_table

# Fit -----------------------------------------------------------------------

behavior <- fit_predator_behavior(trials, handle_age)
summary(behavior$handling$fit)
summary(behavior$post_handling$fit)
summary(behavior$age_model$fit)
summary(behavior$rejection)
behavior$partial

mass_anchors <- aphid_mass_anchors()
ref_age <- reference_prey_age(mass_anchors)
mass_anchors
ref_age

saveRDS(
  list(behavior = behavior, mass_anchors = mass_anchors, ref_age = ref_age),
  "output/params/predator_behavior.rds"
)

# Model vs Kaplan-Meier, by starvation class --------------------------------

starve_classes <- function(d) {
  mutate(d, starve_class = cut(starve, c(0, 8, 16, 24),
                               labels = c("2–8 h", "8–16 h", "16–24 h")))
}

km_by_group <- function(d, time, event) {
  d |>
    group_by(aphid, starve_class) |>
    group_modify(\(g, k) {
      fit <- survival::survfit(survival::Surv(g[[time]], g[[event]]) ~ 1)
      tibble(t = c(0, fit$time), surv = c(1, fit$surv))
    }) |>
    ungroup()
}

model_curves <- function(d, model, tmax) {
  d |>
    summarise(starve = median(starve), .by = c(aphid, starve_class)) |>
    mutate(lp = aft_lp(model, pick(aphid, starve))) |>
    crossing(t = exp(seq(log(0.5), log(tmax), length.out = 200))) |>
    mutate(surv = aft_survival(t, lp, model$sigma, model$dist))
}

plot_duration_fit <- function(d, time, event, model, title, x_lab) {
  d <- starve_classes(d)
  km <- km_by_group(d, time, event)
  fit <- model_curves(d, model, max(d[[time]]))
  ggplot(mapping = aes(t, surv, colour = aphid)) +
    geom_step(data = filter(km, t > 0), linewidth = 0.4, alpha = 0.7) +
    geom_line(data = fit, linewidth = 0.7, linetype = "22") +
    facet_wrap(~starve_class, nrow = 1) +
    scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    scale_colour_manual(values = species_colours, labels = species_labels) +
    labs(
      x = x_lab, y = "Proportion still in behavior", colour = NULL, title = title,
      subtitle = "Steps: Kaplan-Meier from 2008 trials; dashed: fitted model at the median starvation of each class"
    ) +
    theme_model()
}

p_handle <- plot_duration_fit(
  trials, "handle", "handle_end", behavior$handling,
  "Handling time (log-logistic, aphid × starvation)", "Handling time (min, log scale)"
)
ggsave("output/figures/behavior_handling.png", p_handle, width = 8, height = 3.8, dpi = 200)

p_post <- plot_duration_fit(
  mutate(trials, post_handle = pmax(post_handle, 0.05)), "post_handle", "left",
  behavior$post_handling,
  "Post-handling time until leaving the leaf (Weibull, aphid × starvation)",
  "Post-handling time (min, log scale)"
)
ggsave("output/figures/behavior_post_handling.png", p_post, width = 8, height = 3.8, dpi = 200)

# Handling time vs aphid age (2005) and the 2008 reference prey -------------

age_grid <- tibble(aphid_age = seq(1, 12, by = 0.1)) |>
  mutate(lp = aft_lp(behavior$age_model, pick(aphid_age)),
         median = exp(lp)) # log-logistic median = exp(lp)

ref_2008 <- tibble(
  aphid = factor(names(ref_age), levels = c("pea", "bean")),
  aphid_age = ref_age
) |>
  mutate(median = exp(aft_lp(behavior$handling, tibble(aphid = aphid, starve = 2.7))))

p_age <- handle_age |>
  filter(rejected == 0, !is.na(handle)) |>
  ggplot(aes(aphid_age, handle)) +
  geom_line(data = age_grid, aes(y = median), colour = "grey55", linewidth = 0.6) +
  geom_point(aes(colour = aphid, shape = factor(handle_end)),
             position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.8) +
  geom_point(data = ref_2008, aes(y = median, fill = aphid), shape = 23, size = 3.5,
             colour = "white", stroke = 0.8) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = species_colours, labels = species_labels) +
  scale_fill_manual(values = species_colours, labels = species_labels, guide = "none") +
  scale_shape_manual(values = c(`1` = 16, `0` = 1), labels = c(`1` = "Observed", `0` = "Censored")) +
  labs(
    x = "Aphid age (days)", y = "Handling time (min, log scale)", colour = NULL, shape = NULL,
    title = "Handling time increases with aphid age (2005 trials)",
    subtitle = "Grey: fitted median, both species. Diamonds: 2008 fit for size-matched prey (2.7 h starved)"
  ) +
  theme_model()
ggsave("output/figures/behavior_handling_by_age.png", p_age, width = 7, height = 4.5, dpi = 200)

# Provisional aphid mass at age ---------------------------------------------

mass_curve <- crossing(aphid = factor(c("pea", "bean"), levels = c("pea", "bean")),
                       age = seq(0, 12, by = 0.05)) |>
  mutate(mass = purrr::map2_dbl(age, aphid, \(a, sp) aphid_mass(a, filter(mass_anchors, aphid == sp))))

p_mass <- ggplot(mass_curve, aes(age, mass, colour = aphid)) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = tibble(aphid = factor(names(ref_age), levels = c("pea", "bean")),
                  age = ref_age, mass = 0.9),
    shape = 21, fill = "white", size = 2.5, stroke = 1
  ) +
  scale_colour_manual(values = species_colours, labels = species_labels) +
  labs(
    x = "Aphid age (days)", y = "Fresh mass (mg)", colour = NULL,
    title = "Provisional aphid mass at age",
    subtitle = "Assumed exponential growth to adult mass at maturity\nCircles: 0.9 mg prey used in the predator trials"
  ) +
  theme_model()
ggsave("output/figures/behavior_aphid_mass.png", p_mass, width = 6, height = 4, dpi = 200)

# Sampler check ------------------------------------------------------------
#
# Draw handling and post-handling times for the trial covariates (with the
# 2008 reference prey ages, so the age multiplier is 1) and compare with the
# Kaplan-Meier curves of the data, pooled over starvation.

n_rep <- 200
sim <- trials |>
  slice(rep(seq_len(n()), each = n_rep)) |>
  mutate(
    handle = sample_handling_time(aphid, starve, ref_age[as.character(aphid)], behavior, ref_age),
    post_handle = sample_post_handling_time(aphid, starve, behavior)
  )

km_pooled <- function(d, time, event) {
  d |>
    group_by(aphid) |>
    group_modify(\(g, k) {
      fit <- survival::survfit(survival::Surv(g[[time]], g[[event]]) ~ 1)
      tibble(t = c(0, fit$time), surv = c(1, fit$surv))
    }) |>
    ungroup()
}

check <- bind_rows(
  km_pooled(trials, "handle", "handle_end") |> mutate(quantity = "Handling", source = "2008 trials"),
  km_pooled(mutate(sim, handle_end = 1), "handle", "handle_end") |> mutate(quantity = "Handling", source = "Simulated"),
  km_pooled(mutate(trials, post_handle = pmax(post_handle, 0.05)), "post_handle", "left") |>
    mutate(quantity = "Post-handling", source = "2008 trials"),
  km_pooled(mutate(sim, left = 1), "post_handle", "left") |> mutate(quantity = "Post-handling", source = "Simulated")
) |>
  filter(t > 0.1)

p_check <- ggplot(check, aes(t, surv, colour = aphid, linetype = source)) +
  geom_step(linewidth = 0.6) +
  facet_wrap(~quantity, scales = "free_x") +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = species_colours, labels = species_labels) +
  scale_linetype_manual(values = c("solid", "22")) +
  labs(
    x = "Minutes (log scale)", y = "Proportion still in behavior", colour = NULL, linetype = NULL,
    title = "Behavior samplers reproduce the 2008 trials",
    subtitle = paste(n_rep, "draws per trial at each trial's starvation level")
  ) +
  theme_model()
ggsave("output/figures/behavior_sampler_check.png", p_check, width = 8, height = 4, dpi = 200)
