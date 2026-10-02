# Aphid demography in the absence of predators -----------------------------
#
# Fits lifespan (Weibull) and birth (cumulative intensity function) models to
# the clip-cage life tables, checks them against the data, and verifies that
# the event samplers used in the simulation reproduce the life tables.
#
# Outputs:
#   output/params/aphid_demography.rds   fitted parameters (days and minutes)
#   output/figures/aphid_*.png           diagnostic figures

library(dplyr)
library(tidyr)
library(ggplot2)

source("R/data.R")
source("R/aphid_demography.R")
source("R/aphid_events.R")
source("R/plot_theme.R")

dir.create("output/params", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
set.seed(20091215)

clip_cage <- read_clip_cage()
lifespans <- clip_cage_lifespans(clip_cage)

# Fit ---------------------------------------------------------------------

demography <- fit_aphid_demography(clip_cage)
demography$lifespan

# Compare interval-censored fit with the exact-time fit used in 2010.
# 2010 values: bean shape 5.718, scale 29.24 d; pea shape 11.27, scale 23.97 d.
lifespan_exact <- lifespans |>
  nest(.by = aphid) |>
  mutate(fit = purrr::map(data, \(d) fit_lifespan_weibull(d, interval = FALSE))) |>
  unnest(fit) |>
  select(-data)
lifespan_exact

params <- aphid_params_minutes(demography)
saveRDS(
  list(fit = demography, params = params),
  "output/params/aphid_demography.rds"
)

# Lifespan: Kaplan-Meier vs Weibull ----------------------------------------

km <- survival::survfit(survival::Surv(lifespan, died) ~ aphid, data = lifespans)
km_df <- tibble(
  aphid = factor(rep(sub("aphid=", "", names(km$strata)), km$strata),
                 levels = levels(lifespans$aphid)),
  day = km$time,
  surv = km$surv
) |>
  # start each step curve at (0, 1)
  bind_rows(tibble(aphid = factor(levels(lifespans$aphid)), day = 0, surv = 1)) |>
  arrange(aphid, day)

weibull_curves <- bind_rows(
  `Interval-censored Weibull` = demography$lifespan,
  `Exact-time Weibull (2010 method)` = lifespan_exact,
  .id = "model"
) |>
  crossing(day = seq(0, 40, by = 0.25)) |>
  mutate(surv = pweibull(day, shape, scale, lower.tail = FALSE))

p_lifespan <- bind_rows(
  km_df |> mutate(model = "Kaplan-Meier (data)"),
  weibull_curves
) |>
  mutate(model = factor(model, levels = c(
    "Kaplan-Meier (data)", "Interval-censored Weibull", "Exact-time Weibull (2010 method)"
  ))) |>
  ggplot(aes(day, surv, colour = aphid, linetype = model)) +
  geom_step(data = \(d) filter(d, model == "Kaplan-Meier (data)"), linewidth = 0.4, alpha = 0.6) +
  geom_line(data = \(d) filter(d, model != "Kaplan-Meier (data)"), linewidth = 0.6) +
  facet_wrap(~aphid, labeller = labeller(aphid = species_labels)) +
  scale_colour_manual(values = species_colours, guide = "none") +
  scale_linetype_manual(values = c("solid", "22", "11")) +
  guides(linetype = guide_legend(override.aes = list(alpha = c(0.6, 1, 1), linewidth = c(0.4, 0.6, 0.6)))) +
  labs(
    x = "Age (days)", y = "Proportion surviving", linetype = NULL,
    title = "Aphid lifespan without predators"
  ) +
  theme_model()

ggsave("output/figures/aphid_lifespan.png", p_lifespan, width = 7, height = 3.8, dpi = 200)

# Births: cumulative intensity function -----------------------------------

individual_cum <- clip_cage |>
  mutate(cum_offspring = cumsum(coalesce(offspring, 0L)), .by = c(aphid, id))

p_cif <- ggplot(demography$birth_cif, aes(day)) +
  geom_line(
    data = individual_cum,
    aes(day, cum_offspring, group = id),
    colour = "grey82", linewidth = 0.3
  ) +
  geom_line(aes(y = lower, colour = aphid), linetype = "22", linewidth = 0.4) +
  geom_line(aes(y = upper, colour = aphid), linetype = "22", linewidth = 0.4) +
  geom_line(aes(y = cumint, colour = aphid), linewidth = 0.8) +
  facet_wrap(~aphid, labeller = labeller(aphid = species_labels)) +
  scale_colour_manual(values = species_colours, guide = "none") +
  labs(
    x = "Age (days)", y = "Cumulative offspring",
    title = "Birth cumulative intensity function (CIF)",
    subtitle = "Solid: CIF (mean offspring per aphid alive); dashed: 95% CI; grey: individual aphids"
  ) +
  theme_model()

ggsave("output/figures/aphid_birth_cif.png", p_cif, width = 7, height = 3.8, dpi = 200)

# Sampler check: simulate cohorts and re-estimate -------------------------
#
# Simulated cohorts of newborn aphids should reproduce (1) the fitted lifespan
# distribution and (2) the empirical CIF, when the simulated events are
# summarized with the same estimators used on the data.

n_sim <- 5000

sim_events <- params |>
  purrr::map(\(p) simulate_aphid_cohort(n_sim, p)) |>
  bind_rows(.id = "aphid") |>
  mutate(aphid = factor(aphid, levels = levels(lifespans$aphid)),
         age_day = age / minutes_per_day)

sim_lifespans <- sim_events |>
  filter(event == "death") |>
  select(aphid, id, lifespan = age_day)

sim_births <- sim_events |>
  filter(event == "birth") |>
  count(aphid, id, name = "total_offspring") |>
  right_join(select(sim_lifespans, aphid, id), by = c("aphid", "id")) |>
  mutate(total_offspring = coalesce(total_offspring, 0L))

# Simulated CIF: births per aphid alive during each day of age, as in the data.
sim_cif <- sim_events |>
  filter(event == "birth") |>
  mutate(day = ceiling(age_day)) |>
  count(aphid, day, name = "births") |>
  complete(aphid, day = 1:40, fill = list(births = 0L)) |>
  left_join(
    sim_lifespans |>
      mutate(last_day = ceiling(lifespan)) |>
      reframe(day = 1:40, n_at_risk = sapply(1:40, \(d) sum(last_day >= d)), .by = aphid),
    by = c("aphid", "day")
  ) |>
  filter(n_at_risk > 0) |>
  arrange(aphid, day) |>
  mutate(cumint = cumsum(births / n_at_risk), .by = aphid)

sim_surv <- sim_lifespans |>
  reframe(day = seq(0, 40, by = 0.25),
          surv = sapply(day, \(d) mean(lifespan > d)),
          .by = aphid)

p_check <- bind_rows(
  km_df |> mutate(quantity = "Proportion surviving", source = "Clip-cage data", value = surv),
  sim_surv |> mutate(quantity = "Proportion surviving", source = "Simulated cohort", value = surv),
  demography$birth_cif |> mutate(quantity = "Cumulative offspring (CIF)", source = "Clip-cage data", value = cumint),
  sim_cif |> mutate(quantity = "Cumulative offspring (CIF)", source = "Simulated cohort", value = cumint)
) |>
  mutate(quantity = factor(quantity, levels = c("Proportion surviving", "Cumulative offspring (CIF)"))) |>
  ggplot(aes(day, value, colour = aphid, linetype = source)) +
  geom_step(linewidth = 0.6) +
  facet_grid(quantity ~ aphid, scales = "free_y", switch = "y",
             labeller = labeller(aphid = species_labels)) +
  scale_colour_manual(values = species_colours, guide = "none") +
  scale_linetype_manual(values = c("solid", "22")) +
  labs(
    x = "Age (days)", y = NULL, linetype = NULL,
    title = "Event samplers reproduce the clip-cage life tables",
    subtitle = paste(format(n_sim, big.mark = ","), "simulated newborns per species")
  ) +
  theme_model() +
  theme(strip.placement = "outside")

ggsave("output/figures/aphid_sampler_check.png", p_check, width = 7, height = 5.5, dpi = 200)

# Lifetime fecundity: data vs simulation ----------------------------------
#
# Without individual variation, lifetime offspring in the simulation were
# about half as variable as in the data (SD ~7 vs ~14). Each aphid now has a
# gamma fecundity multiplier (fitted below); check that the spread matches.

demography$frailty

fecundity_summary <- bind_rows(
  `Clip-cage data` = lifespans |> filter(died == 1),
  `Simulated cohort` = sim_births,
  .id = "source"
) |>
  summarise(
    n = n(),
    mean = mean(total_offspring),
    sd = sd(total_offspring),
    cv = sd / mean,
    .by = c(aphid, source)
  ) |>
  arrange(aphid, source)
fecundity_summary

p_fecundity <- bind_rows(
  `Clip-cage data` = lifespans |> filter(died == 1),
  `Simulated cohort` = sim_births,
  .id = "source"
) |>
  ggplot(aes(total_offspring, after_stat(density), fill = aphid)) +
  geom_histogram(binwidth = 5, boundary = 0, colour = "white", linewidth = 0.3) +
  facet_grid(source ~ aphid, labeller = labeller(aphid = species_labels)) +
  scale_fill_manual(values = species_colours, guide = "none") +
  labs(
    x = "Lifetime offspring", y = "Density",
    title = "Lifetime fecundity",
    subtitle = "Data: aphids observed until death; simulation includes individual fecundity multipliers"
  ) +
  theme_model()

ggsave("output/figures/aphid_lifetime_fecundity.png", p_fecundity, width = 7, height = 5, dpi = 200)
