# Predator development and survival ----------------------------------------
#
# Fits food thresholds per instar, bean value per kill, larval mortality,
# ad lib intake, and the critical L4 food for pupation from the diet and
# diet-timing experiments; calibrates the bean pupation-failure exponent by
# simulating the vial experiments; and checks the combined model against
# both experiments.
#
# Outputs:
#   output/params/predator_development.rds
#   output/figures/development_*.png

library(dplyr)
library(tidyr)
library(ggplot2)

source("R/data.R")
source("R/predator_development.R")
source("R/vial_sim.R")
source("R/plot_theme.R")

dir.create("output/params", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
set.seed(20080915)

diet_intervals <- read_diet_experiment()
diet_timing <- read_diet_timing()
diet_summary <- read_mac_csv(data_path("diet_summary.csv")) |> janitor::clean_names()

diet_labels <- c(P = "Pea diet", M = "Mixed diet", B = "Bean diet")
diet_colours <- c(P = species_colours[["pea"]], M = "#8c6bb1", B = species_colours[["bean"]])

# Fit -----------------------------------------------------------------------

development <- fit_predator_development(diet_intervals, diet_timing)
development$thresholds[c("threshold", "sigma", "v_bean_alone", "v_bean_mixed")]
exp(development$mortality$log_hazard)
development$mortality$l1_bean_exponent
development$max_intake |> pivot_wider(id_cols = stage, names_from = diet, values_from = kills_per_day)
development$critical$critical_food
coef(development$critical$delay)

# Calibrate pupation failure on bean ----------------------------------------

obs_groups <- bind_rows(
  diet_summary |>
    summarise(obs = mean(pupate), n = n(), .by = diet) |>
    transmute(group = diet, experiment = "Diet", obs, n, obs_ttp = NA_real_),
  diet_timing |>
    filter(trt != "P") |>
    summarise(obs = mean(pupated), n = n(), obs_ttp = mean(trt_to_pupa, na.rm = TRUE), .by = c(trt, timing)) |>
    transmute(group = paste0(trt, timing), experiment = "Diet timing", obs, n, obs_ttp)
)

simulate_groups <- function(gamma, n = 2000) {
  diets <- c("P", "M", "B") |>
    purrr::map(\(d) simulate_vial_larvae(n, d, development, pupation_gamma = gamma) |>
                 mutate(group = d)) |>
    bind_rows()
  switches <- crossing(to = c("B", "S"), day = 1:4) |>
    purrr::pmap(\(to, day) {
      simulate_vial_larvae(n, "P", development, switch = list(day = day, to = to),
                           pupation_gamma = gamma) |>
        filter(switched) |> # experiment included larvae alive at the switch
        mutate(group = paste0(to, day))
    }) |>
    bind_rows()
  bind_rows(diets, switches)
}

group_summary <- function(sims) {
  sims |>
    summarise(sim = mean(pupated), sim_ttp = mean(trt_to_pupa, na.rm = TRUE), .by = group) |>
    left_join(obs_groups, by = "group")
}

binomial_deviance <- function(s) {
  with(s, -2 * sum(n * (obs * log(pmax(sim, 1e-3)) + (1 - obs) * log(pmax(1 - sim, 1e-3)))))
}

calibration <- tibble(gamma = seq(0, 0.3, by = 0.05)) |>
  mutate(deviance = purrr::map_dbl(gamma, \(g) binomial_deviance(group_summary(simulate_groups(g)))))
calibration

pupation_gamma <- calibration$gamma[which.min(calibration$deviance)]
pupation_gamma

saveRDS(
  list(development = development, pupation_gamma = pupation_gamma,
       threshold_sd = 0.2, prepupa_days = 1.5),
  "output/params/predator_development.rds"
)

# Check: simulated vs observed ----------------------------------------------

final <- simulate_groups(pupation_gamma, n = 4000)
check <- group_summary(final)
check

group_order <- c("P", "M", "B", paste0("S", 1:4), paste0("B", 1:4))
group_labels <- c(P = "Pea", M = "Mixed", B = "Bean",
                  S1 = "Starved d1", S2 = "Starved d2", S3 = "Starved d3", S4 = "Starved d4",
                  B1 = "Bean d1", B2 = "Bean d2", B3 = "Bean d3", B4 = "Bean d4")

p_pupation <- check |>
  mutate(
    group = factor(group, levels = group_order),
    experiment = factor(experiment, levels = c("Diet", "Diet timing")),
    se = sqrt(obs * (1 - obs) / n)
  ) |>
  ggplot(aes(y = group)) +
  geom_linerange(aes(xmin = pmax(obs - 2 * se, 0), xmax = pmin(obs + 2 * se, 1)),
                 colour = "grey70", linewidth = 2.5) +
  geom_point(aes(x = obs, shape = "Observed (± 2 SE)"), size = 2.4, colour = "grey30") +
  geom_point(aes(x = sim, shape = "Simulated"), size = 2.6, stroke = 1.1, colour = "#2a78d6") +
  facet_grid(experiment ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = group_labels, limits = rev) +
  scale_shape_manual(values = c(16, 4)) +
  scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  labs(
    x = "Larvae pupating", y = NULL, shape = NULL,
    title = "Pupation: model vs experiments",
    subtitle = "Diet: reared on one diet from hatching. Diet timing: pea-reared, switched on L4 day 1–4."
  ) +
  theme_model()
ggsave("output/figures/development_pupation_check.png", p_pupation, width = 7, height = 5, dpi = 200)

# Stage durations
obs_durations <- diet_intervals |>
  filter(stage %in% larval_stages) |>
  count(id, diet, stage, name = "days") |>
  left_join(select(diet_summary, id, pupate), by = "id") |>
  filter(pupate == 1) |>
  mutate(source = "Observed")

sim_durations <- final |>
  filter(group %in% c("P", "M", "B"), pupated) |>
  select(diet = group, L1:L4) |>
  pivot_longer(L1:L4, names_to = "stage", values_to = "days") |>
  mutate(source = "Simulated")

p_durations <- bind_rows(obs_durations, sim_durations) |>
  mutate(diet = factor(diet, levels = c("P", "M", "B")), stage = factor(stage, levels = larval_stages)) |>
  summarise(mean = mean(days), lo = quantile(days, 0.1), hi = quantile(days, 0.9), .by = c(diet, stage, source)) |>
  ggplot(aes(stage, mean, colour = diet, shape = source)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.4) +
  facet_wrap(~diet, labeller = labeller(diet = diet_labels)) +
  scale_colour_manual(values = diet_colours, guide = "none") +
  scale_shape_manual(values = c(Observed = 16, Simulated = 1)) +
  labs(
    x = NULL, y = "Days in instar", shape = NULL,
    title = "Instar durations of larvae that pupated",
    subtitle = "Points: mean; bars: 10th–90th percentiles"
  ) +
  theme_model()
ggsave("output/figures/development_durations_check.png", p_durations, width = 8, height = 3.8, dpi = 200)

# Food at molt ----------------------------------------------------------------

p_thresholds <- development$thresholds$data |>
  mutate(diet = factor(diet, levels = c("P", "M", "B"))) |>
  ggplot(aes(stage, units, colour = diet)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1),
             size = 1.3, alpha = 0.8) +
  geom_point(
    data = tibble(stage = factor(larval_stages), units = development$thresholds$threshold),
    aes(stage, units), inherit.aes = FALSE, shape = 95, size = 14, colour = "grey40"
  ) +
  scale_y_log10() +
  scale_colour_manual(values = diet_colours, labels = diet_labels) +
  labs(
    x = NULL, y = "Food eaten in instar (units, log scale)", colour = NULL,
    title = "Food eaten per completed instar",
    subtitle = sprintf("1 unit = one size-matched pea aphid; a bean kill counts %.2f (bean diet) or %.2f (mixed). Bars: fitted thresholds",
                       development$thresholds$v_bean_alone, development$thresholds$v_bean_mixed)
  ) +
  theme_model()
ggsave("output/figures/development_thresholds.png", p_thresholds, width = 8, height = 4.5, dpi = 200)
