# Two-plant aphid–coccinellid predation model: design

## Purpose

Explore indirect effects between two aphid species that share a predator but differ in their
suitability as prey. The predator is *Hippodamia convergens* larvae. The prey are the pea aphid
(*Acyrthosiphon pisum*, high quality) and the bean aphid (*Aphis fabae*, low quality / "toxic").
Two pathways are of interest:

- **Short-term (behavioral).** Long handling times, partial consumption, post-meal lethargy and
  quicker plant departure after eating bean aphids change predation pressure on pea aphids.
- **Longer-term (developmental).** A bean-heavy diet slows predator development and lowers its
  survival. That changes how many aphids the predator eats over time, and of which species.

## Decisions so far (2026-10-02)

| Topic | Decision |
|---|---|
| Model type | Continuous-time discrete-event simulation (DES), individual-based, clock in minutes |
| System | **One plant species (fava bean, *Vicia faba*), two aphid species (pea and bean aphid), one predator species (*Hippodamia convergens*).** Two fava plants, identical as hosts; both aphid species can occur on either plant, and the initial composition is a scenario setting. Phrases like "pea plant" mean "the fava plant holding pea aphids". |
| Predators | One larva per run |
| Predator development | Dynamic: L1 → L4 → pupation, driven by what it eats |
| Run end | Predator pupates or dies, then the run continues to a fixed length so aphid dynamics play out |
| Run length | Several weeks (default experiment: 40 days). **Aphid density dependence** added 2026-10-02 (see Aphid processes) |
| Encounter rate | No data. Free parameter; explore by sensitivity analysis or use literature values for *H. convergens* |
| Code style | Modern R: dplyr, tidyr, ggplot2 |

## Data inventory and interpretation

All aphids were reared on fava beans (*Vicia faba*), and all insects were kept at about 24 °C with a
16:8 photoperiod. *H. convergens* adults came from a commercial supplier; larvae were reared
individually from hatching (< 24 h) on excess pea aphids unless assigned to a diet treatment. Sources:
dissertation Ch. 1 (diet) and Ch. 2 (behavior), `~/Dropbox/Research/UNL/Dissertation/latex/`.

**Prey size matching.** Every predator experiment used size-matched prey: large (adult) bean aphids
and pea aphids of similar size, chosen by eye. Apterous adults weigh about 0.9 mg for bean and
about 3.8 mg for pea. So predator handling, consumption and development data all describe
**~0.9 mg prey**, while the simulated aphid populations include every age and size.

### Aphid demography: `ClipCageData.xlsx` (Dec 2009)
- 32 bean and 28 pea aphids in clip cages, followed daily from L1 to death, with offspring counted each day.
- `Status2`: 1 = died, 0 = censored (escaped or lost). `Status` (f/m/o/d) was recorded for pea only.
- Bean start dates recorded as 2010-12-15 are a **typo for 2009-12-15**.
- An extra column notes "escaped while handling" (pea ID 31, which is censored).

### Predator diet performance: `diet_performance_II.xlsx` = `develop.csv`, `feed.csv`, `diet_summary.csv` (Sep 2008)
- Larvae reared from L1 on Bean (B), Pea (P) or Mixed (M) diets.
- `develop`: daily stage and the `Disabled` / `Dead` flags.
- `feed`: daily counts of aphids by status (Ch. 1 methods):
  - `fed` = number supplied, adjusted daily to the previous day's kills;
  - `eaten` = **killed**: the number supplied minus live and dead aphids, i.e. aphids showing evidence of piercing;
  - `dead` = dead with no evidence of piercing (background mortality in the vial). **Excluded from consumption.**
- Development time to adult (Ch. 1, Table 1): bean 26.3 d, mixed 26.3 d, pea 17.0 d. Adult mass: 10.9, 12.1 and 21.5 mg.
- On the mixed diet, larvae killed significantly **more bean than pea** aphids over the larval period (sign test, p = 0.024).
- `diet_summary`: survival to adult (Eclose) is B 13%, M 45%, P 70% (this is `Survival.pdf`). Time to pupation is about 22 days on B and M and about 12 days on P.

### Diet switching: `DietTimingData.xlsx`
- From L4 day 1–4 (`Timing`), larvae were switched to bean (B) or no food (S). Trt P with `Timing` = NA are the controls.
- `First`…`Fourth`, `Pupa` = days in each stage. `T.weight` = mass at the start of treatment (g). `A.weight` = adult mass (g).
- The discontinued second experiment (`DietTimingDataII.csv`, not in `data/`) used B1/B3/S1/S3 = bean or starved for 1 or 3 days, plus Block/Clutch codes.

### Single-encounter behavior: `handle_depart_move.csv`, `intensive.csv`, `inactivity.csv` (2008)
- Predators: **4th-instar** larvae, 1 or 2 days into the instar (`Age`), starved 2–24 h before the trial.
- Arena: a compound fava leaf (2 leaflets, about 48 cm²) on an agar plate (`Agar` = new or reused plate). The larva was moved in on a piece of stem or a brush. Prey was an adult bean aphid or a size-matched pea aphid placed at its mouthparts.
- **Handling** runs from when the larva secures the aphid until it **moves away from the feeding site**. **Post-handling time** (`Depart`) runs from the end of handling until the larva is no longer touching the leaf or stem. Patch residence time = handling + post-handling.
- Partial consumption: pea 0%, bean 85%.
- `Handle` (min) with `HandleEnd` as the event flag. `Partial`, `Carry`, `Vomit`.
- `Depart` (min) = patch (leaf) residence time after handling, with `Left` as the event flag. `PRT.short` = time until first revisiting the stem or departing.
- `Starve` = hours without food before the trial. `Brush` and `Harass` are nuisance covariates left out of the final models.
- `intensive` and `inactivity`: bout durations in **seconds**. Video tracking was truncated at about 90 min (5394 s).

### Handling time vs aphid age: `HandleAge.csv` (2005)
- **Lower trust than the 2008 data.** These were exploratory trials, and predator instar and condition can't be confirmed. Where the two conflict, 2008 wins. 2005 is used only where 2008 has no information: the *relative* effect of aphid age on handling, and rejection probability.
- Handling time in **minutes**, by aphid age in days (1–12). Covers partial consumption, rejection and censoring.

### Aphid vials: `survival.csv` and `Aphid Survival Data.xls`
- Vials of 20 aphids, probably over about 20 h. Most likely meant to show that predators can eat fewer bean aphids because each one induces lethargy. The predator's presence and the exact protocol are uncertain, so this is **not used for fitting**; at most it is a qualitative check.

## Model structure

### Entities and state
- **Plant** (`1`, `2`): the aphids on it.
- **Aphid**: id, species, plant, birth time, next event (birth or death) and its time.
- **Predator**: stage (L1–L4), accumulated food value, gut / time since last meal, current plant, current behavioral state.

### Aphid processes (as in the 2011 NetLogo model, `GrowthDevelopment0.2.nlogo`)
- **Death**: Weibull lifespan, sampled conditional on current age:
  `t = scale * ((age/scale)^shape - log(U))^(1/shape) - age`.
- **Birth**: a non-homogeneous Poisson process from a piecewise-linear cumulative intensity function of age, sampled by inversion (Leemis 2004).
- **Individual fecundity**: each aphid draws a multiplier m ~ Gamma(mean 1, SD about 0.11 for pea and 0.27 for bean) at birth (not inherited). It scales the birth intensity: the next birth is where the CIF has risen by E/m. The SD comes from a negative binomial fit of lifetime offspring, offset by log CIF(lifespan). Without it, simulated lifetime fecundity varied only about half as much as observed.
- Pea aphid 23 (lived 24 d, never reproduced) is excluded from all fits as an abnormal individual.
- Each aphid's next event is whichever of death and birth comes first. Newborns appear on their mother's plant.
- **Density dependence** (after the ALMaSS aphid model, [Thomsen, Duan & Topping 2024](https://fesmj.pensoft.net/article/123747/)):
  - **Delayed crowding mortality**: an extra death hazard ln(1 + k·x) per day, k = 0.012 (pea), 0.02 (bean) per aphid/g (ALMaSS Table 5). x = the species' crowding density `density_lag` days earlier (daily census). The hazard changes only at day boundaries, when every aphid's crowding death time is redrawn; this is exact for a piecewise-constant hazard. ALMaSS's density-independent survivorship (Sa = 0.9) is not used: it represents field predators, and clip-cage lifespans set baseline mortality here.
  - **Winged emigration**: a newborn becomes winged with probability (2.603·x + 0.847·GS − 27.19)% (Carter 1982 via ALMaSS) at the current crowding density, with plant growth stage GS = 3.5. Winged aphids live (and can be eaten) on the plant but don't reproduce there, and leave at maturity. Emigrants are counted in the output.
  - **Crowding density** for a species = (own count + α × other species' count) / plant biomass (g fresh). α = `competition_alpha` (default 0.5): α = 0 is ALMaSS-like (no direct competition) and α = 1 counts all aphids equally. Without predators on a shared plant, α = 1 lets pea aphids exclude bean aphids (bean has the higher k), while α = 0.5 gives coexistence (~530 pea, ~140 bean).
  - **Lag 1 day**, not ALMaSS's 4: on one plant, with growth ~0.38/day, lags of 2–4 days gave boom-bust cycles (pea 40 ↔ 1,700). ALMaSS works on 10 × 10 m field cells with background mortality. With a 1-day lag, colonies level off (pea alone ~600, bean alone ~400 on a 20 g plant) with no late decline.
  - **Plant biomass** 20 g fresh: ~370 cm² leaf area per side ÷ specific leaf area 25.7 mm²/mg ≈ 1.4 g dry leaf, plus stem, at ~90% water.
  - Replaced earlier versions: births-only thinning against a carrying capacity, without and then with a 1-day averaged lag. A colony at capacity replaced losses almost instantly and then cycled once its age structure synchronized.
- No density dependence for now.

### Predator processes
1. **Search** on the current plant (see Simulation engine: encounters are built from 2008 tracking data, body size and plant area, with aphid aggregation).
2. **Capture**: success = species capture (free) × a prey-size factor (small larvae cannot subdue large aphids). Pea aphids escape well (they drop from the plant); bean aphids are easy to catch. This trades off against prey quality: pea is high quality but hard to catch, bean low quality but easy.
3. **Rejection**: a captured aphid is rejected (not eaten) with probability logit⁻¹(−2.48 + 0.160 · age in days), from 2005. Species had no effect.
4. **Handling** (`R/predator_behavior.R`): log-logistic AFT ~ aphid × starvation, fitted to the 2008 trials (size-matched prey). For other prey ages, handling is scaled by (age / ref_age)^1.30, where the exponent comes from the 2005 trials (species had no effect once age was included). ref_age is the age of the 2008 prey: 6.5 d for bean (adult) and 4.2 d for pea (0.9 mg on the provisional mass curve).
5. **Partial consumption**: bean 85%, pea 0% (2008). The fraction eaten in a partial meal is unknown, so it is a free parameter.
6. **Food value** = aphid mass at age × fraction eaten. Mass at age is **provisional** (`R/aphid_mass.R`): exponential growth from a neonate (pea 0.15 mg from the literature, still unverified; bean scaled by the same neonate:adult ratio) to adult mass (pea 3.8 mg, bean 0.9 mg) at maturity (pea 7.5 d, bean 6.5 d).
7. **Post-handling / stay or leave**: Weibull AFT ~ aphid × starvation for the time to leave the plant after a meal, fitted to 2008 `Depart`. Each meal restarts the clock (a renewal process, as in the Cox analyses). On a plant, leave times are scaled up from leaf to plant by area, and the clock pauses while the larva is satiated. Leaving starts travel to the other plant.
8. **Starvation** in the model = hours since the last meal, clamped to the trial range of 2–24 h when predicting handling and post-handling times.
9. **Development** (`R/predator_development.R`; fitted in `analysis/03_predator_development.R`). Food is counted in **units**: one size-matched pea aphid (0.9 mg) fully eaten. An aphid of mass m is worth m/0.9 units if pea and v·m/0.9 if bean. This calibrated v replaces the separate partial-consumption fraction in item 6, so partial consumption isn't counted twice. Fitted from each larva's own observed daily kills in the diet experiment:
   - **Thresholds** per instar: L1 9.0, L2 16.9, L3 31.1, L4 100.4 units. Individual log-threshold SD in the simulation: 0.2.
   - **Bean value per kill**: v = 0.50 when bean is the only food, but **0.15 when pea is also eaten**. Mixed-diet larvae killed many bean aphids and gained little from them. Rule for the simulation: v = 0.15 if pea was eaten in the previous 24 h, otherwise 0.50.
   - **Molting** happens when an instar's food reaches its threshold. L4 then enters a **1.5-day non-feeding pre-pupa** and pupates; the predator leaves the simulation at pupation.
   - **Critical L4 food**: a larva with L4 food ≥ about 25 units (individual value ~ logistic with location 25.3 and scale 8.8) pupates even if it stops eating. The delay to pupation is 4.12 − 0.030·food days. Lower-food larvae that stop eating die.
   - **Gut-limited intake**: maximum kills/day per instar (pea diet: L1 3.7, L2 9.2, L3 16.2, L4 40.7, using the first 2 days of L4) cap intake through a satiation/digestion rule in the full model. Kill rates were similar on bean, so the gut fills per aphid **killed**, not per unit of value.
10. **Mortality** (daily hazard, cloglog model on the check intervals):
   - Baseline: L1 0.035, L2–L3 0.014, L4 0.051, pupa 0.005 per day.
   - **L1 bean toxicity**: hazard × (1 + bean/(1 + pea))^0.45, using kills since hatching. This is "amount eaten, buffered by pea", as chosen.
   - **Pupation failure**: a larva reaching the L4 threshold pupates with probability (1 + bean/(1 + pea))^−0.15. γ = 0.15 was calibrated by simulating the diet and diet-timing experiments. Bean-fed 4th instars often ate well past the threshold and still died ("attempting to pupate").
   - **Starvation**: no time-to-death data ("died as larva" only). Free parameter, still to be set when the full model is built.

#### Development fit: known misfits (see `output/figures/development_*`)
- Mixed-diet L1 is too fast in the simulation (3.2 d vs 5.3 d observed). Mixed-diet L1 larvae needed more units than the common threshold.
- Diet-timing groups switched late (L4 day 3–4) pupated 100% in the experiment vs 82–94% in the simulation. That experiment had lower background mortality than the 2008 diet experiment that the hazards come from.
- Larvae switched to bean late in L4 pupated within 1–1.6 d, as fast as starved larvae, but the simulation has them keep eating bean to the threshold (3 d). Starved from L4 day 2: 61% pupate in the simulation vs 82% observed.
- **All larvae in every experiment were reared on pea aphids** before use, so the 2005 and 2008 behavior trials both used bean-naive larvae (a known criticism of the design: neither prey was novel in the same way). Rearing history therefore does not explain the 2005 vs 2008 bean handling difference. The diet-experiment larvae on B and M diets, however, ate bean from hatching.
- Separately, **bean handling time vs vial kill rates**: bean-diet L4s killed ~25 bean aphids/day. At the 2008 handling times (~70+ min each) that would take more than 24 h a day; the 2005 handling times (~15–20 min) fit easily. Larvae with repeated bean experience may handle bean faster than the pea-reared larvae in the 2008 trials.

#### Known issues with the behavior parameters
- **Bean handling differs between years — mostly a starvation effect.** At matched starvation (≤ 4.5 h, the 2005 range), 2008 bean handling has a median of 52 min (IQR 40–70, n = 7), against 20 min (IQR 11–51, n = 31) in 2005 for adult-sized bean aphids. Across all 2008 trials (2–24 h) the median is 126 min. So most of the apparent conflict came from comparing hungrier 2008 larvae with recently fed 2005 larvae. The 2008 model captures the effect through its aphid × starvation term. The remaining ~2.6× gap rests on 7 low-starvation 2008 trials.
- (Earlier note) **Bean handling differs between years.** The 2005 and 2008 trials defined handling the same way (until the larva moved away from the feeding site). Pea handling agrees between them, but for adult-sized bean aphids 2005 gives about 15–20 min and 2008 about 70–130 min. The model uses 2008, which is larger, documented (4th instar, video) and the source of the lethargy results. Possible causes: predator instar or history in 2005 (undocumented), or a change in the bean aphid culture.
- Behavior data come only from 4th instars. Handling for earlier instars needs a scaling assumption (to be set when development is modeled).

### Calibration targets for predator development and survival
- Development time per stage and survival to pupation and to adult under B, M and P diets (`develop`, `diet_summary`).
- Daily consumption by age and diet (`feed`, eaten only).
- Responses to bean or starvation from L4 (`DietTimingData`).

## Simulation engine (`R/simulation.R`, `R/params.R`, `R/scenarios.R`)

- **Clock and events.** Continuous time in minutes. Aphid next-event times are held in vectors and the earliest is found with `which.min()`. The predator is a state machine with its own event schedule (encounter, rejection end, handling end, gut ready, leave, arrive, death, starvation, starvation check, pupation). Aphid counts are recorded at a census each day.
- **Predator states**: search → (encounter) → rejecting | handling → search or satiated; travel between plants; pre-pupa; done.
- **Search rate** (area searched per minute) = speed × activity × detection width × (body length / L4 length)².
  - speed = 114 mm/min, the median while moving, from 2008 tracking;
  - activity = proportion of time moving, from a 2008 quasi-binomial fit ~ aphid × starvation. It is ~0.4 after a pea meal, and 0.32 → 0.17 after a bean meal as starvation goes 2 → 24 h (**lethargy**). It is set at each meal and applies until the next;
  - detection width = 3.5 mm for an L4 (body width 2.1 mm plus ~0.7 mm detection distance each side);
  - body length L1 1.9, L2 2.95, L3 4.57, L4 7.0 mm, so an L1 searches ~7% of the L4 area per minute.
- **Aphid aggregation**: colony area = colony_min_area (10 cm²) + aphid density (adult-mass equivalents) / colony_density (4 mg/cm²), capped at the plant area. A larva that has found the colony (it hatched beside it, or has eaten on this plant) searches only the colony area. On arriving at a plant it searches the whole plant (2,500 cm²) until it eats.
- **Encounters**: the larva meets aphids at rate (search rate / search space) × number of aphids on its plant, and each encounter is with a random aphid. The rate is redrawn whenever aphid numbers on that plant change (valid because the process is memoryless).
- **Capture per encounter**: pea 0.1, bean 0.8. Adult *H. convergens* consumed 1 of 72 pea aphids encountered in alfalfa field arenas ([Nelson & Rosenheim 2006](https://rosenheim.faculty.ucdavis.edu/wp-content/uploads/sites/137/2014/09/Nelson-and-Rosenheim-EEA-2006.pdf)), and ladybirds trigger dropping by pea aphids >3× as often as damsel bugs (Losey & Denno 1998). There are no data for larvae or for bean aphids, which rarely drop. Dixon (1959): capture improves with larval instar and is better on young aphids; 22 of 50 unfed *Adalia* first instars starved before their first capture.
- **Attack outcome**: success = capture[species] × 1 / (1 + (prey length / (prey_size_ratio × larval length))^4), with aphid length = 2.4 × mass^(1/3) mm. prey_size_ratio = 2 so that an L1 vs 0.9 mg prey gets ~0.9, as diet-experiment L1s killed ~3.7 size-matched aphids/day. A failed attack or a rejection (2005 model) costs 1 min, and the aphid survives.
- **Satiation**: the gut empties exponentially (digestion_rate = 4/day, free). The larva attacks only when the gut has room for one size-matched prey (gut ≤ C − 1). C for each instar is set so that, with ad lib prey, the steady cycle of handling plus digestion pause gives the observed ad lib kill rate: C − 1 = e^(−k(T−h)) / (1 − e^(−kT)), with T = 1/(kills/day) and h = mean pea handling (15 min). Gut fill per kill = aphid mass / 0.9 for both species.
- **Hunger** for the behavior models = hours since the last meal, clamped to 2–24 h.
- **Bean handling with experience**: bean handling multiplier = m + (1 − m)·exp(−bean meals / learn_meals). Naive larvae follow the 2008 trials (bean-naive, pea-reared). m is calibrated so that simulated bean-diet vials match the 2008 L4 kill rate of ~25 bean/day (`analysis/04`). learn_meals = 20 (free; "slow learners").
- **Leaving a plant**: the 2008 post-handling times (Weibull, species × starvation) are for leaving one 48 cm² leaf. Walking off a leaf on a plant leads to another leaf of the same plant, so leaving the plant = giving up on n_leaves = plant_area / trial_leaf_area × leave_scaling (≈ 52) leaves in turn. The leave time is the **sum of n_leaves leaf-level draws**, redrawn after each meal.
  - On arrival at a plant, the pea model at current hunger applies; a hatchling doesn't leave before its first meal.
  - The leave clock pauses while the larva is satiated.
  - Two rejected versions: multiplying one draw by 52 kept the fitted early-departure spike (Weibull shape 0.7), and resetting the clock on every encounter made it worse. In both, larvae walked off dense colonies within hours of a meal and starved on the empty plant.
  - **Consequence**: a larva on a growing colony essentially never leaves (summed times are days). Between-plant effects therefore need conditions that make a larva leave: a depleted or small colony, a poor starting plant, or several larvae.
- **Travel**: the plants are separate, so the larva walks down one, across the soil and up the other (plant_path = 1.2 m). Travel time = path / (speed × activity), ~25–60 min.
- **Starvation**: death at lognormal time since the last meal (median L1 1.5, L2 2, L3 2.5, L4 4 d; sdlog 0.25; free). The clock pauses during handling. An L4 with food ≥ its critical value that has gone 24 h without a meal enters the pre-pupa after the fitted delay.
- **Mortality**: background plus L1 bean hazard, applied as a cumulative-hazard budget that is updated whenever the hazard changes.
- **Vial mode** (calibration and validation): prey are replaced as eaten and held at the size-matched age, there is no aphid demography, the larva doesn't leave, and the arena is 50 cm² with capture = 1.
- **Default plant and scenario**: a medium fava plant (~3–4 weeks, ~30 cm, 5–6 compound leaves): 800 cm² searchable (both leaf surfaces plus stem), ~20 g fresh biomass, plants 0.9 m apart by walking path. Each colony is founded by **10 aphids of mixed ages** drawn from the species' stable age distribution (Euler–Lotka on the fitted life table: λ ≈ 1.46/day pea, 1.44/day bean; ~7% adults), giving ~30–45 aphids by day 3 (as 2 adult founders did). Ages are drawn inside each run with the run's seed. The larva hatches on **day 3** beside the plant-1 colony; runs last 35 days.
  - Mixed-age founders soften, but don't remove, the predator-free decline after the colony reaches capacity (pea at day 35: ~290 vs ~190 with 2 adult founders). The rest comes from births-only density dependence: births nearly stop at capacity, so the colony ages together and cycles with a period near the aphid lifespan. Removing it would need crowding to act on deaths or emigration (alates), which isn't modelled.
  - **Why these values** (exploratory runs, 2026-10-03). On 2,500 cm² plants a larva hatching on day 7 never depleted a colony and never left it: giving up on ~52 leaves takes days. Depletion needs the larva to arrive during early colony growth. Hatching on day 3 (~20 aphids) cut pea aphid-days by ~46% and wiped out the colony in ~40% of runs. Hatching on day 4–5 had almost no effect, and on day 0 hatchlings mostly starved.
  - On an 800 cm² plant (~17 leaves to give up on), ~15% of larvae left for the other plant about a day after the colony died out. On 1,200 cm² none did.
  - The window is narrow, so hatch day is a factor in the sensitivity analysis.
- **Assumed (no data) parameters**: detection width, plant area, colony_min_area, colony_density, prey_size_ratio and size steepness, capture (pea 0.3, bean 0.8), leave_scaling, plant_path, digestion_rate, learn_meals, starvation medians, rejection_time, plant_biomass, competition_alpha, density_lag. They are explored by sensitivity analysis.
- **Not yet represented**: alate production / aphid emigration, instar effects on handling time (behavior data are L4 only), separate colonies per aphid species on a shared plant, and more than one larva.

## Sensitivity analysis (`R/sensitivity.R`, `analysis/06_sensitivity.R`)

- **Method**: Morris elementary-effects screening with `sensitivity::morris()`: 4 levels, grid jump 2. Effects are scaled by factor range (`scale = TRUE`). Factors with wide ranges are sampled on a log10 scale.
- **Factors and ranges** (15; `free_parameter_ranges()`): detection width 2–5 mm; plant area 400–1,600 cm² (log); colony minimum area 2–40 cm² (log); colony packing 1–10 mg/cm² (log); prey size limit 1–4 (log); capture pea 0.02–0.3 (log); capture bean 0.5–1.0; leave-time scaling ×0.25–4 (log); path between plants 0.45–1.8 m (log); digestion 2–8/day (log); bean learning 5–100 meals (log); starvation tolerance ×0.5–2 (log); plant biomass 10–40 g (log); competition α 0–1; density lag 0–3 d (whole days); hatch day 2–5. 10 trajectories → 170 parameter sets.
- **Scenarios** (35 days; 10 mixed-age founders per colony; larva hatches on plant 1 on the hatch-day factor): Pea | pea, Pea | bean, Pea + bean | none. 20 replicate runs each per parameter set, with common random numbers (seeds depend only on scenario and replicate). No-predator baselines depend only on plant biomass, α and density lag, and are run once per combination in the design (20 runs each).
- **Outputs**: suppression of pea aphid-days on plant 1 in each scenario; indirect effects ie_separate = supp(Pea | bean) − supp(Pea | pea) and ie_shared = supp(Pea + bean | none) − supp(Pea | pea), the latter including direct competition; pupation probability in each scenario.

## Repository layout
```
R/          functions: data import/cleaning, fitting helpers, samplers, DES engine, scenarios
analysis/   numbered scripts: fit each parameter set, calibrate, validate, run experiments
docs/       this design document
data/       raw data (read-only)
```

## Open items
- Parameterizing the encounter rate: literature values for *H. convergens*, or scenario ranges.
- Travel time and cost between plants.
- Aphid crowding: plant biomass, competition α and the density lag have no direct data (ALMaSS values are field-scale); all three are in the sensitivity analysis. Colony packing density has no published values.
- How Starve (a lab treatment) maps onto the model's internal hunger state. Currently: hours since the last meal, clamped to 2–24 h.
- Interpret the Morris screening (analysis/06) and decide which free parameters need data, literature values, or a finer (e.g. variance-based) analysis.
- Whether `survival.csv` can be used to check per-encounter lethargy.
