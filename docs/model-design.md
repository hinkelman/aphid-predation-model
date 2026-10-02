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
| Run length | Several weeks eventually. Start with shorter runs **without aphid density dependence**; revisit carrying capacity later |
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
- No density dependence for now.

### Predator processes
1. **Search** on the current plant. Encounters occur at rate `a * N_species` (free parameter `a`, possibly stage-dependent).
2. **Capture**: success is species-specific and a **free parameter**. Pea aphids escape well (they drop from the plant); bean aphids are easy to catch. This trades off against prey quality: pea is high quality but hard to catch, bean low quality but easy.
3. **Rejection**: a captured aphid is rejected (not eaten) with probability logit⁻¹(−2.48 + 0.160 · age in days), from 2005. Species had no effect.
4. **Handling** (`R/predator_behavior.R`): log-logistic AFT ~ aphid × starvation, fitted to the 2008 trials (size-matched prey). For other prey ages, handling is scaled by (age / ref_age)^1.30, where the exponent comes from the 2005 trials (species had no effect once age was included). ref_age is the age of the 2008 prey: 6.5 d for bean (adult) and 4.2 d for pea (0.9 mg on the provisional mass curve).
5. **Partial consumption**: bean 85%, pea 0% (2008). The fraction eaten in a partial meal is unknown, so it is a free parameter.
6. **Food value** = aphid mass at age × fraction eaten. Mass at age is **provisional** (`R/aphid_mass.R`): exponential growth from a neonate (pea 0.15 mg from the literature, still unverified; bean scaled by the same neonate:adult ratio) to adult mass (pea 3.8 mg, bean 0.9 mg) at maturity (pea 7.5 d, bean 6.5 d).
7. **Post-handling / stay or leave**: Weibull AFT ~ aphid × starvation for the time to leave the plant after a meal, fitted to 2008 `Depart`. Each meal restarts the clock (a renewal process, as in the Cox analyses). Leaving starts travel to the other plant (travel time is a free parameter).
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
- **Encounters**: rate = (search_rate / plant_area) × Σ_species capture × N on the current plant. Redrawn whenever aphid numbers on that plant change, which is valid because the process is memoryless. The prey individual is chosen at random within the chosen species.
- **Satiation**: the gut empties exponentially (digestion_rate = 4/day, free). The larva attacks only when the gut has room for one size-matched prey (gut ≤ C − 1). C for each instar is set so that, with ad lib prey, the steady cycle of handling plus digestion pause gives the observed ad lib kill rate: C − 1 = e^(−k(T−h)) / (1 − e^(−kT)), with T = 1/(kills/day) and h = mean pea handling (15 min). Gut fill per kill = aphid mass / 0.9 for both species.
- **Hunger** for the behavior models = hours since the last meal, clamped to 2–24 h.
- **Bean handling with experience**: bean handling multiplier = m + (1 − m)·exp(−bean meals / learn_meals). Naive larvae follow the 2008 trials (bean-naive, pea-reared). m is calibrated so that simulated bean-diet vials match the 2008 L4 kill rate of ~25 bean/day (`analysis/04`). learn_meals = 20 (free; "slow learners").
- **Leaving a plant**: after each meal, the leave time comes from the post-handling Weibull model (species × starvation). On arrival, or with no meal yet, the larva leaves after an exponential giving-up time (mean 120 min, free). Travel takes 60 min (free).
- **Starvation**: death at lognormal time since the last meal (median L1 1.5, L2 2, L3 2.5, L4 4 d; sdlog 0.25; free). The clock pauses during handling. An L4 with food ≥ its critical value that has gone 24 h without a meal enters the pre-pupa after the fitted delay.
- **Mortality**: background plus L1 bean hazard, applied as a cumulative-hazard budget that is updated whenever the hazard changes.
- **Vial mode** (calibration and validation): prey are replaced as eaten and held at the size-matched age, there is no aphid demography, the larva doesn't leave, and the arena is 50 cm² with capture = 1.
- **Free parameters with provisional defaults**: search_rate 2 cm²/min, plant_area 400 cm², capture (pea 0.3, bean 0.8), digestion_rate, learn_meals, giving_up_time, travel_time, starvation medians, rejection_time. These should be explored by sensitivity analysis before drawing conclusions.
- **Not yet represented**: aphid density dependence, larval instar effects on capture of large prey and on handling (behavior data are L4 only), reduced search during post-bean inactivity bouts, and larval predators other than one larva.

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
- Aphid carrying capacity for longer runs.
- How Starve (a lab treatment) maps onto the model's internal hunger state. Currently: hours since the last meal, clamped to 2–24 h.
- Sensitivity analysis of the free parameters, especially search_rate/plant_area, capture, learn_meals, giving_up_time and travel_time.
- Whether `survival.csv` can be used to check per-encounter lethargy.
