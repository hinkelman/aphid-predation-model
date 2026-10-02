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
| Space | Two fava bean plants. Both aphid species can occur on either plant; the initial composition is a scenario setting |
| Predators | One larva per run |
| Predator development | Dynamic: L1 → L4 → pupation, driven by what it eats |
| Run end | Predator pupates or dies, then the run continues to a fixed length so aphid dynamics play out |
| Run length | Several weeks eventually. Start with shorter runs **without aphid density dependence**; revisit carrying capacity later |
| Encounter rate | No data. Free parameter; explore by sensitivity analysis or use literature values for *H. convergens* |
| Code style | Modern R: dplyr, tidyr, ggplot2 |

## Data inventory and interpretation

All aphids were reared on fava beans.

### Aphid demography: `ClipCageData.xlsx` (Dec 2009)
- 32 bean and 28 pea aphids in clip cages, followed daily from L1 to death, with offspring counted each day.
- `Status2`: 1 = died, 0 = censored (escaped or lost). `Status` (f/m/o/d) was recorded for pea only.
- Bean start dates recorded as 2010-12-15 are a **typo for 2009-12-15**.
- An extra column notes "escaped while handling" (pea ID 31, which is censored).

### Predator diet performance: `diet_performance_II.xlsx` = `develop.csv`, `feed.csv`, `diet_summary.csv` (Sep 2008)
- Larvae reared from L1 on Bean (B), Pea (P) or Mixed (M) diets.
- `develop`: daily stage and the `Disabled` / `Dead` flags.
- `feed`: daily counts of aphids by status:
  - `fed` = number offered;
  - `eaten` = consumed (counted as **killed**);
  - `dead` = not clearly eaten or fed upon but dead, either killed and left uneaten or died on its own. **Excluded from consumption.**
- `diet_summary`: survival to adult (Eclose) is B 13%, M 45%, P 70% (this is `Survival.pdf`). Time to pupation is about 22 days on B and M and about 12 days on P.

### Diet switching: `DietTimingData.xlsx`
- From L4 day 1–4 (`Timing`), larvae were switched to bean (B) or no food (S). Trt P with `Timing` = NA are the controls.
- `First`…`Fourth`, `Pupa` = days in each stage. `T.weight` = mass at the start of treatment (g). `A.weight` = adult mass (g).
- The discontinued second experiment (`DietTimingDataII.csv`, not in `data/`) used B1/B3/S1/S3 = bean or starved for 1 or 3 days, plus Block/Clutch codes.

### Single-encounter behavior: `handle_depart_move.csv`, `intensive.csv`, `inactivity.csv` (2008)
- Arena: a fava leaf laid flat on an agar plate (`Agar` = new or reused plate). The larva was brought in on a piece of stem laid next to the leaf, and one aphid was placed in front of it to trigger an attack. Behavior after the attack was then recorded, including search on the leaf.
- `Handle` (min) with `HandleEnd` as the event flag. `Partial`, `Carry`, `Vomit`.
- `Depart` (min) = patch (leaf) residence time after handling, with `Left` as the event flag. `PRT.short` = time until first revisiting the stem or departing.
- `Starve` = hours without food before the trial. `Brush` and `Harass` are nuisance covariates left out of the final models.
- `intensive` and `inactivity`: bout durations in **seconds**. Video tracking was truncated at about 90 min (5394 s).

### Handling time vs aphid age: `HandleAge.csv` (2005)
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
2. **Attack / handling**: Weibull handling time by aphid species, adjusted for hunger (Starve) and aphid age (HandleAge). Bean aphids are more often partially consumed.
3. **Post-meal behavior**: inactivity and intensive search, with bean-induced lethargy.
4. **Stay or leave**: the hazard of leaving the plant depends on Aphid × Starve (Cox / parametric fits to `Depart`). Leaving starts travel to the other plant (travel time is a free parameter).
5. **Development**: each meal adds food value (bean < pea; partial meals count less). The larva molts when it reaches a stage-specific threshold, and the death hazard rises with poor nutrition or starvation.

### Calibration targets for predator development and survival
- Development time per stage and survival to pupation and to adult under B, M and P diets (`develop`, `diet_summary`).
- Daily consumption by age and diet (`feed`, eaten only).
- Responses to bean or starvation from L4 (`DietTimingData`).

## Planned repository layout
```
R/          functions: data import/cleaning, fitting helpers, samplers, DES engine
analysis/   scripts that fit each parameter set and plot fits against data (ggplot2)
sim/        scenario definitions and simulation runs
docs/       this design document
data/       raw data (read-only)
```

## Open items
- Parameterizing the encounter rate: literature values for *H. convergens*, or scenario ranges.
- Travel time and cost between plants.
- Aphid carrying capacity for longer runs.
- How Starve (a lab treatment) maps onto the model's internal hunger state.
- Whether `survival.csv` can be used to check per-encounter lethargy.
