[README.md](https://github.com/user-attachments/files/31686956/README.md)
# Dataset and analysis scripts for:

**"Nymphal habitat drydown prolongs development but preserves adult body size and reproductive strategies in the water strider *Gerris latiabdominis* (Hemiptera: Gerridae)"**

**Author:** Manabu Kishi  
**Journal:** Aquatic Insects (submitted)  
**Contact:** kishi@hotmail.co.jp  
**ORCID:** 0009-0001-8955-5189

---

## Contents

This repository contains three data files and three R analysis scripts supporting the results reported in the manuscript.

| File | Description |
| --- | --- |
| `Nymphal_period_data.csv` | Nymphal period (days) of *G. latiabdominis* reared under four drydown conditions |
| `Fecundity.csv` | Reproductive traits (preoviposition period and total number of eggs laid) of adult *G. latiabdominis* |
| `Body_weight_data.csv` | Fresh body weight measurements of adult *G. latiabdominis* after emergence |
| `GLM_analysis.R` | R script for all statistical analyses reported in Table 1 (GLMs) |
| `LMM_Fig2_analysis.R` | R script for LMM analysis comparing body weight trajectories between reproductive and non-reproductive females (pooled data; asterisks in Fig. 2) |
| `LMM_TableS1.R` | R script for LMM analyses of body weight trajectories in four categories separately by experimental group (Table S1) |

---

## Experimental groups

Four experimental groups differed in the substrate (water surface or wet paper) experienced during different nymphal instars:

| Group | Nymphal substrate |
| --- | --- |
| A (control) | Water surface throughout all instars |
| B | Wet paper during 1st–3rd instars; water surface during 4th–5th instars |
| C | Water surface during 1st–3rd instars; wet paper during 4th–5th instars |
| D | Wet paper throughout all instars |

All nymphs were reared individually in plastic cases (34 × 23.5 × 4.5 cm) under a photoperiod of 15.5L–8.5D at 20 ± 2°C. Four replicates of 45 nymphs per group were established, for a total of 180 nymphs per group.

---

## File descriptions

### Nymphal_period_data.csv

Nymphal period data in frequency format (number of individuals that completed development in each number of days).

| Column | Type | Description |
| --- | --- | --- |
| `Group` | character | Experimental group (A, B, C, or D) |
| `Nymphal_period_days` | integer | Number of days from hatching to adult emergence |
| `Number_of_individuals` | integer | Number of individuals with that nymphal period |

* Rows: 65 (one row per unique Group × Nymphal_period_days combination)
* Total individuals: A=161, B=155, C=146, D=136

### Fecundity.csv

Reproductive trait data for adult *G. latiabdominis* maintained as single pairs for 30 days after adult emergence. Each row represents one individual.

| Column | Type | Description |
| --- | --- | --- |
| `Group` | character | Experimental group (A, B, C, or D) |
| `Individual_No` | integer | Individual identifier (unique within each Group) |
| `Sex` | character | Sex of the individual (female or male) |
| `Reproductive_status_of_female` | character | Whether the female in the pair laid at least one egg during the 30-day observation period (yes or no). For males, this records the reproductive status of their paired female. |
| `Preoviposition_period` | numeric | Number of days from adult emergence to first oviposition (days). NA for non-reproductive females and all males. |
| `Total_number_of_eggs` | numeric | Total number of eggs laid during the 30-day observation period. NA for non-reproductive females and all males. |

* Rows: 320 (40 females + 40 males × 4 groups)
* Reproductive females per group: A=15, B=13, C=15, D=7

### Body_weight_data.csv

Fresh body weight measurements taken every 5 days after adult emergence (day 0 to day 30) for all individuals. Each row represents one measurement occasion for one individual. These data are used to generate Fig. 2 and Fig. S1, and to conduct the LMM analyses reported in Table S1 of the manuscript.

| Column | Type | Description |
| --- | --- | --- |
| `Group` | character | Experimental group (A, B, C, or D) |
| `Individual_No` | integer | Individual identifier (unique within each Group) |
| `Sex` | character | Sex of the individual (female or male) |
| `Reproductive_female_in_pair` | character | Whether the female in the pair laid at least one egg during the 30-day observation period (yes or no). For females, this reflects their own reproductive status; for males, this reflects the reproductive status of their paired female. |
| `Days_after_adult_emergence` | integer | Days after adult emergence at time of measurement (0, 5, 10, 15, 20, 25, or 30) |
| `Body_weight_mg` | numeric | Fresh body weight in milligrams (mg). NA if the individual died before the measurement occasion. |

* Rows: 2240 (40 females + 40 males × 4 groups × 7 time points)
* Non-NA measurements: 2024 (216 NA values due to mortality during the observation period)

### GLM_analysis.R

R script performing all generalized linear model (GLM) analyses reported in Table 1, including:

* Analysis 1: Nymphal period — Poisson GLM with log link; likelihood ratio test; Bonferroni-corrected pairwise comparisons (*emmeans*)
* Analysis 2: Preoviposition period — Gaussian GLM with identity link; F-test
* Analysis 3: Total number of eggs laid — Negative binomial GLM with log link; likelihood ratio test
* Analysis 4: Fresh body weight at adult emergence — Gaussian GLM with identity link (sex as covariate); F-test
* Analysis 5: Survival rate at adult emergence — Binomial GLM with logit link; likelihood ratio test; Bonferroni-corrected pairwise comparisons
* Analysis 6: Proportion of reproductive females — Binomial GLM with logit link; likelihood ratio test

**Required input files:** Nymphal_period_data.csv, Fecundity.csv, Body_weight_data.csv (all must be in the R working directory)

**Required R packages:** tidyverse 2.0.0, MASS 7.3-65, emmeans 2.0.1, multcomp 1.4-30

**R version:** 4.5.1 (R Core Team, 2025)

### LMM_Fig2_analysis.R

R script for the linear mixed model (LMM) analysis comparing fresh body weight trajectories between reproductive and non-reproductive females, using data pooled across all four experimental groups. Results are used to determine the asterisks shown in Fig. 2.

* Model: Body weight ~ Reproductive status × Days + (1|Individual_ID), AR(1) covariance
* Post-hoc Bonferroni-corrected pairwise contrasts at each time point (P < 0.0071; 7 time points tested)

**Required input file:** Body_weight_data.csv (must be in the R working directory)

**Required R packages:** tidyverse 2.0.0, nlme 3.1-169, emmeans 2.0.1

**R version:** 4.5.1 (R Core Team, 2025)

### LMM_TableS1.R

R script for the LMM analyses examining the effects of experimental group, days after adult emergence, and their interaction on fresh body weight, fitted separately for four categories. Results are reported in Table S1.

* Four categories: reproductive females, non-reproductive females, males with reproductive females, males with non-reproductive females
* Likelihood ratio tests for group main effect, days main effect, and group × days interaction

**Required input file:** Body_weight_data.csv (must be in the R working directory)

**Required R packages:** tidyverse 2.0.0, nlme 3.1-169

**R version:** 4.5.1 (R Core Team, 2025)

---

## Software

All analyses were performed in R version 4.5.1 (R Core Team, 2025) with the following packages:

* tidyverse 2.0.0
* MASS 7.3-65
* emmeans 2.0.1
* multcomp 1.4-30
* nlme 3.1-169

---

## Notes

* `Individual_No` values are unique within each Group but not across groups. To create a globally unique identifier, combine Group and Individual_No (e.g., `paste(Group, Individual_No, sep = "_")` in R).
* `Nymphal_period_data.csv` is in frequency format. To expand to individual-level data for GLM analysis, use `uncount(Number_of_individuals)` in R (tidyr package).
* `Body_weight_mg` values of zero are not present; NA indicates death prior to the measurement occasion.

---

## License

These data are made available under the Creative Commons Attribution 4.0 International License (CC BY 4.0). You are free to share and adapt the data, provided appropriate credit is given.

---

## Citation

If you use these data, please cite:

Kishi M (*year*) Nymphal habitat drydown prolongs development but preserves adult body size and reproductive strategies in the water strider *Gerris latiabdominis* (Hemiptera: Gerridae). *Aquatic Insects* doi:XXXXXX
