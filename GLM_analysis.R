# ============================================================================
# Statistical Analysis for:
# "Larval habitat drydown prolongs development but preserves adult body size and
#  reproductive strategies in the water strider Gerris latiabdominis (Hemiptera: Gerridae)"
#
# Author  : Manabu Kishi
# Journal : Aquatic Insects (submitted)
# Date    : 2026-04-29
#
# Description:
#   This script performs the statistical analyses reported in Table 1 of the
#   above manuscript. Generalized linear models (GLMs) with appropriate error
#   distributions are used to compare life-history traits (larval period,
#   preoviposition period, total number of eggs laid, and fresh body weight at
#   adult emergence) among four experimental groups (A-D) differing in larval
#   substrate moisture conditions.
#
# Input files (must be in the working directory):
#   - Larval_period_data.csv
#   - Fecundity.csv
#   - Body_weight_data.csv
#
# Software:
#   R version 4.5.1 (R Core Team, 2025)
#   Packages: tidyverse, MASS 7.3-65, emmeans 2.0.1, multcomp 1.4-30, nlme 3.1-169, tidyverse 2.0.0
# ============================================================================


# ============================================================================
# 1. Package Loading
# ============================================================================

packages <- c("tidyverse", "MASS", "readxl", "emmeans", "multcomp")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


# ============================================================================
# 2. Data Loading
# ============================================================================

cat("Loading data files...\n\n")

larval_data <- read_csv("Larval_period_data.csv")

fecundity_data <- read_csv("Fecundity.csv") %>%
  mutate(Group = factor(Group, levels = c("A", "B", "C", "D")))

weight_data <- read_csv("Body_weight_data.csv") %>%
  mutate(Group = factor(Group, levels = c("A", "B", "C", "D")))

cat("Data files loaded successfully.\n\n")


# ============================================================================
# 3. Data Preparation
# ============================================================================

# Expand frequency-format larval period data to individual-level data
larval_expanded <- larval_data %>%
  uncount(Number_of_individuals) %>%
  mutate(Group = factor(Group, levels = c("A", "B", "C", "D")))

# Subset to reproductive females only (for preoviposition period and egg analyses)
repro_females <- fecundity_data %>%
  filter(Sex == "female", Reproductive_status_of_female == "yes")


# ============================================================================
# ANALYSIS 1: Larval Period
# GLM with Poisson distribution and log link
# ============================================================================

cat("========================================\n")
cat("ANALYSIS 1: Larval Period (Poisson GLM)\n")
cat("========================================\n\n")

larval_summary <- larval_expanded %>%
  group_by(Group) %>%
  summarise(
    n    = n(),
    mean = mean(Larval_period_days),
    SD   = sd(Larval_period_days),
    .groups = "drop"
  )
print(larval_summary)

model_larval <- glm(
  Larval_period_days ~ Group,
  data   = larval_expanded,
  family = poisson(link = "log")
)

cat("\nGLM Summary:\n")
print(summary(model_larval))

model_null_larval <- glm(
  Larval_period_days ~ 1,
  data   = larval_expanded,
  family = poisson(link = "log")
)

lr_larval <- anova(model_null_larval, model_larval, test = "Chisq")
cat("\nLikelihood Ratio Test (Group effect):\n")
print(lr_larval)


# Post-hoc pairwise comparisons (Bonferroni correction)

emm_larval      <- emmeans(model_larval, ~ Group)
pairwise_larval <- contrast(emm_larval, method = "pairwise", adjust = "bonferroni")

cat("\nEstimated Marginal Means:\n")
print(emm_larval, type = "response")

cat("\nPairwise Comparisons (Bonferroni-corrected):\n")
print(pairwise_larval, type = "response")

# Letter assignment based on pairwise results
# All six pairwise comparisons are significant (p < 0.05) after Bonferroni correction;
# therefore each group receives a unique letter.
letters_larval <- c(A = "a", B = "b", C = "c", D = "d")

larval_table <- larval_summary %>%
  mutate(
    letter           = letters_larval[as.character(Group)],
    Mean_with_letter = paste0(sprintf("%.2f", mean), letter)
  )

cat("\nSummary for Table 1 (larval period):\n")
print(larval_table[, c("Group", "n", "mean", "SD", "letter")], row.names = FALSE)


# ============================================================================
# ANALYSIS 2: Preoviposition Period
# GLM with Gaussian distribution and identity link
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS 2: Preoviposition Period (Gaussian GLM)\n")
cat("========================================\n\n")

preov_summary <- repro_females %>%
  drop_na(Preoviposition_period) %>%
  group_by(Group) %>%
  summarise(
    n    = n(),
    mean = mean(Preoviposition_period),
    SD   = sd(Preoviposition_period),
    .groups = "drop"
  )
print(preov_summary)

model_preov <- glm(
  Preoviposition_period ~ Group,
  data      = repro_females,
  family    = gaussian(link = "identity"),
  na.action = na.omit
)

cat("\nGLM Summary:\n")
print(summary(model_preov))

model_null_preov <- glm(
  Preoviposition_period ~ 1,
  data      = repro_females,
  family    = gaussian(link = "identity"),
  na.action = na.omit
)

f_preov <- anova(model_null_preov, model_preov, test = "F")
cat("\nF-test (Group effect):\n")
print(f_preov)


# ============================================================================
# ANALYSIS 3: Total Number of Eggs Laid
# GLM with Negative Binomial distribution (to account for overdispersion)
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS 3: Total Number of Eggs (Negative Binomial GLM)\n")
cat("========================================\n\n")

eggs_summary <- repro_females %>%
  drop_na(Total_number_of_eggs) %>%
  group_by(Group) %>%
  summarise(
    n        = n(),
    mean     = mean(Total_number_of_eggs),
    variance = var(Total_number_of_eggs),
    .groups  = "drop"
  )
print(eggs_summary)

model_eggs <- glm.nb(
  Total_number_of_eggs ~ Group,
  data      = repro_females,
  na.action = na.omit
)

cat("\nGLM Summary:\n")
print(summary(model_eggs))

model_null_eggs <- glm.nb(
  Total_number_of_eggs ~ 1,
  data      = repro_females,
  na.action = na.omit
)

lr_eggs <- anova(model_null_eggs, model_eggs, test = "Chisq")
cat("\nLikelihood Ratio Test (Group effect):\n")
print(lr_eggs)


# ============================================================================
# ANALYSIS 4: Fresh Body Weight at Adult Emergence
# GLM with Gaussian distribution; Sex included as a covariate
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS 4: Fresh Body Weight at Emergence (Gaussian GLM)\n")
cat("========================================\n\n")

weight_emergence <- weight_data %>%
  filter(Days_after_adult_emergence == 0) %>%
  drop_na(Body_weight_mg) %>%
  mutate(
    Group = factor(Group, levels = c("A", "B", "C", "D")),
    Sex   = factor(Sex,   levels = c("female", "male"))
  )

weight_summary <- weight_emergence %>%
  group_by(Group, Sex) %>%
  summarise(
    n    = n(),
    mean = mean(Body_weight_mg),
    SD   = sd(Body_weight_mg),
    .groups = "drop"
  )
print(weight_summary)

model_weight <- glm(
  Body_weight_mg ~ Group + Sex,
  data   = weight_emergence,
  family = gaussian(link = "identity")
)

cat("\nGLM Summary:\n")
print(summary(model_weight))

model_null_weight <- glm(
  Body_weight_mg ~ Sex,
  data   = weight_emergence,
  family = gaussian(link = "identity")
)

f_weight <- anova(model_null_weight, model_weight, test = "F")
cat("\nF-test (Group effect, controlling for Sex):\n")
print(f_weight)


# ============================================================================
# SUMMARY: Statistical Test Results (Table 1)
# ============================================================================

cat("\n========================================\n")
cat("SUMMARY: Statistical Test Results for Table 1\n")
cat("========================================\n\n")

lr_eggs_stat <- 2 * (as.numeric(logLik(model_eggs)) - as.numeric(logLik(model_null_eggs)))
lr_eggs_p    <- 1 - pchisq(lr_eggs_stat, df = 3)

summary_table <- data.frame(
  Trait      = c("Larval period", "Preoviposition period",
                 "Total number of eggs", "Fresh body weight"),
  Family     = c("Poisson", "Gaussian", "Negative Binomial", "Gaussian"),
  Test       = c("LR chi-sq", "F", "LR chi-sq", "F"),
  Statistic  = c(
    sprintf("chi2 = %.2f", lr_larval$Deviance[2]),
    sprintf("F = %.4f",    f_preov$F[2]),
    sprintf("chi2 = %.2f", lr_eggs_stat),
    sprintf("F = %.4f",    f_weight$F[2])
  ),
  df         = c(
    "3, 594",
    sprintf("3, %d", nrow(repro_females %>% drop_na(Preoviposition_period)) - 4),
    sprintf("3, %d", nrow(repro_females %>% drop_na(Total_number_of_eggs))  - 4),
    sprintf("3, %d", nrow(weight_emergence) - 5)  # 5 parameters: intercept + 3 group + sex
  ),
  P_value    = c(
    format(lr_larval$`Pr(>Chi)`[2], digits = 3),
    format(f_preov$`Pr(>F)`[2],    digits = 3),
    format(lr_eggs_p,               digits = 3),
    format(f_weight$`Pr(>F)`[2],   digits = 3)
  ),
  stringsAsFactors = FALSE
)

print(summary_table, row.names = FALSE)

cat("\nAnalysis complete.\n")


# ============================================================================
# ANALYSIS 5: Survival Rate at Adult Emergence
# GLM with Binomial distribution and logit link
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS 5: Survival Rate (Binomial GLM)\n")
cat("========================================\n\n")

# Survival counts per group (from larval period data)
survival_data <- read_csv("Larval_period_data.csv") %>%
  group_by(Group) %>%
  summarise(n_survived = sum(Number_of_individuals), .groups = "drop") %>%
  mutate(
    Group    = factor(Group, levels = c("A","B","C","D")),
    n_total  = 180,
    n_failed = n_total - n_survived
  )

cat("Survival summary:\n")
print(survival_data)

model_surv <- glm(
  cbind(n_survived, n_failed) ~ Group,
  data   = survival_data,
  family = binomial(link = "logit")
)

model_null_surv <- glm(
  cbind(n_survived, n_failed) ~ 1,
  data   = survival_data,
  family = binomial(link = "logit")
)

cat("\nGLM Summary:\n")
print(summary(model_surv))

lr_surv <- anova(model_null_surv, model_surv, test = "Chisq")
cat("\nLikelihood Ratio Test (Group effect):\n")
print(lr_surv)

# Post-hoc pairwise comparisons (Bonferroni)
emm_surv      <- emmeans(model_surv, ~ Group)
pairwise_surv <- contrast(emm_surv, method = "pairwise", adjust = "bonferroni")
cat("\nPairwise Comparisons (Bonferroni-corrected):\n")
print(pairwise_surv, type = "response")


# ============================================================================
# ANALYSIS 6: Proportion of Reproductive Females
# GLM with Binomial distribution and logit link
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS 6: Proportion of Reproductive Females (Binomial GLM)\n")
cat("========================================\n\n")

# Reproductive female counts (n = 40 per group)
repro_data <- read_csv("Body_weight_data.csv") %>%
  rename(Reproductive_status = Reproductive_female_in_pair) %>%
  filter(Sex == "female") %>%
  group_by(Group) %>%
  summarise(
    n_females = n_distinct(Individual_No),
    n_repro   = n_distinct(Individual_No[Reproductive_status == "yes"]),
    .groups   = "drop"
  ) %>%
  mutate(
    Group      = factor(Group, levels = c("A","B","C","D")),
    n_nonrepro = n_females - n_repro
  )

cat("Reproductive female summary:\n")
print(repro_data)

model_repro_prop <- glm(
  cbind(n_repro, n_nonrepro) ~ Group,
  data   = repro_data,
  family = binomial(link = "logit")
)

model_null_repro <- glm(
  cbind(n_repro, n_nonrepro) ~ 1,
  data   = repro_data,
  family = binomial(link = "logit")
)

cat("\nGLM Summary:\n")
print(summary(model_repro_prop))

lr_repro <- anova(model_null_repro, model_repro_prop, test = "Chisq")
cat("\nLikelihood Ratio Test (Group effect):\n")
print(lr_repro)

cat("\nNo post-hoc comparisons performed (overall group effect not significant).\n")

cat("\nAnalysis complete.\n")
