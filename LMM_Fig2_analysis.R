# ============================================================================
# LMM Analysis: Body Weight Trajectories of Reproductive vs Non-reproductive Females
# (Pooled across Groups A-D; used for Fig. 2)
#
# Author  : Manabu Kishi
# Journal : Aquatic Insects (submitted)
# Date    : 2026-04-29
#
# Description:
#   This script compares the trajectory of fresh body weight change after adult
#   emergence between reproductive and non-reproductive females, using data
#   pooled across all four experimental groups (A-D). A linear mixed model
#   (LMM) with AR(1) covariance structure is fitted to account for repeated
#   measurements within individuals. Post-hoc pairwise contrasts with Bonferroni
#   correction are used to identify time points at which the two groups differ
#   significantly (asterisks in Fig. 2).
#
# Input file (must be in the working directory):
#   - Body_weight_data.csv
#
# Software:
#   R version 4.5.1 (R Core Team, 2025)
#   Packages: tidyverse 2.0.0, nlme 3.1-169, emmeans 2.0.1
# ============================================================================


# ============================================================================
# 1. Package Loading
# ============================================================================

packages <- c("tidyverse", "nlme", "emmeans")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


# ============================================================================
# 2. Data Loading and Preparation
# ============================================================================

cat("Loading data...\n\n")

weight_data <- read_csv("Body_weight_data.csv") %>%
  mutate(Group = factor(Group, levels = c("A", "B", "C", "D")))

# Subset to females only, all groups pooled
# Create unique individual ID (Individual_No is unique only within each Group)
females_pooled <- weight_data %>%
  filter(Sex == "female") %>%
  drop_na(Body_weight_mg) %>%
  mutate(
    Reproductive_status = factor(
      Reproductive_female_in_pair,
      levels = c("yes", "no"),
      labels = c("Reproductive", "Non-reproductive")
    ),
    Days          = Days_after_adult_emergence,
    Individual_ID = factor(paste(Group, Individual_No, sep = "_"))
  )

cat("Sample sizes (females, pooled across groups A-D):\n")
females_pooled %>%
  group_by(Reproductive_status, Days) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()


# ============================================================================
# 3. Descriptive Statistics
# ============================================================================

cat("\nDescriptive Statistics (mean and SE by reproductive status and time point):\n")
females_pooled %>%
  group_by(Reproductive_status, Days) %>%
  summarise(
    n    = n(),
    mean = mean(Body_weight_mg),
    SE   = sd(Body_weight_mg) / sqrt(n()),
    .groups = "drop"
  ) %>%
  print()


# ============================================================================
# 4. LMM: Reproductive status × Time interaction
# ============================================================================

cat("\n========================================\n")
cat("LMM: Reproductive Status × Time (AR(1) covariance)\n")
cat("========================================\n\n")

# Full model (with interaction)
model_full <- lme(
  Body_weight_mg ~ Reproductive_status * Days,
  random      = ~ 1 | Individual_ID,
  correlation = corAR1(form = ~ Days | Individual_ID),
  data        = females_pooled,
  method      = "ML"
)

# Reduced model (without interaction)
model_reduced <- lme(
  Body_weight_mg ~ Reproductive_status + Days,
  random      = ~ 1 | Individual_ID,
  correlation = corAR1(form = ~ Days | Individual_ID),
  data        = females_pooled,
  method      = "ML"
)

cat("Full model summary:\n")
print(summary(model_full))

cat("\nLikelihood Ratio Test (Reproductive_status × Time interaction):\n")
lr_interaction <- anova(model_reduced, model_full)
print(lr_interaction)


# ============================================================================
# 5. Post-hoc Pairwise Contrasts at Each Time Point (Bonferroni correction)
# ============================================================================

cat("\n========================================\n")
cat("Post-hoc Contrasts: Reproductive vs Non-reproductive at Each Time Point\n")
cat("(Bonferroni correction; 7 time points tested)\n")
cat("========================================\n\n")

emm_repro <- emmeans(
  model_full,
  ~ Reproductive_status | Days,
  at = list(Days = c(0, 5, 10, 15, 20, 25, 30))
)

contrasts_repro <- contrast(
  emm_repro,
  method = "pairwise",
  adjust = "bonferroni"
)

cat("Estimated Marginal Means at each time point:\n")
print(emm_repro)

cat("\nPairwise contrasts (Reproductive - Non-reproductive) at each time point:\n")
print(contrasts_repro)

cat("\nSignificance summary (Bonferroni-corrected, alpha = 0.05/7 = 0.0071):\n")
contrasts_df <- as.data.frame(contrasts_repro)
contrasts_df$Significant <- ifelse(contrasts_df$p.value < 0.0071, "*", "ns")
print(contrasts_df[, c("Days", "estimate", "SE", "p.value", "Significant")],
      row.names = FALSE)


# ============================================================================
# 6. Summary for Fig. 2 and Methods
# ============================================================================

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("Interaction (Reproductive_status × Days):\n")
cat(sprintf("  chi2 = %.3f, df = %d, P = %.4f\n\n",
    lr_interaction$L.Ratio[2],
    lr_interaction$df[2] - lr_interaction$df[1],
    lr_interaction$`p-value`[2]))

cat("Significant time points (asterisks in Fig. 2):\n")
sig_days <- contrasts_df$Days[contrasts_df$Significant == "*"]
if (length(sig_days) > 0) {
  cat(sprintf("  Days: %s\n", paste(sig_days, collapse = ", ")))
} else {
  cat("  None\n")
}

cat("\nAnalysis complete.\n")
