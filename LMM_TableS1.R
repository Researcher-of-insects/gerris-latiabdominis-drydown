# ============================================================================
# LMM Analysis for Table S1:
# Body weight trajectories in four categories of G. latiabdominis
# (Group effect, Days effect, Group × Days interaction)
#
# Author  : Manabu Kishi
# Journal : Aquatic Insects (submitted)
# Date    : 2026-05-01
#
# Description:
#   This script fits linear mixed models (LMM) with AR(1) covariance structure
#   separately for four categories (reproductive females, non-reproductive
#   females, males with reproductive females, males with non-reproductive
#   females). Likelihood ratio tests assess the significance of (1) the group
#   main effect, (2) the days main effect, and (3) the group × days interaction.
#   Results are reported in Table S1 of the manuscript.
#
# Input file (must be in the working directory):
#   - Body_weight_data.csv
#
# Software:
#   R version 4.5.1 (R Core Team, 2025)
#   Packages: tidyverse 2.0.0, nlme 3.1-169
# ============================================================================


# ============================================================================
# 1. Package Loading
# ============================================================================

if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse"); library(tidyverse)
}
if (!require("nlme", character.only = TRUE)) {
  install.packages("nlme"); library(nlme)
}


# ============================================================================
# 2. Data Loading and Preparation
# ============================================================================

cat("Loading data...\n\n")

weight_data <- read_csv("Body_weight_data.csv") %>%
  rename(
    Reproductive_status = Reproductive_female_in_pair,
    Days                = Days_after_adult_emergence
  ) %>%
  drop_na(Body_weight_mg) %>%
  mutate(
    Group         = factor(Group, levels = c("A", "B", "C", "D")),
    Individual_ID = factor(paste(Group, Individual_No, sep = "_"))
  )

# Define four categories
weight_data <- weight_data %>%
  mutate(
    Category = case_when(
      Sex == "female" & Reproductive_status == "yes" ~ "Reproductive females",
      Sex == "female" & Reproductive_status == "no"  ~ "Non-reproductive females",
      Sex == "male"   & Reproductive_status == "yes" ~ "Males with reproductive females",
      Sex == "male"   & Reproductive_status == "no"  ~ "Males with non-reproductive females"
    ),
    Category = factor(Category, levels = c(
      "Reproductive females",
      "Non-reproductive females",
      "Males with reproductive females",
      "Males with non-reproductive females"
    ))
  )

cat("Sample sizes per category:\n")
weight_data %>%
  group_by(Category) %>%
  summarise(
    n_obs = n(),
    n_ind = n_distinct(Individual_ID),
    .groups = "drop"
  ) %>%
  print()
cat("\n")


# ============================================================================
# 3. LMM Function: likelihood ratio tests for Group, Days, Group×Days
# ============================================================================

run_lmm <- function(data, category_name) {

  cat(strrep("=", 70), "\n")
  cat("Category:", category_name, "\n")
  cat(strrep("=", 70), "\n\n")

  sub <- data %>%
    filter(Category == category_name) %>%
    droplevels()

  cat(sprintf("  n_obs = %d,  n_ind = %d\n\n",
              nrow(sub), n_distinct(sub$Individual_ID)))

  # ---- Model 1: Full model (Group + Days + Group×Days) ----
  m_full <- lme(
    Body_weight_mg ~ Group * Days,
    random      = ~ 1 | Individual_ID,
    correlation = corAR1(form = ~ Days | Individual_ID),
    data        = sub,
    method      = "ML"
  )

  # ---- Model 2: No interaction (Group + Days) ----
  m_nointer <- lme(
    Body_weight_mg ~ Group + Days,
    random      = ~ 1 | Individual_ID,
    correlation = corAR1(form = ~ Days | Individual_ID),
    data        = sub,
    method      = "ML"
  )

  # ---- Model 3: Days only (no group effect) ----
  m_nogroup <- lme(
    Body_weight_mg ~ Days,
    random      = ~ 1 | Individual_ID,
    correlation = corAR1(form = ~ Days | Individual_ID),
    data        = sub,
    method      = "ML"
  )

  # ---- Model 4: Null model (intercept only) ----
  m_null <- lme(
    Body_weight_mg ~ 1,
    random      = ~ 1 | Individual_ID,
    correlation = corAR1(form = ~ Days | Individual_ID),
    data        = sub,
    method      = "ML"
  )

  # ---- Likelihood Ratio Tests ----

  # Group × Days interaction: full vs no-interaction
  lr_inter <- anova(m_nointer, m_full)
  cat("LRT: Group × Days interaction (full vs. no-interaction model):\n")
  print(lr_inter)

  # Group main effect: no-interaction vs days-only
  lr_group <- anova(m_nogroup, m_nointer)
  cat("\nLRT: Group main effect (no-interaction vs. days-only model):\n")
  print(lr_group)

  # Days main effect: days-only vs null
  lr_days <- anova(m_null, m_nogroup)
  cat("\nLRT: Days main effect (days-only vs. null model):\n")
  print(lr_days)

  # ---- Summary ----
  chi2_inter <- lr_inter$L.Ratio[2]
  df_inter   <- lr_inter$df[2] - lr_inter$df[1]
  p_inter    <- lr_inter$`p-value`[2]

  chi2_group <- lr_group$L.Ratio[2]
  df_group   <- lr_group$df[2] - lr_group$df[1]
  p_group    <- lr_group$`p-value`[2]

  chi2_days  <- lr_days$L.Ratio[2]
  df_days    <- lr_days$df[2] - lr_days$df[1]
  p_days     <- lr_days$`p-value`[2]

  cat("\n--- Summary for Table S1 ---\n")
  cat(sprintf("  Group effect:       χ² = %.3f, df = %d, P = %.4f  %s\n",
              chi2_group, df_group, p_group,
              ifelse(p_group < 0.05, "***" , "ns")))
  cat(sprintf("  Days effect:        χ² = %.3f, df = %d, P = %.6f  %s\n",
              chi2_days,  df_days,  p_days,
              ifelse(p_days  < 0.05, "***" , "ns")))
  cat(sprintf("  Group × Days inter: χ² = %.3f, df = %d, P = %.4f  %s\n",
              chi2_inter, df_inter, p_inter,
              ifelse(p_inter < 0.05, "***" , "ns")))
  cat("\n")

  # Return values for final summary table
  list(
    category   = category_name,
    n_obs      = nrow(sub),
    n_ind      = n_distinct(sub$Individual_ID),
    chi2_group = chi2_group, df_group = df_group, p_group = p_group,
    chi2_days  = chi2_days,  df_days  = df_days,  p_days  = p_days,
    chi2_inter = chi2_inter, df_inter = df_inter, p_inter = p_inter
  )
}


# ============================================================================
# 4. Run LMM for All Four Categories
# ============================================================================

categories <- levels(weight_data$Category)

results_list <- lapply(categories, function(cat_name) {
  run_lmm(weight_data, cat_name)
})


# ============================================================================
# 5. Final Summary Table (Table S1 format)
# ============================================================================

cat(strrep("=", 70), "\n")
cat("FINAL SUMMARY TABLE S1\n")
cat(strrep("=", 70), "\n\n")

summary_df <- do.call(rbind, lapply(results_list, function(r) {
  data.frame(
    Category   = r$category,
    n_obs      = r$n_obs,
    n_ind      = r$n_ind,
    Effect     = c("Group", "Days", "Group x Days"),
    chi2       = round(c(r$chi2_group, r$chi2_days,  r$chi2_inter), 3),
    df         = c(r$df_group,   r$df_days,   r$df_inter),
    P_value    = formatC(c(r$p_group, r$p_days, r$p_inter),
                         format = "f", digits = 4),
    Sig        = ifelse(c(r$p_group, r$p_days, r$p_inter) < 0.05, "*", "ns"),
    stringsAsFactors = FALSE
  )
}))

print(summary_df, row.names = FALSE)

cat("\n* P < 0.05; ns = not significant\n")
cat("\nAnalysis complete.\n")
