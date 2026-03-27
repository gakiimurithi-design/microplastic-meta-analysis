# Load required packages
library(meta)
library(metafor)
library(dmetar)
library(ggplot2)
library(dplyr)

# Set working directory
setwd("C:/MetaAnalysis")

# Read data
data <- read.csv("microplastic_meta_data.csv", stringsAsFactors = FALSE)

# Calculate effect sizes (correlation between duration and abundance)
# For studies with multiple time points, we calculate r using the full dataset
# Here we provide code for the meta-analysis of the 12 studies with reported correlations

# Create data frame for meta-analysis
meta_data <- data.frame(
  study = c("Huang et al. 2020", "Li et al. 2020", "Li et al. 2022", 
            "Jia et al. 2025", "Li et al. 2025b", "Liu et al. 2025",
            "Wang et al. 2022", "Zhang et al. 2023", "Feng et al. 2021",
            "Liu et al. 2022", "Rezaei et al. 2022", "Abbasi et al. 2021"),
  n = c(45, 20, 54, 28, 42, 66, 30, 15, 100, 25, 10, 30),
  r = c(0.86, 0.89, 0.92, 0.94, 0.98, 0.88, 0.82, 0.85, 0.79, 0.83, 0.75, 0.78),
  stringsAsFactors = FALSE
)

# Fisher's z transformation
meta_data$z <- 0.5 * log((1 + meta_data$r) / (1 - meta_data$r))
meta_data$se_z <- 1 / sqrt(meta_data$n - 3)
meta_data$var_z <- meta_data$se_z^2

# Random-effects meta-analysis
meta_model <- rma(yi = z, vi = var_z, data = meta_data, method = "REML")
summary(meta_model)

# Back-transform to correlation coefficient
pooled_z <- coef(meta_model)[1]
pooled_r <- (exp(2 * pooled_z) - 1) / (exp(2 * pooled_z) + 1)
ci_lower_z <- pooled_z - 1.96 * sqrt(meta_model$se^2)
ci_upper_z <- pooled_z + 1.96 * sqrt(meta_model$se^2)
ci_lower_r <- (exp(2 * ci_lower_z) - 1) / (exp(2 * ci_lower_z) + 1)
ci_upper_r <- (exp(2 * ci_upper_z) - 1) / (exp(2 * ci_upper_z) + 1)

cat("Pooled correlation (r):", round(pooled_r, 3), "\n")
cat("95% CI:", round(ci_lower_r, 3), "-", round(ci_upper_r, 3), "\n")

# Heterogeneity
cat("I-squared:", round(meta_model$I2, 1), "%\n")
cat("Q-statistic p-value:", format(meta_model$QEp, digits = 4), "\n")

# Forest plot
forest(meta_model, 
       slab = meta_data$study,
       xlab = "Fisher's z",
       main = "Forest Plot: Correlation between Mulching Duration and MP Abundance")

# Meta-regression for accumulation rates (simulated data for illustration)
# Actual data from individual time points would be used here
# This is a simplified version using the quadratic model from the main text

# Create data frame for time series
time_series <- data.frame(
  duration = rep(seq(1, 40, by = 1), each = 3),
  mp_abundance = 10^(3.12 + 0.042 * seq(1, 40, by = 1) + 0.00089 * (seq(1, 40, by = 1))^2) + 
                 rnorm(120, 0, 0.2 * 10^(3.12 + 0.042 * seq(1, 40, by = 1) + 0.00089 * (seq(1, 40, by = 1))^2))
)

# Quadratic meta-regression
quad_model <- lm(log10(mp_abundance) ~ duration + I(duration^2), data = time_series)
summary(quad_model)

# Subgroup analysis (simulated with categorical moderator)
# For actual analysis, we would use the full dataset with moderators

# Publication bias: funnel plot
funnel(meta_model, main = "Funnel Plot for Publication Bias")

# Egger's regression test
egger_test <- regtest(meta_model, model = "lm")
print(egger_test)

# Trim-and-fill analysis
tf_model <- trimfill(meta_model)
summary(tf_model)