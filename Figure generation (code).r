# Load required packages
library(meta)
library(metafor)
library(ggplot2)

# ----------------------------------------------------------------------
# 1. Prepare data (replace with your actual data)
# ----------------------------------------------------------------------
# Example data (replace with your extracted data)
meta_data <- data.frame(
  study = c("Huang et al. 2020", "Li et al. 2020", "Li et al. 2022", 
            "Jia et al. 2025", "Li et al. 2025b", "Liu et al. 2025",
            "Wang et al. 2022", "Zhang et al. 2023", "Feng et al. 2021",
            "Liu et al. 2022", "Rezaei et al. 2022", "Abbasi et al. 2021"),
  year = c(2020, 2020, 2022, 2025, 2025, 2025,
           2022, 2023, 2021, 2022, 2022, 2021),
  r = c(0.86, 0.89, 0.92, 0.94, 0.98, 0.88,
        0.82, 0.85, 0.79, 0.83, 0.75, 0.78),
  n = c(45, 20, 54, 28, 42, 66,
        30, 15, 100, 25, 10, 30),
  polymer = factor(c("PE", "PE", "PE", "PE", "PE", "PE",
                     "PE", "PE", "PE", "PE", "Other", "Other")),
  climate = factor(c("arid", "arid", "humid", "arid", "arid", "arid",
                     "humid", "arid", "humid", "humid", "arid", "arid")),
  soil = factor(c("sandy", "sandy", "loam", "sandy", "sandy", "sandy",
                  "loam", "sandy", "loam", "loam", "clay", "sandy"))
)

# ----------------------------------------------------------------------
# 2. Helper function to run meta-analysis and return funnel data
# ----------------------------------------------------------------------
run_meta <- function(data) {
  # Fisher's z transformation
  data$z <- 0.5 * log((1 + data$r) / (1 - data$r))
  data$var_z <- 1 / (data$n - 3)
  # Random-effects model
  m <- rma(yi = z, vi = var_z, data = data, method = "REML")
  return(m)
}

# ----------------------------------------------------------------------
# 3. S8.1 Funnel plot for overall analysis
# ----------------------------------------------------------------------
overall_meta <- run_meta(meta_data)

# Save as PDF (or png)
pdf("S8.1_funnel_overall.pdf", width = 6, height = 6)
funnel(overall_meta, 
       main = "Funnel Plot for Overall Analysis",
       xlab = "Fisher's z", ylab = "Standard Error",
       back = "white")
# Add pseudo-95% confidence limits (automatically included by funnel())
dev.off()

# ----------------------------------------------------------------------
# 4. S8.2 Subgroup funnel plots
# ----------------------------------------------------------------------
# Define subgroups
subgroups <- list(
  PE = list(data = subset(meta_data, polymer == "PE"),
            title = "Polyethylene Subgroup"),
  Other = list(data = subset(meta_data, polymer == "Other"),
               title = "Biodegradable/Other Subgroup"),
  Arid = list(data = subset(meta_data, climate == "arid"),
              title = "Arid Climate Subgroup"),
  Humid = list(data = subset(meta_data, climate == "humid"),
               title = "Humid Climate Subgroup"),
  Sandy = list(data = subset(meta_data, soil == "sandy"),
               title = "Sandy Soil Subgroup"),
  Loam = list(data = subset(meta_data, soil == "loam"),
              title = "Loam Soil Subgroup")
)

# Create a PDF with multiple pages
pdf("S8.2_subgroup_funnels.pdf", width = 6, height = 6, onefile = TRUE)
for (sg in names(subgroups)) {
  sub_data <- subgroups[[sg]]$data
  if (nrow(sub_data) > 2) {  # need at least 3 studies for a funnel plot
    m_sub <- run_meta(sub_data)
    funnel(m_sub, 
           main = subgroups[[sg]]$title,
           xlab = "Fisher's z", ylab = "Standard Error",
           back = "white")
  } else {
    plot.new()
    text(0.5, 0.5, paste("Insufficient studies (n =", nrow(sub_data), 
                         ") for funnel plot"), cex = 1.2)
  }
}
dev.off()

# ----------------------------------------------------------------------
# 5. S8.3 Cumulative meta-analysis by publication year
# ----------------------------------------------------------------------
# Order by publication year (ascending)
meta_data_ordered <- meta_data[order(meta_data$year), ]

# Prepare data for cumulative analysis
cum_years <- sort(unique(meta_data_ordered$year))
cum_estimates <- data.frame(year = numeric(), r = numeric(), ci_lower = numeric(), ci_upper = numeric())

for (i in seq_along(cum_years)) {
  # Subset studies up to current year
  subset_data <- meta_data_ordered[meta_data_ordered$year <= cum_years[i], ]
  if (nrow(subset_data) >= 2) {
    m_cum <- run_meta(subset_data)
    # Back-transform to correlation
    cum_r <- (exp(2 * m_cum$b) - 1) / (exp(2 * m_cum$b) + 1)
    ci_lower <- (exp(2 * (m_cum$b - 1.96 * sqrt(m_cum$se))) - 1) /
                (exp(2 * (m_cum$b - 1.96 * sqrt(m_cum$se))) + 1)
    ci_upper <- (exp(2 * (m_cum$b + 1.96 * sqrt(m_cum$se))) - 1) /
                (exp(2 * (m_cum$b + 1.96 * sqrt(m_cum$se))) + 1)
    cum_estimates <- rbind(cum_estimates, 
                           data.frame(year = cum_years[i], r = cum_r,
                                      ci_lower = ci_lower, ci_upper = ci_upper))
  }
}

# Plot cumulative meta-analysis
p <- ggplot(cum_estimates, aes(x = year, y = r)) +
  geom_point(size = 3) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_hline(yintercept = 0.82, linetype = "dashed", color = "red") +
  labs(title = "Cumulative Meta-Analysis by Publication Year",
       x = "Publication Year", y = "Pooled Correlation (r)") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.5, 1.0))

ggsave("S8.3_cumulative_meta.pdf", p, width = 7, height = 5)

# ----------------------------------------------------------------------
# 6. Additional: Egger's test and trim-and-fill for overall (if needed)
# ----------------------------------------------------------------------
# Egger's regression test
egger <- regtest(overall_meta, model = "lm")
cat("Egger's test: t =", egger$zval, ", p =", egger$pval, "\n")

# Trim-and-fill analysis
tf <- trimfill(overall_meta)
cat("Trim-and-fill imputed studies:", tf$k0, "\n")
cat("Adjusted pooled r:", (exp(2 * tf$b) - 1) / (exp(2 * tf$b) + 1), "\n")