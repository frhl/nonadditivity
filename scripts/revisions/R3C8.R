library(ggplot2)
library(dplyr)
library(tidyr)
library(Hmisc)

k <- 6   # Successes
n <- 21  # Trials

# 1. Hmisc Implementation
binconf(k, n)

# 2. compare with manual Wilson Score Verification 
p <- k / n
z <- qnorm(0.975) # 95% CI

# Calculate components
denom  <- 1 + z^2/n
center <- (p + z^2/(2*n)) / denom
margin <- z * sqrt((p*(1-p)/n) + (z^2/(4*n^2))) / denom

# Print comparison
cat(sprintf("\nObserved: %.1f%%\nManual CI:  %.1f%% - %.1f%%\n", 
            p * 100, 
            (center - margin) * 100, 
            (center + margin) * 100))



#### analysis of rare variants.. How many Ns?

# Function: Calculate N for 80% power
get_n <- function(maf, beta, model = "Additive", alpha = 5e-8, power = 0.8) {
  
  z_score <- qnorm(1 - alpha/2) + qnorm(power)
  
  # Variance of G depends on model
  if (model == "Additive") {
    var_g <- 2 * maf * (1 - maf)
  } else {
    # Recessive bottleneck: homozygotes only
    p_hom <- maf^2
    var_g <- p_hom * (1 - p_hom)
  }
  
  h2 <- var_g * beta^2
  return((z_score^2) / h2)
}

# Define scenarios and calculate N row-wise
sim_data <- expand.grid(
  MAF = c(0.05, 0.01, 0.001),
  Model = c("Additive", "Recessive"),
  Beta = 2
) %>%
  rowwise() %>%
  mutate(
    N = get_n(MAF, Beta, Model),
    Label = sprintf("MAF: %s", MAF) # Cleaner labelling
  ) %>%
  ungroup()

# Summary Table: Ratio of Recessive to Additive cost
sim_data %>%
  pivot_wider(names_from = Model, values_from = N) %>%
  mutate(
    Ratio = round(Recessive / Additive, 1),
    Additive = scales::comma(Additive),
    Recessive = scales::comma(Recessive)
  ) %>%
  select(MAF, Beta, Additive, Recessive, Ratio) %>%
  print()

# Plotting
ggplot(sim_data, aes(x = as.factor(MAF), y = N, fill = Model)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = scales::comma(round(N))), 
            position = position_dodge(0.9), vjust = -0.5, size = 3.5) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "The 'Homozygote Bottleneck' for Rare Variants",
    subtitle = "Sample size required for 80% Power (Alpha = 5e-8)",
    y = "Sample Size (Log Scale)",
    x = "Minor Allele Frequency (MAF)"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
