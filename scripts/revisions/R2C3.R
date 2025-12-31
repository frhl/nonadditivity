library(ggplot2)
library(dplyr)
library(latex2exp)
library(bravastring)

set.seed(32)

# --- 1. Encodings ---
get_general_orthogonal_encoding <- function(g) {
  n <- length(g)
  r <- sum(g == 0) / n; h <- sum(g == 1) / n; a <- sum(g == 2) / n
  if(h == 0 | a == 0 | r == 0) return(rep(NA, n))
  
  denom_sq <- h * a * r * (a + r - (a - r)^2)
  if(denom_sq <= 0) return(rep(NA, n))
  
  scaling <- 1 / sqrt(denom_sq)
  enc <- numeric(n)
  enc[g == 0] <- (-h * a) * scaling
  enc[g == 1] <- (2 * a * r) * scaling
  enc[g == 2] <- (-h * r) * scaling
  return(enc)
}

get_hwe_encoding <- function(g) {
  p <- mean(g) / 2
  if(p <= 0 | p >= 1) return(rep(NA, length(g)))
  enc <- numeric(length(g))
  enc[g == 0] <- -p / (1 - p)
  enc[g == 1] <- 1
  enc[g == 2] <- -(1 - p) / p
  return(enc)
}

# --- 2. Simulation Loop ---

# We will test 8 distinct levels of inbreeding (0 to 0.7)
f_levels <- seq(0, 0.35, by=0.05)
n_variants <- 200 # Variants per facet
n_samples <- 10000

all_results <- data.frame()

print("Simulating facets...")

for (f in f_levels) {
  
  for(i in 1:n_variants) {
    # 1. Genotypes with specific F
    maf <- runif(1, 0.2, 0.5)
    p <- maf; q <- 1 - p
    
    prob_aa <- q^2 + p*q*f
    prob_Aa <- 2*p*q*(1 - f)
    prob_AA <- p^2 + p*q*f
    
    g <- sample(c(0, 1, 2), n_samples, replace = TRUE, 
                prob = c(prob_aa, prob_Aa, prob_AA))
    if(var(g) == 0) next
    
    # 2. Phenotype: Purely Additive (Truth = Null Dominance)
    y <- 0.3 * g + rnorm(n_samples)
    
    # 3. Encodings
    enc_orth <- get_general_orthogonal_encoding(g)
    enc_hwe  <- get_hwe_encoding(g)
    if(any(is.na(enc_orth))) next
    
    # 4. Univariate Tests (Dominance term only)
    p_hwe <- summary(lm(y ~ enc_hwe))$coefficients["enc_hwe", "Pr(>|t|)"]
    p_orth <- summary(lm(y ~ enc_orth))$coefficients["enc_orth", "Pr(>|t|)"]
    
    # Store
    all_results <- rbind(all_results, data.frame(
      F_Val = f, # Store numeric F for ordering
      F_Level = paste0("F = ", f),
      P_Standard = p_hwe,
      P_Orthogonal = p_orth
    ))
  }
}

# --- 3. Data Preparation for Plotting ---

# Convert to -log10
all_results$logP_Standard <- -log10(all_results$P_Standard)
all_results$logP_Orthogonal <- -log10(all_results$P_Orthogonal)

# Ensure facets are ordered numerically by F, not alphabetically
all_results$F_Level <- factor(all_results$F_Level, 
                              levels = paste0("F = ", sort(unique(all_results$F_Val))))

# Calculate Correlation per Facet for the Label
cor_labels <- all_results %>%
  group_by(F_Level) %>%
  summarize(
    # Calculate Pearson correlation of the log values
    cor_val = cor(logP_Standard, logP_Orthogonal, use = "complete.obs"),
    label = paste0("r = ", sprintf("%.2f", cor_val))
  )

# --- Calculate Lambda GC (Genomic Inflation Factor) ---
lambda_gc_results <- all_results %>%
  group_by(F_Val, F_Level) %>%
  summarize(
    lambda_standard = calc_inflation(P_Standard),
    lambda_orthogonal = calc_inflation(P_Orthogonal),
    n_variants = n(),
    .groups = "drop"
  ) %>%
  arrange(F_Val)

cat("\n=== Genomic Inflation Factors (Lambda GC) ===\n")
print(lambda_gc_results)

cat("\n=== Summary Statistics ===\n")
cat(sprintf("Maximum lambda_GC for Standard (HWE-dependent) encoding: %.3f (at F = %.2f)\n",
            max(lambda_gc_results$lambda_standard),
            lambda_gc_results$F_Val[which.max(lambda_gc_results$lambda_standard)]))
cat(sprintf("Maximum lambda_GC for Orthogonal encoding: %.3f (at F = %.2f)\n",
            max(lambda_gc_results$lambda_orthogonal),
            lambda_gc_results$F_Val[which.max(lambda_gc_results$lambda_orthogonal)]))
cat(sprintf("Median lambda_GC across F levels for Standard encoding: %.3f\n",
            median(lambda_gc_results$lambda_standard)))
cat(sprintf("Median lambda_GC across F levels for Orthogonal encoding: %.3f\n",
            median(lambda_gc_results$lambda_orthogonal)))

# Overall lambda across ALL p-values (pooled across all F levels)
overall_lambda_standard <- calc_inflation(all_results$P_Standard)
overall_lambda_orthogonal <- calc_inflation(all_results$P_Orthogonal)
cat(sprintf("\nOverall lambda_GC (all %d variants pooled) for Standard encoding: %.3f\n",
            nrow(all_results), overall_lambda_standard))
cat(sprintf("Overall lambda_GC (all %d variants pooled) for Orthogonal encoding: %.3f\n",
            nrow(all_results), overall_lambda_orthogonal))

cat(sprintf("\nMinimum lambda_GC for Standard encoding: %.3f (at F = %.2f)\n",
            min(lambda_gc_results$lambda_standard),
            lambda_gc_results$F_Val[which.min(lambda_gc_results$lambda_standard)]))
cat(sprintf("Minimum lambda_GC for Orthogonal encoding: %.3f (at F = %.2f)\n\n",
            min(lambda_gc_results$lambda_orthogonal),
            lambda_gc_results$F_Val[which.min(lambda_gc_results$lambda_orthogonal)]))

# --- 4. Plotting ---

p <- ggplot(all_results, aes(x = logP_Standard, y = logP_Orthogonal)) +
  geom_point(color = "darkblue", size = 1.5, alpha = 0.9) +
  
  # Diagonal reference line (Identity)
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  
  # Faceting
  facet_wrap(~ F_Level, nrow = 2) +
  
  # Add Correlation Labels (Top Right Corner)
  # Inf places text at the edge of the plot area
  geom_text(data = cor_labels, aes(label = label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic", inherit.aes = FALSE) +
  
  # Labels with TeX formatting
  labs(
    title = "Comparison of Nonadditive encodings as F (Inbreeding) increases", 
    subtitle = "200 variants and 5000 samples",
    y = TeX(r"($-\log_{10}(P_{Orthogonal\ Lassen\ 2025})$)"),
    x = TeX(r"($-\log_{10}(P_{Orthogonal\ Others})$)")
  ) +
  
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 12)
  )

print(p)
ggsave("~/Desktop/f_coeffecient_faceted_p_value_comparison_tex.pdf", p, width = 8, height = 4.5)


# single plot (important too!)
single_dt <- all_results[all_results$F_Val == 0,]
ggplot(single_dt, aes(x = logP_Standard, y = logP_Orthogonal)) +
  geom_point(color = "darkblue", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(
    title = "Comparison of Nonadditive encodings under HWE", 
    subtitle = "200 variants and 5000 samples",
    y = TeX(r"($-\log_{10}(P_{Orthogonal\ Lassen\ 2025})$)"),
    x = TeX(r"($-\log_{10}(P_{Orthogonal\ Others})$)")
  ) +
  
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 12)
  )
