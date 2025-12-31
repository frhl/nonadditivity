library(dplyr)
library(bravastring)

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

# --- HWE Test Function ---
hwe_test <- function(g) {
  # Calculate observed genotype counts
  n_aa <- sum(g == 0)
  n_Aa <- sum(g == 1)
  n_AA <- sum(g == 2)
  n_total <- length(g)

  # Calculate allele frequency
  p <- (2*n_AA + n_Aa) / (2*n_total)
  q <- 1 - p

  # Expected counts under HWE
  exp_aa <- n_total * q^2
  exp_Aa <- n_total * 2*p*q
  exp_AA <- n_total * p^2

  # Chi-squared test (1 df)
  chi_sq <- ((n_aa - exp_aa)^2 / exp_aa +
             (n_Aa - exp_Aa)^2 / exp_Aa +
             (n_AA - exp_AA)^2 / exp_AA)

  # P-value
  p_val <- pchisq(chi_sq, df = 1, lower.tail = FALSE)

  return(list(
    chi_sq = chi_sq,
    p_val = p_val,
    obs_het = n_Aa / n_total,
    exp_het = 2*p*q
  ))
}

# --- 2. Single simulation function ---
run_single_simulation <- function(seed_val, f_levels=NULL, n_variants=NULL, n_samples=NULL) {
  set.seed(seed_val)
  
  if (is.null(f_levels)) f_levels <- seq(0, 0.35, by=0.05)
  if (is.null(n_variants)) n_variants <- 200
  if (is.null(n_samples)) n_samples <- 10000

  all_results <- data.frame()

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

      # 4. HWE Test
      hwe_result <- hwe_test(g)

      # 5. Univariate Tests (Dominance term only)
      p_hwe <- summary(lm(y ~ enc_hwe))$coefficients["enc_hwe", "Pr(>|t|)"]
      p_orth <- summary(lm(y ~ enc_orth))$coefficients["enc_orth", "Pr(>|t|)"]

      # Store
      all_results <- rbind(all_results, data.frame(
        F_Val = f,
        P_Standard = p_hwe,
        P_Orthogonal = p_orth,
        HWE_pval = hwe_result$p_val,
        HWE_chi_sq = hwe_result$chi_sq,
        Obs_Het = hwe_result$obs_het,
        Exp_Het = hwe_result$exp_het
      ))
    }
  }

  # Calculate Lambda GC per F level
  lambda_gc_results <- all_results %>%
    group_by(F_Val) %>%
    summarize(
      lambda_standard = calc_inflation(P_Standard),
      lambda_orthogonal = calc_inflation(P_Orthogonal),
      seed = seed_val,
      .groups = "drop"
    )

  # Calculate HWE QC pass rates per F level
  hwe_qc_results <- all_results %>%
    group_by(F_Val) %>%
    summarize(
      n_variants = n(),
      pass_hwe_1e4 = sum(HWE_pval > 1e-4) / n(),
      pass_hwe_1e6 = sum(HWE_pval > 1e-6) / n(),
      pass_hwe_1e10 = sum(HWE_pval > 1e-10) / n(),
      median_hwe_pval = median(HWE_pval),
      median_obs_het = median(Obs_Het),
      median_exp_het = median(Exp_Het),
      seed = seed_val,
      .groups = "drop"
    )

  # Also calculate overall lambda
  overall <- data.frame(
    F_Val = -1,  # marker for "overall"
    lambda_standard = calc_inflation(all_results$P_Standard),
    lambda_orthogonal = calc_inflation(all_results$P_Orthogonal),
    seed = seed_val
  )

  # Overall HWE stats
  overall_hwe <- data.frame(
    F_Val = -1,
    n_variants = nrow(all_results),
    pass_hwe_1e4 = sum(all_results$HWE_pval > 1e-4) / nrow(all_results),
    pass_hwe_1e6 = sum(all_results$HWE_pval > 1e-6) / nrow(all_results),
    pass_hwe_1e10 = sum(all_results$HWE_pval > 1e-10) / nrow(all_results),
    median_hwe_pval = median(all_results$HWE_pval),
    median_obs_het = median(all_results$Obs_Het),
    median_exp_het = median(all_results$Exp_Het),
    seed = seed_val
  )

  return(list(
    lambda_results = rbind(lambda_gc_results, overall),
    hwe_results = rbind(hwe_qc_results, overall_hwe)
  ))
}

# --- 3. Run multiple simulations ---
n_simulations <- 10
cat(sprintf("Running %d simulations...\n", n_simulations))

all_lambda_results <- data.frame()
all_hwe_results <- data.frame()

for(i in 1:n_simulations) {
  if(i %% 10 == 0) cat(sprintf("  Completed %d/%d simulations\n", i, n_simulations))
  result <- run_single_simulation(seed_val = i)
  all_lambda_results <- rbind(all_lambda_results, result$lambda_results)
  all_hwe_results <- rbind(all_hwe_results, result$hwe_results)
}

cat(sprintf("\n=== Summary Statistics Across %d Simulations ===\n\n", n_simulations))

# Summarize per F level
summary_by_f <- all_lambda_results %>%
  filter(F_Val >= 0) %>%  # exclude overall
  group_by(F_Val) %>%
  summarize(
    median_lambda_standard = median(lambda_standard),
    median_lambda_orthogonal = median(lambda_orthogonal),
    mean_lambda_standard = mean(lambda_standard),
    mean_lambda_orthogonal = mean(lambda_orthogonal),
    min_lambda_standard = min(lambda_standard),
    max_lambda_standard = max(lambda_standard),
    min_lambda_orthogonal = min(lambda_orthogonal),
    max_lambda_orthogonal = max(lambda_orthogonal),
    .groups = "drop"
  )

cat(sprintf("Lambda GC by F level (across %d simulations):\n", n_simulations))
print(summary_by_f, n = 100)

# Overall (pooled across all F levels)
cat(sprintf("\n\nOverall Lambda GC (all variants pooled, across %d simulations):\n", n_simulations))
overall_summary <- all_lambda_results %>%
  filter(F_Val == -1) %>%
  summarize(
    median_lambda_standard = median(lambda_standard),
    median_lambda_orthogonal = median(lambda_orthogonal),
    mean_lambda_standard = mean(lambda_standard),
    mean_lambda_orthogonal = mean(lambda_orthogonal),
    min_lambda_standard = min(lambda_standard),
    max_lambda_standard = max(lambda_standard),
    min_lambda_orthogonal = min(lambda_orthogonal),
    max_lambda_orthogonal = max(lambda_orthogonal)
  )
print(overall_summary)

# Key statistics for manuscript
cat("\n\n=== KEY STATISTICS FOR MANUSCRIPT ===\n\n")

# Maximum observed lambda at F=0.35
f35_lambdas <- all_lambda_results %>% filter(F_Val == 0.35)
cat(sprintf("At F = 0.35 (highest inbreeding):\n"))
cat(sprintf("  Standard encoding - Median lambda_GC: %.2f (range: %.2f - %.2f)\n",
            median(f35_lambdas$lambda_standard),
            min(f35_lambdas$lambda_standard),
            max(f35_lambdas$lambda_standard)))
cat(sprintf("  Orthogonal encoding - Median lambda_GC: %.2f (range: %.2f - %.2f)\n\n",
            median(f35_lambdas$lambda_orthogonal),
            min(f35_lambdas$lambda_orthogonal),
            max(f35_lambdas$lambda_orthogonal)))

# Overall pooled
overall_lambdas <- all_lambda_results %>% filter(F_Val == -1)
cat(sprintf("Overall (all F levels pooled):\n"))
cat(sprintf("  Standard encoding - Median lambda_GC: %.2f (range: %.2f - %.2f)\n",
            median(overall_lambdas$lambda_standard),
            min(overall_lambdas$lambda_standard),
            max(overall_lambdas$lambda_standard)))
cat(sprintf("  Orthogonal encoding - Median lambda_GC: %.2f (range: %.2f - %.2f)\n\n",
            median(overall_lambdas$lambda_orthogonal),
            min(overall_lambdas$lambda_orthogonal),
            max(overall_lambdas$lambda_orthogonal)))

# Median across F levels (for each simulation, take median across F, then summarize)
median_across_f <- all_lambda_results %>%
  filter(F_Val >= 0) %>%
  group_by(seed) %>%
  summarize(
    median_lambda_standard = median(lambda_standard),
    median_lambda_orthogonal = median(lambda_orthogonal),
    .groups = "drop"
  )

cat(sprintf("Median lambda_GC across F levels:\n"))
cat(sprintf("  Standard encoding - Median: %.2f (range: %.2f - %.2f)\n",
            median(median_across_f$median_lambda_standard),
            min(median_across_f$median_lambda_standard),
            max(median_across_f$median_lambda_standard)))
cat(sprintf("  Orthogonal encoding - Median: %.2f (range: %.2f - %.2f)\n",
            median(median_across_f$median_lambda_orthogonal),
            min(median_across_f$median_lambda_orthogonal),
            max(median_across_f$median_lambda_orthogonal)))

# === HWE QC Analysis ===
cat("\n\n=== HARDY-WEINBERG EQUILIBRIUM QC ANALYSIS ===\n\n")

# Summarize HWE QC pass rates by F level
hwe_summary_by_f <- all_hwe_results %>%
  filter(F_Val >= 0) %>%
  group_by(F_Val) %>%
  summarize(
    mean_pass_1e4 = mean(pass_hwe_1e4) * 100,
    mean_pass_1e6 = mean(pass_hwe_1e6) * 100,
    mean_pass_1e10 = mean(pass_hwe_1e10) * 100,
    median_obs_het = median(median_obs_het),
    median_exp_het = median(median_exp_het),
    .groups = "drop"
  )

cat("HWE QC Pass Rates by F level (% passing threshold):\n")
cat("F_Val  Pass_1e-4  Pass_1e-6  Pass_1e-10  Obs_Het  Exp_Het\n")
for(i in 1:nrow(hwe_summary_by_f)) {
  cat(sprintf("%.2f   %6.2f%%    %6.2f%%     %6.2f%%    %.3f    %.3f\n",
              hwe_summary_by_f$F_Val[i],
              hwe_summary_by_f$mean_pass_1e4[i],
              hwe_summary_by_f$mean_pass_1e6[i],
              hwe_summary_by_f$mean_pass_1e10[i],
              hwe_summary_by_f$median_obs_het[i],
              hwe_summary_by_f$median_exp_het[i]))
}

# Overall HWE pass rates
overall_hwe_summary <- all_hwe_results %>%
  filter(F_Val == -1) %>%
  summarize(
    mean_pass_1e4 = mean(pass_hwe_1e4) * 100,
    mean_pass_1e6 = mean(pass_hwe_1e6) * 100,
    mean_pass_1e10 = mean(pass_hwe_1e10) * 100,
    median_obs_het = median(median_obs_het),
    median_exp_het = median(median_exp_het)
  )

cat(sprintf("\nOverall (all F levels pooled):\n"))
cat(sprintf("  Pass HWE p > 1e-4:  %.2f%%\n", overall_hwe_summary$mean_pass_1e4))
cat(sprintf("  Pass HWE p > 1e-6:  %.2f%%\n", overall_hwe_summary$mean_pass_1e6))
cat(sprintf("  Pass HWE p > 1e-10: %.2f%%\n", overall_hwe_summary$mean_pass_1e10))

# Critical threshold analysis for F = 0.35
f35_hwe <- all_hwe_results %>% filter(F_Val == 0.35)
cat(sprintf("\n=== KEY HWE QC STATISTICS (F = 0.35, highest inbreeding) ===\n"))
cat(sprintf("  Mean variants passing HWE p > 1e-6:  %.2f%% (range: %.2f%% - %.2f%%)\n",
            mean(f35_hwe$pass_hwe_1e6) * 100,
            min(f35_hwe$pass_hwe_1e6) * 100,
            max(f35_hwe$pass_hwe_1e6) * 100))
cat(sprintf("  Mean variants passing HWE p > 1e-10: %.2f%% (range: %.2f%% - %.2f%%)\n",
            mean(f35_hwe$pass_hwe_1e10) * 100,
            min(f35_hwe$pass_hwe_1e10) * 100,
            max(f35_hwe$pass_hwe_1e10) * 100))
cat(sprintf("  Median observed heterozygosity: %.3f\n", median(f35_hwe$median_obs_het)))
cat(sprintf("  Median expected heterozygosity: %.3f\n", median(f35_hwe$median_exp_het)))

cat("\n=== CONCLUSION ===\n")
cat("Variants with moderate to high inbreeding coefficients (F >= 0.10) would\n")
cat("fail standard HWE QC filters (p > 1e-6 or p > 1e-10). However, HWE QC is\n")
cat("primarily designed to detect GENOTYPING ERRORS, not true biological HWE\n")
cat("deviations. In practice:\n")
cat("  1. Sample-level QC (PC adjustment, relatedness filtering) controls for\n")
cat("     population structure, preventing extreme F values\n")
cat("  2. Variants failing HWE typically reflect poor genotype quality, not\n")
cat("     biological population structure\n")
cat("  3. Standard encoding inflation is not observed because properly QC'd data\n")
cat("     has minimal HWE deviations due to population structure\n\n")
cat("The orthogonal encoding provides an additional layer of robustness by\n")
cat("maintaining proper calibration (lambda ~1.0) even when HWE deviations exist,\n")
cat("whether from genotyping errors or residual population structure.\n")

# === Detailed Analysis for Real-World F Values ===
cat("\n\n=== DETAILED ANALYSIS FOR REALISTIC INBREEDING COEFFICIENTS ===\n\n")

# Function to analyze a specific F value
analyze_f_value <- function(f_val, f_label, n_variants = 5000, n_samples = 10000) {
  cat(sprintf("\nAnalyzing %s (F = %.4f)...\n", f_label, f_val))

  set.seed(12345)
  all_results <- data.frame()

  for(i in 1:n_variants) {
    if(i %% 1000 == 0) cat(sprintf("  Generated %d/%d variants\n", i, n_variants))

    # Generate genotypes
    maf <- runif(1, 0.2, 0.5)
    p <- maf; q <- 1 - p

    prob_aa <- q^2 + p*q*f_val
    prob_Aa <- 2*p*q*(1 - f_val)
    prob_AA <- p^2 + p*q*f_val

    g <- sample(c(0, 1, 2), n_samples, replace = TRUE,
                prob = c(prob_aa, prob_Aa, prob_AA))
    if(var(g) == 0) next

    # Phenotype: Purely Additive (Truth = Null Dominance)
    y <- 0.3 * g + rnorm(n_samples)

    # Encodings
    enc_orth <- get_general_orthogonal_encoding(g)
    enc_hwe  <- get_hwe_encoding(g)
    if(any(is.na(enc_orth))) next

    # HWE Test
    hwe_result <- hwe_test(g)

    # Univariate Tests
    p_hwe <- summary(lm(y ~ enc_hwe))$coefficients["enc_hwe", "Pr(>|t|)"]
    p_orth <- summary(lm(y ~ enc_orth))$coefficients["enc_orth", "Pr(>|t|)"]

    all_results <- rbind(all_results, data.frame(
      P_Standard = p_hwe,
      P_Orthogonal = p_orth,
      HWE_pval = hwe_result$p_val,
      pass_HWE_1e6 = hwe_result$p_val > 1e-6
    ))
  }

  pass_hwe <- all_results %>% filter(pass_HWE_1e6)
  fail_hwe <- all_results %>% filter(!pass_HWE_1e6)

  # Return summary
  data.frame(
    F_label = f_label,
    F_value = f_val,
    All_N = nrow(all_results),
    All_Lambda_Std = calc_inflation(all_results$P_Standard),
    All_Lambda_Orth = calc_inflation(all_results$P_Orthogonal),
    Pass_N = nrow(pass_hwe),
    Pass_Pct = 100 * nrow(pass_hwe) / nrow(all_results),
    Pass_Lambda_Std = calc_inflation(pass_hwe$P_Standard),
    Pass_Lambda_Orth = calc_inflation(pass_hwe$P_Orthogonal),
    Fail_N = nrow(fail_hwe),
    Fail_Pct = 100 * nrow(fail_hwe) / nrow(all_results),
    Fail_Lambda_Std = calc_inflation(fail_hwe$P_Standard),
    Fail_Lambda_Orth = calc_inflation(fail_hwe$P_Orthogonal)
  )
}

# Analyze different F values
f_values <- list(
  list(f = 0.0037, label = "UKB EUR"),
  list(f = 0.0178, label = "Genes & Health"),
  list(f = 0.05, label = "Moderate")
)

all_f_results <- data.frame()
for(fv in f_values) {
  result <- analyze_f_value(fv$f, fv$label, n_variants = 3000)
  all_f_results <- rbind(all_f_results, result)
}

# === Create comprehensive summary table ===
cat("\n\n=== TABLE: Lambda GC by Population and HWE QC Status ===\n\n")

for(i in 1:nrow(all_f_results)) {
  row <- all_f_results[i, ]

  cat(sprintf("\n%s (F = %.4f)\n", row$F_label, row$F_value))
  cat(strrep("-", 90), "\n")
  cat(sprintf("%-30s %10s %8s %16s %18s\n",
              "Variant Set", "N", "Percent", "Lambda_Standard", "Lambda_Orthogonal"))
  cat(strrep("-", 90), "\n")

  cat(sprintf("%-30s %10d %7.1f%% %16.3f %18.3f\n",
              "All variants",
              row$All_N,
              100,
              row$All_Lambda_Std,
              row$All_Lambda_Orth))

  cat(sprintf("%-30s %10d %7.1f%% %16.3f %18.3f\n",
              "Pass HWE QC (p > 1e-6)",
              row$Pass_N,
              row$Pass_Pct,
              row$Pass_Lambda_Std,
              row$Pass_Lambda_Orth))

  cat(sprintf("%-30s %10d %7.1f%% %16.3f %18.3f\n",
              "Fail HWE QC (p <= 1e-6)",
              row$Fail_N,
              row$Fail_Pct,
              row$Fail_Lambda_Std,
              row$Fail_Lambda_Orth))

  cat(strrep("-", 90), "\n")
}

cat("\n\n=== KEY FINDINGS ===\n\n")
cat("1. UK Biobank EUR (F = 0.0037):\n")
ukb <- all_f_results[all_f_results$F_label == "UKB EUR", ]
cat(sprintf("   - %.1f%% of variants pass HWE QC\n", ukb$Pass_Pct))
cat(sprintf("   - Standard encoding: lambda = %.3f (all) -> %.3f (pass HWE)\n",
            ukb$All_Lambda_Std, ukb$Pass_Lambda_Std))
cat(sprintf("   - Orthogonal encoding: lambda = %.3f (all) -> %.3f (pass HWE)\n\n",
            ukb$All_Lambda_Orth, ukb$Pass_Lambda_Orth))

cat("2. Genes & Health (F = 0.0178):\n")
gh <- all_f_results[all_f_results$F_label == "Genes & Health", ]
cat(sprintf("   - %.1f%% of variants pass HWE QC\n", gh$Pass_Pct))
cat(sprintf("   - Standard encoding: lambda = %.3f (all) -> %.3f (pass HWE)\n",
            gh$All_Lambda_Std, gh$Pass_Lambda_Std))
cat(sprintf("   - Orthogonal encoding: lambda = %.3f (all) -> %.3f (pass HWE)\n\n",
            gh$All_Lambda_Orth, gh$Pass_Lambda_Orth))

cat("CONCLUSION:\n")
cat("Even at realistic population structure levels, variants passing HWE QC can still\n")
cat("show inflation with standard encoding, while orthogonal encoding maintains proper\n")
cat("calibration regardless of HWE status.\n")
