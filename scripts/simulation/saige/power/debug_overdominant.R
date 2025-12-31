#!/usr/bin/env Rscript
# Debug overdominant architecture - check if dominance encoding captures signal

library(data.table)
library(dplyr)

# Load metadata to see what effect sizes were simulated
metadata_file <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/simulated_phenotypes/power/UKB.eur.power_phenos.051025.metadata.tsv.gz"

if (!file.exists(metadata_file)) {
  cat("Metadata file not found. Looking for it...\n")
  # Try to find it
  system("find /Users/flassen/Projects/11_wes_ko_ukbb_nexus -name '*power*metadata*' 2>/dev/null | head -5")
  stop("Please provide correct path to metadata file")
}

meta <- fread(metadata_file)

cat("=== Overdominant Architecture Metadata ===\n\n")

# Filter for overdominant, h2=0.01, M=10
od_meta <- meta %>%
  filter(architecture == "overdominant",
         h2_total == 0.01,
         M == 10,
         replicate == 4)

cat("Number of causal variants:", nrow(od_meta), "\n\n")

if (nrow(od_meta) > 0) {
  cat("Effect sizes for overdominant variants:\n")
  print(od_meta %>%
    select(variant_id, maf, beta_A, beta_D, h2_A, h2_D, h2_total) %>%
    head(10))

  cat("\n=== Summary ===\n")
  cat(sprintf("Mean |beta_A|: %.4f\n", mean(abs(od_meta$beta_A))))
  cat(sprintf("Mean |beta_D|: %.4f\n", mean(abs(od_meta$beta_D))))
  cat(sprintf("Mean h2_A: %.4f (%.1f%% of total)\n",
              mean(od_meta$h2_A),
              100 * mean(od_meta$h2_A / od_meta$h2_total)))
  cat(sprintf("Mean h2_D: %.4f (%.1f%% of total)\n",
              mean(od_meta$h2_D),
              100 * mean(od_meta$h2_D / od_meta$h2_total)))
}

# Now load actual results
results_file <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/saige_power_analysis/UKB.auto.arch_overdominant_h2_001_M_10_rep_4.dominance.eur.txt.gz"

if (file.exists(results_file)) {
  cat("\n=== SAIGE Results for Overdominant (Dominance Model) ===\n\n")

  results <- fread(results_file)
  results$p.value <- as.numeric(results$p.value)

  # Get the causal variants
  if (nrow(od_meta) > 0) {
    causal_ids <- od_meta$variant_id
    causal_results <- results %>% filter(MarkerID %in% causal_ids)

    cat("Causal variants results (dominance model):\n")
    print(causal_results %>%
      arrange(p.value) %>%
      select(MarkerID, p.value, BETA, SE, AF_Allele2) %>%
      head(10))

    cat("\n=== Comparison: Expected vs Observed ===\n")
    comparison <- merge(
      od_meta %>% select(variant_id, beta_D, h2_D, maf),
      causal_results %>% select(MarkerID, p.value, BETA, SE, AF_Allele2),
      by.x = "variant_id",
      by.y = "MarkerID"
    )

    comparison$log10p <- -log10(comparison$p.value)
    comparison <- comparison %>% arrange(desc(log10p))

    print(comparison)

    cat("\nNumber of significant causal variants (P < 1e-5):",
        sum(comparison$p.value < 1e-5), "/", nrow(comparison), "\n")
  }
} else {
  cat("\nResults file not found:", results_file, "\n")
}

# Also check additive model results for comparison
additive_results_file <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/saige_power_analysis/UKB.auto.arch_overdominant_h2_001_M_10_rep_4.additive.eur.txt.gz"

if (file.exists(additive_results_file) && nrow(od_meta) > 0) {
  cat("\n=== SAIGE Results for Overdominant (Additive Model) ===\n\n")

  add_results <- fread(additive_results_file)
  add_results$p.value <- as.numeric(add_results$p.value)

  causal_add <- add_results %>% filter(MarkerID %in% od_meta$variant_id)

  cat("Causal variants results (additive model):\n")
  print(causal_add %>%
    arrange(p.value) %>%
    select(MarkerID, p.value, BETA, SE) %>%
    head(10))

  cat("\nNumber of significant causal variants (P < 1e-5):",
      sum(causal_add$p.value < 1e-5), "/", nrow(causal_add), "\n")
}
