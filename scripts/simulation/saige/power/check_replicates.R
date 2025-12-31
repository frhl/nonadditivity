#!/usr/bin/env Rscript
# Check how replicates are combined

library(data.table)
library(dplyr)

# Quick check of what's in the additive data
results_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/saige_power_analysis"

H2 <- 0.01
N_CAUSAL <- 10

# Load just the additive architecture files
files <- list.files(results_dir, pattern = "arch_additive.*h2_001.*M_10.*\\.txt\\.gz$", full.names = TRUE)

cat("=== Files matching additive, h2=0.01, M=10 ===\n")
print(basename(files))

cat("\n=== Breakdown by model and replicate ===\n")
for (f in files) {
  filename <- basename(f)
  parts <- strsplit(filename, "\\.")[[1]]

  # Get replicate
  param_string <- parts[3]
  param_parts <- strsplit(param_string, "_")[[1]]
  rep_idx <- which(param_parts == "rep")
  replicate <- param_parts[rep_idx + 1]

  # Get model
  model <- parts[4]

  cat(sprintf("Replicate %s, Model %s: %s\n", replicate, model, filename))
}

cat("\n=== Loading one replicate to check significant variants ===\n")
# Load replicate 1, additive model
test_file <- files[grep("rep_1.*additive", files)[1]]
cat("File:", basename(test_file), "\n")

df <- fread(test_file)
df$p.value <- as.numeric(df$p.value)

sig_variants <- df %>%
  filter(p.value < 1e-5) %>%
  arrange(p.value) %>%
  select(MarkerID, p.value, BETA, SE)

cat("\nSignificant variants (P < 1e-5):", nrow(sig_variants), "\n")
print(sig_variants)

cat("\n=== If we combine all replicates and models ===\n")
cat("Total files:", length(files), "\n")
cat("If 10 causal variants × ", length(files), " files = ", 10 * length(files), " total tests\n")
cat("But unique variants across replicates should be ~10\n")
