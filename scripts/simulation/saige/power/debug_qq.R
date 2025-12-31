#!/usr/bin/env Rscript
# Quick diagnostic script to check what's going on with the QQ plot

library(data.table)
library(dplyr)

# Load one file as a test
test_file <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/saige_power_analysis/UKB.auto.arch_additive_h2_001_M_10_rep_3.additive.eur.txt.gz"

df <- fread(test_file)

# Check p-value distribution
cat("=== P-value Summary ===\n")
cat("Class:", class(df$p.value), "\n")
print(summary(df$p.value))

cat("\n=== Top 10 most significant p-values ===\n")
print(head(df[order(df$p.value), .(MarkerID, p.value, BETA, SE)], 10))

cat("\n=== Bottom 10 least significant p-values ===\n")
print(head(df[order(-df$p.value), .(MarkerID, p.value, BETA, SE)], 10))

cat("\n=== P-value quantiles ===\n")
print(quantile(df$p.value, probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)))

cat("\n=== Number of significant p-values ===\n")
cat("P < 0.05:", sum(df$p.value < 0.05, na.rm = TRUE), "\n")
cat("P < 0.01:", sum(df$p.value < 0.01, na.rm = TRUE), "\n")
cat("P < 0.001:", sum(df$p.value < 0.001, na.rm = TRUE), "\n")
cat("P < 1e-5:", sum(df$p.value < 1e-5, na.rm = TRUE), "\n")

# Calculate -log10(p) range
df$log10p <- -log10(df$p.value)
cat("\n=== -log10(P) range ===\n")
print(summary(df$log10p))
cat("Max -log10(P):", max(df$log10p, na.rm = TRUE), "\n")
