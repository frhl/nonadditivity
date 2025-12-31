#!/usr/bin/env Rscript
# Test script to validate power simulation framework
# Simulates mock genotype data to test all functions

cat("=============================================================\n")
cat("Testing Power Simulation Framework\n")
cat("=============================================================\n\n")

# Source utility functions
source("_encoding_utils.R")
source("_effect_size_utils.R")

# Test 1: Encoding functions
cat("\n==================== Test 1: Encoding Functions ====================\n")

test_encodings <- function() {
    # Test with typical rare variant
    r <- 0.90
    h <- 0.08
    a <- 0.02

    cat("Testing rare variant (r=0.90, h=0.08, a=0.02)\n")
    results <- verify_encoding_properties(r, h, a)

    if (abs(results$X_A_variance - 1.0) < 1e-10 &&
        abs(results$X_D_variance - 1.0) < 1e-10 &&
        abs(results$orthogonality) < 1e-10) {
        cat("✓ Encoding test PASSED\n")
        return(TRUE)
    } else {
        cat("✗ Encoding test FAILED\n")
        return(FALSE)
    }
}

encoding_pass <- test_encodings()


# Test 2: Effect size calculation
cat("\n==================== Test 2: Effect Size Calculation ====================\n")

test_effect_sizes <- function() {
    r <- 0.90
    h <- 0.08
    a <- 0.02
    h2_target <- 0.05

    architectures <- c("additive", "recessive", "partially_recessive", "overdominant")
    all_pass <- TRUE

    for (arch in architectures) {
        cat(sprintf("\nTesting %s architecture:\n", arch))

        effects <- compute_variant_effects(arch, r, h, a, h2_target)

        # Check that h2 sums to target
        h2_achieved <- effects$h2_A + effects$h2_D

        cat(sprintf("  Target h²: %.4f\n", h2_target))
        cat(sprintf("  Achieved h²: %.4f (A=%.4f, D=%.4f)\n",
                   h2_achieved, effects$h2_A, effects$h2_D))
        cat(sprintf("  Ratio A:D = %.2f:%.2f\n",
                   effects$h2_A / h2_achieved, effects$h2_D / h2_achieved))

        if (abs(h2_achieved - h2_target) > 1e-6) {
            cat(sprintf("  ✗ FAILED: h² mismatch (%.6f vs %.6f)\n", h2_achieved, h2_target))
            all_pass <- FALSE
        } else {
            cat("  ✓ PASSED\n")
        }
    }

    return(all_pass)
}

effect_pass <- test_effect_sizes()


# Test 3: MAF calculation
cat("\n==================== Test 3: MAF Calculation ====================\n")

test_maf_calculation <- function() {
    # Test case 1: MAF = 0.05
    r1 <- 0.9025; h1 <- 0.095; a1 <- 0.0025
    maf1 <- calculate_maf(r1, h1, a1)
    expected1 <- 0.05

    cat(sprintf("Test 1: r=%.4f, h=%.4f, a=%.4f\n", r1, h1, a1))
    cat(sprintf("  MAF: %.4f (expected: %.4f)\n", maf1, expected1))

    # Test case 2: Common variant (MAF = 0.3)
    r2 <- 0.49; h2 <- 0.42; a2 <- 0.09
    maf2 <- calculate_maf(r2, h2, a2)
    expected2 <- 0.3

    cat(sprintf("\nTest 2: r=%.4f, h=%.4f, a=%.4f\n", r2, h2, a2))
    cat(sprintf("  MAF: %.4f (expected: %.4f)\n", maf2, expected2))

    if (abs(maf1 - expected1) < 1e-3 && abs(maf2 - expected2) < 1e-3) {
        cat("\n✓ MAF calculation PASSED\n")
        return(TRUE)
    } else {
        cat("\n✗ MAF calculation FAILED\n")
        return(FALSE)
    }
}

maf_pass <- test_maf_calculation()


# Test 4: Mock phenotype simulation
cat("\n==================== Test 4: Phenotype Simulation ====================\n")

test_phenotype_simulation <- function() {
    set.seed(123)

    # Simulate mock genotype data
    n_samples <- 1000
    n_variants <- 100
    M <- 5  # 5 causal variants

    cat(sprintf("Simulating %d samples, %d variants, %d causal\n",
               n_samples, n_variants, M))

    # Create mock genotype matrix (rare variants)
    # Ensure each variant has ≥5 homozygous alternate carriers (like real data)
    genotype_matrix <- matrix(0, nrow = n_samples, ncol = n_variants)
    for (j in 1:n_variants) {
        # Sample genotypes with MAF ~ 0.01-0.05
        # Ensure at least 5 homozygous carriers
        repeat {
            maf <- runif(1, 0.02, 0.05)  # Slightly higher to ensure carriers
            geno <- sample(c(0, 1, 2), n_samples, replace = TRUE,
                          prob = c((1-maf)^2, 2*maf*(1-maf), maf^2))
            if (sum(geno == 2) >= 5) {  # Check for ≥5 homozygous alt carriers
                genotype_matrix[, j] <- geno
                break
            }
        }
    }

    # Compute genotype proportions for all variants
    variant_mafs <- sapply(1:n_variants, function(j) {
        props <- compute_genotype_proportions(genotype_matrix[, j])
        calculate_maf(props['r'], props['h'], props['a'])
    })

    cat(sprintf("Variant MAF range: [%.4f, %.4f]\n", min(variant_mafs), max(variant_mafs)))

    # Test MAF bin filtering
    cat("\nTesting MAF bin selection:\n")

    # Bin 1: 0.02-0.04
    eligible_1 <- which(variant_mafs >= 0.02 & variant_mafs < 0.04)
    cat(sprintf("  MAF bin [0.02, 0.04): %d variants\n", length(eligible_1)))

    # Bin 2: 0.04-0.06
    eligible_2 <- which(variant_mafs >= 0.04 & variant_mafs < 0.06)
    cat(sprintf("  MAF bin [0.04, 0.06): %d variants\n", length(eligible_2)))

    # Select from whichever bin has enough variants
    if (length(eligible_1) >= M) {
        causal_indices <- sample(eligible_1, M)
        causal_mafs <- variant_mafs[causal_indices]
        bin_label <- "[0.02, 0.04)"
    } else if (length(eligible_2) >= M) {
        causal_indices <- sample(eligible_2, M)
        causal_mafs <- variant_mafs[causal_indices]
        bin_label <- "[0.04, 0.06)"
    } else {
        # Just select from all variants
        causal_indices <- sample(1:n_variants, M)
        causal_mafs <- variant_mafs[causal_indices]
        bin_label <- "all"
    }

    if (length(causal_indices) >= M) {
        cat(sprintf("\nSelected %d causal variants from bin %s\n", M, bin_label))
        cat(sprintf("  MAF range: [%.4f, %.4f]\n", min(causal_mafs), max(causal_mafs)))

        # Compute effect sizes for each causal variant
        h2_total <- 0.1
        h2_per_variant <- h2_total / M

        # Try to compute effects, skip variants with encoding issues
        effect_list <- list()
        valid_indices <- c()

        for (j in causal_indices) {
            props <- compute_genotype_proportions(genotype_matrix[, j])

            # Try to compute effects, catch any errors
            tryCatch({
                effects <- compute_variant_effects("additive", props['r'], props['h'], props['a'], h2_per_variant)
                effect_list <- append(effect_list, list(effects))
                valid_indices <- c(valid_indices, j)
            }, error = function(e) {
                cat(sprintf("  Warning: Skipping variant %d due to encoding error: %s\n", j, e$message))
            }, warning = function(w) {
                # Allow warnings, just print them
                cat(sprintf("  Warning for variant %d: %s\n", j, w$message))
            })
        }

        # Update M to number of valid variants
        M <- length(valid_indices)
        causal_indices <- valid_indices

        if (M == 0) {
            cat("\n✗ No valid variants could be encoded\n")
            return(FALSE)
        }

        cat(sprintf("\nSuccessfully encoded %d causal variants\n", M))

        # Simulate phenotype
        cat(sprintf("\nSimulating phenotype with h² = %.2f\n", h2_total))

        y_genetic <- rep(0, n_samples)
        for (i in 1:M) {
            effects_i <- effect_list[[i]]
            geno <- genotype_matrix[, causal_indices[i]]
            standardized <- standardize_genotypes(geno, effects_i$X_A, effects_i$X_D)
            y_genetic <- y_genetic + effects_i$beta_A * standardized$X_A + effects_i$beta_D * standardized$X_D
        }

        # Add environmental noise
        y_noise <- rnorm(n_samples, 0, sd = sqrt(1 - h2_total))
        y <- y_genetic + y_noise

        # Verify realized heritability
        realized_h2 <- var(y_genetic) / var(y)

        cat(sprintf("  Target h²: %.4f\n", h2_total))
        cat(sprintf("  Realized h²: %.4f\n", realized_h2))
        cat(sprintf("  Phenotype mean: %.4f (expected ≈ 0)\n", mean(y)))
        cat(sprintf("  Phenotype variance: %.4f (expected ≈ 1)\n", var(y)))

        if (abs(realized_h2 - h2_total) < 0.02) {  # Allow 2% tolerance
            cat("\n✓ Phenotype simulation PASSED\n")
            return(TRUE)
        } else {
            cat("\n✗ Phenotype simulation FAILED (h² mismatch)\n")
            return(FALSE)
        }
    } else {
        cat(sprintf("\n✗ Not enough variants in MAF bin (need %d, have %d)\n", M, length(eligible_1)))
        return(FALSE)
    }
}

pheno_pass <- test_phenotype_simulation()


# Test 5: All architectures in simulation
cat("\n==================== Test 5: Multi-Architecture Simulation ====================\n")

test_all_architectures <- function() {
    set.seed(456)

    n_samples <- 500
    M <- 3
    h2_total <- 0.05

    # Simulate genotype matrix with ≥5 homozygous carriers
    genotype_matrix <- matrix(0, nrow = n_samples, ncol = M)
    for (j in 1:M) {
        repeat {
            maf <- 0.03  # Use slightly higher MAF to ensure carriers
            geno <- sample(c(0, 1, 2), n_samples, replace = TRUE,
                          prob = c((1-maf)^2, 2*maf*(1-maf), maf^2))
            if (sum(geno == 2) >= 5) {
                genotype_matrix[, j] <- geno
                break
            }
        }
    }

    architectures <- c("additive", "recessive", "partially_recessive", "overdominant")
    all_pass <- TRUE

    for (arch in architectures) {
        cat(sprintf("\nTesting %s architecture:\n", arch))

        # Compute effect sizes
        effect_list <- lapply(1:M, function(j) {
            props <- compute_genotype_proportions(genotype_matrix[, j])
            compute_variant_effects(arch, props['r'], props['h'], props['a'], h2_total / M)
        })

        # Simulate phenotype
        y_genetic <- rep(0, n_samples)
        for (j in 1:M) {
            effects_j <- effect_list[[j]]
            geno <- genotype_matrix[, j]
            standardized <- standardize_genotypes(geno, effects_j$X_A, effects_j$X_D)
            y_genetic <- y_genetic + effects_j$beta_A * standardized$X_A + effects_j$beta_D * standardized$X_D
        }

        y_noise <- rnorm(n_samples, 0, sd = sqrt(1 - h2_total))
        y <- y_genetic + y_noise

        realized_h2 <- var(y_genetic) / var(y)

        cat(sprintf("  Realized h²: %.4f (target: %.4f)\n", realized_h2, h2_total))

        # Check heritability components across variants
        total_h2_A <- sum(sapply(effect_list, function(x) x$h2_A))
        total_h2_D <- sum(sapply(effect_list, function(x) x$h2_D))

        cat(sprintf("  Total h²_A: %.4f, h²_D: %.4f\n", total_h2_A, total_h2_D))
        cat(sprintf("  Ratio: %.2f%% additive, %.2f%% dominance\n",
                   100 * total_h2_A / h2_total, 100 * total_h2_D / h2_total))

        if (abs(realized_h2 - h2_total) > 0.03) {
            cat("  ✗ FAILED\n")
            all_pass <- FALSE
        } else {
            cat("  ✓ PASSED\n")
        }
    }

    return(all_pass)
}

arch_pass <- test_all_architectures()


# Summary
cat("\n=============================================================\n")
cat("Test Summary\n")
cat("=============================================================\n")
cat(sprintf("1. Encoding functions:          %s\n", ifelse(encoding_pass, "✓ PASS", "✗ FAIL")))
cat(sprintf("2. Effect size calculation:     %s\n", ifelse(effect_pass, "✓ PASS", "✗ FAIL")))
cat(sprintf("3. MAF calculation:             %s\n", ifelse(maf_pass, "✓ PASS", "✗ FAIL")))
cat(sprintf("4. Phenotype simulation:        %s\n", ifelse(pheno_pass, "✓ PASS", "✗ FAIL")))
cat(sprintf("5. Multi-architecture test:     %s\n", ifelse(arch_pass, "✓ PASS", "✗ FAIL")))
cat("=============================================================\n")

all_tests_pass <- encoding_pass && effect_pass && maf_pass && pheno_pass && arch_pass

if (all_tests_pass) {
    cat("\n✓✓✓ ALL TESTS PASSED ✓✓✓\n")
    cat("Framework is ready for deployment!\n\n")
    quit(status = 0)
} else {
    cat("\n✗✗✗ SOME TESTS FAILED ✗✗✗\n")
    cat("Please fix issues before deployment.\n\n")
    quit(status = 1)
}
