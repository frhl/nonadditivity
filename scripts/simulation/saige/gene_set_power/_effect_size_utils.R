#!/usr/bin/env Rscript
# Effect size calculation and scaling utilities
# Implements genetic architecture definitions and heritability scaling

source("_encoding_utils.R")


#' Define genetic architecture
#'
#' Returns expected phenotypic values for each genotype class [0, 1, 2]
#' Based on different genetic architectures
#'
#' @param architecture String: "additive", "partially_recessive_X", "recessive", or "overdominant"
#'   For partially recessive, specify heterozygote value after underscore (e.g., "partially_recessive_0.1")
#' @return Numeric vector of length 3: g_j = [g_0, g_1, g_2]
#' @examples
#' define_genetic_architecture("additive")  # Returns [0, 1, 2]
#' define_genetic_architecture("recessive")  # Returns [0, 0, 2]
#' define_genetic_architecture("partially_recessive_0.1")  # Returns [0, 0.1, 2]
#' define_genetic_architecture("partially_recessive_0.2")  # Returns [0, 0.2, 2]
define_genetic_architecture <- function(architecture) {
    # Check if this is a partially recessive architecture with specified heterozygote value
    if (grepl("^partially_recessive_", architecture)) {
        # Extract heterozygote value from architecture name
        het_value_str <- sub("^partially_recessive_", "", architecture)
        het_value <- as.numeric(het_value_str)

        if (is.na(het_value) || het_value < 0 || het_value > 2) {
            stop(sprintf("Invalid heterozygote value for partially_recessive: %s (must be between 0 and 2)",
                        het_value_str))
        }

        g_j <- c(0, het_value, 2)
        return(g_j)
    }

    # Standard architectures
    valid_architectures <- c("additive", "partially_recessive", "recessive", "overdominant")

    if (!architecture %in% valid_architectures) {
        stop(sprintf("Invalid architecture: %s. Must be one of: %s or 'partially_recessive_X' where X is the heterozygote value",
                    architecture, paste(valid_architectures, collapse = ", ")))
    }

    g_j <- switch(architecture,
        "additive" = c(0, 1, 2),
        "partially_recessive" = c(0, 0.1, 2),  # Default for backward compatibility
        "recessive" = c(0, 0, 2),
        "overdominant" = c(0, 2, 0)
    )

    return(g_j)
}


#' Solve for raw effect sizes from genetic architecture
#'
#' Solves the system: g_j = intercept + beta_A_raw * X^A + beta_D_raw * X^D
#' Uses projections since X^A and X^D are orthogonal
#' The intercept is absorbed (doesn't contribute to variance)
#'
#' @param g_j Genetic architecture vector [g_0, g_1, g_2]
#' @param X_A Additive encoding vector (length 3)
#' @param X_D Dominance encoding vector (length 3)
#' @param r Proportion of homozygous reference genotypes
#' @param h Proportion of heterozygous genotypes
#' @param a Proportion of homozygous alternate genotypes
#' @return Named list with beta_A_raw, beta_D_raw, h2_A_raw, h2_D_raw, h2_total_raw
solve_raw_effects <- function(g_j, X_A, X_D, r, h, a) {
    # Inner product under genotype proportions
    inner_product <- function(X, Y) {
        r * X[1] * Y[1] + h * X[2] * Y[2] + a * X[3] * Y[3]
    }

    # Center g_j to remove intercept (which doesn't contribute to variance)
    X_0 <- c(1, 1, 1)
    intercept <- inner_product(g_j, X_0)
    g_j_centered <- g_j - intercept * X_0

    # Since X^A and X^D are orthogonal with unit variance and zero mean:
    # beta_A_raw = <g_j_centered, X^A>
    # beta_D_raw = <g_j_centered, X^D>
    beta_A_raw <- inner_product(g_j_centered, X_A)
    beta_D_raw <- inner_product(g_j_centered, X_D)

    # Unscaled heritability components (since Var(X^A) = Var(X^D) = 1)
    h2_A_raw <- beta_A_raw^2
    h2_D_raw <- beta_D_raw^2
    h2_total_raw <- h2_A_raw + h2_D_raw

    # Verify reconstruction (including intercept)
    g_j_reconstructed <- intercept * X_0 + beta_A_raw * X_A + beta_D_raw * X_D
    reconstruction_error <- max(abs(g_j - g_j_reconstructed))

    if (reconstruction_error > 1e-6) {
        warning(sprintf("Genetic architecture reconstruction error: %.6e", reconstruction_error))
    }

    return(list(
        beta_A_raw = beta_A_raw,
        beta_D_raw = beta_D_raw,
        h2_A_raw = h2_A_raw,
        h2_D_raw = h2_D_raw,
        h2_total_raw = h2_total_raw,
        intercept = intercept,
        reconstruction_error = reconstruction_error
    ))
}


#' Scale effect sizes to achieve target heritability
#'
#' Applies uniform scaling to both beta_A and beta_D to achieve target h2
#' while preserving the natural additive/dominance ratio
#'
#' @param beta_A_raw Unscaled additive effect size
#' @param beta_D_raw Unscaled dominance effect size
#' @param h2_target Target total heritability for this variant
#' @param h2_total_raw Unscaled total heritability (beta_A_raw^2 + beta_D_raw^2)
#' @return Named list with beta_A, beta_D, h2_A, h2_D, h2_total
scale_to_target_heritability <- function(beta_A_raw, beta_D_raw, h2_target, h2_total_raw = NULL) {
    # Calculate h2_total_raw if not provided
    if (is.null(h2_total_raw)) {
        h2_total_raw <- beta_A_raw^2 + beta_D_raw^2
    }

    # Avoid division by zero
    if (h2_total_raw < 1e-12) {
        warning("Raw heritability is near zero. Returning zero effect sizes.")
        return(list(
            beta_A = 0,
            beta_D = 0,
            h2_A = 0,
            h2_D = 0,
            h2_total = 0,
            scaling_factor = 0
        ))
    }

    # Uniform scaling factor
    scaling_factor <- sqrt(h2_target / h2_total_raw)

    # Apply scaling
    beta_A <- beta_A_raw * scaling_factor
    beta_D <- beta_D_raw * scaling_factor

    # Scaled heritability components
    h2_A <- beta_A^2
    h2_D <- beta_D^2
    h2_total <- h2_A + h2_D

    # Verify target achieved
    if (abs(h2_total - h2_target) > 1e-10) {
        warning(sprintf("Target heritability not achieved: %.6f vs %.6f", h2_total, h2_target))
    }

    return(list(
        beta_A = beta_A,
        beta_D = beta_D,
        h2_A = h2_A,
        h2_D = h2_D,
        h2_total = h2_total,
        scaling_factor = scaling_factor
    ))
}


#' Compute effect sizes for a single variant
#'
#' Complete pipeline: architecture -> raw effects -> scaled effects
#'
#' @param architecture Genetic architecture type
#' @param r Proportion of homozygous reference genotypes
#' @param h Proportion of heterozygous genotypes
#' @param a Proportion of homozygous alternate genotypes
#' @param h2_target Target total heritability for this variant
#' @return Named list with all effect sizes and heritability components
compute_variant_effects <- function(architecture, r, h, a, h2_target) {
    # Compute encodings
    X_A <- compute_additive_encoding(r, h, a)
    X_D <- compute_dominance_encoding(r, h, a, X_A)

    # Define architecture
    g_j <- define_genetic_architecture(architecture)

    # Solve for raw effects
    raw_effects <- solve_raw_effects(g_j, X_A, X_D, r, h, a)

    # Scale to target heritability
    scaled_effects <- scale_to_target_heritability(
        raw_effects$beta_A_raw,
        raw_effects$beta_D_raw,
        h2_target,
        raw_effects$h2_total_raw
    )

    # Combine results
    return(c(
        list(
            architecture = architecture,
            r = r, h = h, a = a,
            h2_target = h2_target,
            X_A = X_A,
            X_D = X_D,
            g_j = g_j
        ),
        raw_effects,
        scaled_effects
    ))
}


#' Simulate phenotypes from causal variants
#'
#' Implements Equation 5 from the manuscript
#'
#' @param genotype_matrix Matrix of genotypes (rows = individuals, cols = variants), coded 0/1/2
#' @param effect_list List of effect size results from compute_variant_effects()
#' @param h2_total Total heritability across all variants
#' @param seed Random seed for reproducibility
#' @return Numeric vector of simulated phenotypes
simulate_phenotypes <- function(genotype_matrix, effect_list, h2_total, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n_individuals <- nrow(genotype_matrix)
    n_variants <- ncol(genotype_matrix)

    if (length(effect_list) != n_variants) {
        stop("Number of variants in effect_list must match genotype_matrix columns")
    }

    # Initialize phenotype with genetic component
    y_genetic <- rep(0, n_individuals)

    # Add contribution from each causal variant
    for (j in 1:n_variants) {
        effects_j <- effect_list[[j]]

        # Standardize genotypes for this variant
        genotypes_j <- genotype_matrix[, j]
        standardized <- standardize_genotypes(genotypes_j, effects_j$X_A, effects_j$X_D)

        # Add genetic contribution: beta_A * X^A + beta_D * X^D
        y_genetic <- y_genetic + effects_j$beta_A * standardized$X_A + effects_j$beta_D * standardized$X_D
    }

    # Add environmental noise
    environmental_noise <- rnorm(n_individuals, mean = 0, sd = sqrt(1 - h2_total))
    y <- y_genetic + environmental_noise

    # Verify realized heritability
    realized_h2 <- var(y_genetic) / var(y)
    if (abs(realized_h2 - h2_total) > 0.01) {  # Allow 1% tolerance
        warning(sprintf("Realized heritability (%.4f) differs from target (%.4f)",
                       realized_h2, h2_total))
    }

    return(y)
}


# Testing and examples
if (!interactive()) {
    cat("Testing effect size utilities...\n\n")

    # Test case: Additive architecture with rare variant
    cat("=== Test 1: Additive architecture ===\n")
    r <- 0.90
    h <- 0.08
    a <- 0.02
    h2_target <- 0.05

    effects <- compute_variant_effects("additive", r, h, a, h2_target)

    cat(sprintf("Architecture: %s\n", effects$architecture))
    cat(sprintf("Target h2: %.4f\n", effects$h2_target))
    cat(sprintf("Raw effects: beta_A = %.4f, beta_D = %.4f\n",
               effects$beta_A_raw, effects$beta_D_raw))
    cat(sprintf("Raw h2: total = %.4f (A = %.4f, D = %.4f)\n",
               effects$h2_total_raw, effects$h2_A_raw, effects$h2_D_raw))
    cat(sprintf("Scaled effects: beta_A = %.4f, beta_D = %.4f\n",
               effects$beta_A, effects$beta_D))
    cat(sprintf("Scaled h2: total = %.4f (A = %.4f, D = %.4f)\n",
               effects$h2_total, effects$h2_A, effects$h2_D))
    cat(sprintf("Scaling factor: %.4f\n", effects$scaling_factor))

    cat("\n=== Test 2: Recessive architecture ===\n")
    effects_rec <- compute_variant_effects("recessive", r, h, a, h2_target)

    cat(sprintf("Architecture: %s\n", effects_rec$architecture))
    cat(sprintf("Raw h2: total = %.4f (A = %.4f, D = %.4f)\n",
               effects_rec$h2_total_raw, effects_rec$h2_A_raw, effects_rec$h2_D_raw))
    cat(sprintf("Scaled h2: total = %.4f (A = %.4f, D = %.4f)\n",
               effects_rec$h2_total, effects_rec$h2_A, effects_rec$h2_D))
    cat(sprintf("Ratio h2_A : h2_D = %.4f : %.4f\n",
               effects_rec$h2_A / effects_rec$h2_total,
               effects_rec$h2_D / effects_rec$h2_total))

    cat("\n=== Test 3: Overdominant architecture ===\n")
    effects_od <- compute_variant_effects("overdominant", r, h, a, h2_target)

    cat(sprintf("Architecture: %s\n", effects_od$architecture))
    cat(sprintf("Raw h2: total = %.4f (A = %.4f, D = %.4f)\n",
               effects_od$h2_total_raw, effects_od$h2_A_raw, effects_od$h2_D_raw))
    cat(sprintf("Scaled h2: total = %.4f (A = %.4f, D = %.4f)\n",
               effects_od$h2_total, effects_od$h2_A, effects_od$h2_D))

    cat("\n=== Test 4: Phenotype simulation ===\n")
    # Simulate 1000 individuals, 2 causal variants
    set.seed(123)
    n <- 1000
    genotype_matrix <- matrix(
        sample(c(0, 1, 2), n * 2, replace = TRUE, prob = c(r, h, a)),
        nrow = n, ncol = 2
    )

    # Each variant contributes h2 = 0.025, total h2 = 0.05
    effects1 <- compute_variant_effects("additive", r, h, a, 0.025)
    effects2 <- compute_variant_effects("recessive", r, h, a, 0.025)

    y <- simulate_phenotypes(genotype_matrix, list(effects1, effects2), h2_total = 0.05, seed = 456)

    cat(sprintf("Phenotype mean: %.4f (expected ≈ 0)\n", mean(y)))
    cat(sprintf("Phenotype variance: %.4f (expected ≈ 1)\n", var(y)))
    cat(sprintf("Phenotype range: [%.4f, %.4f]\n", min(y), max(y)))

    cat("\nAll tests completed!\n")
}
