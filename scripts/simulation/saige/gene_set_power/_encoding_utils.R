#!/usr/bin/env Rscript
# Orthogonal genotype encoding utilities
# Implements additive and dominance encodings without assuming HWE

#' Compute genotype proportions from genotype data
#'
#' @param genotypes Numeric vector of genotypes coded as 0, 1, 2
#' @return Named vector with proportions: r (0/0), h (0/1), a (1/1)
#' @examples
#' genotypes <- c(0, 0, 1, 1, 2)
#' compute_genotype_proportions(genotypes)
#' # Returns: r=0.4, h=0.4, a=0.2
compute_genotype_proportions <- function(genotypes) {
    if (!all(genotypes %in% c(0, 1, 2))) {
        stop("Genotypes must be coded as 0, 1, or 2")
    }

    n <- length(genotypes)
    r <- sum(genotypes == 0) / n  # Homozygous reference (wildtype)
    h <- sum(genotypes == 1) / n  # Heterozygous (monoallelic)
    a <- sum(genotypes == 2) / n  # Homozygous alternate (biallelic)

    # Verify proportions sum to 1
    if (abs(r + h + a - 1.0) > 1e-10) {
        stop("Genotype proportions do not sum to 1")
    }

    return(c(r = r, h = h, a = a))
}


#' Compute orthogonal additive encoding (X^A)
#'
#' Implements Equation 2 from the manuscript
#' Returns a 3-vector representing standardized genotype-class deviations
#'
#' @param r Proportion of homozygous reference genotypes (0/0)
#' @param h Proportion of heterozygous genotypes (0/1)
#' @param a Proportion of homozygous alternate genotypes (1/1)
#' @return Numeric vector of length 3: X^A for genotypes [0, 1, 2]
compute_additive_encoding <- function(r, h, a) {
    # Verify proportions sum to 1
    if (abs(r + h + a - 1.0) > 1e-10) {
        stop("Genotype proportions must sum to 1")
    }

    # X^A' (unscaled additive component)
    X_A_prime <- c(0, 1, 2)

    # X^0 (intercept, already normalized)
    X_0 <- c(1, 1, 1)

    # Inner product under genotype proportions: <X, Y> = r*X[1]*Y[1] + h*X[2]*Y[2] + a*X[3]*Y[3]
    inner_product <- function(X, Y) {
        r * X[1] * Y[1] + h * X[2] * Y[2] + a * X[3] * Y[3]
    }

    # Project out intercept: X^A' - <X^A', X^0> * X^0
    proj_intercept <- inner_product(X_A_prime, X_0)
    X_A_centered <- X_A_prime - proj_intercept * X_0

    # Compute norm: ||X||^2 = <X, X>
    norm_squared <- inner_product(X_A_centered, X_A_centered)

    # Normalize to unit variance
    X_A <- X_A_centered / sqrt(norm_squared)

    # Verify unit variance
    var_X_A <- inner_product(X_A, X_A)
    if (abs(var_X_A - 1.0) > 1e-10) {
        warning("Additive encoding does not have unit variance")
    }

    return(X_A)
}


#' Compute orthogonal dominance encoding (X^D)
#'
#' Implements Equation 3 from the manuscript
#' Returns a 3-vector orthogonal to X^A with unit variance
#'
#' @param r Proportion of homozygous reference genotypes (0/0)
#' @param h Proportion of heterozygous genotypes (0/1)
#' @param a Proportion of homozygous alternate genotypes (1/1)
#' @param X_A Optional pre-computed additive encoding (for efficiency)
#' @return Numeric vector of length 3: X^D for genotypes [0, 1, 2]
compute_dominance_encoding <- function(r, h, a, X_A = NULL) {
    # Verify proportions sum to 1
    if (abs(r + h + a - 1.0) > 1e-10) {
        stop("Genotype proportions must sum to 1")
    }

    # X^D' (unscaled dominance component)
    X_D_prime <- c(0, 1, 0)

    # X^0 (intercept)
    X_0 <- c(1, 1, 1)

    # Compute X^A if not provided
    if (is.null(X_A)) {
        X_A <- compute_additive_encoding(r, h, a)
    }

    # Inner product under genotype proportions
    inner_product <- function(X, Y) {
        r * X[1] * Y[1] + h * X[2] * Y[2] + a * X[3] * Y[3]
    }

    # Project out intercept and additive component (Gram-Schmidt)
    proj_intercept <- inner_product(X_D_prime, X_0)
    proj_additive <- inner_product(X_D_prime, X_A)

    X_D_orthogonal <- X_D_prime - proj_intercept * X_0 - proj_additive * X_A

    # Compute norm
    norm_squared <- inner_product(X_D_orthogonal, X_D_orthogonal)

    # Normalize to unit variance
    X_D <- X_D_orthogonal / sqrt(norm_squared)

    # Check for numerical issues
    if (any(is.na(X_D)) || any(!is.finite(X_D))) {
        stop(sprintf("Dominance encoding contains NA or infinite values. Genotype props: r=%.4f, h=%.4f, a=%.4f",
                    r, h, a))
    }

    # Verify unit variance
    var_X_D <- inner_product(X_D, X_D)
    if (is.na(var_X_D) || !is.finite(var_X_D)) {
        stop(sprintf("Dominance encoding variance is NA or infinite. Genotype props: r=%.4f, h=%.4f, a=%.4f",
                    r, h, a))
    }

    if (abs(var_X_D - 1.0) > 1e-6) {
        warning(sprintf("Dominance encoding variance not unit (%.6f). Genotype props: r=%.4f, h=%.4f, a=%.4f",
                       var_X_D, r, h, a))
    }

    # Verify orthogonality to X^A
    orthogonality <- inner_product(X_A, X_D)
    if (abs(orthogonality) > 1e-6) {
        warning(sprintf("Additive and dominance encodings not orthogonal (%.6e). Genotype props: r=%.4f, h=%.4f, a=%.4f",
                       orthogonality, r, h, a))
    }

    return(X_D)
}


#' Calculate MAF from genotype proportions
#'
#' @param r Proportion of homozygous reference (0/0)
#' @param h Proportion of heterozygous (0/1)
#' @param a Proportion of homozygous alternate (1/1)
#' @return Minor allele frequency
calculate_maf <- function(r, h, a) {
    # Allele frequency = (h + 2*a) / 2
    af <- (h + 2 * a) / 2
    # MAF is the minimum of af and 1-af
    maf <- min(af, 1 - af)
    return(maf)
}


#' Standardize individual genotypes using orthogonal encodings
#'
#' Maps genotype data (0, 1, 2) to standardized additive and dominance scores
#'
#' @param genotypes Numeric vector of genotypes coded as 0, 1, 2
#' @param X_A 3-vector of additive encoding for genotypes [0, 1, 2]
#' @param X_D 3-vector of dominance encoding for genotypes [0, 1, 2]
#' @return Data frame with columns: X_A (additive scores), X_D (dominance scores)
#' @examples
#' genotypes <- c(0, 1, 2, 1, 0)
#' props <- compute_genotype_proportions(genotypes)
#' X_A <- compute_additive_encoding(props['r'], props['h'], props['a'])
#' X_D <- compute_dominance_encoding(props['r'], props['h'], props['a'], X_A)
#' standardize_genotypes(genotypes, X_A, X_D)
standardize_genotypes <- function(genotypes, X_A, X_D) {
    if (!all(genotypes %in% c(0, 1, 2))) {
        stop("Genotypes must be coded as 0, 1, or 2")
    }

    # Map genotypes to encodings
    # genotype = 0 -> X_A[1], X_D[1]
    # genotype = 1 -> X_A[2], X_D[2]
    # genotype = 2 -> X_A[3], X_D[3]

    X_A_individual <- X_A[genotypes + 1]  # +1 for 1-based indexing
    X_D_individual <- X_D[genotypes + 1]

    return(data.frame(
        X_A = X_A_individual,
        X_D = X_D_individual
    ))
}


#' Verify orthogonal encoding properties
#'
#' Diagnostic function to check all mathematical properties
#'
#' @param r Proportion of homozygous reference genotypes
#' @param h Proportion of heterozygous genotypes
#' @param a Proportion of homozygous alternate genotypes
#' @return List of verification results
verify_encoding_properties <- function(r, h, a) {
    X_A <- compute_additive_encoding(r, h, a)
    X_D <- compute_dominance_encoding(r, h, a, X_A)

    # Inner product
    inner_product <- function(X, Y) {
        r * X[1] * Y[1] + h * X[2] * Y[2] + a * X[3] * Y[3]
    }

    # Check all properties
    results <- list(
        genotype_props_sum = r + h + a,
        X_A_variance = inner_product(X_A, X_A),
        X_D_variance = inner_product(X_D, X_D),
        orthogonality = inner_product(X_A, X_D),
        X_A_mean = inner_product(X_A, c(1, 1, 1)),
        X_D_mean = inner_product(X_D, c(1, 1, 1))
    )

    cat("=== Encoding Verification ===\n")
    cat(sprintf("Genotype proportions sum to 1: %.10f\n", results$genotype_props_sum))
    cat(sprintf("X^A has unit variance: %.10f\n", results$X_A_variance))
    cat(sprintf("X^D has unit variance: %.10f\n", results$X_D_variance))
    cat(sprintf("X^A and X^D are orthogonal: %.10e\n", results$orthogonality))
    cat(sprintf("X^A has zero mean: %.10e\n", results$X_A_mean))
    cat(sprintf("X^D has zero mean: %.10e\n", results$X_D_mean))

    all_pass <- abs(results$genotype_props_sum - 1.0) < 1e-10 &&
                abs(results$X_A_variance - 1.0) < 1e-10 &&
                abs(results$X_D_variance - 1.0) < 1e-10 &&
                abs(results$orthogonality) < 1e-10 &&
                abs(results$X_A_mean) < 1e-10 &&
                abs(results$X_D_mean) < 1e-10

    cat(sprintf("\nAll checks passed: %s\n", all_pass))

    return(results)
}


# Example usage and testing
if (!interactive()) {
    cat("Testing encoding utilities...\n\n")

    # Test case 1: Typical rare variant
    cat("Test 1: Rare variant (r=0.90, h=0.08, a=0.02)\n")
    verify_encoding_properties(0.90, 0.08, 0.02)

    cat("\n")

    # Test case 2: Common variant
    cat("Test 2: Common variant (r=0.25, h=0.50, a=0.25)\n")
    verify_encoding_properties(0.25, 0.50, 0.25)

    cat("\n")

    # Test case 3: Very rare variant
    cat("Test 3: Very rare variant (r=0.99, h=0.009, a=0.001)\n")
    verify_encoding_properties(0.99, 0.009, 0.001)
}
