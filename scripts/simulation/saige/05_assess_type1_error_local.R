#!/usr/bin/env Rscript

options(scipen = 999)

# Required packages
packages <- c('data.table', 'optparse', 'dplyr', 'stringr')

for (p in packages) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Installing package: ", p, "\n"))
        install.packages(p, repos = "http://cran.us.r-project.org")
        library(p, character.only = TRUE)
    }
}

# Function to calculate Type I error rate and CI
track_type1 <- function(p_values, alpha) {
    # Remove NAs
    p_values <- p_values[!is.na(p_values)]
    n <- length(p_values)
    
    if (n == 0) {
        return(list(rate = NA, ci = c(NA, NA), n = 0))
    }
    
    # Calculate empirical rate
    empirical_rate <- mean(p_values < alpha)
    
    # Calculate 95% Confidence Interval (Standard Error of proportion)
    se <- sqrt((empirical_rate * (1 - empirical_rate)) / n)
    ci_lower <- empirical_rate - 1.96 * se
    ci_upper <- empirical_rate + 1.96 * se
    
    # Clamp CI to [0, 1]
    ci_lower <- max(0, ci_lower)
    ci_upper <- min(1, ci_upper)
    
    return(list(rate = empirical_rate, ci = c(ci_lower, ci_upper), n = n))
}

main <- function(opt) {
    if (is.null(opt$dx_dir)) {
        print_help(opt_parser)
        stop("DNAnexus directory must be supplied.", call.=FALSE)
    }

    # List files in the directory
    cat(paste0("Listing files in: ", opt$dx_dir, "\n"))
    cmd <- paste0("dx ls '", opt$dx_dir, "'")
    files <- system(cmd, intern = TRUE)
    
    # Filter by pattern if provided
    if (!is.null(opt$pattern)) {
        files <- grep(opt$pattern, files, value = TRUE)
    }
    
    cat(paste0("Found ", length(files), " files to process.\n"))
    
    if (length(files) == 0) {
        stop("No files found matching the criteria.")
    }
    
    alphas <- c(0.05, 0.01, 0.001, 0.0001, 0.00001)
    results_list <- list()
    
    for (i in seq_along(files)) {
        f <- files[i]
        cat(paste0("[", i, "/", length(files), "] Processing: ", f, "\n"))
        
        # Construct full path
        full_path <- file.path(opt$dx_dir, f)
        
        # Stream file content using dx cat
        # We use tryCatch to handle potential errors during streaming
        tryCatch({
            # Read data
            # cmd argument in fread allows executing a shell command
            dt <- fread(cmd = paste0("dx cat '", full_path, "'"))
            
            # Identify P-value columns
            p_cols <- grep("\\.P$|^p\\.value$", colnames(dt), value = TRUE)
            
            if (length(p_cols) > 0) {
                for (col in p_cols) {
                    p_vals <- dt[[col]]
                    
                    for (alpha in alphas) {
                        res <- track_type1(p_vals, alpha)
                        
                        results_list[[length(results_list) + 1]] <- data.table(
                            file = f,
                            model = col,
                            alpha = alpha,
                            rate = res$rate,
                            ci_lower = res$ci[1],
                            ci_upper = res$ci[2],
                            n_variants = res$n
                        )
                    }
                }
            } else {
                cat("  Warning: No P-value columns found.\n")
            }
            
        }, error = function(e) {
            cat(paste0("  Error processing file: ", e$message, "\n"))
        })
    }
    
    # Combine results
    if (length(results_list) > 0) {
        out_dt <- rbindlist(results_list)
        
        # Write to output
        if (!is.null(opt$output)) {
            fwrite(out_dt, file = opt$output, sep = "\t")
            cat(paste0("\nResults written to: ", opt$output, "\n"))
        } else {
            print(out_dt)
        }
    } else {
        cat("\nNo results generated.\n")
    }
}

# Add arguments
option_list = list(
    make_option(c("-d", "--dx_dir"), type="character", default=NULL,
        help="DNAnexus directory containing result files", metavar="character"),
    make_option(c("-p", "--pattern"), type="character", default=NULL,
        help="Pattern to filter files (regex)", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="type1_error_summary.tsv",
        help="Output file for summary stats", metavar="character")
)

opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
