#!/usr/bin/env Rscript

# to deal with transcripts, genes and other mappings
devtools::install_github("frhl/bravastring")
library(bravastring)
library(data.table)
library(argparse)

encode_dominance <- function(add){
  a <- sum(add == 2) / length(add)
  h <- sum(add == 1) / length(add)
  r <- sum(add == 0) / length(add)
  const <- 1 / sqrt(h * a * r * (a + r - (a - r)^2))

  dom_hom_ref <- -h * a
  dom_het <- 2 * a * r
  dom_hom_alt <- -h * r
  
  dom <- add
  dom[which(add == 0)] <- const*dom_hom_ref 
  dom[which(add == 1)] <- const*dom_het
  dom[which(add == 2)] <- const*dom_hom_alt
  return(dom)
}

encode_additive <- function(add){
  a <- sum(add == 2) / length(add)
  h <- sum(add == 1) / length(add)
  r <- sum(add == 0) / length(add)
  const <- 1 / sqrt(a + r - (a - r)^2)
  
  dom_hom_ref <- -(h+2*a)
  dom_het <- 1-(h+2*a)
  dom_hom_alt <- 2-(h+2*a)
  
  dom <- add
  dom[which(add == 0)] <- const*dom_hom_ref 
  dom[which(add == 1)] <- const*dom_het
  dom[which(add == 2)] <- const*dom_hom_alt
  return(dom)
}

# rescale dosage to 0-2
rescale_dom_dosage <- function(dosage){
    return(2*((dosage-min(dosage))/(max(dosage)-min(dosage))))
}

encode_recessive <- function(add) {
  rec <- ifelse(add == 2, 1, 0)
  rec <- as.numeric(scale(rec))
  return(rec)
}

lrt <- function(nested, complex){
  df <- df.residual(nested) - df.residual(complex)
  if (df>0){
      nested <- logLik(nested)
      complex <- logLik(complex)
      chisq <- -2 * (as.numeric(nested)-as.numeric(complex))
      p.val <- pchisq(chisq, df = df, lower.tail = FALSE)
      return(p.val)}
  else {
      warning(paste("DF=",df, "but expected df>0. Returning 'NA'"))
      return(NA)
  }
}

# this is needed for very small numbers!
combine_pvalues_log <- function(p_add, p_dom) {
  
  # Convert p-values to chi-squared values (in log space)
  log_chi_add <- qchisq(log(p_add), df=1, lower.tail=FALSE, log.p = TRUE)
  log_chi_dom <- qchisq(log(p_dom), df=1, lower.tail=FALSE, log.p = TRUE)
  
  # Sum the chi-squared values (using log-sum-exp trick for numerical stability)
  max_log_chi <- pmax(log_chi_add, log_chi_dom)
  min_log_chi <- pmin(log_chi_add, log_chi_dom)
  log_chi_both <- max_log_chi + log1p(exp(min_log_chi - max_log_chi))
  
  # and combine
  log_p_both <- pchisq(exp(log_chi_both), df=2, lower.tail=FALSE, log.p=TRUE)
  
  # Return the log p-value directly
  return(log_p_both)
}

# avoid subscript errors when models fail
get_coef_safely <- function(model_summary, term) {
	# Get beta
	beta <- tryCatch({
		coef(model_summary)[term, "Estimate"]
	}, error = function(e) {
		NA_real_
	})
	
	# Get SE separately
	se <- tryCatch({
		coef(model_summary)[term, "Std. Error"]
	}, error = function(e) {
		NA_real_
	})
	
	return(list(beta = beta, se = se))
}

extract_coef_from_lm <- function(fit){
    coefs <- coef(summary(fit))
    d <- data.table(coefs)
    colnames(d) <- c("est", "stderr", "t", "p")
    rownames(d) <- rownames(coefs)
    return(d)
}

main <- function(args){

    # unrelated samples
    unrelated_samples <- fread(args$samples_path)$s

    # Generate simulated phenotype with simulation seed
    set.seed(args$sim_seed)
    
    # Create a data table with EIDs and simulated phenotype
    phenotype_dt <- data.table(eid = unrelated_samples)
    phenotype_dt$sim_pheno <- rnorm(nrow(phenotype_dt), mean = 0, sd = 1)
    
    print(paste("Simulated phenotype with seed", args$sim_seed))
    print(paste("Using", nrow(phenotype_dt), "samples"))
    
    # Use the simulated phenotype as the phenotype column
    args$pheno_col <- "sim_pheno"
    
    # Split the covariates and categorical_covariates strings into vectors
    all_covariates_vector <- strsplit(args$covariates, ",")[[1]]
    categorical_covariates_vector <- strsplit(args$categorical_covariates, ",")[[1]]
    
    # Add dummy covariates if needed for the simulation
    if (length(all_covariates_vector) > 0) {
        # Add dummy covariates for simulation purposes
        for (cov in all_covariates_vector) {
            if (cov == "sex") {
                phenotype_dt$sex <- sample(c(0, 1), nrow(phenotype_dt), replace = TRUE)
            } else {
                phenotype_dt[[cov]] <- rnorm(nrow(phenotype_dt), mean = 0, sd = 1)
            }
        }
    } else {
        # If no covariates specified, create at least one dummy covariate
        phenotype_dt$dummy_covar <- rnorm(nrow(phenotype_dt), mean = 0, sd = 1)
        args$covariates <- "dummy_covar"
        all_covariates_vector <- c("dummy_covar")
        categorical_covariates_vector <- c()
    }
    
    # Set each categorical covariate as a factor
    for (col in categorical_covariates_vector) {
      if (col %in% colnames(phenotype_dt)) {
        phenotype_dt[, (col) := as.factor(get(col))]
      } else {
        warning(paste("Categorical column", col, "not found in phenotypes data.table"))
      }
    }

    # Set each numeric covariate as numeric
    numeric_covariates_vector <- setdiff(all_covariates_vector, categorical_covariates_vector)
    for (col in numeric_covariates_vector) {
      if (col %in% colnames(phenotype_dt)) {
        phenotype_dt[, (col) := as.numeric(get(col))]
      } else {
        warning(paste("Numeric column", col, "not found in phenotypes data.table"))
      }
    }

    # Get samples 
    samples <- unique(phenotype_dt$eid)

    # Read dosage data
    d <- fread(args$dosages_path)
    d <- d[,c(1,3,5)]
    colnames(d) <- c("eid", "id", "ds")
    print(paste(uniqueN(d$id), "genes were found"))
    d <- d[d$eid %in% samples,]

    # subset to markers with at least five homs/CHs
    counts <- data.table(table(d$id[d$ds==2]))
    counts <- counts[counts$N>=5,]
    markers <- counts$V1
    d <- d[d$id %in% counts$V1,]

    # models we run
    (f_null <- paste0("target~", paste0(all_covariates_vector,collapse="+")))
    (f_add <- paste0("target~add+", paste0(all_covariates_vector,collapse="+")))
    (f_rec <- paste0("target~rec+", paste0(all_covariates_vector,collapse="+")))
    (f_dom <- paste0("target~dom+", paste0(all_covariates_vector,collapse="+")))
    (f_dom_rescaled <- paste0("target~dom_rescaled+", paste0(all_covariates_vector,collapse="+")))
    (f_both <- paste0("target~add+dom+", paste0(all_covariates_vector,collapse="+"))) 

    # Create target variable (no need for inverse normal transform as we already simulated from normal)
    phenotype_dt$target <- phenotype_dt[[args$pheno_col]]

    count <- 1
    lst <- list()
    # run over all markers
    for (marker in markers){

        count <- count + 1
        print(paste(marker, count, length(markers)))
        # get dosages for by samples
        marker_dt <- d[d$id == marker]
        marker_dt$id <- NULL

        # get wildtype samples too
        wildtype_samples <- samples[!samples %in% marker_dt$eid]
        wildtype_dt <- data.table(eid=wildtype_samples, ds=0)
        marker_dt <- rbind(marker_dt, wildtype_dt)

        # check match
        stopifnot(all(marker_dt$eid %in% phenotype_dt$eid))
        stopifnot(all(phenotype_dt$eid %in% marker_dt$eid))

        # set order to phenotype dt
        marker_dt <- marker_dt[match(phenotype_dt$eid, marker_dt$eid),]

        # encode additive/nonadditive/recessive
        marker_dt$add <- marker_dt$ds
        marker_dt$dom <- encode_dominance(marker_dt$ds)
        marker_dt$dom_rescaled <- rescale_dom_dosage(marker_dt$dom)
        marker_dt$rec <- encode_recessive(marker_dt$ds)
        
        # save min max dom dom dosage
        min_dosage <- min(marker_dt$dom)
        max_dosage <- max(marker_dt$dom)

        # get DS
        ds <- marker_dt$ds
        a <- sum(ds == 2) / length(ds)
        h <- sum(ds == 1) / length(ds)
        r <- sum(ds == 0) / length(ds)

        # get some indices for extracting dom/add encoding
        idx_ds0 <- which(ds==0)[1]
        idx_ds1 <- which(ds==1)[2]
        idx_ds2 <- which(ds==2)[3]

        add_0 <- marker_dt$add[[idx_ds0]]
        add_1 <- marker_dt$add[[idx_ds1]]
        add_2 <- marker_dt$add[[idx_ds2]]

        dom_0 <- marker_dt$dom[[idx_ds0]]
        dom_1 <- marker_dt$dom[[idx_ds1]]
        dom_2 <- marker_dt$dom[[idx_ds2]]

        # combine with phenotype_dt
        tmp_phenotype_dt <- cbind(phenotype_dt, marker_dt[,-(1:2)])

        # run models
        add_model <- lm(f_add, data=tmp_phenotype_dt)
        rec_model <- lm(f_rec, data=tmp_phenotype_dt)
        dom_model <- lm(f_dom, data=tmp_phenotype_dt)
        dom_rescaled_model <- lm(f_dom_rescaled, data=tmp_phenotype_dt)

        # get model summaries for extracting betas and SEs
        add_summary <- summary(add_model)
        dom_summary <- summary(dom_model)
        dom_rescaled_summary <- summary(dom_rescaled_model)
        rec_summary <- summary(rec_model)

        # Extract coefficients safely
        add_coefs <- get_coef_safely(add_summary, "add")
        dom_coefs <- get_coef_safely(dom_summary, "dom")
        dom_rescaled_coefs <- get_coef_safely(dom_rescaled_summary, "dom_rescaled")
        rec_coefs <- get_coef_safely(rec_summary, "rec")

        use_lrt=FALSE
        if (!use_lrt){
            p_add <- extract_coef_from_lm(add_model)$p[2]
            p_rec <- extract_coef_from_lm(rec_model)$p[2]
            p_dom <- extract_coef_from_lm(dom_model)$p[2]
            # should be using combine_pvalues_log instead!
            chi_both <- qchisq(1-p_add, df=1) + qchisq(1-p_dom, df=1)
            p_both <- 1-pchisq(chi_both, df=2)
            p_add_vs_both <- NA
        } else {
            # need to fit null models if we use LRT
            null_model <- lm(f_null, data=tmp_phenotype_dt)
            both_model <- lm(f_both, data=tmp_phenotype_dt)
            # get Ps
            p_add <- lrt(null_model, add_model)
            p_dom <- lrt(null_model, dom_model)
            p_rec <- lrt(null_model, rec_model)
            p_both <- lrt(null_model, both_model)
            p_add_vs_both <- lrt(add_model, both_model)
        }

        # return all of them
        out <- data.table(
            marker = marker,
             ac = sum(marker_dt$ds),
            hets = sum(marker_dt$ds==1),
            homs = sum(marker_dt$ds==2),
            p_add = p_add,
            p_dom = p_dom,
            p_rec = p_rec,
            p_both = p_both,
            p_add_vs_both = p_add_vs_both,
            beta_add = add_coefs$beta,
            se_add = add_coefs$se,
            beta_dom = dom_coefs$beta,
            se_dom = dom_coefs$se,
            beta_rescaled_dom = dom_rescaled_coefs$beta,
            se_rescaled_dom = dom_rescaled_coefs$se,
            a = a,
            h = h,
            r = r,
            add_ds_0 = add_0,
            add_ds_1 = add_1,
            add_ds_2 = add_2,
            dom_ds_0 = dom_0,
            dom_ds_1 = dom_1,
            dom_ds_2 = dom_2,
            min_dosage = min_dosage,
            max_dosage = max_dosage,
            sim_seed = args$sim_seed
        )

        lst[[marker]] <- out
    }

    # combine
    final <- do.call(rbind,lst)
    final$gene_id <-  transcript_to_gene_id(final$marker)
    final$hgnc_symbol <- gene_id_to_hgnc_symbol(transcript_to_gene_id(final$marker))

    # write outfile
    fwrite(final, args$out_path, col.names=TRUE, sep="\t")
}

parser <- ArgumentParser()
parser$add_argument("--samples_path", default=NULL, required=TRUE, help="Path to samples file")
parser$add_argument("--dosages_path", default=NULL, required=TRUE, help="Path to dosages file")
parser$add_argument("--covariates", default="", required=FALSE, help="Comma-separated list of covariates")
parser$add_argument("--categorical_covariates", default="", required=FALSE, help="Comma-separated list of categorical covariates")
parser$add_argument("--out_path", default=NULL, required=TRUE, help="Output file path")
parser$add_argument("--sim_seed", default=1, type="integer", help="Seed for simulating phenotypes")
args <- parser$parse_args()

main(args)