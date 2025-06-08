#!/usr/bin/env Rscript

# to deal with transcripts, genes and other mappings
devtools::install_github("frhl/bravastring")
install.packages("GMMAT")
library(bravastring)
library(data.table)
library(argparse)
library(GMMAT)

get_p_1df_glmmkin <- function(model){
    coef_estimates <- model$coefficients
    coef_se <- sqrt(diag(model$cov))
    wald_stats <- (coef_estimates / coef_se)^2
    p_values <- 1 - pchisq(wald_stats, df = 1)
    return(p_values[2])
}

get_p_2df_glmmkin <- function(model){
    # note assumes that the relevant coef are the 2nd and 3rd!
    coef_x_xr <- model$coefficients[2:3]
    cov_matrix <- model$cov[2:3,2:3]
    wald_stat <- t(coef_x_xr) %*% solve(cov_matrix) %*% coef_x_xr
    p_value_r <- pchisq(wald_stat, 2, lower.tail = FALSE)
    return(c("both"=as.numeric(p_value_r)))
}

inverse_normal_transform <- function(dt, phenoCol) {
  # Ensure the column exists in the data.table
  if (!(phenoCol %in% colnames(dt))) {
    stop(paste("Column", phenoCol, "not found in the data.table"))
  }
  
  # Perform the inverse normal transformation
  dt[, target := {
    ranks <- rank(.SD[[1]], na.last="keep", ties.method="average")
    n <- sum(!is.na(.SD[[1]]))
    qnorm((ranks - 0.5) / n)
  }, .SDcols = phenoCol]
  
  return(dt)
}

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

encode_recessive <- function(add) {
  rec <- ifelse(add == 2, 1, 0)
  rec <- as.numeric(scale(rec))
  return(rec)
}

lrt <- function(nested, complex){
  df <- df.residual(nested) - df.residual(complex)
  stopifnot(df>0)
  nested <- logLik(nested)
  complex <- logLik(complex)
  chisq <- -2 * (as.numeric(nested)-as.numeric(complex))
  p.val <- pchisq(chisq, df = df, lower.tail = FALSE)
  #p.val <- p.val / log(10)
  return(p.val)
}


main <- function(args){

    # should be add, rec or both
    mixed_model_vec <- unlist(strsplit(args$mixed_model,split=","))

    # sparse GRM
    grm <- Matrix::readMM(args$sparse_grm_path)
    path_grm_samples <- paste0(args$sparse_grm_path, ".sampleIDs.txt")
    grm_samples <- readLines(path_grm_samples)
    dimnames(grm) <- list(grm_samples, grm_samples)
	
	# unrelated samples
	#unrelated_samples <- fread(args$samples_path)$s

	# read data
	phenotype_dt <- fread(args$phenotypes_path)

	# subset to specific samples in GRM
	phenotype_dt <- phenotype_dt[phenotype_dt$eid %in% grm_samples,]
	samples <- unique(phenotype_dt$eid)
	print(paste("Using",nrow(phenotype_dt), "samples"))

	# Split the covariates and categorical_covariates strings into vectors
	all_covariates_vector <- strsplit(args$covariates, ",")[[1]]
	categorical_covariates_vector <- strsplit(args$categorical_covariates, ",")[[1]]
	stopifnot(all(all_covariates_vector %in% colnames(phenotype_dt)))
	stopifnot(args$pheno_col %in% colnames(phenotype_dt))

	# Identify numeric covariates (those in all_covariates but not in categorical_covariates)
	numeric_covariates_vector <- setdiff(all_covariates_vector, categorical_covariates_vector)

	# Set each categorical covariate as a factor
	for (col in categorical_covariates_vector) {
	  if (col %in% colnames(phenotype_dt)) {
		phenotype_dt[, (col) := as.factor(get(col))]
	  } else {
		warning(paste("Categorical column", col, "not found in phenotypes data.table"))
	  }
	}

	# Set each numeric covariate as numeric
	for (col in numeric_covariates_vector) {
	  if (col %in% colnames(phenotype_dt)) {
		phenotype_dt[, (col) := as.numeric(get(col))]
	  } else {
		warning(paste("Numeric column", col, "not found in phenotypes data.table"))
	  }
	}


	# subset it to samples in GRM
	#V1	V2	V3	V4	V5	V6
	#<int>	<chr>	<chr>	<chr>	<int>	<chr>
	#1000451	chr1	ENST00000037502	het	1	chr1:171638687:C:T
	#1002310	chr1	ENST00000037502	het	1	chr1:171636338:G:A

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
	(f_both <- paste0("target~add+dom+", paste0(all_covariates_vector,collapse="+")))

	# INT for phenotype
    phenotype_dt <- inverse_normal_transform(phenotype_dt, args$pheno_col)

    lst <- list()
    count <- 0
	# run over all markers
	for (marker in markers){ # <------------------------ LOOOK HERE

        count <- count + 1
        print(paste(marker, count, length(markers)))
		# get dosages for by samples
		marker_dt <- d[d$id == marker]
		marker_dt$id <- NULL

		# get wildtype samples too
		wildtype_samples <- samples[!samples %in% marker_dt$eid]
		wildtype_dt <- data.table(eid=wildtype_samples, ds=0)
		marker_dt <- rbind(marker_dt, wildtype_dt)
		#print(table(marker_dt$ds))

		# check match
		stopifnot(all(marker_dt$eid %in% phenotype_dt$eid))
		stopifnot(all(phenotype_dt$eid %in% marker_dt$eid))

		# set order to phenotype dt
		marker_dt <- marker_dt[match(phenotype_dt$eid, marker_dt$eid),]

		# ecode additive/nonadditive/recessive
		marker_dt$add <- encode_additive(marker_dt$ds)
		marker_dt$dom <- encode_dominance(marker_dt$ds)
		marker_dt$rec <- encode_recessive(marker_dt$ds)

		# combine with phenotype_dt
		tmp_phenotype_dt <- cbind(phenotype_dt, marker_dt[,-(1:2)])
		
		# fit both
        if (mixed_model_vec %in% "add"){ 
            add_model <- glmmkin(fixed = f_add,
                              data = tmp_phenotype_dt,
                              kins = grm,
                              id = "eid",
                              family = gaussian(link = "identity"))
		    p_add <- get_p_1df_glmmkin(add_model)
        } else {
            p_add <- NA
        }

        if (mixed_model_vec %in% "rec"){ 
            rec_model <- glmmkin(fixed = f_rec,
                              data = tmp_phenotype_dt,
                              kins = grm,
                              id = "eid",
                              family = gaussian(link = "identity"))
		    p_rec <- get_p_1df_glmmkin(rec_model)
        } else {
            p_rec <- NA
        }

        if (mixed_model_vec %in% "both"){
            both_model <- glmmkin(fixed = f_both,
                      data = tmp_phenotype_dt,
                      kins = grm,
                      id = "eid",
                      family = gaussian(link = "identity"))	
		    p_both <- get_p_2df_glmmkin(both_model)
        } else {
            p_both <- NA
        }

		# additive/recessive/btoh
		p_dom <- NA
		#p_both <- get_p_2df_glmmkin(both_model)

		out <- data.table(marker, p_add, p_dom, p_rec, p_both)
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
parser$add_argument("--sparse_grm_path", default=NULL, required=TRUE, help="")
parser$add_argument("--phenotypes_path", default=NULL, required=TRUE, help="")
parser$add_argument("--dosages_path", default=NULL, required=TRUE, help="")
parser$add_argument("--covariates", default=NULL, required=TRUE, help="")
parser$add_argument("--categorical_covariates", default=NULL, required=TRUE, help="")
parser$add_argument("--pheno_col", default=NULL, required=TRUE, help="")
parser$add_argument("--out_path", default=NULL, required=TRUE, help="")
parser$add_argument("--mixed_model", default=NULL, required=TRUE, help="")
args <- parser$parse_args()

main(args)


