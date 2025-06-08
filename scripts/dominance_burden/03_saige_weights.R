#!/usr/bin/env Rscript

library(argparse)
library(data.table)

null_omit <- function(lst) {
    lst[!vapply(lst, is.null, logical(1))]
}

# this function is taken from SAIGE
beta_weight <-function(MAF,weights.beta=c(1,25)){
    n<-length(MAF)
    weights<-rep(0,n)   
    IDX_0<-which(MAF == 0)
    if(length(IDX_0) == n){
      stop("No polymorphic SNPs")
    } else if( length(IDX_0) == 0){
      weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
    } else {
      weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
    }
    return(weights)
}

convert_variant_format <- function(variants) {
  # Use gsub to perform the replacement
  converted_variants <- gsub("(chr[0-9]+):([0-9]+):([A-Z]+):([A-Z]+)", "\\1:\\2_\\3/\\4", variants)
  return(converted_variants)
}

# apply the inverse of the scaling
inverse_scaling <- function(y, a, b) {
    stopifnot(b>a)
    return ((y * (b - a) + 2 * a) / 2)
}

main <- function(args){

    print(args)
    ac_path <- args$weights_by_ac_path
    min_mac <- as.numeric(args$min_mac)

    d <- fread(args$input_path)
    print(head(d))
    M <- d[,c("GENE", "SNP_ID", "ANNOTATION")]
    colnames(M) <- c('gene','variant','anno')
    M$weight <- NA

    # subset to markers in domiannce file 
    markers <- fread(args$markers_path, header=FALSE)$V1
    print(nrow(M))
    M <- M[M$variant %in% markers,]
    stopifnot(nrow(M)>0)
    print(nrow(M))

    # get variants in input path
    input_variants<- unique(M$variant)

    # create a mapping of weights
    if (!is.null(ac_path)){
        write("Loading AC to generate weights..", stdout())
        dt_AC <- fread(ac_path)
        stopifnot(nrow(dt_AC)>1)
        dt_AC[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
        dt_AC[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
        dt_AC[, MAC := pmin(AC_A1, AC_A2)]
        dt_AC[, MAF := MAC/(AC_A1 + AC_A2)]
        dt_AC$weight <- beta_weight(dt_AC$MAF)

        # remove variants that have zero alternate alleles
        dt_AC <- dt_AC[dt_AC$MAC>0,]
        dt_AC <- dt_AC[dt_AC$MAC>=min_mac,]

        # inverse of dominance scaling applied 
        a <- as.numeric(args$dom_scaling_min)
        b <- as.numeric(args$dom_scaling_max)

        # need to figure out the issue with SAIGE - currently fails with negative weights
        #dt_AC$weight <- inverse_scaling(dt_AC$weight, a, b)

        # We remove variants from input that are not present in AC file
        # as AC file is ancestry specific and contains only QCed variants
        M <- M[M$variant %in% unique(dt_AC$SNP), ]

        # create weight mapping
        snp_weights <- dt_AC$weight
        names(snp_weights) <- dt_AC$SNP
    }

    genes <- unique(M$gene)
    out <- lapply(genes, function(g){
      variants <- M$variant[M$gene %in% g]
      annotations <- M$anno[M$gene %in% g]
      nas <- (is.na(variants) | is.na(annotations))
      variants <- variants[!nas]
      annotations <- annotations[!nas]
      accepted <- annotations %in% c('pLoF','damaging_missense', "other_missense","synonymous", "non_coding")
      if (sum(accepted) > 0){
          variants <- variants[accepted]
          annotations <- annotations[accepted]
          row1 <- paste(c(g,'var', variants), collapse = " ")
          row2 <- paste(c(g, 'anno', annotations), collapse = " ")
          if (!is.null(ac_path)){ 
            row3 <- paste(c(g, 'weight', abs(snp_weights[variants])), collapse = " ")
            return(paste0(c(row1, row2, row3), collapse = '\n'))
          } else {
            return(paste0(c(row1, row2), collapse = '\n'))
          } 
      }
    })
    out <- null_omit(out)
    write(paste("Writing to", args$output_path), stdout())
    writeLines(paste(out, collapse = '\n'), args$output_path)
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", required=TRUE, default=NULL, help = "?")
parser$add_argument("--markers_path", required=TRUE, default=NULL, help = "?")
parser$add_argument("--output_path", required=TRUE, default=NULL, help = "?")
parser$add_argument("--delimiter", default=NULL, help = "?")
parser$add_argument("--dom_scaling_max", default=NULL, help = "?")
parser$add_argument("--dom_scaling_min", default=NULL, help = "?")
parser$add_argument("--weights_by_ac_path", default=NULL, help = "?")
parser$add_argument("--min_mac", default=0, help = "?")
args <- parser$parse_args()

main(args)

