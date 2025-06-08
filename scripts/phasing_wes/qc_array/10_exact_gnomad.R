#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    out_prefix <- args$out_prefix
    gnomad_path <- args$gnomad_path
    array_path <- args$array_path

    # get gnomad
    dt_gnomad <- fread(gnomad_path)
    dt_gnomad <- dt_gnomad[dt_gnomad$AN_nfe>0,]
    dt_gnomad[, WGS_AC_A1 := pmin(AC_nfe, AN_nfe - AC_nfe)]
    dt_gnomad[, WGS_AC_A2 := AN_nfe-WGS_AC_A1]
    keep <- c("SNPID", "WGS_AC_A1", "WGS_AC_A2")
    dt_gnomad <- dt_gnomad[,..keep]
    setkey(dt_gnomad, "SNPID")

    # get aray
    dt_array <- fread(array_path)

    # get array
    dt_array[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt_array[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt_array[, MAC := pmin(AC_A1, AC_A2)]
    dt_array[, n_samples := `C(HOM A1)` + `C(HOM A2)` + `C(HET)` + `C(MISSING)`]

    # get samples
    n_samples <- dt_array$n_samples[1]
    dt_array[, MAF := MAC/(AC_A1 + AC_A2)]
    dt_array[, AN := n_samples*2]
    dt_array[, ARRAY_AC_A1 := MAC ]
    dt_array[, ARRAY_AC_A2 := AN-MAC ]

    # only keep a few things
    keep_cols <- c("SNP", "ARRAY_AC_A1", "ARRAY_AC_A2")
    dt_array <- dt_array[,..keep_cols]

    setkey(dt_array, "SNP")

    # merged
    merged <- merge(dt_array, dt_gnomad, all.x=TRUE, by.x="SNP", by.y="SNPID")
    merged[, ARRAY_MAC := pmin(ARRAY_AC_A1, ARRAY_AC_A2)]
    merged[, ARRAY_MAF := ARRAY_MAC / (ARRAY_AC_A1 + ARRAY_AC_A2)]
    merged[, WGS_MAC := pmin(WGS_AC_A1, WGS_AC_A2)]
    merged[, WGS_MAF := WGS_MAC / (WGS_AC_A1 + WGS_AC_A2)]

    # combine with fisher
    fisher <- do.call(rbind, lapply(1:nrow(merged), function(i){
    #fisher <- do.call(rbind, lapply(1:100, function(i){
        
        # compare overlap
        snpid <- merged$SNP[i]
        gnomad_maf <- merged$WGS_MAF[i]
        array_maf <- merged$ARRAY_MAF[i]
        
        cont_table <- matrix(
            c(merged$ARRAY_AC_A1[i], merged$WGS_AC_A1[i],
              merged$ARRAY_AC_A2[i], merged$WGS_AC_A2[i]),
            nrow = 2, byrow = TRUE)
        cont_table[is.na(cont_table)] <- 0
        
        # do fisher's exact test
        pval <- fisher.test(cont_table)$p.value
        out <- cbind(data.table(snpid, gnomad_maf, array_maf, pval))
        return(out)
    }))

    # get cutoff 
    fisher$pval_fail_10 <- fisher$pval < 1e-10
    fisher$pval_fail_20 <- fisher$pval < 1e-20
    fisher$pval_fail_30 <- fisher$pval < 1e-30
    fisher$pval_fail_40 <- fisher$pval < 1e-40
    fisher$pval_fail_50 <- fisher$pval < 1e-50
    fisher$pval_fail_60 <- fisher$pval < 1e-60
    fisher$pval_fail_70 <- fisher$pval < 1e-70
    fisher$pval_fail_80 <- fisher$pval < 1e-80
    fisher$pval_fail_90 <- fisher$pval < 1e-90
    fisher$pval_fail_100 <- fisher$pval < 1e-100
    fisher$pval_fail_150 <- fisher$pval < 1e-150
    fisher$pval_fail_200 <- fisher$pval < 1e-200

    # write outfile
    outfile <- paste0(out_prefix,".txt")
    write(paste("writing to", outfile), stdout())
    fwrite(fisher, outfile, sep='\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "out prefix")
parser$add_argument("--gnomad_path", default=NULL, help = "out prefix")
parser$add_argument("--array_path", default=NULL, help = "out prefix")
args <- parser$parse_args()

main(args)

