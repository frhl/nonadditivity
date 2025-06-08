#!/usr/bin/env Rscript

library(data.table)
library(argparse)
library(ggplot2)

main <- function(args){

    in_dir <- args$in_dir
    out_prefix <- args$out_prefix

    # get files
    files <- list.files(in_dir, pattern="csv.gz", full.names=TRUE)
   
    # read the files
    d <- rbindlist(lapply(files, function(f){
        dt <- fread(f)
        dt$chr <- stringr::str_extract(basename(f), "chr(([0-9]+)|(X))")
        return(dt)
    }))
    d$chr <- factor(d$chr, levels=paste0("chr",c(1:22,"X")))

    # get those variants that are SNPs
    is_snp <- data.frame(rbind(table(d$chr, d$`is SNP`)))
    colnames(is_snp) <- c("n_indel", "n_snp")
    is_snp$chr <- rownames(is_snp)

    # get those variants that are SNPs with at least one read
    d_snp <- d[d$`is SNP` == 1,]
    is_snp_with_reads <- data.frame(rbind(table(d_snp$chr, d_snp$`Num PIR`>0)))
    colnames(is_snp_with_reads) <- c("n_snps_wo_reads", "n_snps_w_reads")
    is_snp_with_reads$chr <- rownames(is_snp_with_reads)

    # get variants that are SNPs and have an updated phase because of the read
    d_snp <- d[(d$`is SNP` == 1) & (d$`Num PIR`>0),]
    is_snp_with_updated_phase <- data.frame(rbind(table(d_snp$chr, d_snp$Switched)))
    colnames(is_snp_with_updated_phase) <- c("n_snps_not_updated_gt", "n_snps_updated_gt")
    is_snp_with_updated_phase$chr <- rownames(is_snp_with_updated_phase)

    # get singletons now 
    d_snp <- d[(d$`is SNP` == 1) & (d$`Num PIR`>0) & (d$AC == 1),]
    is_snp_with_updated_phase_singleton <- data.frame(rbind(table(d_snp$chr, d_snp$Switched)))
    colnames(is_snp_with_updated_phase_singleton) <- c("n_snps_singleton_not_updated_gt", "n_snps_singleton_updated_gt")
    is_snp_with_updated_phase_singleton$chr <- rownames(is_snp_with_updated_phase_singleton)

    # combine these results
    lst <- list(is_snp, is_snp_with_reads, is_snp_with_updated_phase, is_snp_with_updated_phase_singleton)
    combined <- Reduce(merge, lst)

    out <- paste0(out_prefix, ".txt")
    write(paste("writing to", out), stderr())
    fwrite(combined, out, sep="\t")

    # make some plots
    colnames(d) <- c("eid", "pp", "switched", "num_pir", "is_snp", "ac", "gt", "vcf_line", "chr")
    d$vcf_line <- NULL

    p <- ggplot(d, aes(x=num_pir)) +
        geom_histogram() +
        scale_y_continuous(trans="log10") +
        ggtitle("Histogram of informative reads") +
        xlab("# phase informative reads") +
        ylab("Counts")

    out <- paste0(out_prefix, ".pir.hist.pdf")
    ggsave(out, p, width=6, height=5) 


}

parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Input")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()


main(args)

