#!/usr/bin/env Rscript

library(data.table)
library(argparse)


main <- function(args){

    print(args)
    input_file <- args$input_file
    n_breaks <- args$n_breaks
    out_prefix <- args$out_prefix
    ac_min <- args$ac_min
    ac_max <- args$ac_max

    d <- fread(input_file, header=FALSE, sep=" ")
    print(head(d))
    colnames(d) <- c("sid", "vid", "gt", "pp", "ac", "af")

    if (!is.null(ac_min)){
        d <- d[d$ac >= as.numeric(ac_min),]
    }
    if (!is.null(ac_max)){
        d <- d[d$ac <= as.numeric(ac_max),]
    }

    print(head(d))

    d <- d[!is.na(d$pp),]
    total_sites <- nrow(d)
    pp_breaks <- seq(1, 0, length.out=as.numeric(n_breaks))
    n_sites <- unlist(lapply(pp_breaks, function(pp_break){ nrow(d[d$pp >= pp_break,]) }))
    combi <- data.table(pp_breaks, n_sites, total_sites, ac_min, ac_max)                        
    
    out <- paste0(out_prefix, ".txt")
    write(paste("writing to", out), stderr())
    fwrite(combi, out, sep="\t")


}

parser <- ArgumentParser()
parser$add_argument("--input_file", default=NULL, required = TRUE, help = "Input")
parser$add_argument("--ac_min", default=NULL, required = FALSE, help = "Input")
parser$add_argument("--ac_max", default=NULL, required = FALSE, help = "Input")
parser$add_argument("--n_breaks", default=NULL, required = TRUE, help = "how many breaks for sequence of limits")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()


main(args)

