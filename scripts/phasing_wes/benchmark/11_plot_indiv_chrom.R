#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(argparse)
library(ggrastr)

lun <- function(x) length(unique(x))

fread_cumsum <- function(f, pp=NULL){
    d <- fread(f)
    d <- setDT(aggregate(switches~POS,d,FUN=sum))
    d$cumsum <- cumsum(d$switches) #
    d$cumsum_div_lun <- d$cumsum/lun(d$POS)
    setkeyv(d, c("POS"))
    d$switches <- NULL
    d$contig <- stringr::str_extract(f, pattern = "chr[0-9]+")
    d$PP <- pp
    return(d)
}

main <- function(args){
   
    files_all <- list.files(args$dir_switch_low, pattern=args$pattern_switch_low, full.names=TRUE)
    files_pp90 <- list.files(args$dir_switch_high, pattern=args$pattern_switch_high, full.names=TRUE)

    lst_all <- rbindlist(lapply(files_all, function(f) fread_cumsum(f, "0.5")))
    lst_pp90 <- rbindlist(lapply(files_pp90, function(f) fread_cumsum(f, "0.9")))
    mrg <- rbind(lst_all, lst_pp90)
    mrg$contig <- factor(mrg$contig, levels=paste0("chr",1:22))

    my_colors <- c("#65B8EF", "#2266AC")
    names(my_colors) <- c("0.5", "0.9")
    fill_scale2 <- scale_fill_manual(name = "PP", values = my_colors)
    color_scale2 <- scale_color_manual(name = "PP", values = my_colors)

    p1 <- ggplot(mrg, aes(x=POS,y=cumsum,color=PP)) +
        geom_point_rast(alpha=0.9) +
        ylab("Cumulative Switch Errors") +
        xlab("GRCh38 position") +
        color_scale2 +
        facet_wrap(~contig, scales="free")
    
     p2 <- ggplot(mrg, aes(x=POS,y=cumsum_div_lun,color=PP)) +
        geom_point_rast(alpha=0.7) +
        ylab("Cumulative Switch Errors Rate") +
        xlab("GRCh38 position") +
        color_scale2 +
        facet_wrap(~contig, scales="free")

    out_file1 <- paste0(args$out_prefix,".cumulative_SE.pdf")
    ggsave(out_file1, p1, width=18, height=16)

    out_file2 <- paste0(args$out_prefix,".cumulative_SER.pdf")
    ggsave(out_file2, p2, width=18, height=16)

}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--dir_switch_low", default=NULL, required = TRUE, help = "")
parser$add_argument("--pattern_switch_low", default=NULL, required = TRUE, help = "")
parser$add_argument("--dir_switch_high", default=NULL, required = TRUE, help = "")
parser$add_argument("--pattern_switch_high", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


