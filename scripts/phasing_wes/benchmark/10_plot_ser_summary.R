#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(argparse)

# format infinite values
format_inf_value <- function(x){
    return(paste0(gsub("-inf","",x),"+"))
}

# make mac bins pretty
pretty_mac_bin <- function(mac_bin, split = "-"){    
    lst <- unlist(lapply(strsplit(mac_bin, split = split), function(x){
        n <- length(x)
        if (n==1){
            return(x)
        } else {
            start <- as.numeric(as.character(x[1]))
            end <- as.numeric(as.character(x[2]))
            if (start >= 1000){
                start <- tolower(scales::label_number_si()(start))
                end <- tolower(scales::label_number_si()(end))
            }
            combined <- paste0(start,'-',end)
            if (grepl("inf",combined)) combined <- format_inf_value(combined)
            return(combined)
        }

    }))
    return(lst)
}



main <- function(args){
   

    switch_all <- args$input_switch_low
    switch_pp90 <- args$input_switch_high

    # no pp filtering
    d_all <- fread(switch_all)
    bins_all <- unique(d_all$mac_bin)
    d_all$mac_bin_pretty <- pretty_mac_bin(d_all$mac_bin)
    d_all$mac_bin_pretty <- factor(d_all$mac_bin_pretty, levels = unique(d_all$mac_bin_pretty))
    d_all$mac_bin <- factor(d_all$mac_bin, levels = unique(d_all$mac_bin))
    d_all$PP <- 0.50

    # pp filtering 90
    d_pp90 <- fread(switch_pp90)
    bins_pp90 <- unique(d_pp90$mac_bin)
    d_pp90$mac_bin_pretty <- pretty_mac_bin(d_pp90$mac_bin)
    d_pp90$mac_bin_pretty <- factor(d_pp90$mac_bin_pretty, levels = unique(d_pp90$mac_bin_pretty))
    d_pp90$mac_bin <- factor(d_pp90$mac_bin, levels = unique(d_pp90$mac_bin))
    d_pp90$PP <- 0.90

    d <- rbind(d_all, d_pp90)
    d$PP <- factor(d$PP)

    my_colors <- c("#65B8EF", "#2266AC")
    names(my_colors) <- c("0.5", "0.9")
    fill_scale2 <- scale_fill_manual(name = "PP", values = my_colors)
    color_scale2 <- scale_color_manual(name = "PP", values = my_colors)

    plt <- ggplot(d,
       aes(
           x=mac_bin_pretty,
           y=100*pointest,
           ymax = 100*upper,
           ymin = 100*lower,
           color = PP,
           group = PP,
       )) +
        theme_bw() +
        geom_line(linetype = "dashed", size = 0.8) +
        geom_pointrange(size=0.95) +
        ylab('% Switch Error rate') + 
        xlab('Mac bin') +
        color_scale2 +
        labs(color = "Variant origin") +
        scale_y_continuous(trans = 'log10', breaks=c(0.2, 0.5, 1, 2, 5, 10, 20, 30, 50)) +    
        annotation_logticks(sides='l')  +
        theme(
            legend.position = "right",
            axis.text=element_text(size=15),
            axis.title=element_text(size=15,face="bold"),
            axis.title.x = element_text(margin=ggplot2::margin(t=10)),
            axis.title.y = element_text(margin=ggplot2::margin(r=10)),
            plot.title = element_text(hjust=0.5),
            axis.text.x = element_text(angle = 45, vjust = 0.96, hjust=0.95)
        ) +
        ggtitle("SHAPEIT5 400k Europeans", "~100 trios for evaluation")


    out_file <- paste0(args$out_prefix,".pdf")
    write(paste("writing",out_file), stderr())
    ggsave(out_file, plt, width=10, height=7)

}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_switch_low", default=NULL, required = TRUE, help = "")
parser$add_argument("--input_switch_high", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


