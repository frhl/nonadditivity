#!/usr/bin/env Rscript

library(argparse)
library(data.table)


str_extract_anc <- function(x){
    return(stringr::str_extract(tolower(x), pattern="(eur)|(eas)|(sas)|(nfe)|(afr)"))
}

order_variants <- function(variants) {
  # Separate the chromosome and position
  chr_pos <- strsplit(variants, ":")
  
  # Extract chromosomes and positions
  chromosomes <- sapply(chr_pos, function(x) x[1])
  positions <- sapply(chr_pos, function(x) as.numeric(x[2]))
  
  # Replace 'chrX' with a value that sorts after numeric chromosomes
  chromosomes <- gsub("chr", "", chromosomes)
  chromosomes <- ifelse(chromosomes == "X", "23", chromosomes)
  
  # Convert chromosomes to numeric
  numeric_chromosomes <- as.numeric(chromosomes)
  
  # Use order on numeric chromosomes and positions
  return(order(numeric_chromosomes, positions))
}


main <- function(args){
 
    input_dir <- args$input_dir
    input_pattern <- args$input_pattern 
    columns <- args$columns
    out_prefix <- args$out_prefix

    print(paste("looking for files in", input_dir))
    files <- list.files(input_dir, pattern=input_pattern, full.names=TRUE)
    stopifnot(length(files)>0)

    # get ancestries
    ancestries <- unique(str_extract_anc(basename(files))) 
    stopifnot(length(ancestries)>0)

    # genreally going to read ID,AC,AN,MAF,AC_Het,AC_Hom
    columns <- unlist(strsplit(columns, split=","))

    # read all files and check ancestry column by adding to column name
    l <- lapply(ancestries, function(cur_anc){
        files_anc <- files[grepl(files, pattern=paste0("\\.",cur_anc,"\\."))]
        #files_anc <- c(head(files_anc, 1),tail(files_anc,1))
        # deal with chrX
        files_anc <- gsub("chrX", "chr23", files_anc)
        # order files
        chr_numbers <- as.numeric(gsub(".*_chr(\\d+).*", "\\1", files_anc))
        files_anc <- files_anc[order(chr_numbers)]
        print(files_anc)
        m <- do.call(rbind, lapply(files_anc, function(f){
            d <- fread(gsub("chr23","chrX",f))
            stopifnot(all(columns %in% colnames(d)))
            d <- d[,..columns]
            colnames(d)[-1] <- paste0(colnames(d)[-1],".", toupper(cur_anc))
            return(d)
        }))
        print(head(m))
        return(m)
    })

    # combine all files based on ID as key, and keep
    # non-intersecting items too (full join).
    d <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "ID" ), l)

    # order by chr:postion
    d <- d[order_variants(d$ID),]
    
    # remove variants not in the right format.
    ok <- grepl(d$ID, pattern="^chr[0-9X]+:[0-9]+:[A-Z]+:[A-Z]+$")
    print(sum(!ok))
    d <- d[ok,]

    # write out
    outfile <- paste0(out_prefix, ".txt.gz")
    write(paste("writing to", outfile), stderr())
    fwrite(d, outfile, sep="\t", quote=FALSE, na="NA")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_dir", default=NULL, help = "")
parser$add_argument("--input_pattern", default=NULL, help = "")
parser$add_argument("--columns", default=NULL, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)

