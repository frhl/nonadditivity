#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){
  
  interval_path <- args$interval_path
  current_intervals_path <- args$current_intervals_path
  padding <- as.numeric(args$padding)
  out_path <- args$out_path
  
  # load files
  stopifnot(file.exists(interval_path))
  stopifnot(file.exists(current_intervals_path))
  dt <- fread(interval_path)
  idx <- fread(current_intervals_path, header=FALSE)
  
  #dt <- fread("~/Downloads/intervals_min_1000_chr21.tsv")
  #idx <- fread("~/Downloads/phasing_intervals_u1000_s100000_o25000_chr21.tsv", header = FALSE)
  #overlap_desired <- 50
  chrom <- unique(dt$chrom)
  stopifnot(length(chrom) == 1)
  
  # extracting positions
  idx <- data.table(do.call(rbind, strsplit(idx$V1, split="\\:")))
  idx <- data.table(do.call(rbind, strsplit(idx$V2, split="\\-")))
  colnames(idx) <- c("start", "stop")
  stopifnot(all(idx$start %in% dt$pos))
  stopifnot(all(idx$stop %in% dt$pos))
  
  padded_chunks <- list()
  for (i in (1:(nrow(idx)))){
    
    # get some bools
    first_chunk <- i == 1
    last_chunk <- i == nrow(idx)
    
    # get intervals
    chunk_start <- idx$start[i]
    chunk_stop <- idx$stop[i]
    chunk_start_idx <- which(dt$pos == chunk_start)
    chunk_stop_idx <- which(dt$pos == chunk_stop)
    
    # deal with first, mid and last chunk
    if ((first_chunk) && (last_chunk)){
      chunk_new_start <- chunk_start_idx
      chunk_new_stop <- chunk_stop_idx    
    } else if (first_chunk){
      chunk_new_start <- chunk_start_idx
      chunk_new_stop <- chunk_stop_idx+padding
    } else if (last_chunk){
      chunk_new_start <- chunk_start_idx-padding
      chunk_new_stop <- chunk_stop_idx
    } else {
      chunk_new_start <- chunk_start_idx-padding
      chunk_new_stop <- chunk_stop_idx+padding
    }
    
    # get trimmed chunks
    padded_chunks[[i]] <- data.table(
      idx_chunk_start = chunk_new_start,
      idx_chunk_end = chunk_new_stop,
      pos_chunk_start = dt$pos[chunk_new_start],
      pos_chunk_end = dt$pos[chunk_new_stop]
    )
  }
  
  # finalize output
  out <- rbindlist(padded_chunks)
  out$region <- paste0(chrom,":", out$pos_chunk_start, "-", out$pos_chunk_end)
  outfile <- out_path
  write(paste("writing to", outfile), stdout())
  print(out)
  fwrite(out, outfile, sep='\t')
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--interval_path", default=NULL, help = "Current full interval path to slice through")
parser$add_argument("--current_intervals_path", default=NULL, help = "The intervals that are currently being used")
parser$add_argument("--padding", default=NULL, help = "Desired padding in thousands (e.g 5)")
parser$add_argument("--out_path", default=NULL, help = "out prefix (.tsv)")
args <- parser$parse_args()

main(args)

