#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

  interval_path <- args$interval_path
  current_intervals_path <- args$current_intervals_path
  overlap_desired <- as.numeric(args$overlap_desired)
  out_path <- args$out_path

  # load files
  stopifnot(file.exists(interval_path))
  stopifnot(file.exists(current_intervals_path))
  dt <- fread(interval_path)
  idx <- fread(current_intervals_path, header=FALSE)

  #dt <- fread("~/Downloads/intervals_min_1000_chr21.tsv")
  #idx <- fread("~/Downloads/intervals_u1000_s100000_o25000_chr21.tsv", header = FALSE)
  #overlap_desired <- 5
  print(unique(dt$chrom))
  chrom <- unique(dt$chrom)
  stopifnot(length(chrom) == 1)
  
  # extracting positions
  idx <- data.table(do.call(rbind, strsplit(idx$V1, split="\\:")))
  idx <- data.table(do.call(rbind, strsplit(idx$V2, split="\\-")))
  colnames(idx) <- c("start", "stop")
  stopifnot(all(idx$start %in% dt$pos))
  stopifnot(all(idx$stop %in% dt$pos))
 
  trimmed_chunks <- list()
  for (i in (2:max(2,nrow(idx)))){
    
    # get some bools
    first_chunk <- i == 2
    last_chunk <- i == (nrow(idx)+1)
    
    # intervals A and B
    chunk_a_start <- idx$start[i-1]
    chunk_a_stop <- idx$stop[i-1]
    chunk_b_start <- idx$start[i]
    chunk_b_stop <- idx$stop[i]
    
    # ensure that they are overlapping
    #stopifnot(chunk_a_stop > chunk_b_start)
    
    # get overlap index
    chunk_a_start_idx <- which(dt$pos == chunk_a_start)
    chunk_a_stop_idx <- which(dt$pos == chunk_a_stop)
    chunk_b_start_idx <- which(dt$pos == chunk_b_start)
    chunk_b_stop_idx <- which(dt$pos == chunk_b_stop)
    
    # get overlaps
    overlaps <- chunk_a_stop_idx-chunk_b_start_idx
    stopifnot(overlap_desired < overlaps)
    stopifnot(overlaps %% overlap_desired == 0)
    chunk_idx_shift <- (overlaps - overlap_desired ) / 2
    stopifnot(chunk_idx_shift == round(chunk_idx_shift))
    
    # deal with first, mid and last chunk 
    if (first_chunk & last_chunk){
      chunk_cur_new_start <- chunk_a_start_idx
      chunk_cur_new_stop <- chunk_a_stop_idx
    } else if (first_chunk){
      chunk_cur_new_start <- chunk_a_start_idx
      chunk_cur_new_stop <- chunk_a_stop_idx - chunk_idx_shift
    } else {
      chunk_cur_new_start <- chunk_a_start_idx + chunk_idx_shift 
      chunk_cur_new_stop <- chunk_a_stop_idx - chunk_idx_shift 
    }
    
    # get trimmed chunks
    trimmed_chunks[[i]] <- data.table(
      idx_chunk_start = chunk_cur_new_start, 
      idx_chunk_end = chunk_cur_new_stop,
      pos_chunk_start = dt$pos[chunk_cur_new_start],
      pos_chunk_end = dt$pos[chunk_cur_new_stop]
    )
  }
  
  i <- i + 1
  # add last chunk
  chunk_last_start <- idx$start[i-1]
  chunk_last_stop <- idx$stop[i-1]
  # get indexes
  chunk_last_start_idx <- which(dt$pos == chunk_last_start)
  chunk_last_stop_idx <- which(dt$pos == chunk_last_stop)
  # get new indexes
  chunk_cur_new_start <- chunk_last_start_idx + chunk_idx_shift 
  chunk_cur_new_stop <- chunk_last_stop_idx
  
  trimmed_chunks[[i]] <- data.table(
    idx_chunk_start = chunk_cur_new_start, 
    idx_chunk_end = chunk_cur_new_stop,
    pos_chunk_start = dt$pos[chunk_cur_new_start],
    pos_chunk_end = dt$pos[chunk_cur_new_stop]
  )

  # finalize output
  out <- rbindlist(trimmed_chunks)
  out$region <- paste0(chrom,":", out$pos_chunk_start, "-", out$pos_chunk_end)
  outfile <- out_path
  write(paste("writing to", outfile), stdout())
  fwrite(out, outfile, sep='\t')
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--interval_path", default=NULL, help = "Current full interval path to slice through")
parser$add_argument("--current_intervals_path", default=NULL, help = "The intervals that are currently being used")
parser$add_argument("--overlap_desired", default=NULL, help = "Desired overlap in thousands (e.g 5)")
parser$add_argument("--out_path", default=NULL, help = "out prefix (.tsv)")
args <- parser$parse_args()

main(args)


