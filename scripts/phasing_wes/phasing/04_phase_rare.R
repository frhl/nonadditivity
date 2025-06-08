#!/usr/bin/env Rscript

library(argparse)
library(data.table)


# list duplicates files
dx.list.files.duplicated <- function(dir, pattern, full.names=FALSE){
  files <- dx.list.files(dir, rm.file.id=TRUE, full.names=full.names)
  files <- files[grepl(files, pattern=pattern)]
  files <- files[which(duplicated(files))]
  return(files)
}


# list files
dx.list.files <- function(dir="", full.names=FALSE, pattern=NULL, rm.file.id=FALSE){
  cmd <- paste("dx ls", dir)
  out <- system(cmd, intern=TRUE)
  if (!is.null(pattern)){
    out <- out[grepl(pattern=pattern, out)]
  }
  if (full.names==TRUE) {
    out <- file.path(dir, out)
  }
  if (rm.file.id){
    out <- gsub("\\ :\\ +file-.*", "", out)
  }
  return(out)
}


main <- function(args){

    target_dir <- args$target_dir
    regex_pattern <- args$regex_pattern
    n_expt <- args$files_expt_count

    # check for duplicates
    files <- dx.list.files.duplicated(target_dir, pattern=regex_pattern, full.names=TRUE)
    if (length(files) > 0){
        stop("Duplicated files were found:\n", paste0(files, collapse="\n"))
    }
    write("* No duplicates.. OK!", stderr())

    # do we have the number of expected files
    files <- dx.list.files(target_dir, pattern=regex_pattern)
    files <- files[!grepl(files, pattern="\\.csi")]
    n <- length(files)
    if (n != n_expt){
        msg <- paste("Found",n,"files of", n_expt,"expected files. Resolve this before continuing")
        stop(msg)
    }
    write("* All files present.. OK!", stderr())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--target_dir", default=NULL, help = "output dir")
parser$add_argument("--files_expt_count", default=NULL, help = "how many files are expected?")
parser$add_argument("--regex_pattern", default=NULL, help = "file in output dir using regex pattern")
args <- parser$parse_args()

main(args)


