#!/usr/bin/env Rscript

# to deal with transcripts, genes and other mappings
devtools::install_github("frhl/bravastring")
library(bravastring)
library(data.table)
library(argparse)

main <- function(args){

	 # List and filter files
    lst <- list.files(args$in_dir, full.names=TRUE, pattern=args$pattern)
    print(paste(length(lst), "items found!"))
    lst <- lst[grepl(lst, pattern=".txt.gz")]
    print(paste(length(lst), "items found!"))

    # combine it all	
	mylist <- lapply(lst, function(f){
	  d <- fread(f)
      if ((nrow(d)>0) & (ncol(d)>5)){
          d$annotation <- str_extract_annotation(basename(f))
          d$ancestry <- str_extract_ancestry(basename(f))
          return(d)
      } else {
          warning(paste(f,"had no rows. Skipping.."))
          return(NULL)
      }
	})

	# write outfile
    mylist <- mylist[lengths(mylist) != 0]
    final <- do.call(rbind, mylist)
    fwrite(final, args$out_path, col.names=TRUE, sep="\t")

}

parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required=TRUE, help="")
parser$add_argument("--pattern", default=NULL, required=TRUE, help="")
parser$add_argument("--out_path", default=NULL, required=TRUE, help="")
args <- parser$parse_args()

main(args)


