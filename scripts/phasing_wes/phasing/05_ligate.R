library(argparse)
library(stringr)

# list files
dx.list.files <- function(dir="", full.names=FALSE, pattern=NULL){
  cmd <- paste("dx ls", dir)
  out <- system(cmd, intern=TRUE)
  if (!is.null(pattern)){
    out <- out[grepl(pattern=pattern, out)]
  }
  if (full.names==TRUE) {
    out <- file.path(dir, out)
  }
  return(out)
}

main <- function(args){

    files <- dx.list.files(dirname(args$in_prefix), pattern = basename(args$in_prefix), full.names = TRUE)
    if (length(files) < 2) stop("expected at least two files!")
    files <- files[grepl(".bcf$", files)]
    files <- files[grepl(".tr[0-9]+_[0-9]+", files)]
    files <- files[!grepl(".csi$", files)]
    files <- files[!grepl("mis10", files)]
    
    write(files,stderr())
    #files <- files[grepl("\\.trimmed\\.", files)]
    stopifnot(length(files) > 0)

    chunks <- stringr::str_extract(files,"[0-9]+of[0-9]+")
    chunks <- data.frame(do.call(rbind, strsplit(chunks, "of")))
    chunks$X1 <- as.numeric(chunks$X1)
    chunks$X2 <- as.numeric(chunks$X2)

    stopifnot(length(unique(chunks$X2)) == 1)
    stopifnot(length(unique(chunks$X1)) == nrow(chunks))

    out <- files[order(chunks$X1)]
    if (!is.null(args$out_prefix)){
        out <- paste0(args$out_prefix, out)
    }
    
    #out <- out[c(1,2)]
    write(paste(out, collapse = " "), stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_prefix", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_prefix", default=NULL, required = FALSE, help = "Add a prefix to the output")
args <- parser$parse_args()

main(args)


