#!/usr/bin/env Rscript

# to deal with transcripts, genes and other mappings
devtools::install_github("frhl/bravastring")
library(bravastring)
library(data.table)
library(argparse)


null_omit <- function(lst) {
        lst[!vapply(lst, is.null, logical(1))]
}

# this function is taken from SAIGE
beta_weight <-function(MAF,weights.beta=c(1,25)){
    n<-length(MAF)
    weights<-rep(0,n)
    IDX_0<-which(MAF == 0)
    if(length(IDX_0) == n){
      stop("No polymorphic SNPs")
    } else if( length(IDX_0) == 0){
      weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
    } else {
      weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
    }
    return(weights)
}

main <- function(args){
    print(args)
    input_ac <- args$input_ac
    canonical_transcripts_only <- args$canonical_transcripts_only
    out_prefix <- args$out_prefix

    # get AC file only
    dt <- fread(args$input_ac)
    colnames(dt) <- c("ID", "TRANSCRIPT", "AC", "AN", "AF", "MAF")

    # map back to gene id if a transcript has been given
    # and subset to canonical transcripts if need be
    # note: required 'bravastring' installed from github
    is_transcript <- any(grepl(dt$TRANSCRIPT, pattern="ENST"))
    if (is_transcript){
        if (args$canonical_transcripts_only){
            cat("Subsetting to canonical transcripts..\n")
            dt <- dt[is_canonical(dt$TRANSCRIPT),]
        }
        dt$gene_id <- transcript_to_gene_id(dt$TRANSCRIPT)
    } else {
        dt$gene_id <- dt$TRANSCRIPT
    }

    # get weight
    dt$weight <- beta_weight(dt$MAF)
    gene_weights <- dt$weight
    #names(gene_weights) <- dt$TRANSCRIPT
    names(gene_weights) <- dt$ID

    # get genesets we want to run from github
    system("git clone https://github.com/frhl/genesets_public.git", intern=FALSE)
    genesets <- list.files("genesets_public/exported/", full.names=TRUE, pattern=".txt.gz")
    all_genesets_paths <- genesets[!grepl("saige_group", genesets)]
    print(all_genesets_paths)

    ## now create saige mapping
    out_paths <- paste0(args$out_prefix, ".paths")
    paths_lst <- list()
    for (geneset_path in all_genesets_paths){
        
        #print(geneset_path)
        stopifnot(file.exists(geneset_path))
        # read genesets
        geneset_name <- gsub(".txt.gz", "",basename(geneset_path))
        gs <- fread(geneset_path)

        # deal with gTEX data as column names are different
        if ("tissue" %in% colnames(gs)){
             cols_keep <- c("geneset", "values")
             gs <- gs[,..cols_keep]
             colnames(gs) <- c("geneset", "gene_id")
        } 

        map <- gs$geneset
        names(map) <- gs$gene_id
        
        # create deep copy of d
        d <- dt
        d$geneset <- map[d$gene_id]
        print(nrow(d))

        # only test if there are at least two one gene involved
        if ((nrow(d)>0) & (sum(!is.na(d$geneset))>0)) {
            count <- data.table(table(d$geneset))
            count <- count[count$N>0]
            d <- d[d$geneset %in% count$V1, ]
            d$annotation <- "pseudo_variant_geneset"
          
            if (nrow(d)>0){ 

                # create saige specific group file
                genesets <- unique(d$geneset)
                if (length(genesets) > 0){
                    out <- lapply(genesets, function(g){
                        
                        # note, that we use id which corresponds to 
                        # chr:pos:ref:alt, and not gene_id as saige can't use that
                        # newer note: in newer versions of saige it ONLY works
                        # if we use the ID from the VCF/plink file.
                        # even newer: after a call with Wei, we need to use
                        # chr:pos:ref:alt
                        
                        variants <- d$ID[d$geneset %in% g]
                        #variants <- d$ID[d$geneset %in% g]
                        annotations <- d$annotation[d$geneset %in% g]
                        weights <- gene_weights[variants]
                        stopifnot(length(weights)==length(variants))
                        row1 <- paste(c(g, 'var', variants), collapse = " ")
                        row2 <- paste(c(g, 'anno', annotations), collapse = " ")
                        if (!is.null(input_ac)){
                            row3 <- paste(c(g, 'weight', weights), collapse = " ")
                            return(paste0(c(row1, row2, row3), collapse = '\n'))
                        } else {
                            return(paste0(c(row1, row2), collapse = '\n'))
                        } 
                    })
                    
                    # write to geneset
                    out <- null_omit(out)
                    output_path_geneset <- paste0(args$out_prefix, ".", geneset_name, ".txt") 
                    write(paste("Writing to", output_path_geneset), stdout())
                    writeLines(paste(out, collapse = '\n'), output_path_geneset)
                    
                    # write to a file containing all (basename) paths to the genesets
                    paths_lst[[geneset_name]] <- data.table(geneset_name, output_path_geneset)
                    #write(output_path_geneset, file = out_paths, append = TRUE, sep = "\n")
                }
            }                          
        }
    }

    # get rid of genesets so that they are not uploaded
    write(paste0("Wrote overview file to ", out_paths,".."), stdout())
    system("rm -rf genesets_public/", intern=FALSE)
    
    # collapse geneset paths and write the names
    paths <- do.call(rbind, paths_lst)
    fwrite(paths, out_paths, col.names=FALSE, sep="\t")

}

parser <- ArgumentParser()
parser$add_argument("--input_ac", default=NULL, required=TRUE, help="Input AF file")
parser$add_argument("--canonical_transcripts_only", action='store_true', help="Use canonical transcripts only")
parser$add_argument("--out_prefix", default=NULL, required=TRUE, help="Where should the results be written?")
args <- parser$parse_args()

main(args)


