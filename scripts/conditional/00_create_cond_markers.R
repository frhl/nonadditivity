# code run locally to get regions around (canonical) transcripts 
# that will be uses to subset the imputed data.

library(bravastring) # install via github frhl/bravastring
library(data.table)

# mapping to chromsome
to_chr <- fread("~/Projects/04_genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
to_chr <- to_chr[,c("ensembl_gene_id", "chromosome_name")]
to_chr <- to_chr[!duplicated(to_chr),]
to_chr <- dict(to_chr$ensembl_gene_id, to_chr$chromosome_name)

files <- list(
  "~/Projects/01_2df_git/data/cts_traits_qc_final_pLoF_damaging_missense.txt",
  "~/Projects/01_2df_git/data/olink_qc_combined_pLoF.txt",
  "~/Projects/01_2df_git/data/olink_qc_combined_pLoF_damaging_missense.txt"
)


process_file <- function(infile){
  
  # process file
  d <- fread(infile)
  d$gene_id <- transcript_to_gene_id(d$MarkerID)
  d[, id := paste0(hgnc_symbol, "_", trait)]
  d <- d[!is.na(d$p.value.dom),]
  d <- d[d$AC_Allele2.rec>=10,]
  d$annotation <- str_extract_annotation(infile)
  
  # get FDR <0.1
  #d$FDR.dom <- p.adjust(d$p.value.dom, method="fdr")
  d <- d[d$p.value < 1e-6,]
  
  if (!"title" %in% colnames(d)) d$title <- NA

  # check that we can map to chromosome
  stopifnot(all(d$gene_id %in% names(to_chr)))
  d$chromosome <- to_chr[d$gene_id]
  
  # Create output dataset
  result <- d[, .(
    dataset = basename(infile),
    annotation = annotation,
    ensembl_gene_id = gene_id,
    ensembl_transcript_id = MarkerID,
    hgnc_symbol = hgnc_symbol,
    chromosome = chromosome,
    trait = trait,
    title = title,
    p_value = p.value,
    p_value_dom = p.value.dom
  )]
  
  return(result)
}

# Create simple BED format from your data
create_simple_bed <- function(data, output_file) {
  # Create basic BED format with only required columns
  bed_df <- data.frame(
    chromosome = paste0("chr", data$chromosome),  # Add 'chr' prefix
    start = data$start_position,
    end = data$end_position,
    stringsAsFactors = FALSE
  )
  
  # Sort by chromosome and start position
  bed_df <- bed_df[order(bed_df$chromosome, bed_df$start), ]
  
  # Write to BED file without headers
  write.table(
    bed_df,
    file = output_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  return(bed_df)
}

# get upstream and downstream of each for canonical transcripts
intervals <- fread("~/Projects/04_genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
intervals <- intervals[, .(ensembl_transcript_id, start_position, end_position)]
intervals <- intervals[!duplicated(intervals),]

# we combine with CANONICAL transcripts.
all_results <- do.call(rbind, lapply(files, process_file))
all_results <- merge(all_results, intervals, by="ensembl_transcript_id", all.x=TRUE)
all_results[!duplicated(all_results$hgnc_symbol),]

#all_results <- all_results[all_results$chromosome=="21"]
unique_results <- all_results[!duplicated(all_results$hgnc_symbol), ]
bed_data <- create_simple_bed(unique_results, "~/Desktop/df2_add_cond_intervals_27jan2025.bed")
fwrite(all_results, "~/Desktop/df2_add_cond_intervals_27jan2025.txt", sep="\t", na="NA", quote=FALSE)



