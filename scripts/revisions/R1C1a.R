# Install packages if you don't have them
# install.packages("readxl")
# install.packages("data.table")

library(readxl)
library(data.table)
library(latex2exp)
library(ggrepel)


# resrtict to those with MAF < 5% 
# NOTE: that any assoc with P>=0.01 have been removed.
d_full <- fread(cmd = "zcat < ~/Desktop/dhindsa2023/azphewas-com-470k-exwas-proteomics.csv.xz")
d <- d_full # deep copy
d <- d[d$MAF<0.05,]

# we only focus on HIGH impact variants,
# see https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
csqs_ontology_HIGH <- c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
                        "stop_gained", "frameshift_variant", "stop_lost", "start_lost")
d <- d[`Consequence type` %in% csqs_ontology_HIGH]

# clean it up - define column mapping
col_mapping <- c(
  GENE    = "Gene",
  PROTEIN = "Phenotype",
  VARIANT = "Variant",
  MAF     = "MAF",
  CSQS    = "Consequence type",
  AA      = "No. AA genotypes",
  AB      = "No. AB genotypes",
  BB      = "No. BB genotypes",
  MODEL   = "Model",
  P       = "p-value",
  BETA    = "Effect size",
  SE      = "Effect size standard error"
)

# select and rename columns
d <- d[, ..col_mapping]
setnames(d, names(col_mapping))

# post-processing
d[, GENE := gsub("\\'", "", GENE)]
d[, PROTEIN := tstrsplit(PROTEIN, "#", fixed = TRUE)[[3]]]

# restrict to relevant models and at least two biallelics
d <- d[MODEL != "dominant",]

d_wide <- dcast(d, GENE + PROTEIN + VARIANT + MAF + CSQS + AA + AB + BB ~ MODEL,
                value.var = c("P", "BETA", "SE"))

# Rename ALL columns unique to d_wide with EXWAS prefix
setnames(d_wide, 
         old = c("VARIANT", "MAF", "CSQS", "AA", "AB", "BB",
                 "P_genotypic", "P_recessive", 
                 "BETA_genotypic", "BETA_recessive",
                 "SE_genotypic", "SE_recessive"),
         new = c("EXWAS_VARIANT", "EXWAS_MAF", "EXWAS_CSQS", 
                 "EXWAS_AA", "EXWAS_AB", "EXWAS_BB",
                 "EXWAS_P_ADD", "EXWAS_P_REC",
                 "EXWAS_BETA_ADD", "EXWAS_BETA_REC",
                 "EXWAS_SE_ADD", "EXWAS_SE_REC"))


f <- fread("~/Downloads/olink_plof_synonymous_table_with_betas.txt")

# define column mapping for pLoF columns only
col_mapping <- c(
  TRANSCRIPT    = "MarkerID",
  GENE_ID       = "gene_id",
  GENE          = "hgnc_symbol",
  PROTEIN       = "trait",
  AC            = "AC_Allele2.pLoF",
  BC            = "AC_Allele2.rec.pLoF",
  P_ADD         = "p.value.pLoF",
  BETA_ADD      = "BETA.pLoF",
  SE_ADD        = "SE.pLoF",
  P_REC         = "p.value.rec.pLoF",
  BETA_REC      = "BETA.rec.pLoF",
  SE_REC        = "SE.rec.pLoF",
  P_DOM         = "p.value.dom.pLoF",
  BETA_DOM      = "BETA.dom.pLoF",
  SE_DOM        = "SE.dom.pLoF",
  P_COND_ADD    = "p.value.cond.pLoF",
  P_COND_DOM    = "p.value.dom.cond.pLoF"
)

# select and rename columns
f <- f[, ..col_mapping]
setnames(f, names(col_mapping))

# Scale recessive beta and SE from [0,0,2] to [0,0,1] coding
# why? EXWAS used [0,0,1] and CH pLoF used [0,0,2], so we need to multiply CH pLoF with 2
# to ensure that they are on the same scale
f[, BETA_REC := BETA_REC * 2]
f[, SE_REC := SE_REC * 2]

# Now merge with f
combined <- merge(f, d_wide, 
                  by = c("GENE", "PROTEIN"), 
                  all.x = TRUE)

# restrict to overlapping only
combined <- combined[!is.na(EXWAS_VARIANT),]

# restrict to those where non-additive effects are present
combined <- combined[P_DOM < 0.05,]

# restrict to those where we have at least 3 homozygotes
combined <- combined[EXWAS_BB >= 4, ]

combined[combined$PROTEIN=="LGALS3"]

# Recessive model comparison
ggplot(combined, aes(x = EXWAS_BETA_REC, y = BETA_REC)) +
  geom_errorbar(aes(ymin = BETA_REC - 1.96*SE_REC, 
                    ymax = BETA_REC + 1.96*SE_REC), 
                width = 0) +
  geom_errorbarh(aes(xmin = EXWAS_BETA_REC - 1.96*EXWAS_SE_REC, 
                     xmax = EXWAS_BETA_REC + 1.96*EXWAS_SE_REC), 
                 height = 0) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "xRecessive Models [0,0,2]",
    x = "Dhindsa 2023, Beta (Recessive)",
    y = "Lassen 2025, Beta (Recessive)",
    subtitle = sprintf("n = %d, r = %.3f", 
                       sum(!is.na(combined$BETA_REC) & !is.na(combined$EXWAS_BETA_REC)),
                       cor(combined$EXWAS_BETA_REC, combined$BETA_REC, use = "complete.obs"))
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# Additive model comparison
ggplot(combined, aes(x = EXWAS_BETA_ADD, y = BETA_ADD)) +
  geom_errorbar(aes(ymin = BETA_ADD - 1.96*SE_ADD, 
                    ymax = BETA_ADD + 1.96*SE_ADD), 
                alpha = 1, width = 0) +
  geom_errorbarh(aes(xmin = EXWAS_BETA_ADD - 1.96*EXWAS_SE_ADD, 
                     xmax = EXWAS_BETA_ADD + 1.96*EXWAS_SE_ADD), 
                 alpha = 1, height = 0) +
  geom_point(alpha = 1, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Additive Models [0,1,2]",
    x = "Dhindsa 2023, Beta (additive)",
    y = "Lassen 2025, Beta (Additive)",
    subtitle = sprintf("n = %d, r = %.3f", 
                       sum(!is.na(combined$BETA_ADD) & !is.na(combined$EXWAS_BETA_ADD)),
                       cor(combined$EXWAS_BETA_ADD, combined$BETA_ADD, use = "complete.obs"))
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


# Convert natural log to log10: log10(OR) = ln(OR) / ln(10)
combined[, `:=`(
  EXWAS_LOG10OR_ADD = EXWAS_BETA_ADD / log(10),
  EXWAS_LOG10OR_REC = EXWAS_BETA_REC / log(10),
  LOG10OR_ADD = BETA_ADD / log(10),
  LOG10OR_REC = BETA_REC / log(10),
  # SE also needs to be converted
  EXWAS_LOG10OR_SE_ADD = EXWAS_SE_ADD / log(10),
  EXWAS_LOG10OR_SE_REC = EXWAS_SE_REC / log(10),
  LOG10OR_SE_ADD = SE_ADD / log(10),
  LOG10OR_SE_REC = SE_REC / log(10)
)]


# Create labels for points (you can customize which points to label)
combined[, LABEL := paste0(GENE, "-", PROTEIN)]

# only those where we have betas
combined <- combined[!is.na(combined$BETA_ADD) & !is.na(combined$BETA_REC) & 
  !is.na(combined$EXWAS_BETA_ADD) & !is.na(combined$EXWAS_BETA_REC),]

# Recessive model comparison
ggplot(combined, aes(x = EXWAS_LOG10OR_REC, y = LOG10OR_REC)) +
  geom_errorbar(aes(ymin = LOG10OR_REC - 1.96*LOG10OR_SE_REC, 
                    ymax = LOG10OR_REC + 1.96*LOG10OR_SE_REC), 
                alpha = 0.3, width = 0) +
  geom_errorbarh(aes(xmin = EXWAS_LOG10OR_REC - 1.96*EXWAS_LOG10OR_SE_REC, 
                     xmax = EXWAS_LOG10OR_REC + 1.96*EXWAS_LOG10OR_SE_REC), 
                 alpha = 0.3, height = 0) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_text_repel(aes(label = LABEL),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Recessive Models [0,0,2]",
    x = TeX("Dhindsa 2023, $\\log_{10}(OR)$ (Recessive)"),
    y = TeX("Lassen 2025, $\\log_{10}(OR)$ (Recessive)"),
    subtitle = sprintf("n = %d, r = %.3f", 
                       sum(!is.na(combined$LOG10OR_REC) & !is.na(combined$EXWAS_LOG10OR_REC)),
                       cor(combined$EXWAS_LOG10OR_REC, combined$LOG10OR_REC, use = "complete.obs"))
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("~/Downloads/dhinsa_lassen_rec.png", width=4.5, height=4.5)

# Additive model comparison
ggplot(combined, aes(x = EXWAS_LOG10OR_ADD, y = LOG10OR_ADD)) +
  geom_errorbar(aes(ymin = LOG10OR_ADD - 1.96*LOG10OR_SE_ADD, 
                    ymax = LOG10OR_ADD + 1.96*LOG10OR_SE_ADD), 
                alpha = 0.3, width = 0) +
  geom_errorbarh(aes(xmin = EXWAS_LOG10OR_ADD - 1.96*EXWAS_LOG10OR_SE_ADD, 
                     xmax = EXWAS_LOG10OR_ADD + 1.96*EXWAS_LOG10OR_SE_ADD), 
                 alpha = 0.3, height = 0) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_text_repel(aes(label = LABEL),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Additive Models [0,1,2]",
    x = TeX("Dhindsa 2023, $\\log_{10}(OR)$ (Additive)"),
    y = TeX("Lassen 2025, $\\log_{10}(OR)$ (Additive)"),
    subtitle = sprintf("n = %d, r = %.3f", 
                       sum(!is.na(combined$LOG10OR_ADD) & !is.na(combined$EXWAS_LOG10OR_ADD)),
                       cor(combined$EXWAS_LOG10OR_ADD, combined$LOG10OR_ADD, use = "complete.obs"))
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("~/Downloads/dhinsa_lassen_add.png", width=4.5, height=4.5)