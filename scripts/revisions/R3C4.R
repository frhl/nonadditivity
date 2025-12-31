
library(data.table)
library(bravastring)
library(ggplot2)
library(latex2exp)

# olink traits from SAIGE
d <- fread("~/Projects/01_2df_git/data/olink_qc_combined_pLoF.txt")
d$conditional <- FALSE
d$conditional.dom <- FALSE
d$annotation <- "pLoF"
d$gene_id <- transcript_to_gene_id(d$MarkerID)
d[, id := paste0(hgnc_symbol, "_", trait)]
d <- d[!is.na(d$p.value.dom) & AC_Allele2.rec>=10]

# bonferroni
n_genes <- uniqueN(d$MarkerID)
n_traits <- uniqueN(d$trait)
pval_T <- 0.05/(n_genes*n_traits)

# get sig pairs (correct for additive first, then dom)
sig_d <- d[p.value<pval_T]
n_to_correct_for <- nrow(sig_d)
sig_pairs <- sig_d[p.value.dom<(0.05/n_to_correct_for)]

# load regenie results
regenie_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/regenie/olink"
olink_traits <- unique(d$trait)

regenie_olink_list <- list()
for (trait in olink_traits) {
  add_file <- file.path(regenie_dir, paste0("UKB.auto.", trait, ".additive.eur.af05.pp0.90.pLoF.regenie.gz"))
  if (file.exists(add_file)) {
    tmp <- fread(add_file)
    tmp$mode <- "additive"
    tmp$trait <- trait
    regenie_olink_list[[paste0(trait, "_additive")]] <- tmp
  }
  dom_file <- file.path(regenie_dir, paste0("UKB.auto.", trait, ".dominance.eur.af05.pp0.90.pLoF.regenie.gz"))
  if (file.exists(dom_file)) {
    tmp <- fread(dom_file)
    tmp$mode <- "dominance"
    tmp$trait <- trait
    regenie_olink_list[[paste0(trait, "_dominance")]] <- tmp
  }
}

regenie_olink <- rbindlist(regenie_olink_list, fill = TRUE)
regenie_olink$gene_id <- transcript_to_gene_id(regenie_olink$ID)
regenie_olink$hgnc_symbol <- gene_id_to_hgnc_symbol(regenie_olink$gene_id)
regenie_olink$MarkerID <- regenie_olink$ID
regenie_olink[, id := paste0(hgnc_symbol, "_", trait)]
regenie_olink$p.value <- 10^(-regenie_olink$LOG10P)
regenie_olink$Z <- regenie_olink$BETA / regenie_olink$SE

# merge saige/regenie for additive
d_add <- d[, .(MarkerID, trait, id, p.value, BETA, SE)]
d_add[, z := BETA/SE]
setnames(d_add, c("p.value", "BETA", "SE", "z"), c("saige_p", "saige_beta", "saige_se", "saige_z"))

r_add <- regenie_olink[mode == "additive", .(MarkerID, trait, id, p.value, BETA, SE, Z)]
setnames(r_add, c("p.value", "BETA", "SE", "Z"), c("regenie_p", "regenie_beta", "regenie_se", "regenie_z"))

m_add <- merge(d_add, r_add, by = c("MarkerID", "trait", "id"), all.x = TRUE)

# merge saige/regenie for dominance
d_dom <- d[, .(MarkerID, trait, id, p.value.dom, BETA.dom, SE.dom)]
d_dom[, z := BETA.dom/SE.dom]
setnames(d_dom, c("p.value.dom", "BETA.dom", "SE.dom", "z"), c("saige_p", "saige_beta", "saige_se", "saige_z"))

r_dom <- regenie_olink[mode == "dominance", .(MarkerID, trait, id, p.value, BETA, SE, Z)]
setnames(r_dom, c("p.value", "BETA", "SE", "Z"), c("regenie_p", "regenie_beta", "regenie_se", "regenie_z"))

m_dom <- merge(d_dom, r_dom, by = c("MarkerID", "trait", "id"), all.x = TRUE)

# regenie allele is flipped relative to SAIGE
m_add$regenie_z <- -m_add$regenie_z
m_dom$regenie_z <- -m_dom$regenie_z

# get significant gene-traits
m_dom$significant <- m_dom$id %in% sig_pairs$id
m_add$significant <- m_add$id %in% sig_pairs$id

# correlations for sig pairs
m_add_sig <- m_add[m_add$significant==TRUE,]
m_dom_sig <- m_dom[m_dom$significant==TRUE,]

# calculate correlation
dom_cor <- cor(m_dom$saige_z, m_dom$regenie_z, use="complete")

p_olink <- ggplot(m_dom, aes(x=saige_z, y=regenie_z, color=significant)) +
  geom_point(size=2) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey20") +
  scale_color_brewer(palette=3, type="qual") +
  annotate("text", x=min(m_dom$saige_z, na.rm=TRUE),
           y=max(m_dom$regenie_z, na.rm=TRUE),
           label=sprintf("italic(r) == %.3f", dom_cor),
           parse=TRUE, hjust=0, vjust=0.75, size=5) +
  labs(x=TeX("SAIGE $Z$-score"),
       y=TeX("REGENIE $Z$-score"),
       title="Nonadditive Encoding (Olink traits)") +
  theme_bw(base_size=10) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = "none")

p_olink

cat("\n=== OLINK ===\n")
cat("Significant pairs (SAIGE):", nrow(sig_pairs), "\n")
cat("With REGENIE data:", sum(!is.na(m_dom_sig$regenie_z)), "\n")
cat("Correlation (all):", dom_cor, "\n\n")


# medical/cts traits from SAIGE
d_cts <- fread("~/Projects/01_2df_git/data/cts_traits_qc_final_pLoF_damaging_missense.txt")
d_cts$annotation <- "pLoF_damaging_missense"
d_cts$gene_id <- transcript_to_gene_id(d_cts$MarkerID)
d_cts[, id := paste0(hgnc_symbol, "_", trait)]
d_cts <- d_cts[!is.na(p.value.dom) & AC_Allele2.rec>=10 & trait != "p47"]

# bonferroni
n_genes_cts <- uniqueN(d_cts$MarkerID)
n_traits_cts <- uniqueN(d_cts$trait)
pval_T_cts <- 0.05/(n_genes_cts*n_traits_cts)

# get sig pairs (correct for additive first, then dom)
sig_d_cts <- d_cts[p.value<pval_T_cts]
n_to_correct_for_cts <- nrow(sig_d_cts)
sig_pairs_cts <- sig_d_cts[p.value.dom<(0.05/n_to_correct_for_cts)]

# load regenie results (trait names have full descriptions, need to match prefix)
regenie_cts_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/regenie/medical"
cts_traits <- unique(d_cts$trait)

regenie_cts_list <- list()
all_files <- list.files(regenie_cts_dir, pattern = "pLoF_damaging_missense.regenie.gz$")
for (trait in cts_traits) {
  # find files that match this trait
  add_pattern <- paste0("^UKB.auto.", trait, "\\.")
  add_files <- grep(add_pattern, all_files, value=TRUE)
  add_files <- grep("additive", add_files, value=TRUE)
  if (length(add_files) > 0) {
    tmp <- fread(file.path(regenie_cts_dir, add_files[1]))
    tmp$mode <- "additive"
    tmp$trait <- trait
    regenie_cts_list[[paste0(trait, "_additive")]] <- tmp
  }

  dom_files <- grep(add_pattern, all_files, value=TRUE)
  dom_files <- grep("dominance", dom_files, value=TRUE)
  if (length(dom_files) > 0) {
    tmp <- fread(file.path(regenie_cts_dir, dom_files[1]))
    tmp$mode <- "dominance"
    tmp$trait <- trait
    regenie_cts_list[[paste0(trait, "_dominance")]] <- tmp
  }
}

if (length(regenie_cts_list) == 0) {
  stop("No regenie CTS files found! Check file paths and trait names.")
}

regenie_cts <- rbindlist(regenie_cts_list, fill = TRUE)
regenie_cts$gene_id <- transcript_to_gene_id(regenie_cts$ID)
regenie_cts$hgnc_symbol <- gene_id_to_hgnc_symbol(regenie_cts$gene_id)
regenie_cts$MarkerID <- regenie_cts$ID
regenie_cts[, id := paste0(hgnc_symbol, "_", trait)]
regenie_cts$p.value <- 10^(-regenie_cts$LOG10P)
regenie_cts$Z <- regenie_cts$BETA / regenie_cts$SE

# merge saige/regenie for additive
dc_add <- d_cts[, .(MarkerID, trait, id, p.value, BETA, SE)]
dc_add[, z := BETA/SE]
setnames(dc_add, c("p.value", "BETA", "SE", "z"), c("saige_p", "saige_beta", "saige_se", "saige_z"))

rc_add <- regenie_cts[mode == "additive", .(MarkerID, trait, id, p.value, BETA, SE, Z)]
setnames(rc_add, c("p.value", "BETA", "SE", "Z"), c("regenie_p", "regenie_beta", "regenie_se", "regenie_z"))

mc_add <- merge(dc_add, rc_add, by = c("MarkerID", "trait", "id"), all.x = TRUE)

# merge saige/regenie for dominance
dc_dom <- d_cts[, .(MarkerID, trait, id, p.value.dom, BETA.dom, SE.dom)]
dc_dom[, z := BETA.dom/SE.dom]
setnames(dc_dom, c("p.value.dom", "BETA.dom", "SE.dom", "z"), c("saige_p", "saige_beta", "saige_se", "saige_z"))

rc_dom <- regenie_cts[mode == "dominance", .(MarkerID, trait, id, p.value, BETA, SE, Z)]
setnames(rc_dom, c("p.value", "BETA", "SE", "Z"), c("regenie_p", "regenie_beta", "regenie_se", "regenie_z"))

mc_dom <- merge(dc_dom, rc_dom, by = c("MarkerID", "trait", "id"), all.x = TRUE)

# regenie allele is flipped relative to SAIGE
mc_add$regenie_z <- -mc_add$regenie_z
mc_dom$regenie_z <- -mc_dom$regenie_z

# get significant gene-traits
mc_dom$significant <- mc_dom$id %in% sig_pairs_cts$id
mc_add$significant <- mc_add$id %in% sig_pairs_cts$id

# correlations for sig pairs
mc_add_sig <- mc_add[mc_add$significant==TRUE,]
mc_dom_sig <- mc_dom[mc_dom$significant==TRUE,]

# calculate correlation
dom_cor_cts <- cor(mc_dom$saige_z, mc_dom$regenie_z, use="complete")

p_cts <- ggplot(mc_dom, aes(x=saige_z, y=regenie_z, color=significant)) +
  geom_point(size=2, data=mc_dom[mc_dom$significant==FALSE]) +
  geom_point(size=2, data=mc_dom[mc_dom$significant==TRUE]) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey20") +
  scale_color_brewer(palette=3, type="qual") +
  annotate("text", x=min(mc_dom$saige_z, na.rm=TRUE),
           y=max(mc_dom$regenie_z, na.rm=TRUE),
           label=sprintf("italic(r) == %.3f", dom_cor_cts),
           parse=TRUE, hjust=0, vjust=0.75, size=5) +
  labs(x=TeX("SAIGE $Z$-score"),
       y=TeX("REGENIE $Z$-score"),
       title="Nonadditive Encoding (Medical traits)") +
  theme_bw(base_size=10) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = "none")

p_cts

cat("\n=== MEDICAL/CTS ===\n")
cat("Significant pairs (SAIGE):", nrow(sig_pairs_cts), "\n")
cat("With REGENIE data:", sum(!is.na(mc_dom_sig$regenie_z)), "\n")
cat("Correlation (all):", dom_cor_cts, "\n\n")

# save plots to desktop
desktop_path <- "~/Desktop"
ggsave(file.path(desktop_path, "R3C4_olink_dominance.pdf"), p_olink, width=6, height=5)
#ggsave(file.path(desktop_path, "R3C4_medical_dominance.pdf"), p_cts, width=6, height=5)

cat("\nPlots saved to desktop:\n")
cat("  - R3C4_olink_dominance.pdf\n")
cat("  - R3C4_medical_dominance.pdf\n")
