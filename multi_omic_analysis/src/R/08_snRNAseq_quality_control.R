# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("08_snRNAseq_quality_control.R")
message("######################################")

if(!exists("snRNA_so")) {
  snRNA_so <- readRDS(paste(dir_list$rds, "snRNA_so.rds", sep = "/"))
}

snRNA_so_qc <- snRNA_so

# Adding QC parameters
snRNA_so_qc@meta.data <-  snRNA_so_qc@meta.data %>% 
  mutate(percent.mt = PercentageFeatureSet(snRNA_so_qc, pattern = "^MT-"),
         percent.rb = PercentageFeatureSet(snRNA_so_qc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"),
         log10GenesPerUMI = log10(nFeature_RNA) / log10(nCount_RNA))

# Genes with non-zero counts in at least 10 cells
keep_genes <- rownames(snRNA_so_qc)[Matrix::rowSums(JoinLayers(snRNA_so_qc)[["RNA"]]$counts>0) >= 10]

# Filtering cells and genes based on quality criteria
snRNA_so_filt <- snRNA_so_qc %>% 
  subset(nCount_RNA > 3000 & nCount_RNA < 60000 & nFeature_RNA > 1500 & nFeature_RNA < 10000 & percent.mt < 5, features = keep_genes)

saveRDS(snRNA_so_qc, file = paste(dir_list$rds, "snRNA_so_qc.rds", sep = "/"))
saveRDS(snRNA_so_filt, file = paste(dir_list$rds, "snRNA_so_filt.rds", sep = "/"))

message("Quality control completed. Filtered data has been saved.")

rm(keep_genes)
gc()