# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("11_snRNAseq_doublet_filt.R")
message("######################################")

if(!exists("snRNA_so_clust")) {
  snRNA_so_clust <- readRDS(paste(dir_list$rds, "snRNA_so_clust.rds", sep = "/"))
}

set.seed(seed.use)
bp <- BiocParallel::SnowParam(workers = 3, RNGseed = seed.use, type = "SOCK")

message("Predicting doublets with scDblFinder.")

scdbl_res <- snRNA_so_clust %>% 
  as.SingleCellExperiment() %>% 
  scDblFinder::scDblFinder(samples = "sample_id", clusters = "cluster_ann", BPPARAM = bp) %>% 
  colData(.) %>% 
  as.data.frame() %>% 
  select(starts_with("scDblFinder")) %>% 
  arrange(match(rownames(.), colnames(snRNA_so_clust)))

message("Doublets predicted.")

snRNA_so_clust_doubletfilt <- AddMetaData(snRNA_so_clust, metadata = scdbl_res) %>% 
  subset(scDblFinder.class == "singlet")

message("Doublets removed.")

saveRDS(snRNA_so_clust_doubletfilt, file = paste(dir_list$rds, "snRNA_so_clust_doubletfilt.rds", sep = "/"))
rm(bp, scdbl_res)
gc()

message("Data saved.")
