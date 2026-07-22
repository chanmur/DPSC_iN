# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("10_snRNAseq_clustering.R")
message("######################################")

if(!exists("snRNA_so_norm_pca_int")) {
  snRNA_so_norm_pca_int <- readRDS(paste(dir_list$rds, "snRNA_so_norm_pca_int.rds", sep = "/"))
}

set.seed(seed.use)
snRNA_so_clust <- snRNA_so_norm_pca_int %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(random.seed = seed.use, resolution = 0.6) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umap", seed.use = seed.use)

message("Clustering complete.")

snRNA_so_clust <- RenameIdents(snRNA_so_clust,
                               "3" = "Fibroblast-Like_1",
                               "8" = "Fibroblast-Like_2",
                               "6" = "Fibroblast-Like_3",
                               "4" = "Intermediary_1",
                               "2" = "Intermediary_2",
                               "1" = "Maturing_Neuron_1",
                               "0" = "Maturing_Neuron_2",
                               "7" = "GABAergic_Neurons",
                               "5" = "Unknown_Fate")

message("Clusters annotated.")

snRNA_so_clust$cluster_ann <- Idents(snRNA_so_clust)
saveRDS(snRNA_so_clust, file = paste(dir_list$rds, "snRNA_so_clust.rds", sep = "/"))
gc()

message("Data saved.")