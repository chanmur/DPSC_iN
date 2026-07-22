# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("09_snRNAseq_norm_scaling_pca_harmonyInt.R")
message("######################################")

if(!exists("snRNA_so_filt")) {
  snRNA_so_filt <- readRDS(paste(dir_list$rds, "snRNA_so_filt.rds", sep = "/"))
}

message("Normalizing and scaling the data")
## Data Transformation
set.seed(seed.use)
snRNA_so_norm_pca_int <- snRNA_so_filt %>% 
  NormalizeData(verbose = F) %>% # Log Normalization 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, 
                       verbose = F) %>% 
  ScaleData(verbose = F) %>%
  RunPCA(features = VariableFeatures(.),
         seed.use = seed.use,
         verbose = F)

message("Harmony integration")

## Data Integration with Harmony
set.seed(seed.use)
snRNA_so_norm_pca_int <- IntegrateLayers(object = snRNA_so_norm_pca_int, 
                                         method = HarmonyIntegration,
                                         orig.reduction = "pca",
                                         new.reduction = "harmony", 
                                         verbose = FALSE,
                                         harmony.params = list(group.by.vars = "donor_id", seed.use = seed.use))

# Adding standard deviation information for the Harmony-corrected Principal Components
snRNA_so_norm_pca_int[["harmony"]]@stdev <- apply(Embeddings(snRNA_so_norm_pca_int, "harmony"), 2, sd)

# Joining RNA count layers for further data analysis
snRNA_so_norm_pca_int <- JoinLayers(snRNA_so_norm_pca_int)

saveRDS(snRNA_so_norm_pca_int, file = paste(dir_list$rds, "snRNA_so_norm_pca_int.rds", sep = "/"))
gc()
message("Data successfully transformed and saved.")