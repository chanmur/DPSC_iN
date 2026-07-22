# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("12_snRNAseq_clust_characterisation.R")
message("######################################")

if(!exists("snRNA_so_clust_doubletfilt")) {
  snRNA_so_clust_doubletfilt <- readRDS(paste(dir_list$rds, "snRNA_so_clust_doubletfilt.rds", sep = "/"))
}

## Filtering genes
# Extracting list of human protein-coding genes from org.Hs.eg.db
protein_coding_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "SYMBOL"), columns = c("GENETYPE"), keytype = "SYMBOL")  %>% 
  filter(GENETYPE == "protein-coding") %>% 
  dplyr::pull(SYMBOL) %>% 
  unique()

# Identifying the protein-coding genes in the dataset
genes_to_keep <- intersect(rownames(snRNA_so_clust_doubletfilt), protein_coding_genes)

# Excluding any mitochondrial genes
genes_to_keep <- setdiff(genes_to_keep, grep("^MT-", genes_to_keep, value = TRUE))

# Subsetting the data to only include the filtered genes
snRNA_so_clust_prot <- subset(snRNA_so_clust_doubletfilt, features = genes_to_keep)

## Identifying marker genes
cluster_markers <- FindAllMarkers(snRNA_so_clust_prot)

cluster_markers_filt <- cluster_markers %>% 
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5 & abs(pct.1 - pct.2) > 0.2) %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE)

cluster_markers_filt %>% 
  select(cluster, avg_log2FC, pct.1, pct.2, gene, p_val_adj) %>% 
  fwrite(paste(dir_list$docs, "Cluster_Markers.csv", sep = "/"))

## GO Analysis - ORA
cluster_markers_go <- purrr::map(levels(snRNA_so_clust_prot), function(cluster_ann) {
  
  go_res <- list()
  
  genes <- cluster_markers_filt %>% 
    filter(cluster == cluster_ann & avg_log2FC > 0) %>%
    pull("gene")
  
  go_res[["Upregulated"]] <- clusterProfiler::enrichGO(gene = genes, universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)
  
  go_res[["Upregulated"]]@result %>% 
    fwrite(paste(dir_list$docs, paste(paste("Cluster", cluster_ann, "GO", sep = "_"), "csv", sep = "."), sep = "/"))
  
  
  genes <- cluster_markers_filt %>% 
    filter(cluster == cluster_ann & avg_log2FC < 0) %>%
    pull("gene")
  
  go_res[["Downregulated"]] <- clusterProfiler::enrichGO(gene = genes, universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)
  
  go_res[["Downregulated"]]@result %>% 
    fwrite(paste(dir_list$docs, paste(paste("Cluster", cluster_ann, "GO", sep = "_"), "csv", sep = "."), sep = "/"))
  
  return(go_res)
})

names(cluster_markers_go) <- as.character(levels(snRNA_so_clust_prot))

saveRDS(cluster_markers, paste(dir_list$rds, "cluster_markers.rds", sep = "/"))
saveRDS(cluster_markers_go, paste(dir_list$rds, "cluster_markers_go.rds", sep = "/"))

rm(protein_coding_genes, genes_to_keep)
gc()