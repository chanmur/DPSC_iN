# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("14_snRNAseq_trajectory_inf.R")
message("######################################")

##
# 1 | Computing pseudotime and predicting lineage trajectories
##

message("Computing pseudotime and predicting lineages.")

set.seed(seed.use)
lineages <- as.SlingshotDataSet(getLineages(
  data           = snRNA_so_clust_doubletfilt_module@reductions$umap@cell.embeddings,
  clusterLabels  = snRNA_so_clust_doubletfilt_module$cluster_ann,
  dist.method    = "slingshot",
  start.clus = "Fibroblast-Like_1",
  end.clus = c("GABAergic_Neurons", "Unknown_Fate")))

lineage_traj_curves <- as.SlingshotDataSet(getCurves(
  data          = lineages,
  thresh        = 1e-1,
  stretch       = 1e-1,
  allow.breaks  = F, approx_points = FALSE))

saveRDS(lineage_traj_curves, file = paste(dir_list$rds, "lineage_traj_curves.rds", sep = "/"))

pseudotime <- slingPseudotime(lineage_traj_curves, na = FALSE)
cell_weights <- slingCurveWeights(lineage_traj_curves)

weighted_pseudotime <- rowSums(pseudotime * cell_weights, na.rm = TRUE) / rowSums(cell_weights * !is.na(pseudotime), na.rm = TRUE)

# Adding pseudotime and weights into metadata
snRNA_so_final <- AddMetaData(snRNA_so_clust_doubletfilt_module, 
                              metadata = as.data.frame(pseudotime) %>%  setNames(paste0("pseudotime_", colnames(pseudotime))))
snRNA_so_final <- AddMetaData(snRNA_so_final,
                              metadata = as.data.frame(cell_weights) %>% setNames(paste0("cellWeight_", colnames(cell_weights))))
snRNA_so_final <- AddMetaData(snRNA_so_final,
                              metadata = as.data.frame(weighted_pseudotime) %>% setNames("weighted_pseudotime"))

saveRDS(snRNA_so_final, file = paste(dir_list$rds, "snRNA_so_final.rds", sep = "/"))
## Saving the data as .h5ad file
scCustomize::as.anndata(x = snRNA_so_final, file_path = dir_list$h5ad, file_name = "snRNA_so_final.h5ad")

message("Final seurat object saved as .rds and .h5ad")

##
# 2 | Genes differentially expressed
##

## Identifying trajectory influencing genes

message("Fitting GAM model on highly variable genes. Takes time.")

# Fitting GAM model on highly variable genes
sceGAM <- tradeSeq::fitGAM(
   counts = snRNA_so_final@assays$RNA@data[VariableFeatures(snRNA_so_final),],
   pseudotime = pseudotime,
   cellWeights = cellWeights,
   nknots = 5, 
   verbose = T)

 saveRDS(sceGAM, file = paste(dir_list$rds, "sceGAM.rds", sep = "/"))

message("GAM model successfully generated and saved.")
message("Identifying genes diff exp in different lineages.") 

traj_gene_list <- list()

## Genes diff exp Start vs End in Lineages 1 and 2
traj_gene_list$StartEnd <- tradeSeq::startVsEndTest(models = sceGAM) %>% 
  mutate(adj_pvalue = p.adjust(p = pvalue, method = "BH"))

### lineage 1
traj_gene_list$StartEnd_lin1_filt <- traj_gene_list$StartEnd %>% 
  filter(adj_pvalue < 0.05) %>% 
  filter(waldStat > mean(waldStat)) %>% 
  filter(rownames(.) %in% protein_coding_genes & abs(logFClineage1) > 0.5) %>% 
  select(-logFClineage2) %>% 
  arrange(desc(logFClineage1), desc(waldStat))

### lineage 2
traj_gene_list$StartEnd_lin2_filt <- traj_gene_list$StartEnd %>% 
  filter(adj_pvalue < 0.05) %>% 
  filter(waldStat > mean(waldStat)) %>% 
  filter(rownames(.) %in% protein_coding_genes & abs(logFClineage2) > 0.5) %>% 
  select(-logFClineage1) %>% 
  arrange(desc(logFClineage2), desc(waldStat))

## Genes diff exp between end points of Lin 1 and 2
traj_gene_list$DiffEnd <- tradeSeq::diffEndTest(models = sceGAM) %>% 
  mutate(adj_pvalue = p.adjust(p = pvalue, method = "BH"))

traj_gene_list$DiffEnd_filt <- traj_gene_list$DiffEnd %>% 
  filter(adj_pvalue < 0.05) %>% 
  filter(waldStat > mean(waldStat)) %>% 
  filter(rownames(.) %in% protein_coding_genes & abs(logFC1_2) > 0.5) %>% 
  arrange(desc(logFC1_2), desc(waldStat))

## Saving gene lists
traj_gene_list$StartEnd_lin1_filt %>% 
  rownames_to_column("Gene") %>% 
  fwrite(paste(dir_list$docs, "snRNA_Lin1_StartvsEnd.csv", sep = "/"))

traj_gene_list$StartEnd_lin2_filt %>% 
  rownames_to_column("Gene") %>% 
  fwrite(paste(dir_list$docs, "snRNA_Lin2_StartvsEnd.csv", sep = "/"))

traj_gene_list$DiffEnd_filt %>% 
  rownames_to_column("Gene") %>% 
  fwrite(paste(dir_list$docs, "snRNA_Lin1vs2_DiffEnd.csv", sep = "/"))
# traj_gene_list$BranchingPoint <- tradeSeq::earlyDETest(models = sceGAM) %>% 
#   mutate(adj_pvalue = p.adjust(p = pvalue, method = "BH"))

# traj_gene_list$BranchingPoint %>% filter(!is.na(waldStat)) %>% filter(pvalue < 0.05 & (waldStat > mean(waldStat)) & rownames(.) %in% protein_coding_genes & fcMedian > 0.5) %>% select(waldStat, fcMedian) %>% arrange(desc(fcMedian), desc(waldStat)) %>% View()

saveRDS(traj_gene_list, file = paste(dir_list$rds, "traj_gene_list.rds", sep = "/"))

##
# 3 | GO - ORA
##

message("Gene over-representation on genes diff exp in the trajectories.")

traj_gene_ora <- list()

traj_gene_ora$StartEnd_lin1_up <- traj_gene_list$StartEnd_lin1_filt %>% 
  filter(logFClineage1 > 0) %>%
  rownames(.) %>% 
  clusterProfiler::enrichGO(gene = ., universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)

traj_gene_ora$StartEnd_lin1_up <- traj_gene_ora$StartEnd_lin1_up@result %>% 
  filter(p.adjust < 0.05)

traj_gene_ora$StartEnd_lin1_up %>% 
  fwrite(paste(dir_list$docs, "snRNA_Lin1_StartvsEnd_Up_GO_ORA.csv", sep = "/"))

traj_gene_ora$StartEnd_lin1_down <- traj_gene_list$StartEnd_lin1_filt %>% 
  filter(logFClineage1 < 0) %>%
  rownames(.) %>% 
  clusterProfiler::enrichGO(gene = ., universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)

traj_gene_ora$StartEnd_lin1_down <- traj_gene_ora$StartEnd_lin1_down@result %>% 
  filter(p.adjust < 0.05)

traj_gene_ora$StartEnd_lin1_down %>% 
  fwrite(paste(dir_list$docs, "snRNA_Lin1_StartvsEnd_Down_GO_ORA.csv", sep = "/"))

traj_gene_ora$StartEnd_lin2_up <- traj_gene_list$StartEnd_lin2_filt %>% 
  filter(logFClineage2 > 0) %>%
  rownames(.) %>% 
  clusterProfiler::enrichGO(gene = ., universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)

traj_gene_ora$StartEnd_lin2_up <- traj_gene_ora$StartEnd_lin2_up@result %>% 
  filter(p.adjust < 0.05)

traj_gene_ora$StartEnd_lin2_up %>% 
  fwrite(paste(dir_list$docs, "snRNA_Lin2_StartvsEnd_Up_GO_ORA.csv", sep = "/"))

traj_gene_ora$StartEnd_lin2_down <- traj_gene_list$StartEnd_lin2_filt %>% 
  filter(logFClineage2 < 0) %>%
  rownames(.) %>% 
  clusterProfiler::enrichGO(gene = ., universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)

traj_gene_ora$StartEnd_lin2_down <- traj_gene_ora$StartEnd_lin2_down@result %>% 
  filter(p.adjust < 0.05)

traj_gene_ora$StartEnd_lin2_down %>% 
  fwrite(paste(dir_list$docs, "snRNA_Lin2_StartvsEnd_Down_GO_ORA.csv", sep = "/"))

traj_gene_ora$DiffEnd_up <- traj_gene_list$DiffEnd_filt %>% 
  filter(logFC1_2 > 0) %>%
  rownames(.) %>% 
  clusterProfiler::enrichGO(gene = ., universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)

traj_gene_ora$DiffEnd_up <- traj_gene_ora$DiffEnd_up@result %>% 
  filter(p.adjust < 0.05)

traj_gene_ora$DiffEnd_up %>% 
  fwrite(paste(dir_list$docs, "snRNA_DiffEnd_Up_GO_ORA.csv", sep = "/"))

traj_gene_ora$DiffEnd_down <- traj_gene_list$DiffEnd_filt %>% 
  filter(logFC1_2 < 0) %>%
  rownames(.) %>% 
  clusterProfiler::enrichGO(gene = ., universe = rownames(snRNA_so_clust_prot), ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", pAdjustMethod = "BH", readable = TRUE)

traj_gene_ora$DiffEnd_down <- traj_gene_ora$DiffEnd_down@result %>% 
  filter(p.adjust < 0.05)

traj_gene_ora$DiffEnd_down %>% 
  fwrite(paste(dir_list$docs, "snRNA_DiffEnd_Down_GO_ORA.csv", sep = "/"))

saveRDS(traj_gene_ora, file = paste(dir_list$rds, "traj_gene_ora.rds", sep = "/"))

message("End of scripts.")