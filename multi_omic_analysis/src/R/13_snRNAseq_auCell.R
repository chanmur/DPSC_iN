# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("13_snRNAseq_auCell.R")
message("######################################")

## Gene sets
kegg <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2")
react <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
syngo <- fread(paste(dir_list$data, "syngo1.3_ontologies.csv", sep = "/"), sep = ",", data.table = FALSE)

gene_list <- list()

gene_list[["Oxidative_phosphorylation"]] <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2") %>% filter(gs_name == "KEGG_OXIDATIVE_PHOSPHORYLATION") %>% pull("gene_symbol") %>% as.character()

gene_list[["Glycolysis"]] <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2") %>% filter(gs_name == "KEGG_GLYCOLYSIS_GLUCONEOGENESIS") %>% pull("gene_symbol") %>% as.character()

gene_list[["Lipid_metabolism"]] <- react %>% filter(gs_name == "REACTOME_METABOLISM_OF_LIPIDS") %>% pull("gene_symbol") %>% as.character()

gene_list[["WNT_signaling"]] <- react %>% filter(gs_name == "REACTOME_SIGNALING_BY_WNT") %>% pull("gene_symbol") %>% as.character()

gene_list[["NOTCH_signaling"]] <- react %>% filter(gs_name == "REACTOME_SIGNALING_BY_NOTCH") %>% pull("gene_symbol") %>% as.character()

gene_list[["ECM_KEGG"]] <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2") %>% filter(gs_name == "KEGG_ECM_RECEPTOR_INTERACTION") %>% pull("gene_symbol") %>% as.character()

gene_list[["ECM_React"]] <- react %>% filter(gs_name == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION") %>% pull("gene_symbol") %>% as.character()

gene_list[["BMP_signaling"]] <- react %>% filter(gs_name == "REACTOME_SIGNALING_BY_BMP") %>% pull("gene_symbol") %>% as.character()

# gene_list[["presynapse"]] <- get_syngo_genes(syngo, "presynapse$", "CC")
# gene_list[["postsynapse"]] <- get_syngo_genes(syngo, "postsynapse$", "CC")
gene_list[["presynaptic_signaling"]] <- get_syngo_genes(syngo, "SYNGO:presynprocess", "BP")
gene_list[["presynaptic_signaling"]] <- get_syngo_genes(syngo, "SYNGO:postsynprocess", "BP")
gene_list[["synapse_organisation"]] <- get_syngo_genes(syngo, "GO:0050808", "BP")
gene_list[["synaptic_signaling"]] <- get_syngo_genes(syngo, "GO:0099536", "BP")
gene_list[["synaptic_transport"]] <- get_syngo_genes(syngo, "SYNGO:transport", "BP")
gene_list[["synaptic_metabolism"]] <- get_syngo_genes(syngo, "SYNGO:metabolism", "BP")
gene_list[["trans-synaptic_signaling"]] <- get_syngo_genes(syngo, "GO:0099537", "BP")

gene_list[["sasp"]] <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS") %>% filter(gs_name == "WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP") %>% pull("gene_symbol") %>% as.character()

gene_list[["sen_mayo"]] <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2") %>% filter(gs_name == "SAUL_SEN_MAYO") %>% pull("gene_symbol") %>% as.character()

## AU Cell
set.seed(seed.use)
cell_rankings <- AUCell::AUCell_buildRankings(exprMat = GetAssayData(object = snRNA_so_clust_prot, layer = "counts"))

cells_AUC <- AUCell::AUCell_calcAUC(gene_list, cell_rankings)

auc_matrix <- AUCell::getAUC(cells_AUC)

snRNA_so_clust_doubletfilt_module <- snRNA_so_clust_doubletfilt

for(gs in rownames(auc_matrix)) {
  snRNA_so_clust_doubletfilt_module <- AddMetaData(snRNA_so_clust_doubletfilt_module, 
                                                   metadata = auc_matrix[gs,],
                                                   col.name = gs)
}

saveRDS(snRNA_so_clust_doubletfilt_module, paste(dir_list$rds, "snRNA_so_clust_doubletfilt_module.rds", sep = "/"))
