# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("06_msProt_analysis.R")
message("######################################")

##
# 1 | Data Import
##

ms_prot_data <- fread(paste(dir_list$data, "Supplementary_Data_3.csv", sep = "/"), data.table = FALSE) %>% 
  select(`Protein FDR Confidence: Combined`, Accession, `Gene Symbol`, `Entrez Gene ID`, `Ensembl Gene ID`, `Abundances (by Bio. Rep.): iN,DPSC120214`, `Abundances (by Bio. Rep.): iN,DPSC14227`, `Abundances (by Bio. Rep.): iN,DPSC160609`, `Abundances (by Bio. Rep.): iN,DPSC190709`, `Abundances (by Bio. Rep.): native,DPSC120214`, `Abundances (by Bio. Rep.): native,DPSC14227`, `Abundances (by Bio. Rep.): native,DPSC160609`, `Abundances (by Bio. Rep.): native,DPSC190709`) %>% 
  rename("prot_conf" = "Protein FDR Confidence: Combined", # Renaming columns
         "accession" = "Accession",
         "gene_name" = "Gene Symbol",
         "gene_id" = "Ensembl Gene ID",
         "gene_id_entrez" = "Entrez Gene ID",
         "DPSC120214_iN" = "Abundances (by Bio. Rep.): iN,DPSC120214",
         "DPSC14227_iN" = "Abundances (by Bio. Rep.): iN,DPSC14227",
         "DPSC160609_iN" = "Abundances (by Bio. Rep.): iN,DPSC160609",
         "DPSC190709_iN" = "Abundances (by Bio. Rep.): iN,DPSC190709",
         "DPSC120214_DPSC" = "Abundances (by Bio. Rep.): native,DPSC120214",
         "DPSC14227_DPSC" = "Abundances (by Bio. Rep.): native,DPSC14227",
         "DPSC160609_DPSC" = "Abundances (by Bio. Rep.): native,DPSC160609",
         "DPSC190709_DPSC" = "Abundances (by Bio. Rep.): native,DPSC190709")

ms_prot_data <- ms_prot_data %>% 
  mutate(prot_conf = factor(prot_conf, levels = unique(prot_conf)))

message("MS data imported.")

##
# 2 | Filtering and Log-transforming the data
##

# Filtering the data to only include proteins identified with High confidence (FDR < 0.05) and have measurements (not NA)
ms_prot_data_filt <- ms_prot_data %>%
  filter(prot_conf == "High") %>%
  filter(if_all(matches("^DPSC[0-9]+_(DPSC|iN)"), ~ !is.na(.))) %>%
  rename_with(.cols = starts_with("DPSC"), ~ as.character(sample_name_map[.]))

# Log2 transforming the data
ms_prot_data_filt_log2 <- ms_prot_data_filt %>%
  mutate(across(matches("^DPSC(-iN)?[0-4]"), ~ log2(. + 1)))

message("MS data filtered and log-transformed.")

saveRDS(ms_prot_data_filt_log2, paste(dir_list$rds, "ms_prot_data_filt_log2.rds", sep = "/"))

##
# 3 | Differential Abundance Analysis - Limma + BH correction
##

message("Diff abundance testing with limma")

design <- model.matrix(~0 + metadata$donor_id + metadata$cell_line)
colnames(design) <- c(levels(metadata$donor_id), "iN_vs_DPSC")

ms_prot_data_exp_mat <- ms_prot_data_filt_log2 %>% 
  column_to_rownames("accession") %>% 
  select(metadata$sample_id) %>% 
  as.matrix()

fit <- limma::lmFit(ms_prot_data_exp_mat, design)
fit_cont <- limma::eBayes(fit) 
ms_prot_DPSC_vs_iN <- limma::topTable(fit_cont, coef = "iN_vs_DPSC", adjust.method = "BH", confint = 0.95, number = Inf) %>% 
  mutate(status = case_when(adj.P.Val < 0.05 & logFC > 0.5 ~ "Upregulated", 
                            adj.P.Val < 0.05 & logFC < -0.5 ~ "Downregulated",
                            TRUE ~ "Insignificant"),
         status = factor(status, levels = c("Upregulated", "Downregulated", "Insignificant"))) %>% 
  rownames_to_column("accession") %>% 
  inner_join(ms_prot_data_filt_log2 %>% select(accession, gene_name, gene_id, gene_id_entrez) %>% unique()) %>% 
  arrange(status, -logFC)

ms_prot_DPSC_vs_iN_filt <- ms_prot_DPSC_vs_iN %>% 
  filter(status != "Insignificant")

ms_prot_DPSC_vs_iN_filt %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_proteins_limma_fc0.5_padj0.05.csv", sep = "/"))

message("Diff abundant proteins identified and saved.")

saveRDS(ms_prot_DPSC_vs_iN, paste(dir_list$rds, "ms_prot_DPSC_vs_iN.rds", sep = "/"))
saveRDS(ms_prot_DPSC_vs_iN_filt, paste(dir_list$rds, "ms_prot_DPSC_vs_iN_filt.rds", sep = "/"))

##
# 4 | GO Analysis - Reactome Pathways
##

message("Identifying reactome pathways enriched.")

# Extracting ENTREZ gene IDs of upregulated proteins
genes <- ms_prot_DPSC_vs_iN_filt %>% 
  filter(status == "Upregulated") %>% 
  pull("gene_id_entrez")

ms_prot_DPSC_vs_iN_react_up <- ReactomePA::enrichPathway(gene = genes, organism = "human", universe = ms_prot_DPSC_vs_iN$gene_id_entrez, pAdjustMethod = "BH", readable = TRUE)

# Extracting ENTREZ gene IDs of downregulated proteins
genes <- ms_prot_DPSC_vs_iN_filt %>% 
  filter(status == "Downregulated") %>% 
  pull("gene_id_entrez")

ms_prot_DPSC_vs_iN_react_down <- ReactomePA::enrichPathway(gene = genes, organism = "human", universe = ms_prot_DPSC_vs_iN$gene_id_entrez, pAdjustMethod = "BH", readable = TRUE)

ms_prot_DPSC_vs_iN_react_up@result %>%
  filter(p.adjust < 0.05) %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_proteins_up_reactome_pathways.csv", sep = "/"))

ms_prot_DPSC_vs_iN_react_down@result %>%
  filter(p.adjust < 0.05) %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_proteins_down_reactome_pathways.csv", sep = "/"))

message("Over-represented reactome pathways identified and saved.")

saveRDS(ms_prot_DPSC_vs_iN_react_up, paste(dir_list$rds, "ms_prot_DPSC_vs_iN_react_up.rds", sep = "/"))
saveRDS(ms_prot_DPSC_vs_iN_react_down, paste(dir_list$rds, "ms_prot_DPSC_vs_iN_react_down.rds", sep = "/"))

# Removing objects not required any further
rm(genes, ms_prot_data, ms_prot_data_filt, design, ms_prot_data_exp_mat, fit, fit_cont)
gc()