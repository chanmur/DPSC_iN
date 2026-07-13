# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("04_bulkRNAseq_analysis.R")
message("######################################")

##
# 1 | Data Import 
##

# Importing gene annotations as data table for faster operations
gene_ann <- fread(paste(dir_list$data, "gencode.v47.annotation.gtf", sep = "/"), data.table = TRUE, skip = 5)
gene_ann <- gene_ann[V3 == "gene"]
gene_ann[,`:=` (gene_id = sapply(V9, extract_gene_id), gene_name = sapply(V9, extract_gene_name), gene_type = sapply(V9, extract_gene_type))][, paste("V", 1:9, sep = "") := NULL]

message("Imported Gencode v47 annotations.")
# Importing raw gene counts as data table for faster operations
rna_counts_dt <- fread(file = paste(dir_list$data, "DPSC_iN_bulkRNAseq_raw_gene_counts_batch_merged.csv", sep = "/"))

message("Imported bulk RNAseq counts.")
# Removing columns not required
rna_counts_dt[, grep("Geneid|DPSC[0-9]+_.*", names(rna_counts_dt), invert = TRUE, value = TRUE) := NULL]

# Merging with the annotation table to include gene symbols
rna_counts_dt <- rna_counts_dt[gene_ann, on = .(Geneid = gene_id)]

# Selecting and renaming the count columns
count_cols <- grep("DPSC[0-9]+_.*", names(rna_counts_dt), value = TRUE)
setnames(rna_counts_dt, count_cols, sapply(count_cols, gsub, pattern = "_Aligned.sortedByCoord.out.bam", replacement = ""))

# Converting back to data frame
rna_counts_raw <- as.data.frame(rna_counts_dt)

rna_counts_raw <- rna_counts_raw %>% 
  rename("gene_id" = "Geneid") %>% # Renaming gene_id column for consistency
  mutate(gene_id = gsub(pattern = "\\.\\d+", replacement = "", x = gene_id, perl = TRUE),
         gene_name_unique = make.unique(gene_name, sep = "_")) %>% # Making duplicate gene symbols unique
  filter(rowSums(across(matches("_(DPSC|iN)$"))) > 0 & gene_type == "protein_coding") %>% # Filtering to only include protein-coding genes with non-zero counts
  rename_with(.cols = starts_with("DPSC"), ~ as.character(sample_name_map[.])) # Renaming counts columns by sample IDs.

saveRDS(rna_counts_raw, paste(dir_list$rds, "rna_counts_raw.rds", sep = "/"))

##
# 2 | Differential Gene Expression Analysis - DESeq2
##
sample_info <- metadata %>% 
  column_to_rownames("sample_id")

# Extracting count columns
rna_exp <- rna_counts_raw %>% 
  # Setting gene_name as the rownames
  column_to_rownames("gene_name_unique")

# Setting columns in the same order as the rows of sample_info
rna_exp <- rna_exp[row.names(sample_info)[row.names(sample_info) %in% colnames(rna_exp)]]

# Removing rows with 0 counts across all samples
rna_exp <- rna_exp %>% filter(rowSums(.) > 0)

# Running DESeq2
rna_deseq_obj <- DESeq(DESeqDataSetFromMatrix(countData = rna_exp, colData = sample_info[row.names(sample_info) %in% colnames(rna_exp),], design = ~donor_id + cell_line))

saveRDS(rna_deseq_obj, paste(dir_list$rds, "rna_deseq_obj.rds", sep = "/"))

# Running pairwise comparisons between cell lines
rna_DPSC_vs_iN <- DESeq2::lfcShrink(rna_deseq_obj, contrast = c("cell_line", "iN", "DPSC"), type = "ashr", parallel = TRUE)

message("Completed DESeq2 differential gene expression analysis in bulk RNAseq data.")

rna_DPSC_vs_iN <- as.data.frame(rna_DPSC_vs_iN) %>% 
  rownames_to_column("gene_name_unique") %>% 
  inner_join(rna_counts_raw %>% select(gene_id, gene_name, gene_type, gene_name_unique), by = "gene_name_unique") %>% # Merging results with gene info
  filter(!is.na(padj)) %>% # Removing genes with adjusted p-value of NA
  mutate(status = case_when(padj < 0.05 & log2FoldChange > 0.5 ~ "Upregulated", # Setting criteria for sig. upregulated or downregulated genes
                            padj < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
                            TRUE ~ "Insignificant"),
         status = factor(status, levels = c("Upregulated", "Downregulated", "Insignificant")))

rna_DPSC_vs_iN_filt <- rna_DPSC_vs_iN %>%
  filter(status != "Insignificant")

rna_DPSC_vs_iN_filt %>%
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_genes_deseq2_fc0.5_padj0.05.csv", sep = "/"))

saveRDS(rna_DPSC_vs_iN, paste(dir_list$rds, "rna_DPSC_vs_iN.rds", sep = "/"))
saveRDS(rna_DPSC_vs_iN_filt, paste(dir_list$rds, "rna_DPSC_vs_iN_filt.rds", sep = "/"))

##
# 3 | Extracting transformed data
##

# Extracting normalized counts for visualisation
rna_exp_norm <- as.data.frame(counts(rna_deseq_obj, normalized = TRUE))

rna_exp_norm <- rna_exp_norm %>% 
  rownames_to_column(var = "gene_name_unique") %>%  
  inner_join(rna_counts_raw %>% select(gene_id, gene_name, gene_type, gene_name_unique), by = "gene_name_unique") # Merging with gene info

saveRDS(rna_exp_norm, paste(dir_list$rds, "rna_exp_norm.rds", sep = "/"))

rna_exp_norm %>% 
  select(-gene_name_unique) %>% 
  select(gene_name, gene_id, gene_type, starts_with("DPSC")) %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_bulkRNAseq_protein_coding_normalized_gene_counts.csv", sep = "/"), sep = ",")

# Extracting regularized logarithmic values for PCA
rna_counts_rlog <- rlog(rna_deseq_obj, blind = TRUE)

saveRDS(rna_counts_rlog, paste(dir_list$rds, "rna_counts_rlog.rds", sep = "/"))

message("Normalized bulk RNAseq counts generated.")

##
# 4 | GO Analysis - ORA with ClusterProfiler and Semantic Similarity with rrvgo
##

message("Performing gene over-representation test with ClusterProfiler.")

## ORA
# Extracting ENSEMBL gene IDs of upregulated genes
genes <- rna_DPSC_vs_iN_filt %>% 
  filter(status == "Upregulated") %>% 
  pull("gene_id")

rna_DPSC_vs_iN_ora_up <- clusterProfiler::enrichGO(gene = genes, universe = rna_DPSC_vs_iN$gene_id, ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENSEMBL", pAdjustMethod = "BH", readable = TRUE)

# Extracting ENSEMBL gene IDs of downregulated genes
genes <- rna_DPSC_vs_iN_filt %>% 
  filter(status == "Downregulated") %>% 
  pull("gene_id")

rna_DPSC_vs_iN_ora_down <- clusterProfiler::enrichGO(gene = genes, universe = rna_DPSC_vs_iN$gene_id, ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENSEMBL", pAdjustMethod = "BH", readable = TRUE)

message("Gene over-representation test completed.")

saveRDS(rna_DPSC_vs_iN_ora_up, paste(dir_list$rds, "rna_DPSC_vs_iN_ora_up.rds", sep = "/"))
saveRDS(rna_DPSC_vs_iN_ora_down, paste(dir_list$rds, "rna_DPSC_vs_iN_ora_down.rds", sep = "/"))

## Semantic Similarity Analysis

# RNA Up GO terms filtered
rna_ora_up_filt <- rna_DPSC_vs_iN_ora_up@result %>% filter(p.adjust < 0.05)
# Creating similarity matrix
rna_up_semsim_mat <- rrvgo::calculateSimMatrix(x = rna_ora_up_filt$ID, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
# Grouping similar terms
rna_up_semsim_grouped <- rrvgo::reduceSimMatrix(rna_up_semsim_mat, scores = setNames(-log10(rna_ora_up_filt$p.adjust), rna_ora_up_filt$ID), threshold = 0.7, orgdb = "org.Hs.eg.db")

# RNA Down GO terms filtered
rna_ora_down_filt <- rna_DPSC_vs_iN_ora_down@result %>% filter(p.adjust < 0.05)
# Creating similarity matrix
rna_down_semsim_mat <- rrvgo::calculateSimMatrix(x = rna_ora_down_filt$ID, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
# Grouping similar terms
rna_down_semsim_grouped <- rrvgo::reduceSimMatrix(rna_down_semsim_mat, scores = setNames(-log10(rna_ora_down_filt$p.adjust), rna_ora_down_filt$ID), threshold = 0.7, orgdb = "org.Hs.eg.db")

message("Semantic similarity grouping completed.")

rna_ora_up_filt %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_genes_up_ora.csv", sep = "/"))

rna_ora_down_filt %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_genes_down_ora.csv", sep = "/"))

rna_up_semsim_grouped %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_genes_up_ora_rrvgo_grouped.csv", sep = "/"))

rna_down_semsim_grouped %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_exp_genes_down_ora_rrvgo_grouped.csv", sep = "/"))

saveRDS(rna_up_semsim_grouped, paste(dir_list$rds, "rna_up_semsim_grouped.rds", sep = "/"))
saveRDS(rna_down_semsim_grouped, paste(dir_list$rds, "rna_down_semsim_grouped.rds", sep = "/"))

# Removing objects not required any further
rm(genes, rna_ora_up_filt, rna_ora_down_filt, rna_up_semsim_mat, rna_down_semsim_mat, rna_counts_raw, gene_ann, rna_counts_dt, count_cols, sample_info, rna_exp, rna_deseq_obj)
gc()