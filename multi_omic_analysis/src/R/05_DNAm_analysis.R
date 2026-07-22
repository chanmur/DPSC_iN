# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("05_DNAm_analysis.R")
message("######################################")

##
# 1 | Data Import
##

## IDAT Files
# Path with IDAT files
idat_dir <- paste(dir_list$data, "dnam_idat_files", sep = "/")

# Importing call rate outliers
call_rate_outliers <- fread(paste(dir_list$data, "AppendixB_Probes_CallRates_outliers.txt", sep = "/"), data.table = FALSE)

# Importing bead number outliers
bead_num_outliers <- fread(paste(dir_list$data, "AppendixC_Probes_BeadNumber_outliers.txt", sep = "/"), data.table = FALSE)

# Importing IDAT path files
idat_files <- sesame::searchIDATprefixes(idat_dir)

message("Importing IDAT files as beta values.")

# Importing data as a matrix of beta values
beta_values <- sesame::openSesame(idat_files)

# Replacing IDAT IDs with sample IDs in beta value matrix
colnames(beta_values) <- metadata %>% 
  arrange(match(DNAm_idat_id, colnames(beta_values))) %>% 
  pull("sample_id")

message("Beta values computed.")

saveRDS(beta_values, paste(dir_list$rds, "beta_values.rds", sep = "/"))

## EPICv2 annotations
# Importing EPICv2 hg38 Gencode v41 annotations
epicv2_gencode_ann <- fread(paste(dir_list$data, "EPICv2.hg38.manifest.gencode.v41.tsv", sep = "/"), data.table = FALSE)

epicv2_gencode_ann <- epicv2_gencode_ann %>% 
  rename("probe_id" = "probeID",
         "gene_name" = "genesUniq") %>% # Changing column names
  filter(grepl("^cg", probe_id) & !is.na(CpG_chrm)) # Only including CpG IDs, starting with cg, and excluding sites marked NA

## MANE annotations
# Importing MANE GrCh38 v1.4 annotations
mane_gtf <- rtracklayer::import(paste(dir_list$data, "MANE.GRCh38.v1.4.gtf", sep = "/"), format = "gtf")

# Extract only transcripts - one transcript per row
mane_gtf <- mane_gtf[mane_gtf$type == "transcript" & mane_gtf$tag == "MANE_Select"]

# Setting annotation style based on UCSC
GenomeInfoDb::seqlevelsStyle(mane_gtf) <- "UCSC"

# Converting the transcript coordinates into promoter coordinates:
# TSS is position 1 of the transcript, upstream value is 1500bp in this case, while downstream is 500bp
# For + strand: [TSS - upstream, TSS + downstream]
# For - strand: [TSS - downstream, TSS + upstream]
# This would help identify the CpGs that fall within this promoter interval
promoters_mane <- GenomicRanges::promoters(mane_gtf, upstream = 1500, downstream = 500)

##
# 2 | Quality Control
##

message("Removing low-quality probes.")

# Gathering indices of low quality flagged probes to be filtered out
facility_flagged_probes <- c(call_rate_outliers$Probe_Id, bead_num_outliers$Probe_Id) %>% unique()
probe_prefix <- stringi::stri_replace_first_regex(rownames(beta_values), pattern = "_.*", replacement = "")
exclude_probes_index <- unique(c(which(probe_prefix %in% facility_flagged_probes), which(rownames(beta_values) %in% facility_flagged_probes)))

# Removing SNP and cross-hybrid probes
beta_values_filt <- DMRcate::rmSNPandCH(beta_values)

# Filtering to only include CpGs
beta_values_filt <- beta_values_filt[grep("^cg[0-9]+.*", rownames(beta_values_filt)),]

# Filtering to exclude flagged probes
beta_values_filt <- beta_values_filt[-exclude_probes_index,]

# Filtering out probes with NA beta in all samples
beta_values_filt <- beta_values_filt[rowSums(!is.na(beta_values_filt)) > 0,]

# Filtering to only include probes in the EPICv2 annotation file
beta_values_filt <- beta_values_filt[rownames(beta_values_filt) %in% epicv2_gencode_ann$probe_id,]

message("Low-quality probes removed.")

saveRDS(beta_values_filt, paste(dir_list$rds, "beta_values_filt.rds", sep = "/"))

beta_values_filt %>% 
  as.data.frame() %>% 
  rownames_to_column("Probe_ID") %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_DNAm_betaValues_processed.csv", sep = "/"), col.names = TRUE, sep = ",")

##
# 3 | Aggregating Probes by IDs
##

message("Converting beta values to m-values")

## M-values
# Converting beta values to M values for downstream analysis
m_values_filt <- sesame::BetaValueToMValue(beta_values_filt)

# Merging with Sesame's EPIC v2 manifest annotations
m_values_filt <- m_values_filt %>% 
  as.data.frame() %>% 
  rownames_to_column("probe_id") %>% 
  inner_join(epicv2_gencode_ann, by = "probe_id")

message("Aggregating m-values across duplicate probes")
## Aggregation
# Converting into data table for faster operation and removing suffixes from Probe IDs to later merge probe replicates
m_values_filt_aggr <- as.data.table(m_values_filt)[, probe_id := stringi::stri_replace_first_regex(probe_id, pattern = "_.*", replacement = "")]

setkey(m_values_filt_aggr, probe_id, CpG_chrm, CpG_beg, CpG_end)
val_cols <- grep("^DPSC(-iN)?[0-4]", names(m_values_filt_aggr), value = TRUE)

# Setting value columns as numeric
m_values_filt_aggr[, (val_cols) := lapply(.SD, as.numeric), .SDcols = val_cols]

# Aggregating by average of M values across probe replicates
m_values_filt_aggr <- m_values_filt_aggr[, c(lapply(.SD, mean, na.rm = TRUE), .(n_probes = .N)), by = .(probe_id, CpG_chrm, CpG_beg, CpG_end), .SDcols = val_cols]

# Converting back into data frame
setDF(m_values_filt_aggr)

# Converting into genomic ranges object for downstream operations on genomic coordinates
m_values_filt_aggr_gr <- GenomicRanges::GRanges(seqnames = m_values_filt_aggr$CpG_chrm, 
                                                    ranges = IRanges(start = m_values_filt_aggr$CpG_beg, end = m_values_filt_aggr$CpG_end),
                                                    probe_id = m_values_filt_aggr$probe_id)

saveRDS(m_values_filt_aggr, paste(dir_list$rds, "m_values_filt_aggr.rds", sep = "/"))


##
# 4 | Promoter Methylation
##

message("Aggregating m-values (methylation signals) over promoter regions")

## Methylation per promoter annotation
# Finding CpGs in the mane promoter coordinates. And, ignoring strand since CpG methylation is strand symmetric
cpg_promoters <- findOverlaps(m_values_filt_aggr_gr, promoters_mane, ignore.strand = TRUE)

# Filtering probes mapping to multiple promoters - only selecting those mappings that are closest to a TSS
# Extracting tss coordinates
tss_coords <- ifelse(strand(mane_gtf) == "+", start(mane_gtf), end(mane_gtf))
tss_gr <- GRanges(seqnames(mane_gtf),
                  IRanges(start = tss_coords, width = 1),
                  strand = strand(mane_gtf))
mcols(tss_gr) <- mcols(mane_gtf)[, c("gene_name", "gene_id", "transcript_id")]

# Mapping CpGs to genes based on MANE promoter annotation
# Calculating the distance to TSS for every promoter and filtering
cpg_promoter_gene_map_unique <- data.frame(probe_id = m_values_filt_aggr_gr$probe_id[queryHits(cpg_promoters)],
                                           promoter_idx = subjectHits(cpg_promoters),
                                           distance_to_tss = distance(m_values_filt_aggr_gr[queryHits(cpg_promoters)],
                                                                      tss_gr[subjectHits(cpg_promoters)],
                                                                      ignore.strand = FALSE)) %>% 
  group_by(probe_id) %>%
  slice_min(distance_to_tss, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(gene_name = tss_gr$gene_name[promoter_idx],
         gene_id = tss_gr$gene_id[promoter_idx],
         transcript_id = tss_gr$transcript_id[promoter_idx]) %>% 
  select(-promoter_idx, -distance_to_tss) %>% 
  as.data.frame()

# Filtering m values data frame to only include promoter CpGs
m_values_filt_aggr_promoters <- m_values_filt_aggr %>% 
  inner_join(cpg_promoter_gene_map_unique, by = "probe_id")

## Aggregation by gene
# Mean m values per gene (and transcript - to keep methylation per promoter distinct)
m_values_filt_aggr_promoters_per_gene <- m_values_filt_aggr_promoters %>%
  group_by(gene_name, transcript_id, gene_id) %>%
  summarise(across(matches("^DPSC(-iN)?[0-4]"), mean, na.rm = TRUE),
            n_cpgs = n_distinct(probe_id),
            probe_ids = paste(probe_id, collapse = ",")) %>%
  as.data.frame() %>% 
  filter(n_cpgs >= 3) # promoters with 1-2 cpgs can be quite noisy, so filtering to exclude these before diff meth analysis

saveRDS(m_values_filt_aggr_promoters_per_gene, paste(dir_list$rds, "m_values_filt_aggr_promoters_per_gene.rds", sep = "/"))


##
# 5 | Differential Methylation Analysis on Promoters - limma + BH
##

message("Differential methylation analysis over promoter regions.")

design <- model.matrix(~0 + metadata$donor_id + metadata$cell_line)
colnames(design) <- c(levels(metadata$donor_id), "iN_vs_DPSC")

m_values_exp <- m_values_filt_aggr_promoters_per_gene %>%
  column_to_rownames("gene_id") %>%
  select(matches("^DPSC(-iN)?[0-4]"))

fit <- limma::lmFit(m_values_exp, design)
fit_cont <- limma::eBayes(fit)

dnam_dmr <- limma::topTable(fit_cont, coef = "iN_vs_DPSC", adjust.method = "BH", confint = 0.95, number = Inf) %>%
  mutate(status = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Hypermethylated",
                            adj.P.Val < 0.05 & logFC < 0 ~ "Hypomethylated",
                            TRUE ~ "Insignificant"),
         status = factor(status, levels = c("Hypermethylated", "Hypomethylated", "Insignificant"))) %>%
  rownames_to_column("gene_id") %>%
  arrange(status, -logFC)

dnam_dmr_ann <- inner_join(dnam_dmr, m_values_filt_aggr_promoters_per_gene, by = "gene_id") %>% 
  mutate(across(matches("^DPSC(-iN)?[0-4]"), sesame::MValueToBetaValue))

dnam_dmr_ann_filt <- dnam_dmr_ann %>% filter(status != "Insignificant")

dnam_dmr_ann_filt %>% 
  fwrite(paste(dir_list$docs, "DPSC_iN_diff_meth_genes_limma_padj0.05.csv", sep = "/"))

message("Differential methylated promoters identifies and saved.")

saveRDS(dnam_dmr_ann, paste(dir_list$rds, "dnam_dmr_ann.rds", sep = "/"))
saveRDS(dnam_dmr_ann_filt, paste(dir_list$rds, "dnam_dmr_ann_filt.rds", sep = "/"))

rm(idat_dir, call_rate_outliers, bead_num_outliers, idat_files, beta_values, epicv2_gencode_ann, mane_gtf, promoters_mane, facility_flagged_probes, probe_prefix, exclude_probes_index, beta_values_filt, m_values_filt, val_cols, m_values_filt_aggr, m_values_filt_aggr_gr, cpg_promoters, tss_coords, tss_gr, cpg_promoter_gene_map_unique, m_values_filt_aggr_promoters, m_values_filt_aggr_promoters_per_gene, design, fit, fit_cont, dnam_dmr)
gc()