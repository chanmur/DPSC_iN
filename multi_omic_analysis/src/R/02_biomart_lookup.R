# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("02_biomart_lookup.R")
message("######################################")


biomart_rds_path <- paste(dir_list$rds, "biomart_mapping_ensembl_113.rds", sep = "/")

### Checking if the biomart object already exists ###
if (file.exists(biomart_rds_path)) {
  
  message("Pre-saved BioMart mapping found.")
  message("Loading BioMart mapping from: ", biomart_rds_path)
  
  biomart_mapping_ensembl_113 <- readRDS(biomart_rds_path)
  
  message("BioMart mapping successfully loaded.")
  
} else {
  ### Generating a biomart lookup table for fetch gene ids ###
  # Connecting to Ensembl to fetch the data.
  tryCatch({
    message("No pre-saved BioMart mapping found.")
    message("Attempting to retrieve the mapping from Ensembl v113 ...")
    
    biomart_obj <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 113)
    
    # Exporting two mapping tables since the getBM function limits the number of external columns that can be queried at once
    
    message("Retreiving tables ...")
    
    biomart_mapping_uniprot_1 <- biomaRt::getBM(attributes = c("uniprot_gn_id", "hgnc_symbol", "ensembl_gene_id", "entrezgene_id"), mart = biomart_obj)
    
    biomart_mapping_uniprot_2 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "uniprotswissprot", "uniprotsptrembl"), mart = biomart_obj)
    
    # Merging the tables by ensembl ids
    biomart_mapping_ensembl_113 <- merge(biomart_mapping_uniprot_1, biomart_mapping_uniprot_2, by = "ensembl_gene_id", all = TRUE)
    
    message("Tables successfully retreived.")
    
    # Converting into a data table for faster operations
    setDT(biomart_mapping_ensembl_113)
    
    # Converting into long format for lookup
    biomart_mapping_ensembl_113 <- melt(biomart_mapping_ensembl_113, id.vars = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), measure.vars = c("uniprot_gn_id", "uniprotswissprot", "uniprotsptrembl"), value.name = "accession", na.rm = TRUE)
    
    # Removing rows with NA accession IDs
    biomart_mapping_ensembl_113 <- biomart_mapping_ensembl_113[!is.na(accession), .(ensembl_gene_id, hgnc_symbol, entrezgene_id, accession)] 
    
    # Setting keys for faster lookup
    setkey(biomart_mapping_ensembl_113, ensembl_gene_id, entrezgene_id, hgnc_symbol, accession)
    
    message("BioMart mapping saved to: ", biomart_rds_path)
    
    saveRDS(biomart_mapping_ensembl_113, biomart_rds_path)
    
    rm(biomart_mapping_uniprot_1, biomart_mapping_uniprot_2, biomart_rds_path)
    gc()
    
  }, error = function(e) {
      stop(paste0("Failed to retrieve the BioMart mapping from Ensembl, and no saved mapping file is available at:\n", biomart_rds_path, "\n\nOriginal error: ", conditionMessage(e)), call. = FALSE)
  })
}

