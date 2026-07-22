# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

##
# Functions to extract gene info from gtf file
##
extract_gene_name <- function(gene_info) {
  info_vec <- unlist(strsplit(gene_info, ";"))
  
  gene_name_ext <- grep("gene_name", info_vec, value = TRUE)
  gene_name <- gsub("( |gene_name|\")", "", gene_name_ext)
  
  return(gene_name)
}

extract_gene_id <- function(gene_info) {
  info_vec <- unlist(strsplit(gene_info, ";"))
  
  gene_id_ext <- grep("gene_id", info_vec, value = TRUE)
  gene_id <- gsub("( |gene_id|\")", "", gene_id_ext) 
  
  return(gene_id)
}

extract_gene_type <- function(gene_info) {
  info_vec <- unlist(strsplit(gene_info, ";"))
  
  gene_type_ext <- grep("gene_type", info_vec, value = TRUE)
  gene_type <- gsub("( |gene_type|\")", "", gene_type_ext) 
  
  return(gene_type)
}

get_ensembl_id <- function(gene) {
  return(gene_ann[gene_name == gene, gene_id])
}

##
# Functions to extract gene info from biomart lookup table
##

get_ensembl_from_accn <- function(accn_id) {
  id <- biomart_mapping_ensembl_113[accession %in% accn_id, unique(ensembl_gene_id)]
  id <- id[!is.na(id) & id != ""]
  
  if (length(id) == 0) {return(NULL)}
  
  return(id)
}

get_accn_from_symbol <- function(gene_name) {
  id <- biomart_mapping_ensembl_113[hgnc_symbol %in% gene_name, unique(accession)]
  id <- id[!is.na(id) & id != ""]
  
  if (length(id) == 0) {return(NULL)}
  
  return(id)
}

get_entrez_from_ensembl <- function(gene_id) {
  id <- biomart_mapping_ensembl_113[ensembl_gene_id %in% gene_id, unique(entrezgene_id)]
  id <- id[!is.na(id) & id != ""]
  
  if (length(id) == 0) {return(NULL)}
  
  return(id)
}

##
# Functions to extract sample info
##

fetch_sample_info <- function(sample, var) {
  return(metadata %>% filter(sample_name == sample) %>% pull(var) %>% as.character())
}

##
# Function to access SynGO data
##

get_syngo_genes <- function(syngo_df, ID, domain_type = "BP") {
  genes <- syngo_df %>%
    filter(domain == domain_type & id == ID) %>% 
    pull("hgnc_symbol") %>%
    as.character() %>% 
    str_split_1(pattern = ", ") %>% 
    unique()
  return(genes)
}

add_named_module_scores <- function(seurat_obj, gene_sets, ctrl = 100, seed) {
  for (nm in names(gene_sets)) {
    gene_set_list <- gene_sets[[nm]]
    gene_set_list <- gene_set_list[gene_set_list %in% rownames(seurat_obj)]
    
    if(length(gene_set_list) > 0) {
      seurat_obj <- AddModuleScore(seurat_obj,
                                   features = gene_set_list,
                                   name     = paste0(nm, "_"),
                                   ctrl     = ctrl,
                                   seed     = seed)
      
      colnames(seurat_obj@meta.data)[
        colnames(seurat_obj@meta.data) == paste0(nm, "_1")
      ] <- nm
      
    }
  }
  return(seurat_obj)
}
