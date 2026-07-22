# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("07_snRNAseq_import_data.R")
message("######################################")

message("Importing snRNAseq raw count matrix files and creating a seurat object.")
# Important cell ranger output files
mat_files <- list.files(paste(dir_list$snRNA_mat), pattern = "\\.tar\\.gz$", full.names = TRUE)

snRNA_data_list <- list()

snRNA_data_list$mat_data <- list()
snRNA_data_list$seurat_objects <- list()

for (file in mat_files) {
  
  extract_path <- gsub("\\.tar\\.gz$", "", file)
  
  # Splitting the file name to extract sample name
  folder_name_split <- str_split_1(tail(str_split_1(extract_path, "/"), 1), "_")
  
  # Sample name from file name
  sample_name <- paste(folder_name_split[3], folder_name_split[5], sep ="_")
  # Retrieving other sample metadata from sample name
  sample_id <- fetch_sample_info(sample_name, "sample_id")
  donor <- fetch_sample_info(sample_name, "donor")
  donor_id <- fetch_sample_info(sample_name, "donor_id")
  age <- fetch_sample_info(sample_name, "age")
  sex <- fetch_sample_info(sample_name, "sex")
  
  # Extracting compressed files
  untar(tarfile = file, exdir = dir_list$snRNA_mat)
  
  # Creating dgc matrix from raw data and adding them to the raw data list 
  snRNA_data_list$mat_data[[sample_id]] <- Read10X(extract_path, strip.suffix = TRUE)
  
  colnames(snRNA_data_list$mat_data[[sample_id]]) <- paste(sample_id, 
                                                           colnames(snRNA_data_list$mat_data[[sample_id]]), 
                                                           sep = "-")

    
  # Creating Seurat Object from the dgc matrices
  snRNA_data_list$seurat_objects[[sample_id]] <- CreateSeuratObject(counts = snRNA_data_list$mat_data[[sample_id]],
                                                                    project = "DPSC_iN_snRNAseq",
                                                                    names.field = 1,
                                                                    names.delim = "-")
  
  # Adding metadata to the Seurat Object
  snRNA_data_list$seurat_objects[[sample_id]]@meta.data  <- snRNA_data_list$seurat_objects[[sample_id]]@meta.data %>%
    mutate("sample_name" = rep(sample_name),
           "sample_id" = rep(sample_id),
           "donor" = rep(donor),
           "donor_id" = rep(donor_id),
           "age" = rep(age),
           "sex" = rep(sex))
}

snRNA_so <- merge(snRNA_data_list$seurat_objects[[1]], y = snRNA_data_list$seurat_objects[c(2:length(snRNA_data_list$seurat_objects))])

message("Files imported and seurat object created.")

saveRDS(snRNA_data_list, file = paste(dir_list$rds, "snRNA_data_list.rds", sep = "/"))
saveRDS(snRNA_so, file = paste(dir_list$rds, "snRNA_so.rds", sep = "/"))

rm(snRNA_data_list, sample_name, sample_id, donor, donor_id, age, sex, mat_files, mat_data, file, seurat_objects, extract_path, folder_name_split)
gc()
