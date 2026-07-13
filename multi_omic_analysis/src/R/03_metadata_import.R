# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

message("######################################")
message("03_metadata_import.R")
message("######################################")

### Importing metadata ###

metadata <- fread(paste(dir_list$data, "sample_metadata_v6.csv", sep = "/"), data.table = FALSE)

# Setting metadata levels
metadata <- metadata %>% 
  mutate(cell_line = factor(cell_line, levels = c("DPSC", "iN")),
         donor_id = factor(donor_id, levels = c("D12", "D14", "D16", "D19")),
         sample_id = factor(sample_id, levels = c("DPSC1", "DPSC2", "DPSC3", "DPSC4", "DPSC-iN1", "DPSC-iN2", "DPSC-iN3", "DPSC-iN4")))

# Creating a sample name map to replace the sample_names with the shorter sample_ids
sample_name_map <- setNames(metadata$sample_id, metadata$sample_name)
message("Metadata imported.")
saveRDS(metadata, paste(dir_list$rds, "metadata.rds", sep = "/"))