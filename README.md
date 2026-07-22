# Rapid direct neuronal reprogramming of human dental pulp stem cells reveals heterogeneous neuronal states

This repository contains the code and supporting data files for the preprint:

> **Rapid direct neuronal reprogramming of human dental pulp stem cells reveals heterogeneous neuronal states**
> Abbas, Muralidharan, et al. (2026)

The repository provides:

* Nextflow pipelines for preprocessing bulk RNA-seq and single-nucleus RNA-seq data
* R scripts for multi-omic data analysis
* Code for reproducing the analyses, tables, and plots presented in the manuscript

> **Data availability:** GEO accession numbers and selected supplementary files will be added when the manuscript is published.

---

## Table of Contents

* [Repository Structure](#repository-structure)
* [Workflow Overview](#workflow-overview)
* [Preprocessing Pipelines](#preprocessing-pipelines)

  * [Requirements](#preprocessing-requirements)
  * [Running the Pipelines](#running-the-preprocessing-pipelines)
  * [Pipeline Outputs](#preprocessing-pipeline-outputs)
* [Multi-Omic Data Analysis](#multi-omic-data-analysis)

  * [Required Input Files](#required-input-files)
  * [Software Requirements](#multi-omic-software-requirements)
  * [Running the Analysis](#running-the-multi-omic-analysis)
  * [Analysis Outputs](#multi-omic-analysis-outputs)
* [Session Information](#session-information)

---

## Repository Structure

```text
.
├── pre_processing/
│   ├── bulkRNAseq_nf_pipeline/
│   │   ├── main.nf
│   │   ├── nextflow.config
│   │   ├── params.json
│   │   └── modules/
│   │
│   └── snRNAseq_nf_pipeline/
│       ├── main.nf
│       ├── nextflow.config
│       ├── params.json
│       ├── modules/
│       ├── src/
│       └── bin/
│
├── multi_omic_analysis/
│   ├── Main.Rmd
│   ├── data/
│   ├── rds/
│   └── src/
│       └── R/
│           ├── helper_functions.R
│           ├── 01_setup.R
│           ├── 02_biomart_lookup.R
│           ├── 03_metadata_import.R
│           ├── 04_bulkRNAseq_analysis.R
│           ├── 05_DNAm_analysis.R
│           ├── 06_msProt_analysis.R
│           ├── 07_snRNAseq_import_data.R
│           ├── 08_snRNAseq_quality_control.R
│           ├── 09_snRNAseq_norm_scaling_pca_harmonyInt.R
│           ├── 10_snRNAseq_clustering.R
│           ├── 11_snRNAseq_doublet_filt.R
│           ├── 12_snRNAseq_clust_characterisation.R
│           ├── 13_snRNAseq_auCell.R
│           └── 14_snRNAseq_trajectory_inf.R
│
├── Dockerfile
├── docker-compose.yaml
├── renv.lock
└── README.md
```

### Key Files and Directories

| Path                           | Description                                                                  |
| ------------------------------ | ---------------------------------------------------------------------------- |
| `pre_processing/`              | Nextflow pipelines for preprocessing raw sequencing data                     |
| `bulkRNAseq_nf_pipeline/`      | Bulk RNA-seq preprocessing pipeline                                          |
| `snRNAseq_nf_pipeline/`        | Single-nucleus RNA-seq preprocessing pipeline                                |
| `multi_omic_analysis/Main.Rmd` | Main analysis notebook that runs the downstream analysis and generates plots |
| `multi_omic_analysis/src/R/`   | R scripts for individual analysis stages                                     |
| `multi_omic_analysis/data/`    | Input data for downstream analysis                                           |
| `multi_omic_analysis/rds/`     | Saved R objects generated during analysis                                    |
| `docker-compose.yaml`          | Docker Compose configuration for the RStudio environment                     |
| `renv.lock`                    | R package environment lock file                                              |

---

## Workflow Overview

The repository contains two major workflow components.

### 1. Raw-data preprocessing

Nextflow pipelines are used to preprocess:

* Bulk RNA-seq data
* Single-nucleus RNA-seq data

The two pipelines are independent and should preferably be run in a high-performance computing environment.

### 2. Multi-omic downstream analysis

The downstream analysis integrates:

* Bulk RNA-seq
* Single-nucleus RNA-seq
* DNA methylation
* Mass-spectrometry proteomics

The analysis is run from `multi_omic_analysis/Main.Rmd` using an RStudio Docker environment.

---

# Preprocessing Pipelines

The preprocessing pipelines are located in:

```text
pre_processing/
├── bulkRNAseq_nf_pipeline/
└── snRNAseq_nf_pipeline/
```

Each pipeline is run independently.

## Preprocessing Requirements

### General requirements

1. **FASTQ files**

   Download the FASTQ files from the GEO repository.

   The GEO accession number will be provided when the manuscript is published.

   Store the datasets in separate directories:

   ```text
   fastq/
   ├── bulkRNAseq/
   └── snRNAseq/
   ```

2. **Nextflow**

   Nextflow version 23 or later is required.

3. **Software dependencies**

   Install the software versions specified in each pipeline's `nextflow.config` file.

   As an alternative, and the recommended approach, download the corresponding Singularity or Apptainer containers from:

   https://depot.galaxyproject.org/singularity/

4. **Storage**

   At least **200 GB of available storage** is recommended.

5. **Compute environment**

   An HPC environment is strongly recommended.

---

### snRNA-seq-specific requirements

The single-nucleus RNA-seq pipeline requires:

* 10x Genomics GRCh38 2024 reference index

---

### Bulk RNA-seq-specific requirements

The bulk RNA-seq pipeline requires:

* GENCODE v47 human STAR index
* GENCODE v47 human transcriptome GTF file
* GENCODE v47 human transcriptome BED file

---

## Running the Preprocessing Pipelines

### 1. Clone the repository

```bash
git clone https://github.com/chanmur/DPSC_iN.git
cd DPSC_iN
```
---

### 2. Configure the pipeline parameters

Open the relevant `params.json` file:

```text
pre_processing/bulkRNAseq_nf_pipeline/params.json
```

or:

```text
pre_processing/snRNAseq_nf_pipeline/params.json
```

Update only the file and directory paths.

Use absolute paths for all required files and directories. Placeholders are included in each file to indicate the expected values.

The `out_dir` parameter specifies where the pipeline results will be saved.

---

### 3. Configure the execution environment

Open the corresponding `nextflow.config` file and update the settings for your computing environment.

The configuration file includes:

* Resource requirements
* Executor settings
* Software or container requirements
* Process-specific resource definitions

Adjust these values as needed for your HPC environment.

---

### 4. Run the pipeline

Navigate to the relevant pipeline directory.

For bulk RNA-seq:

```bash
cd pre_processing/bulkRNAseq_nf_pipeline
```

For single-nucleus RNA-seq:

```bash
cd pre_processing/snRNAseq_nf_pipeline
```

Run the pipeline:

```bash
nextflow run main.nf \
  -c nextflow.config \
  -params-file params.json
```

---

### Resume an interrupted pipeline

To resume a previously started pipeline, add the `-resume` option:

```bash
nextflow run main.nf \
  -resume \
  -c nextflow.config \
  -params-file params.json
```

Nextflow will reuse completed tasks from the previous run when possible.

---

## Preprocessing Pipeline Outputs

### Bulk RNA-seq pipeline

The bulk RNA-seq pipeline generates the following directories:

#### `01_qc/`

Contains quality-control reports and results from:

* FastQC
* MultiQC
* RSeQC

#### `02_mapping/`

Contains one directory for each sample.

Each sample directory includes:

* Aligned, sorted, and indexed BAM file
* Normalized BigWig file

#### `03_gene_counts/`

Contains a CSV file with raw gene counts across all samples.

---

### Single-nucleus RNA-seq pipeline

The single-nucleus RNA-seq pipeline generates:

#### `01_count_output/`

Contains the Cell Ranger output directory for each processed sample.

---

# Multi-Omic Data Analysis

The downstream analysis code is located in:

```text
multi_omic_analysis/
```

The main analysis notebook is:

```text
multi_omic_analysis/Main.Rmd
```

This notebook runs the analysis scripts in `multi_omic_analysis/src/R/` and generates the analysis results and figures.

---

## Required Input Files

Add the following files to `multi_omic_analysis/data/`.

### Bulk RNA-seq

```text
DPSC_iN_bulkRNAseq_raw_gene_counts_batch_merged.csv
```

This file will be available from the GEO repository when the manuscript is published.

---

### Proteomics

```text
Supplementary_Data_3.csv
```

This preprocessed proteomics dataset will be available from the supplementary materials associated with the published article.

---

### DNA methylation

Add the DNA methylation IDAT files from GEO to:

```text
multi_omic_analysis/data/dnam_idat_files/
```

---

### Single-nucleus RNA-seq

Add the single-nucleus RNA-seq matrix files from GEO to:

```text
multi_omic_analysis/data/snRNA_mat_files/
```

---

### Reference and annotation files

The following reference files are also required:

* EPICv2 hg38 GENCODE v41 manifest TSV file
* MANE GRCh38 v14 GTF file
* GENCODE v47 human transcriptome GTF file

Update the relevant paths in the analysis configuration or setup scripts after adding these files.

---

## Multi-Omic Software Requirements

The downstream analysis requires:

* R version 4.3.1 or later
* `renv` version 1.0.3 or later
* Docker
* Docker Compose
* A modern web browser
* RStudio, strongly recommended

The session information below records an analysis environment using R 4.5.1 and `renv` 1.1.5.

---

## Running the Multi-Omic Analysis

### 1. Clone the repository

```bash
git clone <repository-url>
cd <repository-directory>
```

---

### 2. Configure the Docker volume

Open:

```text
docker-compose.yaml
```

Under `volumes`, change the path on the left side of the colon to the absolute path of the cloned repository.

For example:

```yaml
volumes:
  - /absolute/path/to/cloned/repository:/home/rstudio/project
```

Use the existing container-side path from `docker-compose.yaml` if it differs from this example.

---

### 3. Start the Docker environment

From the repository root, run:

```bash
docker compose up -d
```

---

### 4. Open RStudio

Open the RStudio container in a web browser.

Use the following credentials:

```text
Username: rstudio
Password: dpsc
```

---

### 5. Create the R project

In RStudio, create a new R project using:

```text
multi_omic_analysis/
```

as the project directory.

---

### 6. Activate the project environment

In the R console, run:

```r
renv::activate()
```

Restart the R session after activation.

---

### 7. Restore the package environment

Run:

```r
renv::restore()
```

Restart the R session again after the package restoration is complete.

---

### 8. Run the analysis

Open:

```text
multi_omic_analysis/Main.Rmd
```

Run all code chunks in order.

---

## Multi-Omic Analysis Outputs

The analysis generates the following directories.

### `docs/`

Contains CSV files with analysis results, including:

* Differentially expressed genes
* Gene Ontology enrichment results
* Other intermediate and final result tables

### `plots/`

Contains all plots generated by the analysis.

Plot filenames correspond to their placement in the manuscript figures.

### `rds/`

Contains RDS files with R objects saved at different stages of the analysis.

These files can be used to:

* Resume analysis from intermediate stages
* Inspect processed objects
* Perform additional exploratory analyses

### `h5ad/`

Contains the final processed single-nucleus RNA-seq dataset in AnnData format:

```text
*.h5ad
```

---

# Session Information

The following session information records the environment used for the multi-omic analysis.

<details>
<summary><strong>View complete R session information</strong></summary>

```text
R version 4.5.1 (2025-06-13)
Platform: aarch64-unknown-linux-gnu
Running under: Ubuntu 24.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/aarch64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/aarch64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so
LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8
 [2] LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8
 [4] LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8
 [6] LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8
 [8] LC_NAME=C
 [9] LC_ADDRESS=C
[10] LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8
[12] LC_IDENTIFICATION=C

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4 graphics grDevices utils datasets methods base

other attached packages:
 [1] sesameData_1.28.0
 [2] ExperimentHub_3.0.0
 [3] AnnotationHub_4.0.0
 [4] BiocFileCache_3.0.0
 [5] dbplyr_2.5.2
 [6] org.Hs.eg.db_3.22.0
 [7] AnnotationDbi_1.72.0
 [8] slingshot_2.18.0
 [9] TrajectoryUtils_1.18.0
[10] SingleCellExperiment_1.32.0
[11] princurve_2.1.6
[12] Seurat_5.5.0
[13] SeuratObject_5.4.0
[14] sp_2.2-1
[15] DESeq2_1.50.2
[16] SummarizedExperiment_1.40.0
[17] Biobase_2.70.0
[18] MatrixGenerics_1.22.0
[19] matrixStats_1.5.0
[20] GenomicRanges_1.62.1
[21] Seqinfo_1.0.0
[22] IRanges_2.44.0
[23] S4Vectors_0.48.1
[24] BiocGenerics_0.56.0
[25] generics_0.1.4
[26] patchwork_1.3.2
[27] ggplot2_4.0.3
[28] stringr_1.6.0
[29] tibble_3.3.1
[30] tidyr_1.3.2
[31] dplyr_1.2.1
[32] data.table_1.18.4

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.23
  [2] splines_4.5.1
  [3] later_1.4.8
  [4] filelock_1.0.3
  [5] polyclip_1.10-7
  [6] preprocessCore_1.72.0
  [7] fastDummies_1.7.6
  [8] lifecycle_1.0.5
  [9] httr2_1.2.2
 [10] globals_0.19.1
 [11] lattice_0.22-7
 [12] MASS_7.3-65
 [13] magrittr_2.0.5
 [14] plotly_4.12.0
 [15] rmarkdown_2.31
 [16] yaml_2.3.12
 [17] httpuv_1.6.17
 [18] otel_0.2.0
 [19] sctransform_0.4.3
 [20] spam_2.11-3
 [21] spatstat.sparse_3.1-0
 [22] reticulate_1.46.0
 [23] cowplot_1.2.0
 [24] pbapply_1.7-4
 [25] DBI_1.3.0
 [26] RColorBrewer_1.1-3
 [27] abind_1.4-8
 [28] Rtsne_0.17
 [29] purrr_1.2.2
 [30] rappdirs_0.3.4
 [31] ggrepel_0.9.8
 [32] irlba_2.3.7
 [33] listenv_0.10.1
 [34] spatstat.utils_3.2-3
 [35] goftest_1.2-3
 [36] RSpectra_0.16-2
 [37] spatstat.random_3.4-5
 [38] fitdistrplus_1.2-6
 [39] parallelly_1.47.0
 [40] svglite_2.2.2
 [41] codetools_0.2-20
 [42] DelayedArray_0.36.1
 [43] tidyselect_1.2.1
 [44] farver_2.1.2
 [45] spatstat.explore_3.8-0
 [46] jsonlite_2.0.0
 [47] wheatmap_0.2.0
 [48] progressr_0.19.0
 [49] ggridges_0.5.7
 [50] survival_3.8-3
 [51] systemfonts_1.3.2
 [52] tools_4.5.1
 [53] ragg_1.5.2
 [54] ica_1.0-3
 [55] Rcpp_1.1.1-1.1
 [56] glue_1.8.1
 [57] gridExtra_2.3
 [58] SparseArray_1.10.10
 [59] xfun_0.57
 [60] withr_3.0.2
 [61] BiocManager_1.30.27
 [62] fastmap_1.2.0
 [63] digest_0.6.39
 [64] R6_2.6.1
 [65] mime_0.13
 [66] textshaping_1.0.5
 [67] colorspace_2.1-2
 [68] scattermore_1.2
 [69] tensor_1.5.1
 [70] dichromat_2.0-0.1
 [71] spatstat.data_3.1-9
 [72] RSQLite_3.52.0
 [73] renv_1.1.5
 [74] httr_1.4.8
 [75] htmlwidgets_1.6.4
 [76] S4Arrays_1.10.1
 [77] uwot_0.2.4
 [78] pkgconfig_2.0.3
 [79] gtable_0.3.6
 [80] blob_1.3.0
 [81] lmtest_0.9-40
 [82] S7_0.2.2
 [83] XVector_0.50.0
 [84] htmltools_0.5.9
 [85] dotCall64_1.2
 [86] scales_1.4.0
 [87] png_0.1-9
 [88] spatstat.univar_3.1-7
 [89] knitr_1.51
 [90] rstudioapi_0.18.0
 [91] tzdb_0.5.0
 [92] reshape2_1.4.5
 [93] nlme_3.1-168
 [94] curl_7.1.0
 [95] cachem_1.1.0
 [96] zoo_1.8-15
 [97] BiocVersion_3.22.0
 [98] KernSmooth_2.23-26
 [99] vipor_0.4.7
[100] parallel_4.5.1
[101] miniUI_0.1.2
[102] ggrastr_1.0.2
[103] pillar_1.11.1
[104] grid_4.5.1
[105] vctrs_0.7.3
[106] RANN_2.6.2
[107] promises_1.5.0
[108] xtable_1.8-8
[109] cluster_2.1.8.2
[110] beeswarm_0.4.0
[111] evaluate_1.0.5
[112] readr_2.2.0
[113] cli_3.6.6
[114] locfit_1.5-9.12
[115] compiler_4.5.1
[116] rlang_1.2.0
[117] crayon_1.5.3
[118] future.apply_1.20.2
[119] labeling_0.4.3
[120] ggbeeswarm_0.7.3
[121] plyr_1.8.9
[122] stringi_1.8.7
[123] viridisLite_0.4.3
[124] deldir_2.0-4
[125] BiocParallel_1.44.0
[126] Biostrings_2.78.0
[127] lazyeval_0.2.3
[128] spatstat.geom_3.7-3
[129] Matrix_1.7-3
[130] RcppHNSW_0.6.0
[131] sesame_1.28.1
[132] hms_1.1.4
[133] bit64_4.8.0
[134] future_1.70.0
[135] conflicted_1.2.0
[136] KEGGREST_1.50.0
[137] shiny_1.13.0
[138] ROCR_1.0-12
[139] igraph_2.3.1
[140] memoise_2.0.1
[141] bit_4.6.0
```

</details>
