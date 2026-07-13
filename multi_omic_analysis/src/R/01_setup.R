# Author: Chandramouli Muralidharan
# Project: DPSC iN multi-omics analysis
# Date: 2026-07-10

### Additional environment setup ### 

###
# Python Environment
###
env = paste(dir_list$home, "venvs/h5ad_env", sep = "/")
if (!dir.exists(env)) {
  reticulate::install_python(version = "3.12")
  py <- reticulate::virtualenv_starter()[[1]]
  check <- file.exists(py)
  system(paste(py, "-m venv --copies", env))
  reticulate::py_install(packages = c("anndata==0.12.10"), envname = env, method = "virtualenv", pip = TRUE)
  system(paste(file.path(env, "bin/pip"), "install anndata==0.12.10"))
}
Sys.setenv(RETICULATE_PYTHON = file.path(env, "bin/python"))
reticulate::use_python(file.path(env, "bin/python"), required = TRUE)
# reticulate::py_config()
# reticulate::py_module_available("anndata")

###
# Environment Cache
###
ExperimentHub::setExperimentHubOption("ASK", FALSE)
sesameData::sesameDataCache("idatSignature")

###
# Extras
###
# The following packages are listed to be considered by renv. Otherwise, the following does nothing.
ashr::ash
harmony::RunHarmony
glmGamPoi::glm_gp
AUCell::AUCell_calcAUC
presto::pseudobulk_deseq2
slingshot::predict
tradeSeq::fitGAM