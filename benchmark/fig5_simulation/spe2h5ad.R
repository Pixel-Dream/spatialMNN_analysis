#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

library(Matrix, lib.loc = "/jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library")

library(magrittr)
library(tidyverse)

source("/users/hzhou1/benchmark/benchmark_helper_func.R")
source("/users/hzhou1/benchmark/fig5_simulation/Simulation_func.R")
source("/users/hzhou1/benchmark/benchmarking_Rwrapper_func.R")

it = 1
save_path = "/users/hzhou1/benchmark/fig 5_simulation/datasets"

for(sample_num in c(2,4,8,16,32,48,64)){
  spe <- readRDS(file.path(save_path, paste0("sim_spe_",it,"_",sample_num,".RDS")))
  use_python("/users/hzhou1/.conda/envs/MENDER/bin")
  dir.create(file.path(save_path, paste0("sim",it,"_",sample_num)), showWarnings = FALSE)
  create_annFiles(spe, "layer", "sample_id", file.path(save_path, paste0("sim",it,"_",sample_num)))
}
