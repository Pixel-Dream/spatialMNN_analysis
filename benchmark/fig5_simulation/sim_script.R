#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(Matrix,lib.loc="/users/hzhou1/R/4.3")
library(Seurat)
library(magrittr)
library(tidyverse)
library(peakRAM)
library(harmony)
library(aricode)
options(future.globals.maxSize = 4000 * 1024^2)
library(atlasClustering)


source("/users/hzhou1/benchmark/benchmark_helper_func.R")
source("/users/hzhou1/benchmark/fig5_simulation/Simulation_func.R")
source("/users/hzhou1/benchmark/benchmarking_Rwrapper_func.R")

it = as.numeric(args[1])
sample_num = as.numeric(args[2])
save_path = args[3]

map_mat <- matrix(c(0.9,0.1,0,0,
                    0,0.9,0.05,0.05,
                    0,0.05,0.95,0,
                    0.1,0,0.1,0.8),nrow = 4, byrow = T)

message(paste("It:", it, "Sample Num:", sample_num))
# Generate Datasets
seu_ls <- my_sim(it=it, n_gene=200, noise = 0.2, ig_ratio=0.9, top_pcs = 4, map_mat = map_mat,
                 n_sample = sample_num, cell_max = 7, segmentation = F)

library(SpatialExperiment)
library(SingleCellExperiment)

sce_ls <- lapply(seu_ls, function(x) {as.SingleCellExperiment(x)})

combined_sce <- do.call(cbind, sce_ls)

# Convert the combined SCE object to a SpatialExperiment object
spe <- SpatialExperiment(
  assays = list(counts = assay(combined_sce, "counts")),
  rowData = rowData(combined_sce),
  colData = colData(combined_sce),
  reducedDims = reducedDims(combined_sce),
  spatialCoords = matrix(c(combined_sce@colData@listData[["coord_x"]],
                          combined_sce@colData@listData[["coord_y"]]),
                        ncol = 2, 
                        dimnames = list(combined_sce@colData@rownames,c("x","y"))) # Adjust this if you have spatial coordinates
)

spe@colData@listData$sample_id <- spe$batch


saveRDS(spe, file.path(save_path, paste0("sim_spe_",it,"_",sample_num,".RDS")))
saveRDS(seu_ls, file.path(save_path, paste0("sim_seurat_ls_",it,"_",sample_num,".RDS")))

#use_python("/users/hzhou1/.conda/envs/MENDER/bin")
#create_annFiles(spe, "layer", "sample_id", "/users/hzhou1/benchmark/fig5_simulation/datasets")