#!/usr/bin/env Rscript
##### For SpatialLIBD Benchmarking #####

args = commandArgs(trailingOnly=TRUE)

library(Matrix,lib.loc="/users/hzhou1/R/4.3")
#library(atlasClustering)
#library(spatialLIBD)
library(Seurat)
library(magrittr)
library(tidyverse)
library(peakRAM)
#library(PRECAST)
#library(BASS)
#library(BayesSpace)
library(harmony)
library(aricode)
options(future.globals.maxSize = 4000 * 1024^2)

source("/users/hzhou1/benchmark/benchmark_helper_func.R")


run_spatialMNN = F
run_spatialMNN_par = F
run_Seurat = F
run_BayesSpace = F
run_BASS = F
run_PRECAST = F
run_SLAT = F
run_MENDER = F
run_BANKSY = F

batch_no = args[1]
method_ls <- str_split(args[2],",") %>% unlist()
save_path = args[3]

if("spatialMNN" %in% method_ls) run_spatialMNN = T
if("spatialMNN_par" %in% method_ls) run_spatialMNN_par = T
if("Seurat" %in% method_ls) run_Seurat = T
if("BayesSpace" %in% method_ls) run_BayesSpace = T
if("BASS" %in% method_ls) run_BASS = T
if("PRECAST" %in% method_ls) run_PRECAST = T
if("SLAT" %in% method_ls) run_SLAT = T
if("MENDER" %in% method_ls) run_MENDER = T
if("BANKSY" %in% method_ls) run_BANKSY = T

message(paste("Running:",args[2]))

## Bench result

bench_res <- data.frame(method = c("spatialMNN","spatialMNN_par","Seurat","BayesSpace",
                                   "BASS","PRECAST","SLAT","MENDER","BANKSY"),
                        dataset = "spatialLIBD",
                        tic = 0,
                        toc = 0,
                        time = 0,
                        mem = -1,
                        ari = -1,
                        ari_vec = "-",
                        nmi = -1,
                        nmi_vec = "-")

fig_ls <- list()

##### Load Data #####
load("/users/hzhou1/benchmark/datasets/starmap_mpfc.RData")

seurat_ls <- list()

for(i in names(starmap_cnts)){
  seurat_ls[[i]] <- CreateSeuratObject(counts = starmap_cnts[[i]],
                                      meta.data = starmap_info[[i]],project = i)
  seurat_ls[[i]]$coord_x <- seurat_ls[[i]]$x
  seurat_ls[[i]]$coord_y <- seurat_ls[[i]]$y
  seurat_ls[[i]]$row <- seurat_ls[[i]]$x
  seurat_ls[[i]]$col <- seurat_ls[[i]]$y
}


# convert to seurat list
if(run_spatialMNN | run_spatialMNN_par | run_Seurat | run_BASS | run_PRECAST){
  library(atlasClustering)
  gc()
}

if(run_BayesSpace | run_BANKSY){
    source("/users/hzhou1/benchmark/benchmarking_Rwrapper_func.R")
    library(SpatialExperiment)
    library(SingleCellExperiment)

    sce_ls <- lapply(seurat_ls, function(x) {as.SingleCellExperiment(x)})

    combined_sce <- do.call(cbind, sce_ls)

    # Convert the combined SCE object to a SpatialExperiment object
    spe <- SpatialExperiment(
      assays = list(counts = assay(combined_sce, "counts")),
      rowData = rowData(combined_sce),
      colData = colData(combined_sce),
      reducedDims = reducedDims(combined_sce),
      spatialCoords = matrix(c(combined_sce@colData@listData[["x"]],
                              combined_sce@colData@listData[["y"]]),
                            ncol = 2, 
                            dimnames = list(combined_sce@colData@rownames,c("x","y"))) # Adjust this if you have spatial coordinates
    )
}

# load h5ad
#library(Matrix, lib.loc = "/jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library")
#spe <- readRDS("~/benchmark/datasets/fig2_STARmap_spe.RDS")
#use_python("/users/hzhou1/.conda/envs/MENDER/bin")
#create_annFiles(spe, "z", "ident", "~/benchmark/datasets/")

if(run_SLAT | run_MENDER){
    source("/users/hzhou1/benchmark/benchmarking_Rwrapper_func.R")
    samplePaths <- c("/users/hzhou1/benchmark/datasets/sample_20180417_BZ5_control.h5ad",
                    "/users/hzhou1/benchmark/datasets/sample_20180419_BZ9_control.h5ad",
                    "/users/hzhou1/benchmark/datasets/sample_20180424_BZ14_control.h5ad")

    sampleNames <- c("sample_20180417_BZ5_control", 
                     "sample_20180419_BZ9_control", 
                     "sample_20180424_BZ14_control")
    if(run_SLAT){
      ref_sample = as.numeric(args[4])
    }

}



if(run_spatialMNN){
    mem_usage <- peakRAM(
      {
        tic <- Sys.time()
        seurat_ls <- stage_1(seurat_ls, cor_threshold = 0.7, nn = 6, nn_2 = 20, cl_resolution = 10,
                             top_pcs = 6, cl_min = 5, find_HVG = F, cor_met = "PC",
                             edge_smoothing = T, use_glmpca = T, verbose = T, num_core = 1)
        
        rtn_ls <- stage_2(seurat_ls, cl_key = "merged_cluster",
                          rtn_seurat = T, nn_2 = 10, method = "MNN",
                          top_pcs = 10, use_glmpca = T, rare_ct = "m", resolution = 0.7)
        
        seurat_ls <- assign_label(seurat_ls, rtn_ls$cl_df, "MNN", 0.7, cl_key = "merged_cluster")
        toc <- Sys.time()
      }
    )

    bench_res <- rkd_res(bench_res,tic,toc,mem_usage,"spatialMNN",seurat_ls,"z","sec_cluster_MNN")

    for(i in names(seurat_ls)){
      fig_ls[[paste0("SM_",i)]] <- draw_slide_bench(seurat_ls[[i]]$coord_x,
                                                    seurat_ls[[i]]$coord_y,
                                                    seurat_ls[[i]]$sec_cluster_MNN,"atlasClustering",i,flip = T)
    }
}

if(run_spatialMNN_par){
    mem_usage <- peakRAM(
      {
        tic <- Sys.time()
        seurat_ls <- stage_1(seurat_ls, cor_threshold = 0.7, nn = 6, nn_2 = 20, cl_resolution = 10,
                             top_pcs = 6, cl_min=5, find_HVG = F, cor_met = "PC",
                             edge_smoothing = T, use_glmpca = T, verbose = T, num_core = 8)
        
        rtn_ls <- stage_2(seurat_ls, cl_key = "merged_cluster",
                          rtn_seurat = T, nn_2 = 10, method = "MNN",
                          top_pcs = 10, use_glmpca = T, rare_ct = "m", resolution = 0.7)
        
        seurat_ls <- assign_label(seurat_ls, rtn_ls$cl_df, "MNN", 0.7, cl_key = "merged_cluster")
        toc <- Sys.time()
      }
    )

    bench_res <- rkd_res(bench_res,tic,toc,mem_usage,"spatialMNN_par",seurat_ls,"z","sec_cluster_MNN")

    for(i in names(seurat_ls)){
      fig_ls[[paste0("SP_",i)]] <- draw_slide_bench(seurat_ls[[i]]$coord_x,
                                                    seurat_ls[[i]]$coord_y,
                                                    seurat_ls[[i]]$sec_cluster_MNN,"atlasClustering",i,flip = T)
    }
}

if(run_Seurat){
    mem_usage <- peakRAM(
      {
         tic <- Sys.time()
         seurat_ls <- lapply(X = seurat_ls, FUN = function(x) {
           x <- NormalizeData(x)
           x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 160)
           })
         features <- SelectIntegrationFeatures(object.list = seurat_ls)
         seurat_ls <- lapply(X = seurat_ls, FUN = function(x) {
           x <- ScaleData(x, features = features, verbose = FALSE)
           x <- RunPCA(x, features = features, verbose = FALSE)
           })
         anchors <- FindIntegrationAnchors(object.list = seurat_ls, anchor.features = features, reduction = "rpca", verbose = F)
         seu_combined <- IntegrateData(anchorset = anchors, verbose = F)
         DefaultAssay(seu_combined) <- "integrated"

             # Run the standard workflow for visualization and clustering
         seu_combined <- ScaleData(seu_combined, verbose = FALSE)
         seu_combined <- RunPCA(seu_combined, npcs = 10, verbose = FALSE)
         #seu_combined <- RunUMAP(seu_combined, reduction = "pca", dims = 1:10)
         seu_combined <- FindNeighbors(seu_combined, reduction = "pca", dims = 1:6)
         seu_combined <- FindClusters(seu_combined, resolution = 0.2)
         toc <- Sys.time()
      }
    )

    bench_res <- rkd_res(bench_res,tic,toc,mem_usage,"Seurat",label_vec = seu_combined@meta.data[["z"]], 
                         result_vec = seu_combined@meta.data[["seurat_clusters"]],
                         sample_vec = lapply(names(seurat_ls),function(i)rep(i,ncol(seurat_ls[[i]])))%>%unlist)

    for(i in names(seurat_ls)){
      idx = seu_combined@meta.data[["orig.ident"]] == i
      fig_ls[[paste0("SR_",i)]] <- 
        draw_slide_bench(seu_combined@meta.data[["x"]][idx],
                         seu_combined@meta.data[["y"]][idx],
                         seu_combined@meta.data[["seurat_clusters"]][idx],"Seurat",i,flip = F)
    }
    rm(seu_combined)
}

if(run_BASS){
    library(BASS)
    gc()
    mem_usage <- peakRAM(
      {
        tic <- Sys.time()
        set.seed(0)
        # Set up BASS object
        BASS <- createBASSObject(lapply(seurat_ls, 
                                    function(seu_obj){
                                      seu_obj@assays[["RNA"]]@layers[["counts"]] %>% 
                                        `colnames<-`(colnames(seu_obj)) %>% `row.names<-`(row.names(seu_obj))
                                      }), 
                             lapply(seurat_ls, 
                                    function(seu_obj){
                                      data.frame(x=seu_obj$x,
                                                 y=seu_obj$y,
                                                 row.names=colnames(seu_obj))
                                      }), 
                             C = 15, R = 4,
      beta_method = "SW")
    
    BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE, doPCA = TRUE, scaleFeature = T, nPC = 20)

        # Run BASS algorithm
        BASS <- BASS.run(BASS)
        BASS <- BASS.postprocess(BASS)
        toc <- Sys.time()
      }
    )

    bench_res <- rkd_res(bench_res,tic,toc,mem_usage,"BASS",
                     label_vec  = lapply(names(seurat_ls),function(i){seurat_ls[[i]]@meta.data$z}) %>% unlist, 
                     result_vec = unlist(BASS@results$z),
                     sample_vec = lapply(names(seurat_ls),function(i)rep(i,ncol(seurat_ls[[i]]))) %>% unlist)

    for(i in names(seurat_ls)){
      fig_ls[[paste0("BA_",i)]] <- draw_slide_bench(seurat_ls[[i]]$x,
                                                    seurat_ls[[i]]$y,
                                                    as.factor(BASS@results$z[[which(names(seurat_ls) == i)]]),
                                                    "BASS",i,flip = F)
    }
    rm(BASS)
}

if(run_PRECAST){
    library(PRECAST)
    gc()

    mem_usage <- peakRAM(
        {
            tic <- Sys.time()
            seuInt <- suppressMessages(PRECAST_test(seurat_ls, k = 4, gene_num = 160))
            toc <- Sys.time()
        }
    )

    layer_vec <- sapply(seq_along(seurat_ls),function(i){
      tmp_df <- data.frame(barcode = row.names(seuInt@meta.data) %>% 
                            str_sub(start = 1, end = 18) %>% .[seuInt@meta.data[["batch"]]==i])
      tmp_df <- left_join(tmp_df,
                          data.frame(barcode = colnames(seurat_ls[[names(seurat_ls)[i]]]),
                                    z = seurat_ls[[names(seurat_ls)[i]]]@meta.data[,"z"]), 
                          by="barcode")
      tmp_df$z
      
    }) %>% unlist()

    bench_res <- rkd_res(bench_res,tic,toc,mem_usage,"PRECAST",
                         label_vec = layer_vec, 
                         result_vec = seuInt@meta.data[["cluster"]], 
                         sample_vec = seuInt@meta.data[["batch"]])

    for(i in names(seurat_ls)){
      idx = seuInt@meta.data[["batch"]] == which(names(seurat_ls)==i)
      fig_ls[[paste0("PC_",i)]] <- draw_slide_bench(seuInt@reductions[["position"]]@cell.embeddings[idx,2],
                                                    seuInt@reductions[["position"]]@cell.embeddings[idx,1],
                                                    seuInt@meta.data[["cluster"]][idx],"PRECAST",i,flip = T)
    }
    rm(seuInt)
}

if(run_BayesSpace){
    library(BayesSpace)
    gc()
    sce <- SingleCellExperiment(assays=list(counts=spe@assays@data@listData[["counts"]]),
                            rowData=NULL,
                            colData=spe@colData)

    sample_names <- unique(sce$ident)

    sce$row <- round(sce$x*10)
    sce$col <- round(sce$y*10)


    for(j in seq_along(sample_names)){
      sce$row[sce$ident == sample_names[j]] <- sce$row[sce$ident == sample_names[j]] + 150000*((j-1)%%3)
      sce$col[sce$ident == sample_names[j]] <- sce$col[sce$ident == sample_names[j]] + 70000*floor((j-1)/3)
    }

    mem_usage <- peakRAM(
      {
        tic <- Sys.time()
        sce <- spatialPreprocess(sce, platform="ST", n.PCs = 6)
        sce = RunHarmony(sce, "ident", verbose = F)
        #sce <- qTune(sce, qs=seq(2, 10), platform="ST", d=7)
        #qPlot(sce)

        sce <- spatialCluster(sce, use.dimred = "HARMONY", q = 4, nrep = 10000)
        toc <- Sys.time()
      }
    )


    bench_res <- rkd_res(bench_res,tic,toc,mem_usage,"BayesSpace",
                         label_vec  = sce@colData@listData[["z"]], 
                         result_vec = sce@colData@listData[["spatial.cluster"]],
                         sample_vec = sce@colData@listData[["ident"]])

    for(i in sample_names){
      idx = sce@colData@listData[["ident"]] == i
      fig_ls[[paste0("BS_",i)]] <- draw_slide_bench(sce@colData@listData[["row"]][idx],
                                                    sce@colData@listData[["col"]][idx],
                                                    sce@colData@listData[["spatial.cluster"]][idx],"BayesSpace",i,flip = F)
    }
    rm(sce)
}

if(run_BANKSY){
    require(Banksy)
    annots_label = "z"
    sample_label = "ident"
    batch = TRUE
    sample_info = data.frame(sample_id = spe$ident, 
                             subject = spe$sample_id) %>% 
                             distinct()
    k_geom = 30
    lambda = 0.2
    res = 0.55 
    npcs = 20
    SEED = 1000
    use_agf = TRUE 
    compute_agf = TRUE
    print(paste("Preprocessing started.", Sys.time()))
    if (!is.null(annots_label)) {
        # replacing 'NA' annotated labels with 'Unknown'
        annots = as.data.frame(as.character(spe[[annots_label]]))
        annots[, 1] = replace_na(annots[, 1], "Unknown")
        spe[[annots_label]] = factor(as.character(annots[, 1]))
        print(paste("Annotated data checked - NA replaced with Unknown.", Sys.time()))
    }
    # Trimming dataset
    imgData(spe) <- NULL
    assay(spe, "logcounts") <- NULL
    reducedDims(spe) <- NULL
    rowData(spe) <- NULL
    print(paste("Trimmed spe object.", Sys.time()))
    # Grouping samples by source
    spe$subject = factor(as.character(lapply(spe[[sample_label]], function(x) {
      sample_info[sample_info[, 1] == x, 2]})))
    print(paste("Added sample group information.", Sys.time())) 
    print(paste("Multisample run with batch correction.", Sys.time()))
    colnames(spe) <- paste0(colnames(spe), "_", spe[[sample_label]])
    mem_usage <- peakRAM({
        tic <- Sys.time()
        # Staggering spatial coordinates
        locs = spatialCoords(spe)
        locs = cbind(locs, sample = factor(spe[[sample_label]]))
        locs_dt = data.table(locs)
        colnames(locs_dt) <- c("sdimx", "sdimy", "group")
        locs_dt[, sdimx := sdimx - min(sdimx), by = group]
        global_max = max(locs_dt$sdimx) * 1.5
        locs_dt[, sdimx := sdimx + group * global_max]
        locs = as.matrix(locs_dt[, 1:2])
        rownames(locs) <- colnames(spe)
        spatialCoords(spe) <- locs
        print(paste("Spatial coordinates of samples staggered.", Sys.time()))
        # Seurat  
        # Identifying HVGs
        seu = as.Seurat(spe, data = NULL)
        seu = FindVariableFeatures(seu, nfeatures = 2000)
        # Normalizing data
        scale_factor = median(colSums(assay(spe, "counts")))
        seu = NormalizeData(seu, scale.factor = scale_factor,
                            normalization.method = "RC")
        # Adding data to spe object and subsetting to HVGs
        assay(spe, "normcounts") <- GetAssayData(seu)
        spe = spe[VariableFeatures(seu), ]
        print(paste("Seurat feature selection and normalisation complete.", Sys.time()))  
        #### Running BANKSY
        print(paste("BANKSY run started.", Sys.time()))
        spe <- computeBanksy(spe, assay_name = "normcounts",
                            compute_agf = compute_agf,
                            k_geom = k_geom)
        spe <- runBanksyPCA(spe, use_agf = use_agf, lambda = lambda,
                            npcs = npcs, seed = SEED)
        # Harmony batch correction
        PCA_label = paste0("PCA_M", as.numeric(use_agf), "_lam", lambda)
        set.seed(SEED)
        harmony_embedding = RunHarmony(data_mat = reducedDim(spe, PCA_label),
                                        meta_data = colData(spe),
                                        vars_use = c(sample_label),
                                        #vars_use = sample_label,
                                        max.iter = 20,
                                        verbose = FALSE)
        reducedDim(spe, "PCA_harmony") <- harmony_embedding
        print(paste("Batch correction completed.", Sys.time()))
        # Banksy clustering
        spe = clusterBanksy(spe, dimred = "PCA_harmony", use_agf = use_agf,
                            lambda = lambda, resolution = res, seed = SEED)
        print(paste("BANKSY clustering completed.", Sys.time()))
        clust_label = names(colData(spe))[startsWith(names(colData(spe)), "clust")]
        toc <- Sys.time()
    })
    
    bench_res <- rkd_res(bench_res,tic,toc,mem_usage,"BANKSY",
                         label_vec  = spe$z, 
                         result_vec = spe[[clust_label[1]]],
                         sample_vec = spe$ident)

    for(i in unique(spe$ident)){
        idx = spe$ident == i
        fig_ls[[paste0("BA_",i)]] <- draw_slide_bench(spe[["x"]][idx],
                                                      spe[["y"]][idx],
                                                      spe[[clust_label[1]]][idx],"BANKSY",i,flip = F)
    }
    
}

if(run_SLAT){
    gc()
    #python_path = "/users/hzhou1/.conda/envs/scSLAT/bin"
    #pyscript_path = "/users/hzhou1/benchmark"
    res <- runSLAT("/users/hzhou1/.conda/envs/scSLAT/bin", 
                   "/users/hzhou1/benchmark",
                   samplePaths, sampleNames,
                   sampleIDs = 'sampleID',
                   domains = 'domainAnnotations',
                   cos = 0.1,
                   ref_sample = ref_sample)
                   
    bench_res <- rkd_res(bench_res,res[["stats"]]$tic,res[["stats"]]$toc,NULL,"SLAT",
                         label_vec  = sapply(names(res)[1:length(sampleNames)],
                                             function(x){
                                                res[[x]]$obs$domainAnnotations
                                             }) %>% unlist, 
                         result_vec = sapply(names(res)[1:length(sampleNames)],
                                             function(x){
                                                res[[x]]$obs$clust_ref
                                             }) %>% unlist,
                         sample_vec = sapply(names(res)[1:length(sampleNames)],
                                             function(x){
                                                res[[x]]$obs$sampleID
                                             }) %>% unlist, mem = res[["stats"]]$mem)

    for(i in names(res)[1:length(sampleNames)]){
      fig_ls[[paste0("SL_",i)]] <- draw_slide_bench(res[[i]]$obsm$spatial[,2],
                                                    res[[i]]$obsm$spatial[,1],
                                                    res[[i]]$obs$clust_ref,"scSLAT",i,flip = F)
    }
    rm(res)
}

if(run_MENDER){
    gc()
    python_path = "/users/hzhou1/.conda/envs/MENDER/bin"
    pyscript_path = "/users/hzhou1/benchmark"
    #sampleIDs = 'sampleID'
    #domains = 'domainAnnotations'
    #scale = 6; mode = 'ring'; radius = 6
    #seed = 101; batch = 'True'; msm_res = 8
    res <- runMENDER(python_path, pyscript_path,
                     samplePaths, sampleNames,
                     sampleIDs = 'sampleID',
                     domains = 'domainAnnotations',
                     scale = 6, mode = 'radius', radius = 150,
                     seed = 101, batch = 'False', msm_res = -0.5)

    bench_res <- rkd_res(bench_res,res[[2]]$tic,res[[2]]$toc,NULL,"SLAT",
                         label_vec  = res[[1]]$obs$domainAnnotations, 
                         result_vec = res[[1]]$obs$MENDER,
                         sample_vec = res[[1]]$obs$sampleID, 
                         mem = res[[2]]$mem)

    for(i in sampleNames){
        idx = res[[1]]$obs$sampleID == i
        fig_ls[[paste0("ME_",i)]] <- draw_slide_bench(res[[1]]$obsm$spatial[idx,2],
                                                      res[[1]]$obsm$spatial[idx,1],
                                                      res[[1]]$obs$MENDER[idx],"MENDER",i,flip = F)
    }
    rm(res)
    
}

# Save result
saveRDS(fig_ls, file.path(save_path,paste(Sys.Date(),batch_no,"STARmap_bench_fig.rds", sep = "_")))
write_tsv(bench_res[bench_res$ari_vec != "-",],
          file.path(save_path,paste(Sys.Date(),batch_no,"STARmap_bench_res.tsv", sep = "_")))
