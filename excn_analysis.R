## Loading libraries

library(tidyverse)
library(Seurat)
library(ggplot2)
library(here)
library(future)
library(RColorBrewer)
library(scCustomize)
library(pheatmap)
library(ggrepel)
options(future.globals.maxSize = 25000 * 1024 ^ 2)

path <- here()
fig_path <- file.path(str_remove(path,"/xenium_epilepsy_repo"),"figures/transcripts_location")

# read the megadata
megadata = read_rds(file.path(str_remove(path,"xenium_epilepsy_repo"),"REALfull_QC_RCTD_clust_processed.rds"))

# subsetting only excn
excn <- megadata %>% 
  subset(predicted.celltype == "ExcN")

#SCtransfrom individual layers again
excn <- SCTransform(excn, assay = "Xenium")

# Run PCA and integrate using harmony
excn <- RunPCA(excn,npcs = 100)
excn <- IntegrateLayers(object = excn, method = HarmonyIntegration,
                      orig.reduction = "pca", new.reduction = "harmony",
                      verbose = TRUE, normalization.method = "SCT")

# If you want you can check how many PCs are significant but it will only serve as 
# a guideline since you'll run clustering on harmony embedding/dimensions, could still be useful
ElbowPlot(excn,ndims = 70)+geom_hline(yintercept = 2.5)

#function to cluster and run UMAP
Xenium_SeuratLK<- function(x=data, DIMs=1:15, cluster.resolution=0.5, red='harmony'){
  x <- FindNeighbors(x, dims = DIMs, reduction = red)
  x <- FindClusters(x, resolution = cluster.resolution)
  x <- RunUMAP(x, dims = DIMs, reduction = red)
  x
}
#excn_test <- Xenium_SeuratLK(x=excn,DIMs=1:15, cluster.resolution=0.5) 
#excn_test <- Xenium_SeuratLK(x=excn,DIMs=1:15, cluster.resolution=0.4) 
excn_test <- Xenium_SeuratLK(x=excn,DIMs=1:15, cluster.resolution=0.3) #maybe this one?

#try the niche approach to define the different subregions within the hippocampus, doesn't work well
df_niche <- BuildNicheAssay(object = megadata, fov = "L14", group.by = "predicted.celltype", niches.k = 6, neighbors.k = 30)

# plotting the niches,not amazing
ImageDimPlot(df_niche, fov = "L14", group.by = "niches", size = 2, cols = "polychrome", dark.background = F)

# plotting excn for H1 to see how they locate
ImageDimPlot(megadata, fov = "H1", group.by= "predicted.celltype",cells = Cells(excn))

