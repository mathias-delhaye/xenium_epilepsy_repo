## Using nnSVG to find the spatially variable genes in the xenium data
## Here try using xenium data directly, but it led to all genes being significantly spatially variable

library(SpatialExperiment)
library(SpatialExperimentIO)
library(nnSVG)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(here)
library(future)
library(RColorBrewer)
library(scCustomize)
library(pheatmap)
library(ggrepel)
library(scran)
options(future.globals.maxSize = 25000 * 1024 ^ 2)

path <- here()

## First in one single file

L15 <- read_rds(file.path(str_remove(path,"xenium_epilepsy_repo"),"L15_Seurat_obj.rds"))

## Initially tried with xenium data directly, putting it in a spatial experiment format
L15_spe <- readXeniumSXE("E:/Human epilepsy project/Xenium/2024-12-06/output-XETG00098__0045129__L15__20241206__212646")

dim(L15_spe)

# keep spots over tissue
L15_spe <- L15_spe[, colData(L15_spe)$transcript_counts >= 1]
dim(L15_spe)

L15_spe <- filter_genes(L15_spe)

# calculate logcounts (log-transformed normalized counts) using scran package
# using library size factors
L15_spe <- computeLibraryFactors(L15_spe)
L15_spe <- logNormCounts(L15_spe)
assayNames(L15_spe)

# select small set of random genes and several known SVGs for 
# faster runtime in this example
set.seed(123)
#ix_random <- sample(seq_len(nrow(spe)), 10)

known_genes <- c("MOBP", "PCP4", "SNAP25")
ix_known <- which(rowData(L15_spe)$gene_name %in% known_genes)

#ix <- c(ix_known, ix_random)

L15_spe_subset <- L15_spe[ix_known, ]
dim(L15_spe_subset)

# set seed for reproducibility
# run nnSVG using a single thread for this example workflow
set.seed(123)
L15_spe_nnSVG <- nnSVG(L15_spe, n_threads = 1)

table(rowData(L15_spe_nnSVG)$padj <= 0.05)
