## Looking at the tot transcript from xenium runs and the ones detected in the segmented areas

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

## First in one single file

L15 <- read_rds(file.path(str_remove(path,"xenium_epilepsy_repo"),"L15_Seurat_obj.rds"))

# seurat matrix of the xenium transcript (in the segmented regions)
L15_xenium_mat <- L15[["Xenium"]]

# list of all transcripts detected (panel + controls + unrecognized)
L15_mol_list <- L15@images$L15@molecules$molecules

# list of all targets from the panel (366 genes)
gene_list <- Features(L15)

# make a tibble from the vector contain the targets name, it will be use to add the tot transcripts and the ones out
transcript_df <- tibble(Features(L15))

#initializing an empty vector for the tot transcripts
tot_transcript <- c()

for (molname in gene_list) { #going through the genes
  tot_transcript <- c(tot_transcript,length(L15_mol_list[[molname]])) #adding the tot transcript for each gene
}

transcript_df$tot_transcript <- tot_transcript

# adding to the tibble the transcripts in the assay xenium (in the segmented areas)
transcript_df$insegmentation_transcript <- as.numeric(rowSums(L15_xenium_mat$counts))
