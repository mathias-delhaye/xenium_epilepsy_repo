## Looking at the tot transcript from xenium runs and the ones detected in the segmented areas

## Loading libraries

library(tidyverse)
library(Seurat)
library(ggplot2)
library(here)
library(future)
library(RColorBrewer)
library(scCustomize)

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

## Now in the megadata

# read the megadata
megadata = read_rds(file.path(str_remove(path,"xenium_epilepsy_repo"),"REALfull_QC_RCTD_clust_processed.rds"))

# Looking into the megadata

megadata_matrix <- megadata@assays$Xenium@layers

# list of all targets from the panel (366 genes)
gene_list <- Features(megadata)

# make a tibble from the vector contain the targets name, it will be use to add the tot transcripts and the ones out
transcript_df <- tibble(gene_list)
samples <- levels(megadata@meta.data$sample_type)
for (i in 1:16){
  tot_transcript = c()
  mol_list <- megadata@images[[i]]@molecules$molecules # get the list of mol from the image in the megadata
  for (molname in gene_list) { #going through the genes
    tot_transcript <- c(tot_transcript,length(mol_list[[molname]])) #adding the tot transcript for each gene
  }
  transcript_df[[paste("tot_transcript_",samples[i],sep = "")]] <- tot_transcript # adding the tot as a new column
  transcript_df[[paste("in_segmentation_transcript_",samples[i],sep = "")]] <- as.numeric(rowSums(megadata@assays$Xenium@layers[[i]])) 
}

transcript_df <- gather(transcript_df,key = 'origin', value = 'transcript_number',!gene_list)

transcript_df <- transcript_df %>%
  extract(origin, into = c("transcript_origin", "sample"), regex = "(.*)_([A-Z]\\d+_[ep])")
