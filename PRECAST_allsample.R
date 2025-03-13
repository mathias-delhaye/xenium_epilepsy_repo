## Trying to apply precast for one sample

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
library(PRECAST)

options(future.globals.maxSize = 25000 * 1024 ^ 2)

path <- here()
fig_path <- file.path(str_remove(path,"/xenium_epilepsy_repo"),"figures/transcripts_location")

# read the megadata
megadata = read_rds(file.path(str_remove(path,"xenium_epilepsy_repo"),"REALfull_QC_RCTD_clust_processed.rds"))
# name of all samples in megadata in order assigned in megadata
sample_names <- names(megadata@images)
# all the genes analyzed with xenium (aka 266 panel + 100 epil specific)
gene_list <- Features(megadata)

#initiating seurat list
seulist <- list()

for (im in sample_names){
  # field of view for L15, where all the genes are stored individually, containing the coordinates of all the transcripts
  fov <- megadata@images[[im]]@molecules$molecules
  
  # create a tibble to store the genes + their coordinates
  transcript_coord <- tibble(gene = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)
  
  for (n in (1:length(gene_list))){
    gene <- gene_list[n] #going through the list of genes
    gene_matrix <- fov[[gene]]@coords #selecting in the fov the gene
    temp_df <- tibble(gene = gene, x = gene_matrix[,1], y = gene_matrix[,2]) # temporary df to add all transcripts with their coordinates (transcript name being just the gene name)
    transcript_coord <- rbind(transcript_coord,temp_df) # joining the temp df with the tibble containing all genes
    rm(temp_df)
  }
  
  # only serves to re-select only the coordinates and the genes to re-initiate the table to analyze with new n_bin
  transcript_coord <- transcript_coord %>% 
    select(gene, x, y)
  
  # Define the dimensions of the space (just take the extreme coordinates for x and y)
  x_size <- ceiling(max(transcript_coord$x))
  y_size <- ceiling(max(transcript_coord$y))
  
  # Number of bins 
  n_bin <- 5000
  
  # Calculate the number of spots in x and y directions
  number_bin_x <- ceiling(sqrt(n_bin * (x_size / y_size)))
  number_bin_y <- ceiling(sqrt(n_bin * (y_size / x_size)))
  
  # Calculate the width and height of each spot
  bin_width <- x_size / number_bin_x
  bin_height <- y_size / number_bin_y
  
  # Assign each gene to a spot
  transcript_coord$bin_x <- ceiling(transcript_coord$x / bin_width) # x coordinate / width of bin
  transcript_coord$bin_y <- ceiling(transcript_coord$y / bin_height) # y coordinate / height of bin
  transcript_coord$bin_number <- (transcript_coord$bin_y-1)*number_bin_x+transcript_coord$bin_x # (bin_y coordinate - 1) * number bin X + coordinate bin_x
  
  # counting the number of transcript for each bin and each gene
  transcript_genexbin <- transcript_coord %>% 
    group_by(bin_number, gene) %>% 
    tally() %>% 
    ungroup()
  
  # make the tibble wider to have nb of row = nb of genes and nb of column = nb of bins
  transcript_genexbin <- transcript_genexbin %>% 
    pivot_wider(names_from = bin_number, values_from = n) %>% 
    arrange(gene) %>% 
    column_to_rownames(var = "gene") %>% 
    rename_with(~ paste0("bin_", .)) #change the name of the bin to add "bin_" and not just having numbers
  
  # find the bins that we kept through previous operation, which are bins with at least 1 transcript
  bin_w_transcript <- colnames(transcript_genexbin)
  # find the missing bins and add empty columns with these missing bins
  missing_bin = setdiff(paste0("bin_",as.character(1:(number_bin_x*number_bin_y))),bin_w_transcript)
  transcript_genexbin[,missing_bin] <- NA
  
  #convert tibble into matrix
  transcript_genexbin_mat <- as.matrix(transcript_genexbin)
  
  # replace NA by 0
  transcript_genexbin_mat[is.na(transcript_genexbin_mat)] <- 0
  
  # re-order the column to have them in the order of bins
  transcript_genexbin_mat <- transcript_genexbin_mat[,paste0("bin_",as.character(1:(number_bin_x*number_bin_y)))]
  
  # Get the coordinates for each bin
  bin_grid <- expand_grid(row = 1:number_bin_y, col = 1:number_bin_x) %>% #create all combination possible with the value given (here being the number of bins in X and Y), with the second argument being the one used first to increase value (1-1, 1-2, 1-3, ... 1-n, 2-1, 2-2, ... 2-n, etc...)
    select(col,row) %>% # re-order the columns to have x first then y
    as.data.frame() # convert into a df
  # name each row with the corresponding bin
  row.names(bin_grid) <- paste0("bin_",as.character(1:(number_bin_x*number_bin_y)))
  
  
  seu <- CreateSeuratObject(
    counts = transcript_genexbin_mat,
    assay = "RNA",
    meta.data = bin_grid
  )
  seulist[[im]] <- seu
}

#saving seurat list

#saveRDS(seulist,file.path(str_remove(path,"/xenium_epilepsy_repo"),"seulist_samplesbinned_prePRECAST.rds"))

# testing different k value

precast_list <- list()

preobj <- CreatePRECASTObject(seuList = seulist)

## check the number of genes/features after filtering step
preobj@seulist

for (k in (15:25)){
  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
  PRECASTObj <- AddAdjList(preobj, platform = "Other_SRT")
  
  ## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
  ## information in the algorithm.
  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)
  
  ### Given K
  PRECASTObj <- PRECAST(PRECASTObj, K = k)
  
  precast_list[[paste0("k-", as.character(k))]] <- PRECASTObj
}

saveRDS(precast_list, file.path(str_remove(path,"/xenium_epilepsy_repo"),"precast_list_k15-25.rds"))

for (k in (5:11)){
  seuInt <- L5_k_15_25[[k]]@seulist[[1]]
  seuInt@meta.data$cluster <- factor(unlist(L5_k_15_25[[k]]@resList$cluster))
  seuInt@meta.data$batch <- 1
  seuInt <- Add_embed(L5_k_15_25[[k]]@resList$hZ[[1]], seuInt, embed_name = "PRECAST")
  posList <- lapply(L5_k_15_25[[k]]@seulist, function(x) cbind(x$row, x$col))
  seuInt <- Add_embed(posList[[1]], seuInt, embed_name = "position")
  Idents(seuInt) <- factor(seuInt@meta.data$cluster)
  
  seuInt
  
  ## The low-dimensional embeddings obtained by PRECAST are saved in PRECAST reduction slot.
  png(paste0("L5_",names(L5_k_15_25[k]),".png"),width = 600, height = 450)
  SpaPlot(seuInt, item = "cluster", point_size = 4, combine = F)[[1]] + cowplot::theme_cowplot() +
    ggplot2::xlab("row") + ggplot2::ylab("col")
  dev.off()
}

seuInt <- precast_list[[1]]@seulist[[1]]
seuInt@meta.data$cluster <- factor(unlist(precast_list[[1]]@resList$cluster[1,1]))
seuInt@meta.data$batch <- 1
seuInt <- Add_embed(precast_list[[1]]@resList$hZ[[1]], seuInt, embed_name = "PRECAST")
posList <- lapply(precast_list[[1]]@seulist, function(x) cbind(x$row, x$col))
seuInt <- Add_embed(posList[[1]], seuInt, embed_name = "position")
Idents(seuInt) <- factor(seuInt@meta.data$cluster)

seuInt

## The low-dimensional embeddings obtained by PRECAST are saved in PRECAST reduction slot.
png(paste0("L5_",names(L5_k_15_25[11]),".png"),width = 600, height = 450)
SpaPlot(seuInt, item = "cluster", point_size = 4, combine = F)[[1]] + cowplot::theme_cowplot() +
  ggplot2::xlab("row") + ggplot2::ylab("col")
dev.off()

resList <- precast_list[[1]]@resList
precast_list[[1]] <- SelectModel(precast_list[[1]])

seuInt <- precast_list[[1]]@seulist
