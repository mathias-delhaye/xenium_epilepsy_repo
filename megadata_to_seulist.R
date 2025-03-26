## Building seurat list containing all samples binned with the all the transcripts detected by xenium contained in each bin

## Loading libraries

library(tidyverse)
library(Seurat)
library(ggplot2)
library(here)
library(ggrepel)

options(future.globals.maxSize = 25000 * 1024 ^ 2)

path <- here()
fig_path <- file.path(str_remove(path,"/xenium_epilepsy_repo"),"figures/precast - k cluster")

# read the megadata
megadata = read_rds(file.path(str_remove(path,"xenium_epilepsy_repo"),"MegaData_v3.rds"))
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
  n_bin <- x_size*y_size/20000 # 20000 arbitrary surface (px^2) chosen
  
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
  # adding the position windows for each bin in X & Y so that we can place segmented area in each bin based on the centroids position
  bin_grid$x_window <- bin_grid$col*bin_width
  bin_grid$y_window <- bin_grid$row*bin_height
  
  seu <- CreateSeuratObject(
    counts = transcript_genexbin_mat,
    assay = "RNA",
    meta.data = bin_grid
  )
  seulist[[im]] <- seu
}

#saving seurat list

saveRDS(seulist,file.path(str_remove(path,"/xenium_epilepsy_repo"),"seulist_transcripts_binned.rds"))
