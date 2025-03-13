## Using nnSVG to find the spatially variable genes in the xenium data
## Here get all transcripts for all genes with coordinates and binned the expression in zones (like for visium)
## Loop it through different nb of bins to see how many SVG obtained with different nb of bins

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

# field of view for L15, where all the genes are stored individually, containing the coordinates of all the transcripts
L15_fov <- L15@images$L15@molecules$molecules

# all the genes analyzed with xenium (aka 266 panel + 100 epil specific)
gene_list <- Features(L15)

# create a tibble to store the genes + their coordinates
transcript_coord <- tibble(gene = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)

for (n in (1:length(gene_list))){
  gene <- gene_list[n] #going through the list of genes
  gene_matrix <- L15_fov[[gene]]@coords #selecting in the fov the gene
  temp_df <- tibble(gene = gene, x = gene_matrix[,1], y = gene_matrix[,2]) # temporary df to add all transcripts with their coordinates (transcript name being just the gene name)
  transcript_coord <- rbind(transcript_coord,temp_df) # joining the temp df with the tibble containing all genes
  rm(temp_df)
}

# initialize the list which will contain all the SVG
list_svg = list()
# Define the dimensions of the space (constant for one sample)
x_size <- ceiling(max(transcript_coord$x))
y_size <- ceiling(max(transcript_coord$y))

for (i in seq(from = 800, to = 5000, by = 100)){
  transcript_coord <- transcript_coord %>% 
    select(gene, x, y)
  
  # Number of spots
  n_bin <- i
  
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
  bin_grid <- expand_grid(y = 1:number_bin_y, x = 1:number_bin_x) %>% #create all combination possible with the value given (here being the number of bins in X and Y), with the second argument being the one used first to increase value (1-1, 1-2, 1-3, ... 1-n, 2-1, 2-2, ... 2-n, etc...)
    select(x,y) %>% # re-order the columns to have x first then y
    as.data.frame() # convert into a df
  # name each row with the corresponding bin
  row.names(bin_grid) <- paste0("bin_",as.character(1:(number_bin_x*number_bin_y)))
  
  # create a basic spatial experiment object containing just enough information for normalization + log transform + nnsvg
  spe <- SpatialExperiment(
    assay = list(counts = transcript_genexbin_mat),
    colData = bin_grid,
    rowData = list(gene_name = gene_list),
    spatialCoordsNames = c("x","y")
  )
  
  # only keeps the bins with transcripts
  spe_subset <- spe[,bin_w_transcript]
  
  # why adding bins w/o transcripts and then removing them? In the coordinate of the bins, it contains all bins, not just the ones with transcripts. As this is the colData argument, the number of bins in colData needs to match the one of the assay.
  
  #filtering to remove lowly expressed genes
  spe_filtered <- filter_genes(spe_subset)
  
  # normalize + log transform
  spe_norm <- computeLibraryFactors(spe_filtered)
  spe_norm <- logNormCounts(spe_norm)
  
  assayNames(spe_norm)
  
  # nn SVG
  spe_nnSVG <- nnSVG(spe_norm)
  
  # creat the svg name pasting SVG_ to number of bin used
  svg_name <- paste0("SVG_", as.character(i))
  
  # adding svg to the list
  list_svg[[svg_name]] <- spe_nnSVG
}

# creating data frame containing all the SVG
SVG_nb <- data.frame(table(rowData(list_svg[["SVG_500"]])$padj<=0.05))[,2]
for (i in svg_name){
  SVG_nb <- rbind(SVG_nb,data.frame(table(rowData(list_svg[[i]])$padj<=0.05))[,2])
}
rownames(SVG_nb) <- c("SVG_500",svg_name)

# I try from 500 bins to 5000 here, all having almost all genes with padj <0.05