library(ggplot2)
library(Seurat)

# Create a list to store plots
plots_list <- list()

precast_obj <- precast_list[[1]]  # Get the PRECAST object for the current k
k_value <- 15  # Assuming k ranges from 15 to 25

# Extract Seurat object
seurat_obj <- precast_obj@seulist

# Ensure the cluster assignments are in the metadata
seurat_obj$precast_cluster <- as.factor(precast_obj$domain$cluster)

# Plot using DimPlot (Seurat) or ggplot2
p <- DimPlot(seurat_obj, group.by = "precast_cluster", cols = rainbow(length(unique(seurat_obj$precast_cluster)))) +
  ggtitle(paste("Spatial Domains for k =", k_value))

# Loop over each k value in precast_list
for (i in seq_along(precast_list)) {
  precast_obj <- precast_list[[i]]  # Get the PRECAST object for the current k
  k_value <- 15 + (i - 1)  # Assuming k ranges from 15 to 25
  
  # Extract Seurat object
  seurat_obj <- precast_obj$SeuratObj
  
  # Ensure the cluster assignments are in the metadata
  seurat_obj$precast_cluster <- as.factor(precast_obj$domain$cluster)
  
  # Plot using DimPlot (Seurat) or ggplot2
  p <- DimPlot(seurat_obj, group.by = "precast_cluster", cols = rainbow(length(unique(seurat_obj$precast_cluster)))) +
    ggtitle(paste("Spatial Domains for k =", k_value))
  
  # Store the plot
  plots_list[[i]] <- p
}

# Display all plots
library(patchwork)  # Allows easy layout of multiple plots
wrap_plots(plots_list, ncol = 3)
