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
fig_path <- file.path(str_remove(path,"/xenium_epilepsy_repo"),"figures/transcripts_location")

# read the megadata
megadata = read_rds(file.path(str_remove(path,"xenium_epilepsy_repo"),"REALfull_QC_RCTD_clust_processed.rds"))

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

# tidying the data
transcript_df_long <- gather(transcript_df,key = 'origin', value = 'transcript_number',!gene_list)

# separating the origin into transcript location (in segmented and tot) and the sample
transcript_df_long <- transcript_df_long %>%
  extract(origin, into = c("transcript_location", "sample"), regex = "(.*)_([A-Z]\\d+_[ep])") %>% 
  mutate(across(where(is.character),as.factor)) #change the character variables into factor ones

# compute the ratio in_segmentation_transcript/tot_/tot_transcript
ratio_transcript <- transcript_df_long %>% 
  group_by(gene_list,sample) %>% 
  summarise(ratio=min(transcript_number)/max(transcript_number)) %>% #min(transcript_number)==in_segmentation & max(transcript_number)==tot_trancript
  ungroup()

# Create a matrix with genes as rows and samples as columns
gene_matrix <- ratio_transcript %>%
  pivot_wider(names_from = sample, values_from = ratio) %>%
  column_to_rownames("gene_list") %>%
  as.matrix()

# Plot heatmap to look at the distribution of the ratio for each gene and each sample
heatmap_0to1scale <- function(x = mat, show_rownames = TRUE, fontsize_row = 7){
  pheatmap(x, 
           cluster_rows = TRUE,  # Cluster similar genes
           cluster_cols = TRUE,
           show_rownames = show_rownames,
           fontsize_row = fontsize_row,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "Heatmap of Transcript Ratios")
}
hm_tot <- heatmap_0to1scale(gene_matrix, show_rownames = FALSE)
ggsave(file.path(fig_path,"hm_tot.png"),hm_tot, width = 9, height = 6, units = "in", dpi = 150, bg="white")

# list of genes from human panel
gene_list_loc <- read.csv(file.path(str_remove(path,"xenium_epilepsy_repo"),"gene_panel.csv"))

# all panel genes
gene_matrix_panel <- gene_matrix[gene_list_loc$Gene,]

# synaptic genes only
syn_gene <- gene_list_loc %>% 
  filter(Location == "synaptic")
gene_matrix_syn <- gene_matrix[syn_gene$Gene,]                         

hm_panel <- heatmap_0to1(gene_matrix_panel, fontsize_row = 5)
ggsave(paste(path,"/figures/transcripts_location/hm_panel.png", sep = ""),hm_panel, width = 9, height = 6, units = "in", dpi = 150, bg="white")
hm_syn <- heatmap_0to1(gene_matrix_syn)
ggsave(paste(path,"/figures/transcripts_location/hm_syn.png", sep = ""),hm_syn, width = 9, height = 6, units = "in", dpi = 150, bg="white")

# summarise the data
summary_stats_sans_l10 <- ratio_transcript %>%
  filter(sample != "L10_p") %>% 
  group_by(gene_list) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    min_ratio = min(ratio, na.rm = TRUE),
    max_ratio = max(ratio, na.rm = TRUE)
  ) %>%
  arrange(mean_ratio) %>% 
  mutate(rank = row_number()) %>% 
  ungroup()

# scatter plot of the gene rank by the mean ratio of transcript in_segmentation/tot_transcripts
# synaptic genes are highlighted

p <- ggplot(summary_stats_sans_l10, aes(x = rank, y = mean_ratio, color = location, size = point_size)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = label), vjust = -1, size = 3, color = "black",max.overlaps = 25, box.padding = 0.5) +  # Add gene names and prevent overlap
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("synaptic" = "red", "other" = "blue")) + 
  scale_size_continuous(range = c(1, 3), guide = "none") +  # Ensures small vs large dot contrast
  theme_minimal() +
  scale_y_continuous(limits = c(0,1))+
  labs(title = "Ranked Scatter Plot of Transcript Ratios",
       x = "Gene Rank",
       y = "Mean Ratio (Soma / Total)",
       color = "Gene Type") # Legend title

ggsave(paste(path,"/figures/transcripts_location/ranked_scattered_plot.png", sep = ""),p, width = 12, height = 6, units = "in", dpi = 150, bg="white")
