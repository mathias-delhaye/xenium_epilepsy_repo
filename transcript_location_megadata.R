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

# compute the ratio in_segmentation_transcript/tot_transcript
ratio_transcript <- transcript_df_long %>% 
  group_by(gene_list,sample) %>% 
  summarise(ratio=min(transcript_number)/max(transcript_number)) %>% #min(transcript_number)==in_segmentation & max(transcript_number)==tot_trancript
  ungroup()

# Create a matrix with genes as rows and samples as columns
gene_matrix <- ratio_transcript %>%
  pivot_wider(names_from = sample, values_from = ratio) %>%
  column_to_rownames("gene_list") %>%
  as.matrix()

# Extract the category for each sample
sample_to_ILAE <- unique(megadata@meta.data[, c("sample_type", "ILAE_score")]) 
sample_to_ILAE <- sample_to_ILAE %>% 
  mutate(sample_type = case_when(
    sample_type == "E008_efr"~ "E008_e",
    sample_type == "E015_efr" ~ "E015_e",
    TRUE ~ sample_type #keep all other values unchanged
  ),
  ILAE_score = case_when(
    sample_type == "E015_e" ~ "C",
    TRUE ~ ILAE_score #keep all other values unchanged
  ))
rownames(sample_to_ILAE) <- NULL

# Plot heatmap to look at the distribution of the ratio for each gene and each sample
heatmap_0to1scale <- function(x = mat, show_rownames = TRUE, fontsize_row = 7, col_names = colnames(x)){
  pheatmap(x, 
           cluster_rows = TRUE,  # Cluster similar genes
           cluster_cols = TRUE,
           show_rownames = show_rownames,
           fontsize_row = fontsize_row,
           labels_col = col_names,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "Heatmap of Transcript Ratios")
}
hm_tot <- heatmap_0to1scale(gene_matrix, show_rownames = FALSE)
hm_tot_ILAE <- heatmap_0to1scale(gene_matrix, show_rownames = FALSE, col_names = sample_to_ILAE$ILAE_score)
ggsave(file.path(fig_path,"hm_tot.png"),hm_tot, width = 9, height = 6, units = "in", dpi = 150, bg="white")
ggsave(file.path(fig_path,"hm_tot_ILAE.png"),hm_tot_ILAE, width = 9, height = 6, units = "in", dpi = 150, bg="white")

# list of genes from human panel
gene_list_loc <- read.csv(file.path(str_remove(path,"xenium_epilepsy_repo"),"gene_panel.csv"))

# all panel genes
gene_matrix_panel <- gene_matrix[gene_list_loc$Gene,]

# synaptic genes only
syn_gene <- gene_list_loc %>% 
  filter(Location == "synaptic")
gene_matrix_syn <- gene_matrix[syn_gene$Gene,]                         

hm_panel <- heatmap_0to1scale(gene_matrix_panel, fontsize_row = 5)
ggsave(paste(fig_path,"/hm_panel.png", sep = ""),hm_panel, width = 9, height = 6, units = "in", dpi = 150, bg="white")
hm_syn <- heatmap_0to1scale(gene_matrix_syn)
ggsave(paste(fig_path,"/hm_syn.png", sep = ""),hm_syn, width = 9, height = 6, units = "in", dpi = 150, bg="white")

# summarise the data
summary_stats <- ratio_transcript %>%
  #filter(sample_type != "L10_p") %>% 
  group_by(gene_list) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    min_ratio = min(ratio, na.rm = TRUE),
    max_ratio = max(ratio, na.rm = TRUE)
  ) %>%
  arrange(mean_ratio) %>% 
  mutate(rank = row_number()) %>% 
  ungroup() %>% 
  mutate(location = if_else(gene_list%in%syn_gene$Gene, "synaptic", "non_synaptic"), gene_list = as.character(gene_list))
# if_else used to assign synaptic if genes in gene_list are in syn_gene$Gene, else non_synaptic, as change gene_list in character as needed for plotting

# scatter plot of the gene rank by the mean ratio of transcript in_segmentation/tot_transcripts
# synaptic genes are highlighted

p <- ggplot(summary_stats, aes(x = rank, y = mean_ratio, color = location)) +
  geom_point(aes(size = ifelse(location == "synaptic", 3, 1)),alpha = 0.7) + # if synaptic gene, size 3, else 1
  geom_text_repel(aes(label = ifelse(location == "synaptic", gene_list, NA)), # if synaptic gene, label = gene_list, else NA
                  vjust = -1, size = 3, color = "black", max.overlaps = 25, box.padding = 0.5) +  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  scale_color_manual(values = c("synaptic" = "red", "other" = "blue")) + 
  scale_size_identity(guide = "none") +  # Directly use mapped sizes without scale modification
  theme_minimal() +
  scale_y_continuous(limits = c(0,1))+
  labs(title = "Ranked Scatter Plot of Transcript Ratios",
       x = "Gene Rank",
       y = "Mean Ratio (Soma / Total)",
       color = "Gene Type") # Legend title

ggsave(paste(fig_path,"/ranked_scattered_plot.png", sep = ""),p, width = 12, height = 6, units = "in", dpi = 150, bg="white")


## Test quality sections in relation with ratio

df_transcript_quality <- megadata@meta.data %>% # table with all the variable that could impact the quality of the sample
  group_by(sample_type) %>% 
  summarize(median_cell_count = median(nCount_Xenium),
            block_age = unique(age_block),
            ILAE = unique(ILAE_score),
            PMI = unique(PMI))

ratio <- c()
for (i in 1:16 ){
  tot_transcript = c()
  mol_list <- megadata@images[[i]]@molecules$molecules # get the list of mol from the image in the megadata
  for (molname in gene_list) { #going through the genes
    tot_transcript <- c(tot_transcript,length(mol_list[[molname]])) #adding the tot transcript for each gene
  }
  ratio <- c(ratio, sum(megadata@assays$Xenium@layers[[i]])/sum(tot_transcript)) # ratio of the tot nb of transcript detected in segmented area for each sample/tot nb of transcript for the sample
}

df_transcript_quality$ratio_transcript = ratio #adding a column for the ratio

df_transcript_quality$ILAE[sample_type=="E008_efr"] = NA

df_transcript_quality$ILAE[sample_type=="E015_efr"] = "C"

p_block_age <- df_transcript_quality %>% 
  #filter(is.na(ILAE)|ILAE != "P") %>% 
  ggplot(aes(x = ratio_transcript, y = block_age, label = sample_type, colour = ILAE))+
  geom_point()+
  geom_text_repel(size = 3,  box.padding = 0.5, max.overlaps = 20)+
  #scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,10))+
  geom_smooth(method = "lm", color = "blue", se = FALSE)
print(p_block_age)

ggsave(paste(fig_path,"/ratio_transcriptxblock_age.png", sep = ""),p_block_age, units = "in", dpi = 150, bg="white")

p_median_transcript <- df_transcript_quality %>% 
  #filter(is.na(ILAE)|ILAE != "P") %>% 
  ggplot(aes(x = ratio_transcript, y = median_cell_count, label = sample_type, colour = ILAE))+
  geom_point()+
  geom_text_repel(size = 3,  box.padding = 0.5, max.overlaps = 20)+
  #scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,10))+
  geom_smooth(method = "lm", color = "blue", se = FALSE)
print(p_median_transcript)

ggsave(paste(fig_path,"/ratio_transcriptxmedian_transcript.png", sep = ""),p_median_transcript, units = "in", dpi = 150, bg="white")
