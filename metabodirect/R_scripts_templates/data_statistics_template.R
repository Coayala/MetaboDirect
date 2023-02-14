# ##############################################################
# 
# MetaboDirect
# Data Statistics step
# version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
# 
# ##############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(factoextra)
  library(ggpubr)
  library(ggnewscale)
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

setwd('%currentdir%')

my_matrix.file <- file.path('%outdir%', '1_preprocessing_output', 'matrix_features.csv')
my_report.file <- file.path('%outdir%', '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_classcomp.file <- file.path('%outdir%', '1_preprocessing_output', 'class_composition.csv')
my_metadata.file <- file.path('%metadata%')
my_outdir <- file.path('%outdir%', '5_statistics')

norm_method <- '%norm_method%'

#### Import data ####

matrix <-  read_csv(my_matrix.file, col_types = cols()) %>% 
  column_to_rownames(var = 'Mass')%>% 
  t()
class_comp <- read_csv(my_classcomp.file, col_types = cols()) %>% 
  column_to_rownames(var = 'SampleID')
metadata <- read_csv(my_metadata.file, col_types = cols())
df <- read_csv(my_report.file)


#### Reformat data files ####

### Change all metadata columns to factors to avoid problems at plotting

for(i in 2:ncol(metadata)){
  metadata[,i] <- factor(as_vector(metadata[,i])) 
}

# Intensity data file

matrix <- matrix[metadata$SampleID,]

class_comp[is.na(class_comp)] <- 0

#### NMDS ####

if(norm_method %in% c('max', 'sum', 'none')){
  dm.method = 'bray'
}else if(norm_method %in% c('mean', 'median', 'z_score', 'minmax')){
  dm.method = 'euclidean'
}else {
  dm.method = 'jaccard'
}

dist.matrix <- vegdist(matrix , method = dm.method) 

set.seed(123)
nmds.log <- capture.output(nmds <- metaMDS(dist.matrix, 
                                           k = 2,
                                           distance = dm.method,
                                           maxit = 999,
                                           trymax = 500,
                                           wascores = TRUE))
filename <- file.path(my_outdir, 'nmds.log')
write_csv(as.data.frame(nmds.log), filename)


filename <- file.path(my_outdir, 'NMDS_stressplot.pdf')

pdf(filename)
capture.output(stressplot(nmds))
dev.off()

nmds.scores <- as.data.frame(scores(nmds, display = 'sites')) %>% 
  rownames_to_column(var = 'SampleID') %>% 
  left_join(metadata, by = 'SampleID')

## Saving data table with NMDS scores

filename <- file.path(my_outdir, 'nmds_scores.csv')
write_csv(nmds.scores, filename)

## Set samples color names

my_colors <- get_palette('Dark2', length(unique(metadata$%group1%)))
names(my_colors) <- unique(metadata$%group1%) 

## NMDS plot

nmds_plot <- ggplot(nmds.scores,
                    aes(x = NMDS1,
                        y = NMDS2,
                        color = %group1%,
                        shape = %group2%)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  labs(title = 'NMDS plot') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'NMDS_plot.png')
ggsave(filename, nmds_plot, dpi = 300)

#### PERMANOVA ####

set.seed(456)

dist.matrix <- vegdist(matrix , method = dm.method) # dist matrix

permanova <- adonis2(dist.matrix ~ %group1% + %group2%, 
                    data=metadata, 
                    permutations=999, 
                    method=dm.method)

filename <- file.path(my_outdir, 'permanova.csv')
write.csv(permanova$aov.tab, filename, row.names = TRUE)

#### PCA of classes ####

pca <- prcomp(class_comp, center = TRUE)
eigen <- get_eigenvalue(pca)

## Saving data table with PCA eigenvalues

filename <- file.path(my_outdir, 'eigen_values_by_compound_class.csv')
write_csv(eigen, filename)

scree_plot <- fviz_eig(pca, addlabels = TRUE) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'scree_plot_by_class.png')
ggsave(filename, scree_plot, dpi = 300)

# Extract sample coordinates for PC1 and PC2
pca_coordinates <- as_tibble(pca$x)
pca_coordinates$SampleID <- rownames(pca$x)
pca_coordinates <- left_join(pca_coordinates, metadata, by ='SampleID')

## Saving data table with PCA coordinates

filename <- file.path(my_outdir, 'pca_coordinates_by_compound_class.csv')
write_csv(pca_coordinates, filename)

pc1 <- paste0('PC1 (', round(eigen$variance.percent[1], digits = 1), '%)')
pc2 <- paste0('PC2 (', round(eigen$variance.percent[2], digits = 1), '%)')

var_loadings <- facto_summarize(pca, element = "var", result=c("coord","contrib","cos2"), axes=c(1,2))
class_color <- get_palette('d3', nrow(var_loadings))

pca_plot <- ggplot(pca_coordinates,
                   aes(x = PC1,
                       y = PC2,
                       color = %group1%,
                       shape = %group2%)) +
  geom_point(size = 2.5, position = position_jitter()) +
  scale_color_manual(values = my_colors) +
  new_scale_color() +
  geom_segment(data = var_loadings,
               aes(x = 0,
                   y = 0,
                   xend = Dim.1 * 2,
                   yend = Dim.2 * 2, 
                   color = name), 
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.03, 'npc')),
               linewidth = 0.8) +
  scale_color_manual(values = class_color) +
  geom_text(data = var_loadings,
            aes(x = ifelse(Dim.1 * 2 > 0, Dim.1 * 2 + 1, Dim.1 * 2 - 1),
                y = ifelse(Dim.2 * 2 > 0, Dim.2 * 2 + 1, Dim.2 * 2 - 1),
                label = name,
                color = name),
            inherit.aes = FALSE,
            hjust = 'inward') +
  guides(color = 'none') +
  theme_bw() +
  labs(title = 'PCA plot by compound classes',
       x = pc1,
       y = pc2) +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'PCA_plot_by_compound_class.png')
ggsave(filename, pca_plot, dpi = 300)

#### PCA by molecular characteristics ####

# Report data file

df_longer <- df %>%
  pivot_longer(metadata$SampleID, names_to = 'SampleID', values_to = 'NormIntensity') %>% 
  filter(NormIntensity != 0) %>% 
  select(NormIntensity, OC, HC, NOSC, GFE, DBE, AI_mod, SampleID) %>% 
  left_join(metadata, by = 'SampleID')

mol_char <- tibble(SampleID = metadata$SampleID)

for(char in c('OC', 'HC', 'NOSC', 'GFE', 'DBE', 'AI_mod')){
  col_name <- paste0(char, '_w')
  char <- syms(char)
  w_df <- df_longer %>% 
    mutate(w_value = (!!! char) * NormIntensity) %>% 
    group_by(SampleID) %>% 
    add_tally(w_value, name = 'sum_w_value') %>% 
    add_tally(NormIntensity, name = 'sum_intensity') %>% 
    ungroup() %>% 
    mutate(w_avg = sum_w_value/sum_intensity) %>%
    select(SampleID, w_avg) %>% 
    distinct()
  colnames(w_df) <- c('SampleID', col_name)
  mol_char <- left_join(mol_char, w_df, by = 'SampleID')
  
}

mol_char <- mol_char %>% 
  column_to_rownames(var = 'SampleID')

pca <- prcomp(mol_char, center = TRUE)
eigen <- get_eigenvalue(pca)

## Saving data table with PCA eigenvalues

filename <- file.path(my_outdir, 'eigen_values_by_molecular_characteristics.csv')
write_csv(eigen, filename)

scree_plot <- fviz_eig(pca, addlabels = TRUE) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'scree_plot_by_molecular_characteristics.png')
ggsave(filename, scree_plot, dpi = 300)

# Extract sample coordinates for PC1 and PC2
pca_coordinates <- as_tibble(pca$x)
pca_coordinates$SampleID <- rownames(pca$x)
pca_coordinates <- left_join(pca_coordinates, metadata, by ='SampleID')

## Saving data table with PCA coordinates

filename <- file.path(my_outdir, 'pca_coordinates_by_molecular_characteristics.csv')
write_csv(pca_coordinates, filename)

pc1 <- paste0('PC1 (', round(eigen$variance.percent[1], digits = 1), '%)')
pc2 <- paste0('PC2 (', round(eigen$variance.percent[2], digits = 1), '%)')

var_loadings <- facto_summarize(pca, element = "var", result=c("coord","contrib","cos2"), axes=c(1,2))
class_color <- get_palette('d3', nrow(var_loadings))

pca_plot <- ggplot(pca_coordinates,
                   aes(x = PC1,
                       y = PC2,
                       color = %group1%,
                       shape = %group2%)) +
  geom_point(size = 2.5, position = position_jitter()) +
  scale_color_manual(values = my_colors) +
  new_scale_color() +
  geom_segment(data = var_loadings,
               aes(x = 0,
                   y = 0,
                   xend = Dim.1 * 1,
                   yend = Dim.2 * 1, 
                   color = name), 
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.03, 'npc')),
               linewidth = 0.8) +
  scale_color_manual(values = class_color) +
  geom_text(data = var_loadings,
            aes(x = ifelse(Dim.1 * 1 > 0, Dim.1 * 1 + 1, Dim.1 * 1 - 1),
                y = ifelse(Dim.2 * 1 > 0, Dim.2 * 1 + 1, Dim.2 * 1 - 1),
                label = name,
                color = name),
            inherit.aes = FALSE,
            hjust = 'inward') +
  guides(color = FALSE) +
  theme_bw() +
  labs(title = 'PCA plot by molecular characteristics',
       x = pc1,
       y = pc2) +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'PCA_plot_by_molecular_characteristics.png')
ggsave(filename, pca_plot, dpi = 300)