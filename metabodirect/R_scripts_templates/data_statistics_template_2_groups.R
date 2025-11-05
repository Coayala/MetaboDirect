# ***************************************************************
# 
# MetaboDirect
# Data Statistics step
# MetaboDirect version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
# 
# ***************************************************************

library(tidyverse)
library(vegan)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggnewscale)


# Values between two '%' are to be replaced by the correct values during the python script

# Defining variables ----

current_dir <- '%currentdir%'
output_dir <- '%outdir%'
group1 <- '%group1%'
group2 <- '%group2%'

setwd(current_dir)

# Loading custom functions ----
source(file.path(output_dir, 'custom_functions.R'))

# Import data ----

## Defining paths
my_matrix.file <- file.path(output_dir, '1_preprocessing_output', 'matrix_features.csv')
my_report.file <- file.path(output_dir, '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_classcomp.file <- file.path(output_dir, '1_preprocessing_output', 'class_composition.csv')
my_metadata.file <- file.path('data/processed/fticr_metadata.csv')
my_outdir <- file.path(output_dir, '5_statistics')

norm_method <- 'sum'

## Loading tables
matrix <-  read_csv(my_matrix.file) %>% 
  column_to_rownames(var = 'Mass')%>% 
  t()
class_comp <- read_csv(my_classcomp.file) %>% 
  column_to_rownames(var = 'SampleID')
metadata <- read_csv(my_metadata.file)
df <- read_csv(my_report.file)
figure_metadata <- read_json(file.path(current_dir, 'index.json'))

# Reformat data files ----

## Change all metadata columns to factors to avoid problems at plotting

metadata <- metadata %>% 
  mutate(across(!SampleID, as.factor))

## Intensity data file

matrix <- matrix[metadata$SampleID,]

# Set samples color names ----

my_colors <- set_names(get_palette('Dark2', length(unique(metadata[[group1]]))),
                       nm = unique(metadata[[group1]]))

# NMDS ----

nmds_obj <- calculate_nmds(matrix, normalized_with = norm_method)

filename <- file.path(my_outdir, '1.3_nmds_scores.csv')
write_csv(nmds_scores$scores, filename)

# PERMANOVA ----

permanova <- calculate_permanova(matrix, normalized_with = norm_method,
                                 variables = c(group1, group2))

permanova_table <- data.frame(permanova) %>% 
  rownames_to_column(var = 'Variable')

filename <- file.path(my_outdir, '2.1_permanova.csv')
write.csv(permanova_table, filename, row.names = TRUE)

nmds_plot <- plot_ordination(nmds_obj$scores, x = NMDS1, y = NMDS2,
                             color_by = group1, col_vec = my_colors,
                             shape_by = group2,
                             add_info = list(permanova = permanova, 
                                             nmds = nmds_obj$nmds),
                             title = 'NMDS Ordination') 

filename <- file.path(my_outdir, '1.4_NMDS_plot.png')
ggsave(filename, nmds_plot$plot, dpi = 300, width = 6, height = 4)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = nmds_plot$plot$labels$title,
                    figure_type = 'Ordination plot',
                    caption = "NMDS ordination of metabolite abundances",
                    description = nmds_plot$meta$description,
                    insights = nmds_plot$meta$insight,
                    x_axis_label = 'NMDS1',
                    y_axis_label = 'NMDS2',
                    grouping_variables = c(group1, group2),
                    data_source = list(
                      input_file = my_matrix.file,
                      columns_used = c(metadata$SampleID, group1, group2)
                    ),
                    r_script_path = file.path(output_dir, '5_statistics', 'data_statistics.R'),
                    analysis_category = 'Statistics',
                    figure_info = list(
                      figure_file = file,
                      figure_width = 6,
                      figure_height = 4,
                      dpi = 300
                    ))

# PCA of Compound classes ----

pca_class <- calculate_pca(class_comp)

## Saving data table with PCA eigenvalues

filename <- file.path(my_outdir, '3.1_eigen_values_by_compound_class.csv')
write_csv(pca_class$eigenvalues, filename)

## Saving data table with PCA screeplot

filename <- file.path(my_outdir, '3.2_scree_plot_by_compound_class.png')
ggsave(filename, pca_class$scree_plot, dpi = 300, width = 6, height = 4.5)

## Saving data table with PCA coordinates

filename <- file.path(my_outdir, '3.3_pca_coordinates_by_compound_class.csv')
write_csv(pca_class$coordinates, filename)

## PCA_plot Compound classes ----

## Variables with axis names

pc1 <- paste0('PC1 (', 
              round(pca_class$eigenvalues$variance.percent[1], digits = 1), '%)')
pc2 <- paste0('PC2 (', 
              round(pca_class$eigenvalues$variance.percent[2], digits = 1), '%)')

pca_plot <- plot_ordination(pca_class$coordinates, x = PC1, y = PC2,
                            color_by = group1, col_vec = my_colors,
                            shape_by = group2,
                            add_info = list(pca = pca_class),
                            title = 'PCA Compound Classes') 

pca_plot$plot <- pca_plot$plot +
  labs(x = pc1,
       y = pc2) +
  new_scale_color() +
  geom_segment(data = pca_class$loadings,
               aes(x = 0,
                   y = 0,
                   xend = Dim.1 * 2,
                   yend = Dim.2 * 2,
                   color = name), 
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.03, 'npc')),
               linewidth = 0.8,
               show.legend = FALSE) +
  geom_text_repel(data = pca_class$loadings,
                  aes(x = ifelse(Dim.1 * 2 > 0, Dim.1 * 2 + 0.2*Dim.1, Dim.1 * 2 - 0.2*Dim.1),
                      y = ifelse(Dim.2 * 2 > 0, Dim.2 * 2 + 0.2*Dim.2, Dim.2 * 2 - 0.2*Dim.2),
                      label = name,
                      color = name),
                  inherit.aes = FALSE,
                  hjust = 'inward',
                  show.legend = FALSE)

filename <- file.path(my_outdir, '3.4_PCA_plot_by_compound_class.png')
ggsave(filename, pca_plot, dpi = 300, width = 6, height = 4)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = pca_plot$plot$labels$title,
                    figure_type = 'Ordination plot',
                    caption = "PCA ordination plot based of compound classes",
                    description = pca_plot$meta$description,
                    insights = pca_plot$meta$insight,
                    x_axis_label = 'PC1',
                    y_axis_label = 'PC2',
                    grouping_variables = c(group1, group2),
                    data_source = list(
                      input_file = my_classcomp.file,
                      columns_used = c(metadata$SampleID, group1, group2)
                    ),
                    r_script_path = file.path(output_dir, '5_statistics', 'data_statistics.R'),
                    analysis_category = 'Statistics',
                    figure_info = list(
                      figure_file = file,
                      figure_width = 6,
                      figure_height = 4,
                      dpi = 300
                    ))

# PCA by molecular characteristics ----

## Getting dataframe with molecular characteristics ----

mol_char <- df %>%
  pivot_longer(metadata$SampleID, names_to = 'SampleID', values_to = 'NormIntensity') %>% 
  filter(NormIntensity != 0) %>% 
  select(NormIntensity, OC, HC, NOSC, GFE, DBE, AI_mod, SampleID) %>% 
  pivot_longer(!c(SampleID, NormIntensity), names_to = 'index', values_to = 'value') %>% 
  mutate(w_value = value * NormIntensity) %>% 
  group_by(SampleID, index) %>% 
  summarise(w_avg = sum(w_value, na.rm = TRUE) / sum(NormIntensity, na.rm = TRUE)) %>% 
  pivot_wider(names_from = 'index', values_from = 'w_avg') %>% 
  column_to_rownames(var = 'SampleID')

pca_molchar <- calculate_pca(mol_char)

## Saving data table with PCA eigenvalues

filename <- file.path(my_outdir, '4.1_eigen_values_by_molecular_characteristics.csv')
write_csv(pca_molchar$eigenvalues, filename)

## Saving data table with PCA screeplot

filename <- file.path(my_outdir, '4.2_scree_plot_by_molecular_characteristics.png')
ggsave(filename, pca_molchar$scree_plot, dpi = 300, width = 6, height = 4.4)

## Saving data table with PCA coordinates

filename <- file.path(my_outdir, '4.3_pca_coordinates_by_molecular_characteristics.csv')
write_csv(pca_molchar$coordinates, filename)

## PCA_plot Compound classes ----

## Variables with axis names

pc1 <- paste0('PC1 (', 
              round(pca_molchar$eigenvalues$variance.percent[1], digits = 1), '%)')
pc2 <- paste0('PC2 (', 
              round(pca_molchar$eigenvalues$variance.percent[2], digits = 1), '%)')

pca_plot <- plot_ordination(pca_molchar$coordinates, x = PC1, y = PC2,
                            color_by = group1, col_vec = my_colors,
                            shape_by = group2,
                            add_info = list(pca = pca_molchar),
                            title = 'PCA Molecular Characteristics') 

pca_plot$plot <- pca_plot$plot +
  labs(x = pc1,
       y = pc2) +
  new_scale_color() +
  geom_segment(data = pca_molchar$loadings,
               aes(x = 0,
                   y = 0,
                   xend = Dim.1 * 2,
                   yend = Dim.2 * 2,
                   color = name), 
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.03, 'npc')),
               linewidth = 0.8,
               show.legend = FALSE) +
  geom_text_repel(data = pca_molchar$loadings,
                  aes(x = ifelse(Dim.1 * 2 > 0, Dim.1 * 2 + 0.2*Dim.1, Dim.1 * 2 - 0.2*Dim.1),
                      y = ifelse(Dim.2 * 2 > 0, Dim.2 * 2 + 0.2*Dim.2, Dim.2 * 2 - 0.2*Dim.2),
                      label = name,
                      color = name),
                  inherit.aes = FALSE,
                  hjust = 'inward',
                  show.legend = FALSE)

filename <- file.path(my_outdir, '4.4_PCA_plot_by_compound_class.png')
ggsave(filename, pca_plot, dpi = 300, width = 6, height = 4)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = pca_plot$plot$labels$title,
                    figure_type = 'Ordination plot',
                    caption = "PCA ordination plot based of molecular characteristics",
                    description = pca_plot$meta$description,
                    insights = pca_plot$meta$insight,
                    x_axis_label = 'PC1',
                    y_axis_label = 'PC2',
                    grouping_variables = c(group1, group2),
                    data_source = list(
                      input_file = my_classcomp.file,
                      columns_used = c(metadata$SampleID, group1, group2)
                    ),
                    r_script_path = file.path(output_dir, '5_statistics', 'data_statistics.R'),
                    analysis_category = 'Statistics',
                    figure_info = list(
                      figure_file = file,
                      figure_width = 6,
                      figure_height = 4,
                      dpi = 300
                    ))