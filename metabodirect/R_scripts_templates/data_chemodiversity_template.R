# ***************************************************************
# 
# MetaboDirect
# Data Exploration step
# version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
# 
# ***************************************************************

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(SYNCSA)

# Values between two '%' are to be replaced by the correct values during the python script

# Defining variables ----

current_dir <- '%currentdir%'
output_dir <- '%outdir%'
group1 <- '%group1%'
group2 <- '%group2%'

setwd(current_dir)

# Loading custom functions ----
source(file.path(output_dir, 'custom_functions.R'))

my_data.file <- file.path(output_dir, '1_preprocessing_output', 'Report_processed_noNorm.csv')
my_metadata.file <- file.path('data/processed/fticr_metadata.csv')
my_outdir <- file.path(output_dir, '4_chemodiversity')

#### Import data ----

df <-  read_csv(my_data.file)
metadata <- read_csv(my_metadata.file)
figure_metadata <- read_json(file.path(current_dir, 'index.json'))

# Reformat data files ----

## Change all metadata columns to factors to avoid problems at plotting

metadata <- metadata %>% 
  mutate(across(!SampleID, as.factor))

## Intensity matrix

intensity_matrix <- df %>%
  select(Mass, all_of(metadata$SampleID)) %>%
  column_to_rownames(var = 'Mass') %>% 
  t()

## Sum normalize intensities

norm_intensity_matrix <- decostand(intensity_matrix, method = 'total')

# Set samples color names ----

my_colors <- set_names(get_palette('Dark2', length(unique(metadata[[group1]]))),
                       nm = unique(metadata[[group1]]))

# Abundance-based diversity ----

## Richness

richness <- specaccum(norm_intensity_matrix, 
                      method = 'random', permutations = 100)

richness_plot <- plot_richness(richness)

#richness_plot <- 

filename <- file.path(my_outdir, '1.1_Diversity_plot_richness.png')
ggsave(filename, richness_plot$plot, dpi = 300, width = 8, height = 8)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = richness_plot$plot$labels$title,
                    figure_type = 'Richness plot',
                    caption = 'Richness accumulation plot',
                    description = richness_plot$meta$description,
                    insights = richness_plot$meta$insight,
                    x_axis_label = 'Sites',
                    y_axis_label = 'Richness',
                    grouping_variables = 'none',
                    data_source = list(
                      input_file = my_data.file,
                      columns_used = metadata$SampleID
                    ),
                    r_script_path = file.path(output_dir, '4_chemodiversity', 'data_chemodiversity.R'),
                    analysis_category = 'Chemodiversity',
                    figure_info = list(
                      figure_file = file,
                      figure_width = 8,
                      figure_height = 8,
                      dpi = 300
                    ))

## Rank abundance ----

rank_plot <- plot_rank_abundance(norm_intensity_matrix)

filename <- file.path(my_outdir, '1.2_Diversity_plot_rank_abundance.png')
ggsave(filename, rank_plot$plot, dpi = 300, width = 8, height = 8)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = rank_plot$plot$labels$title,
                    figure_type = 'Rank abundance plot',
                    caption = 'Rank abundance plot showing the distribution of presence vs relative abundance of metabolites',
                    description = rank_plot$meta$description,
                    insights = rank_plot$meta$insight,
                    x_axis_label = 'Molecular Rank',
                    y_axis_label = 'Total relative abundance',
                    grouping_variables = 'none',
                    data_source = list(
                      input_file = my_data.file,
                      columns_used = metadata$SampleID
                    ),
                    r_script_path = file.path(output_dir, '4_chemodiversity', 'data_chemodiversity.R'),
                    analysis_category = 'Chemodiversity',
                    figure_info = list(
                      figure_file = file,
                      figure_width = 8,
                      figure_height = 8,
                      dpi = 300
                    ))

## Diversity Index ----

diversity_table <- calculate_diversity_index(intensity_matrix, norm_intensity_matrix)

diversity_plot <- plot_diversity_index(diversity_table, group_by = group1,
                                       title = 'Abundance-based Diversity')

filename <- file.path(my_outdir, '2.2_abundance_diversity_plot.png')
ggsave(filename, diversity_plot$plot, dpi = 300, width = 8, height = 8)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = diversity_plot$plot$labels$title,
                    figure_type = 'Diversity plot',
                    caption = 'Abundance based diversity metrics',
                    description = diversity_plot$meta$description,
                    insights = diversity_plot$meta$insight,
                    x_axis_label = group1,
                    y_axis_label = 'Diversity index',
                    grouping_variables = group1,
                    data_source = list(
                      input_file = my_data.file,
                      columns_used = c(metadata$SampleID, group1)
                    ),
                    r_script_path = file.path(output_dir, '4_chemodiversity', 'data_chemodiversity.R'),
                    analysis_category = 'Chemodiversity',
                    figure_info = list(
                      figure_file = file,
                      figure_width = 8,
                      figure_height = 8,
                      dpi = 300
                    ))

# Functional diversity ----

functional_diversity <- calculate_functional_div_index(df, norm_intensity_matrix)

functional_div_plot <- plot_diversity_index(functional_diversity, group_by = group1,
                                            title = 'Functional-based Diversity') 

functional_div_plot$plot <- functional_div_plot$plot +
  labs(subtitle = "Rao's quadratic entropy")

filename <- file.path(my_outdir, '3.2_functional_diversity_plot.png')
ggsave(filename, functional_div_plot$plot, dpi = 300, width = 8, height = 8)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = functional_div_plot$plot$labels$title,
                    figure_type = 'Diversity plot',
                    caption = "Functional based diversity metrics (Rao's quadratic entropy)",
                    description = functional_div_plot$meta$description,
                    insights = functional_div_plot$meta$insight,
                    x_axis_label = group1,
                    y_axis_label = 'Diversity index',
                    grouping_variables = group1,
                    data_source = list(
                      input_file = my_data.file,
                      columns_used = c(metadata$SampleID, group1)
                    ),
                    r_script_path = file.path(output_dir, '4_chemodiversity', 'data_chemodiversity.R'),
                    analysis_category = 'Chemodiversity',
                    figure_info = list(
                      figure_file = file,
                      figure_width = 8,
                      figure_height = 8,
                      dpi = 300
                    ))