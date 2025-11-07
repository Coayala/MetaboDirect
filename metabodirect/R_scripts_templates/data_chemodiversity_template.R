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

library(RColorBrewer)
library(ggpubr)
library(vegan)
library(SYNCSA)
library(jsonlite)
library(tidyverse)

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
my_metadata.file <- file.path('%metadata%')
my_outdir <- file.path(output_dir, '4_chemodiversity')

#### Import data ----

df <-  read_csv(my_data.file)
metadata <- read_csv(my_metadata.file)
figure_metadata <- read_json(file.path(output_dir, 'index.json'))

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

# Initializing JSON index of figures

figure_list <- list()

# Set samples color names ----

my_colors <- set_names(get_palette('Dark2', length(unique(metadata[[group1]]))),
                       nm = unique(metadata[[group1]]))

# Abundance-based diversity ----

## Richness

richness <- specaccum(norm_intensity_matrix, 
                      method = 'random', permutations = 100)

richness_plot <- plot_richness(richness)

#richness_plot <- 

filename <- file.path(my_outdir, '4.1_Diversity_plot_richness.png')
ggsave(filename, richness_plot$plot, dpi = 300, width = 12, height = 6)

update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                   figure_title = richness_plot$plot$labels$title)

save_figure_metadata(
  plot_obj = richness_plot,  
  figure_id = str_remove(basename(filename), '.png'),
  analysis_module = '4. Chemodiversity',
  figure_type = 'Richness plot',
  caption = paste0('Species (richness) accumulation plot showing how many metabolites ',
                   'are detected based on samplimg effort'),
  figure_file = filename,
  group_aes = 'none',
  modifiable_aesthetics = 'none',
  data_source = list(data = my_data.file,
                     sample_metadata = my_metadata.file),
  r_script_path = file.path(my_outdir, 'data_chemodiversity.R'),
  functions_used = list(plot = 'plot_richness()'),
  resolution = 300,
  height = 12,
  width = 6
)

## Rank abundance ----

rank_plot <- plot_rank_abundance(norm_intensity_matrix)

filename <- file.path(my_outdir, '4.2_Diversity_plot_rank_abundance.png')
ggsave(filename, rank_plot$plot, dpi = 300, width = 12, height = 6)

update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                   figure_title = richness_plot$plot$labels$title)

save_figure_metadata(
  plot_obj = rank_plot,  
  figure_id = str_remove(basename(filename), '.png'),
  analysis_module = '4. Chemodiversity',
  figure_type = 'Rank abundance plot',
  caption = paste0('Rank abundance plot showing the distribution of presence vs ',
                   'relative abundance of metabolites'),
  figure_file = filename,
  group_aes = 'none',
  modifiable_aesthetics = 'none',
  data_source = list(data = my_data.file,
                     sample_metadata = my_metadata.file),
  r_script_path = file.path(my_outdir, 'data_chemodiversity.R'),
  functions_used = list(plot = 'plot_rank_abundance()'),
  resolution = 300,
  height = 12,
  width = 6
)

## Diversity Index ----

diversity_table <- calculate_diversity_index(intensity_matrix, norm_intensity_matrix)

diversity_plot <- plot_diversity_index(diversity_table, group_by = group1,
                                       title = 'Abundance-based Diversity')

filename <- file.path(my_outdir, '4.3.1_abundance_diversity_plot.png')
ggsave(filename, diversity_plot$plot, dpi = 300, width = 12, height = 6)

update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                   figure_title = diversity_plot$plot$labels$title)

save_figure_metadata(
  plot_obj = diversity_plot,  
  figure_id = str_remove(basename(filename), '.png'),
  analysis_module = '4. Chemodiversity',
  figure_type = 'Diversity plot',
  caption = paste0('Abundance based diversity metrics including: ',
                   'Chao1, Gini-Simpson and Shannon index'),
  figure_file = filename,
  group_aes = c('fill', 'x'),
  modifiable_aesthetics = c('fill', 'x'),
  data_source = list(data = my_data.file,
                     sample_metadata = my_metadata.file),
  r_script_path = file.path(my_outdir, 'data_chemodiversity.R'),
  functions_used = list(calculate = 'calculate_diversity_index()',
                        plot = 'plot_diversity_index()'),
  resolution = 300,
  height = 12,
  width = 6
)

# Functional diversity ----

functional_diversity <- calculate_functional_div_index(df, norm_intensity_matrix)

functional_div_plot <- plot_diversity_index(functional_diversity, group_by = group1,
                                            title = 'Functional-based Diversity') 


filename <- file.path(my_outdir, '4.3.2_functional_diversity_plot.png')
ggsave(filename, functional_div_plot$plot, dpi = 300, width = 12, height = 6)

update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                   figure_title = functional_div_plot$plot$labels$title)

save_figure_metadata(
  plot_obj = functional_div_plot,  
  figure_id = str_remove(basename(filename), '.png'),
  analysis_module = '4. Chemodiversity',
  figure_type = 'Diversity plot',
  caption = paste0("Functional based diversity metrics (Rao's quadratic entropy) based on: ",
                   'Elemental composition, NOSC (reactivity) and ',
                   'DBE/AI_mod (insaturation_and_aromaticity'),
  figure_file = filename,
  group_aes = c('fill', 'x'),
  modifiable_aesthetics = c('fill', 'x'),
  data_source = list(data = my_data.file,
                     sample_metadata = my_metadata.file),
  r_script_path = file.path(my_outdir, 'data_chemodiversity.R'),
  functions_used = list(calculate = 'calculate_functional_div_index()',
                        plot = 'plot_diversity_index()'),
  resolution = 300,
  height = 12,
  width = 6
)

index_figures <- read_json(file.path(output_dir, 'index.json'))

index_figures[[3]] <- list(
  analysis_module = '4. Chemodiversity',
  script_name = file.path(my_outdir, 'data_chemodiversity.R'),
  last_run = Sys.time(),
  figures = figure_list
)

write_json(index_figures, 
           file.path(output_dir, 'index.json'), 
           auto_unbox = TRUE, 
           pretty = TRUE)
