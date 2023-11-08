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
my_metadata.file <- file.path('%metadata%')
my_outdir <- file.path(output_dir, '4_chemodiversity')

#### Import data ----

df <-  read_csv(my_data.file)
metadata <- read_csv(my_metadata.file)

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

richness <- calculate_richness(norm_intensity_matrix)

richness_plot <- ggplot(richness,
                        aes(x = Sites,
                            y = Richness,
                            group = Sites)) +
  geom_boxplot(fill = 'yellow') +
  geom_hline(yintercept = max(richness$Richness) * 0.85, color = 'red') +
  theme_bw() +
  labs(title = 'Richness plot') +
  custom_theme()

filename <- file.path(my_outdir, '1.1_Diversity_plot_richness.png')
ggsave(filename, richness_plot, dpi = 300, width = 8, height = 8)

## Rank abundance ----

rank_sums <- calculate_rank_abundance(norm_intensity_matrix)

rank_plot <- ggplot(rank_sums,
                    aes(x = position,
                        y = intensity_sums,
                        color = presence_perc)) +
  geom_point() +
  scale_color_distiller(palette = 'RdYlBu', limits = c(0, 100)) +
  theme_bw() +
  labs(title = 'Rank abundance plot',
       y = 'Total relative abundance',
       x = 'Molecular rank',
       color = 'Presence') +
  custom_theme()

rank_plot

filename <- file.path(my_outdir, '1.2_Diversity_plot_rank_abundance.png')
ggsave(filename, rank_plot, dpi = 300, width = 8, height = 8)

## Diversity Index ----

diversity_table <- calculate_diversity_index(intensity_matrix, norm_intensity_matrix)

diversity_plot <- plot_diversity_index(diversity_table, group_by = group1,
                                       title = 'Abundance-based Diversity')

filename <- file.path(my_outdir, '2.2_abundance_diversity_plot.png')
ggsave(filename, diversity_plot, dpi = 300, width = 8, height = 8)

# Functional diversity ----

functional_diversity <- calculate_functional_div_index(df, norm_intensity_matrix)

functional_div_plot <- plot_diversity_index(functional_diversity, group_by = group1,
                                            title = 'Functional-based Diversity') +
  labs(subtitle = "Rao's quadratic entropy")

filename <- file.path(my_outdir, '3.2_functional_diversity_plot.png')
ggsave(filename, functional_div_plot, dpi = 300, width = 8, height = 8)
