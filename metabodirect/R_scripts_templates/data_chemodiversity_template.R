# #############################################################
# 
# MetaboDirect
# Data Exploration step
# version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
# 
# #############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ggpubr)
  library(vegan)
  library(SYNCSA)
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

setwd('%currentdir%')

my_data.file <- file.path('%outdir%', '1_preprocessing_output', 'Report_processed_noNorm.csv')
my_metadata.file <- file.path('%metadata%')
my_outdir <- file.path('%outdir%', '4_chemodiversity')

#### Import data ####

df <-  read_csv(my_data.file, col_types = cols())
metadata <- read_csv(my_metadata.file, col_types = cols())

#### Get intensity matrix for diversity analysis ####

intensity_matrix <- df %>%
  select(Mass, all_of(metadata$SampleID)) %>%
  pivot_longer(!Mass, names_to = 'SampleID', values_to = 'intensity') %>%
  pivot_wider(names_from = 'Mass', values_from = 'intensity') %>%
  column_to_rownames(var = 'SampleID')

# Sum normalize intensities

sample_sum = rowSums(intensity_matrix)
norm_intensity_matrix <- intensity_matrix / sample_sum


# Get color palettes

my_colors <- get_palette('Dark2', length(unique(metadata$'%group1%')))
names(my_colors) <- unique(metadata$'%group1%') 

#### Abundance-based diversity ####

# Richness

richness <- specaccum(norm_intensity_matrix, method = 'random', permutations = 100)

richness_long <- as_tibble(richness$perm, rownames = 'Sites') %>%
  pivot_longer(!Sites, names_to = 'permutation', values_to = 'Richness') %>%
  select(-permutation)

richness_long$Sites <- as.numeric(richness_long$Sites)

richness_plot <- ggplot(richness_long,
                        aes(x = Sites,
                            y = Richness,
                            group = Sites)) +
  geom_boxplot(fill = 'yellow') +
  geom_hline(yintercept = max(richness_long$Richness) * 0.85, color = 'red') +
  theme_bw() +
  labs(title = 'Richness plot') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'Diversity_plot_richness.png')
ggsave(filename, richness_plot, dpi = 300, width = 8, height = 8)

# Rank abundance

rank_sums <- tibble(peak = colnames(intensity_matrix),
                    intensity_sums = colSums(intensity_matrix),
                    presence = colSums(intensity_matrix != 0)) %>%
  arrange(intensity_sums) %>%
  mutate(position = n():1,
         presence_perc = presence/nrow(intensity_matrix)*100)

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
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'Diversity_plot_rank_abundance.png')
ggsave(filename, rank_plot, dpi = 300, width = 8, height = 8)

# Diversity Index

diversity_table <- tibble(SampleID = rownames(norm_intensity_matrix),
                          Shannon = diversity(norm_intensity_matrix, index = 'shannon'),
                          Gini_Simpson = diversity(norm_intensity_matrix, index = 'simpson'),
                          Chao1 = estimateR(round(intensity_matrix))['S.chao1',])

filename <- file.path(my_outdir, 'abundance_diversity.csv')
write_csv(diversity_table, filename)


diversity_plot <- diversity_table %>% 
  pivot_longer(!SampleID, names_to = 'index', values_to = 'values') %>% 
  left_join(metadata, by = 'SampleID') %>% 
  ggplot() +
  geom_boxplot(aes(x = %group1%,
                   y = values,
                   fill = %group1%)) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~index, scales = 'free_y') +
  labs(title = 'Abundance-based Diversity') +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'abundance_diversity_plot.png')
ggsave(filename, diversity_plot, dpi = 300, width = 8, height = 8)

#### Functional diversity ####

# Elemental composition
el_traits <- df %>% 
  select(Mass, C, H, O, N, S, P) %>% 
  column_to_rownames(var = 'Mass')

# Reactivity
rx_traits <- df %>% 
  select(Mass, NOSC)%>% 
  column_to_rownames(var = 'Mass')

# Insaturation/aromaticity
ia_traits <- df %>% 
  select(Mass, DBE, AI_mod)%>% 
  column_to_rownames(var = 'Mass')

functional_diversity <- tibble(SampleID = rownames(norm_intensity_matrix),
                               Elemental_composition = rao.diversity(norm_intensity_matrix, el_traits)$FunRao,
                               Reactivity = rao.diversity(norm_intensity_matrix, rx_traits)$FunRao,
                               Insaturation_and_aromaticity = rao.diversity(norm_intensity_matrix, ia_traits)$FunRao)

filename <- file.path(my_outdir, 'functional_diversity.csv')
write_csv(functional_diversity, filename)

functional_plot <- functional_diversity %>% 
  pivot_longer(!SampleID, names_to = 'index', values_to = 'values') %>% 
  left_join(metadata, by = 'SampleID') %>% 
  ggplot() +
  geom_boxplot(aes(x = %group1%,
                   y = values,
                   fill = %group1%)) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~index, scales = 'free_y') +
  labs(title = 'Functional-based Diversity',
       subtitle = "Rao's quadratic entropy") +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(face = 'italic', hjust = 0.5))

filename <- file.path(my_outdir, 'functional_diversity_plot.png')
ggsave(filename, functional_plot, dpi = 300, width = 8, height = 8)