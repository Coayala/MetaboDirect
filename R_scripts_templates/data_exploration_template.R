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
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

workdir <- file.path(getwd(), '%outdir%', '3_exploratory')
setwd(workdir)

my_data.file <- file.path('..', '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_metadata.file <- file.path('..', '..', '%metadata%')
my_elcomp.file <- file.path('..', '1_preprocessing_output', 'elemental_composition.csv')
my_classcomp.file <- file.path('..', '1_preprocessing_output', 'class_composition.csv')
my_outdir <- getwd()

#### Import data ####

df <-  read_csv(my_data.file, col_types = cols())
metadata <- read_csv(my_metadata.file, col_types = cols())
el_comp <- read_csv(my_elcomp.file, col_types = cols())
class_comp <- read_csv(my_classcomp.file, col_types = cols())

#### Reformat data files ####

# Intensity data file
df_longer <- df %>%
  pivot_longer(metadata$SampleID, names_to = 'SampleID', values_to = 'NormIntensity') %>% 
  filter(NormIntensity > 0)
df_longer <- left_join(df_longer, metadata, by = 'SampleID')

# Class composition
class_comp <- class_comp %>%
  pivot_longer(!SampleID, names_to = 'Class', values_to = 'Count')
class_comp <- left_join(class_comp, metadata, by = 'SampleID')

# Elemental composition
el_comp <- el_comp %>%
  pivot_longer(!SampleID, names_to = 'Element', values_to = 'Count')
el_comp <- left_join(el_comp, metadata, by = 'SampleID')


#### Plot - Van Krevelen Diagram ####

vk_plot <- ggplot(df_longer,
                  aes(x = OC,
                      y = HC,
                      color = Class)) +
  geom_point() +
  theme_bw() +
  scale_fill_brewer(palette = 'Set3') +
  labs(x = 'O:C',
       y = 'H:C',
       title = 'Van Krevelen Diagram') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))+
  facet_grid(rows = vars(%group1%),
             cols = vars(%group2%));

filename <- file.path(my_outdir, 'vk_plot.png')
ggsave(filename, vk_plot, dpi = 300, width = 8);

#### Plot - Density Diagram ####

gfe_density <- ggplot(df_longer,
                      aes(x = GFE,
                          fill = %group1%)) +
  geom_density(alpha = 0.4) +
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw() +
  labs(title = 'GFE Density plot') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
  facet_grid(rows = vars(%group1%),
             cols = vars(%group2%))

filename <- file.path(my_outdir, 'GFE_density.png')
ggsave(filename, gfe_density, dpi = 300)

#### Plot - GFE violin ####

gfe_violin <- ggplot(df_longer,
                     aes(x = %group1%,
                         y = GFE,
                         fill = %group1%)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, show.legend = F) +
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw() +
  labs(title = 'GFE Violin plot') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))+
  facet_grid(cols = vars(%group2%))
  
filename <- file.path(my_outdir, 'GFE_violin.png')
ggsave(filename, gfe_density, dpi = 300)

#### Plot - Class comp bar####

class_bar <- class_comp %>% 
  group_by(%group1%, %group2%, Class)  %>% 
  summarise(Count = mean(Count), .groups = 'drop') %>% 
  ggplot(aes(x = %group1%,
             y = Count,
             fill = Class)) +
  geom_col() +
  theme_bw() +
  scale_fill_brewer(palette = 'Set3') +
  labs(title = 'Class Composition') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
  facet_grid(cols = vars(%group2%))

filename <- file.path(my_outdir, 'Class_composition.png')
ggsave(filename, class_bar, dpi = 300)

#### Plot - Elemental comp bar####

el_bar <- el_comp %>% 
  group_by(%group1%, %group2%, Element)  %>% 
  summarise(Count = mean(Count), .groups = 'drop') %>% 
  ggplot(aes(x = %group1%,
             y = Count,
             fill = Element)) +
  geom_col() +
  theme_bw() +
  scale_fill_brewer(palette = 'Accent') +
  labs(title = 'Elemental Composition') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
  facet_grid(cols = vars('%group2%'))

filename <- file.path(my_outdir, 'Elemental_composition.png')
ggsave(filename, el_bar, dpi = 300)
