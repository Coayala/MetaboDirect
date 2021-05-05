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
  library(ggnewscale)
  require(ggpubr)
  require(rstatix)
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

setwd('%currentdir%')

my_data.file <- file.path('%outdir%', '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_metadata.file <- file.path('%metadata%')
my_elcomp.file <- file.path('%outdir%', '1_preprocessing_output', 'elemental_composition.csv')
my_classcomp.file <- file.path('%outdir%', '1_preprocessing_output', 'class_composition.csv')
classification.file <- file.path('%Metabo_HOME%', 'data', 'compound_class_table.csv')
my_outdir <- file.path('%outdir%', '3_exploratory')

#### Import data ####

df <-  read_csv(my_data.file, col_types = cols())
metadata <- read_csv(my_metadata.file, col_types = cols())
el_comp <- read_csv(my_elcomp.file, col_types = cols())
class_comp <- read_csv(my_classcomp.file, col_types = cols())
classification <- read_csv(classification.file, col_types = cols())

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

## Compound class rectangles (for plotting)

class_rect <-  geom_rect(data = classification,
                         aes(xmin = OC_low,
                             xmax = OC_high,
                             ymin = HC_low,
                             ymax = HC_high,
                             color = Class),
                         fill = NA,
                         size = 1.5,
                         inherit.aes = FALSE, 
                         linetype = 'dashed') 

for(var in c('GFE', 'AI_mod', 'DBE')){
  var_s = syms(var)
  vk_plot <- ggplot(df_longer,
                    aes(x = OC,
                        y = HC,
                        color = (!!! var_s))) +
    geom_point() +
    theme_bw() +
    scale_color_viridis_c(direction = -1) +
    labs(x = 'O:C',
         y = 'H:C',
         title = paste0('Van Krevelen Diagram by ', var),
         color = var) +
    new_scale_color() +
    class_rect  +
    scale_color_brewer(palette = 'Dark2') +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
    facet_wrap(%group1% ~ %group2%)
  
  filename <- file.path(my_outdir, paste0('vk_', var, '.png'))
  ggsave(filename, vk_plot, dpi = 300, width = 8)
}

## Group comparisons

# Grouping variable 1

df_w <- df_longer %>% 
  pivot_wider(values_from = '%group1%', names_from = '%group1%', names_prefix = 'GRP_') %>% 
  select(Mass, HC, OC, contains('GRP_')) %>% 
  group_by(Mass) %>% 
  fill(contains('GRP_'), .direction = 'downup')  %>% 
  distinct() %>% 
  unite(Presence, contains('GRP_'), sep = ', ', na.rm = TRUE)

vk_plot <- ggplot(df_w,
                  aes(x = OC,
                      y = HC,
                      color = Presence)) +
  geom_point() +
  theme_bw() +
  scale_color_brewer(palette = 'Paired') +
  new_scale_color() +
  class_rect  +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'O:C',
       y = 'H:C',
       title = paste0('Van Krevelen Diagram by %group1%')) +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(my_outdir, 'comparison_vk_%group1%.png')
ggsave(filename, vk_plot, dpi = 300, width = 8)

# Grouping variable 2

group2 <- '%group2%'

if(!is.null(group2)){
  df_w <- df_longer %>% 
    pivot_wider(values_from = '%group2%', names_from = '%group2%', names_prefix = 'GRP_') %>% 
    select(Mass, HC, OC, contains('GRP_')) %>% 
    group_by(Mass) %>% 
    fill(contains('GRP_'), .direction = 'downup')  %>% 
    distinct() %>% 
    unite(Presence, contains('GRP_'), sep = ', ', na.rm = TRUE)
  
  vk_plot <- ggplot(df_w,
                    aes(x = OC,
                        y = HC,
                        color = Presence)) +
    geom_point() +
    theme_bw() +
    scale_color_brewer(palette = 'Paired') +
    new_scale_color() +
    class_rect  +
    scale_color_brewer(palette = 'Dark2') +
    labs(x = 'O:C',
         y = 'H:C',
         title = paste0('Van Krevelen Diagram by %group2%')) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
  
  filename <- file.path(my_outdir, 'comparison_vk_%group2%.png')
  ggsave(filename, vk_plot, dpi = 300, width = 8)
}


#### Plot - Density Diagram ####

for(var in c('GFE', 'AI_mod', 'DBE')){
  var_s = syms(var)
  density_plot <- ggplot(df_longer,
                         aes(x = (!!! var_s),
                             fill = %group1%)) +
    geom_density(alpha = 0.4) +
    scale_fill_brewer(palette = 'Dark2') +
    theme_bw() +
    labs(title = paste(var, 'density plot'),
         x = var) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
    facet_wrap(%group1% ~ %group2%)
  
  filename <- file.path(my_outdir, paste0('density_plot_', var, '.png'))
  ggsave(filename, density_plot, dpi = 300)
}
#### Plot - Violin ####

for(var in c('GFE', 'AI_mod', 'NOSC')){
  
  var_s = syms(var)
  formula <- formula(paste0(var, ' ~ Habitat'))
  stat.test <- df_longer %>% 
    group_by(%group2%) %>% 
    tukey_hsd(formula) %>% 
    add_y_position(step.increase = 1)
  
  violin_plot <- ggplot(df_longer,
                        aes(x = %group1%,
                            y = (!!! var_s),
                            fill = %group1%)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.size = 0.2, show.legend = F) +
    scale_fill_brewer(palette = 'Dark2') +
    theme_bw() +
    stat_pvalue_manual(stat.test,
                       label = 'p.adj.signif', inherit.aes = FALSE,
                       hide.ns = TRUE) +
    labs(title = paste(var, 'violin plot'),
         y = var) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
    facet_grid(cols = vars(%group2%))
  
  filename <- file.path(my_outdir, paste0('violin_plot_', var, '.png'))
  ggsave(filename, violin_plot, dpi = 300, width = 8)
}

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
  facet_grid(cols = vars(%group2%))

filename <- file.path(my_outdir, 'Elemental_composition.png')
ggsave(filename, el_bar, dpi = 300)
