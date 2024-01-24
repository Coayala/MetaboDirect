# ***************************************************************
# 
# MetaboDirect
# Data Exploration step
# MetaboDirect version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
# 
# ***************************************************************

# Loading libraries ----

library(ggvenn)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(vegan)
library(UpSetR)
library(tidyverse)

# Defining variables ----

current_dir <- '%currentdir%'
output_dir <- '%outdir%'
group1 <- '%group1%'
group2 <- '%group2%'

setwd(current_dir)

# Loading custom functions ----
source(file.path(output_dir, 'custom_functions.R'))

class_colors <- get_palette(palette = 'Set3', k = 9)
names(class_colors) <- c(classification$Class, 'Other')


# Import data ----

## Defining paths
my_data.file <- file.path(output_dir, '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_metadata.file <- file.path('%metadata%')
my_elcomp.file <- file.path(output_dir, '1_preprocessing_output', 'elemental_composition.csv')
my_classcomp.file <- file.path(output_dir, '1_preprocessing_output', 'class_composition.csv')
my_outdir <- file.path(output_dir, '3_exploratory')

## Loading tables
df <- read_csv(my_data.file)
metadata <- read_csv(my_metadata.file)
el_comp <- read_csv(my_elcomp.file)
class_comp <- read_csv(my_classcomp.file)

# Reformat data files ----

## Change all metadata columns to factors to avoid problems at plotting

metadata <- metadata %>% 
  mutate(across(!SampleID, as.factor))

## Intensity data file
df_longer <- df %>%
  pivot_longer(metadata$SampleID, names_to = 'SampleID', values_to = 'NormIntensity') %>% 
  filter(NormIntensity != 0) %>% 
  left_join(metadata, by = 'SampleID')

## Class composition
class_comp <- class_comp %>%
  pivot_longer(!SampleID, names_to = 'Class', values_to = 'Count') %>% 
  left_join(metadata, by = 'SampleID')

## Elemental composition
el_comp <- el_comp %>%
  pivot_longer(!SampleID, names_to = 'Element', values_to = 'Count') %>% 
  left_join(metadata, by = 'SampleID')

# Plotting ----

## setting colors for grouping variable ----

my_colors <- set_names(get_palette('Dark2', length(unique(metadata[[group1]]))),
                       nm = unique(metadata[[group1]]))

## Setting list of indexes to plot ----
thermo_idx <- set_names(c('GFE', 'AI_mod', 'DBE', 'NOSC'),
                        nm = c('GFE', 'AI_mod', 'DBE', 'NOSC'))

## Van Krevelen plot ----
vk_plots <- map(thermo_idx, ~ plot_van_krevelen(df_longer, color_by = .x, facet_by = group1) + 
                  scale_color_viridis_c(direction = -1))

plot_filenames <- file.path(my_outdir, paste0('1.', 1:4, '_van_krevelen_', thermo_idx,'.png'))
walk2(vk_plots, plot_filenames, ~ggsave(.y, .x, dpi = 300, width = 12, height = 6))

## Density plot ----
density_plots <- map(thermo_idx, ~ plot_density(df_longer,
                                                index = .x,
                                                color_by = group1,
                                                facet_col = group1,
                                                facet_row = group2) +
                       scale_fill_manual(values = my_colors))


plot_filenames <- file.path(my_outdir, paste0('2.', 1:4, '_density_plot_', thermo_idx,'.png'))
walk2(density_plots, plot_filenames, ~ggsave(.y, .x, dpi = 300, width = 12, height = 6))

## Violin plot ----

violin_plots <- map(thermo_idx, ~ plot_violin(df_longer,
                                              index = .x,
                                              color_by = group1,
                                              facet_by = group2,
                                              title = .x) +
                      scale_fill_manual(values = my_colors))

plot_filenames <- file.path(my_outdir, paste0('3.', 1:4, '_violin_plot_', thermo_idx,'.png'))
walk2(violin_plots, plot_filenames, ~ggsave(.y, .x, dpi = 300, width = 12, height = 6))

## Weighted Thermodynamical indices ----

# Magnitude-weighted values

wt_df <- map(thermo_idx, ~ calculate_weighted(df_longer, index = .x))
wt_violins <- imap(wt_df, ~ plot_violin(df = .x,
                                        index = 'value',
                                        color_by = group1,
                                        facet_by = 'weighted_idx',
                                        title = .y) +
                     scale_fill_manual(values = my_colors))

plot_filenames <- file.path(my_outdir, paste0('4.', 1:4, '_weighted_index_', thermo_idx,'.png'))
walk2(wt_violins, plot_filenames, ~ggsave(.y, .x, dpi = 300, width = 12, height = 6))

## Plot - Class comp bar ----

class_bar <- plot_comp_bar(class_comp, composition = 'Class',
                           group = c(group1, group2), title = 'Molecular Class')

filename <- file.path(my_outdir, '5_Composition_by_class.png')
ggsave(filename, class_bar, dpi = 300, width = 8, height = 8)

## Plot - Elemental comp bar ----

el_bar <- plot_comp_bar(el_comp, composition = 'Element',
                        group = c(group1, group2), title = 'Elemental Composition')

filename <- file.path(my_outdir, '6_Composition_by_element.png')
ggsave(filename, el_bar, dpi = 300, width = 8, height = 8)

# Comparisons ----

groups <- c(group1, group2)

walk(groups, function(g){
  
  comparison_dir <- file.path(my_outdir, paste0('Comparisons-', g))
  if(!dir.exists(comparison_dir)) dir.create(comparison_dir)
  
  df_group <- df_longer %>% 
    select(Mass, HC, OC, Presence = all_of(g)) %>% 
    distinct() %>% 
    group_by(Mass, HC, OC) %>% 
    arrange(Presence) %>% 
    summarise(Presence = paste0(Presence, collapse = ','))
  
  vk_plot <- plot_van_krevelen(df_group, color_by = 'Presence')
  
  filename <- file.path(comparison_dir, paste0('1_vk_', g, '_all_data.png'))
  ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)
  
  ## Upset plot all features group 1
  
  group_values <- unique(df_longer[[g]])
  mass_list <- map(group_values, ~select_masses(df_longer, group = g, value = .x))
  names(mass_list) <- group_values
  
  filename <- file.path(comparison_dir, paste0('2_upset_plot_', g, '_all_data.png'))
  png(filename, width = 2400, height = 2400, res = 300)
  print(upset(fromList(mass_list), order.by = "freq"))
  dev.off()
  
  ## Pairwise comparisons
  
  comb <- combn(unique(as.character(df_longer[[g]])), 2)
  comb_colors <- set_names(c(get_palette('Dark2', length(unique(df_longer[[g]]))), 'gray80'),
                         nm = c(unique(as.character(df_longer[[g]])), 'Shared'))
  
  walk2(comb[1,], comb[2,], function(val1, val2){
    filt_colors <- comb_colors[c(val1, val2, 'Shared')]
    
    df_comparisons <- df_longer %>% 
      rename(Presence = all_of(g)) %>% 
      filter(Presence %in% c(val1, val2)) %>% 
      group_by(Mass, HC, OC, GFE, Class) %>% 
      arrange(Presence) %>% 
      summarise(Presence = paste0(unique(Presence), collapse = ',')) %>% 
      mutate(Presence = ifelse(str_detect(Presence, ','), 'Shared', Presence),
             Presence = factor(Presence, levels = c(val1, val2, 'Shared'))) %>% 
      ungroup()
    
    vk_comb <- plot_van_krevelen(df_comparisons, color_by = 'Presence') +
      scale_color_manual(values = filt_colors)
    
    violin_comb <- plot_violin(df_comparisons, index = 'GFE', color_by = 'Presence', 
                               title = 'Comparing GFE') +
      scale_fill_manual(values = filt_colors)
    
    class_comb <- df_comparisons %>% 
      group_by(Presence) %>% 
      count(Class, name = 'Count') %>% 
      plot_comp_bar(df = ., composition = 'Class', group = 'Presence', 
                    title = 'Molecular Class')
    
    combined <- ggarrange(vk_comb, ggarrange(violin_comb, class_comb, nrow = 2, align = 'v'),
                          ncol = 2, widths = c(1.5, 1)) %>% 
      annotate_figure(., top = text_grob(paste0(val1, ' vs ', val2), face = 'bold', size = 20))
    
    filename <- file.path(comparison_dir, paste0('combined_graph_', val1, '_vs_', val2, '.png'))
    ggsave(filename, combined, dpi = 300, width = 18, height = 8)
    
    # Getting masses
    mass_list <- map(c(val1, val2), ~select_masses(df_longer, group = g, value = .x))
    names(mass_list) <- c(val1, val2)
    
    # Upset
    filename <- file.path(comparison_dir, paste0('upset_plot_', val1, '_vs_', val2, '.png'))
    png(filename, width = 2400, height = 2400, res = 300)
    print(upset(fromList(mass_list), order.by = "freq"))
    dev.off()
    
    venn_plot <- ggvenn(mass_list, fill_color = as.character(filt_colors[1:2])) +
      labs(title = 'Number of metabolites (assigned molecular formula)') +
      theme(plot.title = element_text(face = 'bold', hjust = 0.5))
    
    filename <- file.path(comparison_dir, paste0('venn_diagram_', val1, '_vs_', val2, '.png'))
    ggsave(filename, venn_plot, dpi = 300, width = 7, height = 5, bg = 'white')
    
  })
  
})


