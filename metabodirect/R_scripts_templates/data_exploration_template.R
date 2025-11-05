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
library(jsonlite)
library(tidyverse)

# Defining variables ----

current_dir <- '%currentdir%'
output_dir <- '%outdir%'
group1 <- '%group1%'


setwd(current_dir)

# Loading custom functions ----
source(file.path(output_dir, 'custom_functions.R'))

class_colors <- get_palette(palette = 'Set3', k = 9)
names(class_colors) <- c(classification$Class, 'Other')

# Import data ----

## Defining paths
my_data.file <- file.path(output_dir, '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_metadata.file <- file.path('data/processed/fticr_metadata.csv')
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

# Initializing JSON with figure metadata

figure_metadata <- list()

# Plotting ----

## setting colors for grouping variable ----

my_colors <- set_names(get_palette('Dark2', length(unique(metadata[[group1]]))),
                       nm = unique(metadata[[group1]]))

## Setting list of indexes to plot ----
thermo_idx <- set_names(c('GFE', 'AI_mod', 'DBE', 'NOSC'),
                        nm = c('GFE', 'AI_mod', 'DBE', 'NOSC'))

## Van Krevenlen plot ----
vk_plots <- map(thermo_idx, function(idx){
  res <- plot_van_krevelen(df_longer, color_by = idx, facet_col = group1) 
  
  res$plot <- res$plot + 
    scale_color_viridis_c(direction = -1)
  
  return(res)
})

plot_filenames <- file.path(my_outdir, paste0('1.', 1:4, '_van_krevelen_', thermo_idx,'.png'))
walk2(vk_plots, plot_filenames, function(vk, file){
  ggsave(file, vk$plot, dpi = 300, width = 12, height = 6)
  
  add_figure_metadata(figure_id = str_remove(basename(file), '.png'),
                      figure_title = vk$plot$labels$title,
                      figure_type = 'van Krevelen Diagram',
                      caption = 'Van Krevelen diagram showing H/C vs O/C for all detected compounds',
                      description = vk$meta$description,
                      insights = vk$meta$insight,
                      x_axis_label = 'O:C',
                      y_axis_label = 'H:C',
                      grouping_variables = c(group1),
                      data_source = list(
                        input_file = my_data.file,
                        columns_used = c('OC',
                                         'HC',
                                         str_extract(vk$plot$labels$title, 'GFE|DBE|AI_mod|NOSC'), 
                                         group1)
                      ),
                      r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                      analysis_category = 'Exploratory',
                      figure_info = list(
                        figure_file = file,
                        figure_width = 12,
                        figure_height = 6,
                        dpi = 300
                      ))
  
})

## Density plot ----
density_plots <- map(thermo_idx, function(idx){
  res <- plot_density(df_longer,
                      index = idx,
                      color_by = group1,
                      facet_col = group1) 
  
  res$plot <- res$plot +
    scale_fill_manual(values = my_colors)
  
  return(res)
})


plot_filenames <- file.path(my_outdir, paste0('2.', 1:4, '_density_plot_', thermo_idx,'.png'))
walk2(density_plots, plot_filenames, function(dens, file){
  ggsave(file, dens$plot, dpi = 300, width = 12, height = 6)
  
  add_figure_metadata(figure_id = str_remove(basename(file), '.png'),
                      figure_title = dens$plot$labels$title,
                      figure_type = 'Density plot',
                      caption = 'Density plot of thermodynamic index values',
                      description = dens$meta$description,
                      insights = dens$meta$insight,
                      x_axis_label = str_extract(dens$plot$labels$title, 'GFE|DBE|AI_mod|NOSC'),
                      y_axis_label = 'density',
                      grouping_variables = c(group1),
                      data_source = list(
                        input_file = my_data.file,
                        columns_used = c(str_extract(dens$plot$labels$title, 'GFE|DBE|AI_mod|NOSC'), 
                                         group1)
                      ),
                      r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                      analysis_category = 'Exploratory',
                      figure_info = list(
                        figure_file = file,
                        figure_width = 12,
                        figure_height = 6,
                        dpi = 300
                      ))
  
})

## Violin plot ----

violin_plots <- map(thermo_idx, function(idx){
  res <- plot_violin(df_longer,
                     index = idx,
                     color_by = group1,
                     title = idx) 
  
  res$plot <- res$plot +
    scale_fill_manual(values = my_colors)
  
  return(res)
})

plot_filenames <- file.path(my_outdir, paste0('3.', 1:4, '_violin_plot_', thermo_idx,'.png'))
walk2(violin_plots, plot_filenames, function(viol, file){
  ggsave(file, viol$plot, dpi = 300, width = 12, height = 6)
  
  add_figure_metadata(figure_id = str_remove(basename(file), '.png'),
                      figure_title = viol$plot$labels$title,
                      figure_type = 'Violin and boxplot',
                      caption = 'Boxplot and violin plots showing differences between thermodynamic indexes among sample groups',
                      description = viol$meta$description,
                      insights = viol$meta$insight,
                      x_axis_label = group1,
                      y_axis_label = str_extract(viol$plot$labels$title, 'GFE|DBE|AI_mod|NOSC'),
                      grouping_variables = c(group1),
                      data_source = list(
                        input_file = my_data.file,
                        columns_used = c(str_extract(viol$plot$labels$title, 'GFE|DBE|AI_mod|NOSC'), 
                                         group1)
                      ),
                      r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                      analysis_category = 'Exploratory',
                      figure_info = list(
                        figure_file = file,
                        figure_width = 12,
                        figure_height = 6,
                        dpi = 300
                      ))
})

## Weighted Thermodynamical indices ----

# Magnitude-weighted values

wt_df <- map(thermo_idx, ~ calculate_weighted(df_longer, index = .x))
wt_violins <- imap(wt_df, function(idx_df, idx){
  res <- plot_violin(df = idx_df,
                     index = 'value',
                     color_by = group1,
                     facet_by = 'weighted_idx',
                     title = idx,
                     calculate_stat_signif = FALSE) 
  
  res$plot <- res$plot +
    scale_fill_manual(values = my_colors)
  
  return(res)
})

plot_filenames <- file.path(my_outdir, paste0('4.', 1:4, '_weighted_index_', thermo_idx,'.png'))
walk2(wt_violins, plot_filenames, function(viol, file){
  ggsave(file, viol$plot, dpi = 300, width = 12, height = 6)
  
  add_figure_metadata(figure_id = str_remove(basename(file), '.png'),
                      figure_title = viol$plot$labels$title,
                      figure_type = 'Violin and boxplot',
                      caption = 'Boxplot and violin plots showing weighted and magnitude-averaged thermodynamic indexes',
                      description = viol$meta$description,
                      insights = viol$meta$insight,
                      x_axis_label = group1,
                      y_axis_label = 'Thermodynamic index value',
                      grouping_variables = c(group1, 'weighted_index'),
                      data_source = list(
                        input_file = my_data.file,
                        columns_used = c(str_extract(viol$plot$labels$title, 'GFE|DBE|AI_mod|NOSC'), 
                                         group1)
                      ),
                      r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                      analysis_category = 'Exploratory',
                      figure_info = list(
                        figure_file = file,
                        figure_width = 12,
                        figure_height = 6,
                        dpi = 300
                      ))
})

## Plot - Class comp bar ----

class_bar <- plot_comp_bar(class_comp, composition = 'Class',
                           group = group1, title = 'Molecular Class')

filename <- file.path(my_outdir, '5_Composition_by_class.png')
ggsave(filename, class_bar$plot, dpi = 300, width = 8, height = 8)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = class_bar$plot$labels$title,
                    figure_type = 'Composition bar',
                    caption = 'Stacked bar plot showing the molecular composition of the detected masses',
                    description = class_bar$meta$description,
                    insights = class_bar$meta$insight,
                    x_axis_label = group1,
                    y_axis_label = 'Relative abundance (percentage)',
                    grouping_variables = group1,
                    data_source = list(
                      input_file = my_classcomp.file,
                      columns_used = c('Class',
                                       'Perc_count',
                                       group1)
                    ),
                    r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                    analysis_category = 'Exploratory',
                    figure_info = list(
                      figure_file = filename,
                      figure_width = 8,
                      figure_height = 8,
                      dpi = 300
                    ))


## Plot - Elemental comp bar ----

el_bar <- plot_comp_bar(el_comp, composition = 'Element',
                        group = group1, title = 'Elemental Composition')

filename <- file.path(my_outdir, '6_Composition_by_element.png')
ggsave(filename, el_bar$plot, dpi = 300, width = 8, height = 8)

add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                    figure_title = el_bar$plot$labels$title,
                    figure_type = 'Composition bar',
                    caption = 'Stacked bar plot showing the elemental composition of the detected masses',
                    description = el_bar$meta$description,
                    insights = el_bar$meta$insight,
                    x_axis_label = group1,
                    y_axis_label = 'Relative abundance (percentage)',
                    grouping_variables = group1,
                    data_source = list(
                      input_file = my_elcomp.file,
                      columns_used = c('Element',
                                       'Perc_count',
                                       group1)
                    ),
                    r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                    analysis_category = 'Exploratory',
                    figure_info = list(
                      figure_file = filename,
                      figure_width = 8,
                      figure_height = 8,
                      dpi = 300
                    ))

# Comparisons ----

groups <- group1

walk(groups, function(g){
  
  comparison_dir <- file.path(my_outdir, paste0('Comparisons-', g))
  if(!dir.exists(comparison_dir)) dir.create(comparison_dir)
  
  df_group <- df_longer %>% 
    select(Mass, HC, OC, Presence = all_of(g), Class) %>% 
    distinct() %>% 
    group_by(Mass, HC, OC, Class) %>% 
    arrange(Presence) %>% 
    summarise(Presence = paste0(Presence, collapse = ','))
  
  vk_plot <- plot_van_krevelen(df_group, color_by = 'Presence')
  
  filename <- file.path(comparison_dir, paste0('1_vk_', g, '_all_data.png'))
  ggsave(filename, vk_plot$plot, dpi = 300, width = 8, height = 8)
  
  add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                      figure_title = vk_plot$plot$labels$title,
                      figure_type = 'van Krevelen Diagram',
                      caption = 'Van Krevelen diagram showing H/C vs O/C for all detected compounds',
                      description = vk_plot$meta$description,
                      insights = vk_plot$meta$insight,
                      x_axis_label = 'O:C',
                      y_axis_label = 'H:C',
                      grouping_variables = 'none',
                      data_source = list(
                        input_file = my_data.file,
                        columns_used = c('OC',
                                         'HC',
                                         'Presence')
                      ),
                      r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                      analysis_category = 'Exploratory',
                      figure_info = list(
                        figure_file = file,
                        figure_width = 8,
                        figure_height = 8,
                        dpi = 300
                      ))
  
  ## Upset plot all features group 1
  
  group_values <- unique(df_longer[[g]])
  mass_list <- map(group_values, ~select_masses(df_longer, group = g, value = .x))
  names(mass_list) <- group_values
  
  upset_p <- plot_upset(mass_list, g)
  
  filename <- file.path(comparison_dir, paste0('2_upset_plot_', g, '_all_data.png'))
  png(filename, width = 10, height = 8, res = 300, unit = 'in')
  print(upset_p$plot)
  dev.off()
  
  add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                      figure_title = paste0('Upset plot - ', g),
                      figure_type = 'Upset plot',
                      caption = 'Upset plot shows intersections between sample sets',
                      description = upset_p$meta$description,
                      insights = upset_p$meta$insight,
                      x_axis_label = 'Intersection Size',
                      y_axis_label = 'Sample set',
                      grouping_variables = 'none',
                      data_source = list(
                        input_file = my_data.file,
                        columns_used = c('Mass')
                      ),
                      r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                      analysis_category = 'Exploratory',
                      figure_info = list(
                        figure_file = file,
                        figure_width = 10,
                        figure_height = 8,
                        dpi = 300
                      ))
  
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
    
    vk_comb <- plot_van_krevelen(df_comparisons, color_by = 'Presence') 
    
    vk_comb$plot <- vk_comb$plot +
      scale_color_manual(values = filt_colors)
    
    violin_comb <- plot_violin(df_comparisons, index = 'GFE', color_by = 'Presence', 
                               title = 'Comparing GFE') 
    
    violin_comb$plot <- violin_comb$plot +
      scale_fill_manual(values = filt_colors)
    
    class_comb <- df_comparisons %>% 
      group_by(Presence) %>% 
      count(Class, name = 'Count') %>% 
      plot_comp_bar(df = ., composition = 'Class', group = 'Presence', 
                    title = 'Molecular Class')
    
    combined <- ggarrange(vk_comb$plot, ggarrange(violin_comb$plot, class_comb$plot, nrow = 2, align = 'v'),
                          ncol = 2, widths = c(1.5, 1)) %>% 
      annotate_figure(., top = text_grob(paste0(val1, ' vs ', val2), face = 'bold', size = 20))
    
    filename <- file.path(comparison_dir, paste0('combined_graph_', val1, '_vs_', val2, '.png'))
    ggsave(filename, combined, dpi = 300, width = 18, height = 8)
    
    add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                        figure_title = paste(c(val1, 'vs', val2), collapse = ' '),
                        figure_type = 'Combined figure',
                        caption = 'Upset plot shows intersections between sample sets',
                        description = list(vk_comb$meta$description,
                                           violin_comb$meta$description,
                                           class_comb$meta$description),
                        insights = list(vk_comb$meta$insight,
                                        violin_comb$meta$insight,
                                        class_comb$meta$insight),
                        x_axis_label = list('O:C', 'Presence', 'Presence'),
                        y_axis_label = list('H:C', 'GFE', 'Relative abundance (percentage)'),
                        grouping_variables = g,
                        data_source = list(
                          input_file = c(my_data.file, my_classcomp.file),
                          columns_used = c('Mass', 'OC', 'HC', 'GFE', g)
                        ),
                        r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                        analysis_category = 'Exploratory',
                        figure_info = list(
                          figure_file = file,
                          figure_width = 18,
                          figure_height = 8,
                          dpi = 300
                        ))
    
    # Getting masses
    mass_list <- map(c(val1, val2), ~select_masses(df_longer, group = g, value = .x))
    names(mass_list) <- c(val1, val2)
    
    # Upset
    upset_comb <- plot_upset(mass_list, g)
    
    filename <- file.path(comparison_dir, paste0('upset_plot_', val1, '_vs_', val2, '.png'))
    png(filename, width = 10, height = 8, res = 300, unit = 'in')
    print(upset_comb$plot)
    dev.off()
    
    add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                        figure_title = paste0('Upset plot - ', val1, ' and ', val2),
                        figure_type = 'Upset plot',
                        caption = 'Upset plot shows intersections between sample sets',
                        description = upset_comb$meta$description,
                        insights = upset_comb$meta$insight,
                        x_axis_label = 'Intersection Size',
                        y_axis_label = 'Sample set',
                        grouping_variables = 'none',
                        data_source = list(
                          input_file = my_data.file,
                          columns_used = c('Mass')
                        ),
                        r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                        analysis_category = 'Exploratory',
                        figure_info = list(
                          figure_file = file,
                          figure_width = 10,
                          figure_height = 8,
                          dpi = 300
                        ))
    
    venn_plot <- plot_venn(mass_list, filt_colors, val1, val2)
    
    filename <- file.path(comparison_dir, paste0('venn_diagram_', val1, '_vs_', val2, '.png'))
    ggsave(filename, venn_plot$plot, dpi = 300, width = 7, height = 5, bg = 'white')
    
    add_figure_metadata(figure_id = str_remove(basename(filename), '.png'),
                        figure_title = paste0('Venn diagram - ', val1, ' and ', val2),
                        figure_type = 'Venn diagram',
                        caption = 'Venn diagram showing metabolites shared between two samples',
                        description = venn_plot$meta$description,
                        insights = venn_plot$meta$insight,
                        x_axis_label = '',
                        y_axis_label = '',
                        grouping_variables = 'none',
                        data_source = list(
                          input_file = my_data.file,
                          columns_used = c('Mass',
                                           g)
                        ),
                        r_script_path = file.path(output_dir, '3_exploratory', 'data_exploration.R'),
                        analysis_category = 'Exploratory',
                        figure_info = list(
                          figure_file = file,
                          figure_width = 7,
                          figure_height = 5,
                          dpi = 300
                        ))
    
  })
  
})

write_json(figure_metadata, file.path(output_dir, 'index.json'))