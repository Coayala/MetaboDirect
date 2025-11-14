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

# Initializing JSON index of figures

figure_list <- list()

# Plotting ----

## setting colors for grouping variable ----

my_colors <- set_names(get_palette('Dark2', length(unique(metadata[[group1]]))),
                       nm = unique(metadata[[group1]]))

## Setting list of indexes to plot ----
thermo_idx <- set_names(c('GFE', 'AI_mod', 'DBE', 'NOSC'),
                        nm = c('GFE', 'AI_mod', 'DBE', 'NOSC'))

## Van Krevenlen plot ----
vk_plots <- map(thermo_idx, function(idx){
  res <- plot_van_krevelen(df_longer, color_by = idx, 
                           facet_col = group1, facet_row = group2) 
  
  res$plot <- res$plot + 
    scale_color_viridis_c(direction = -1)
  
  return(res)
})

plot_filenames <- file.path(my_outdir, paste0('3.1.', 1:4, '_van_krevelen_', thermo_idx,'.png'))
walk2(vk_plots, plot_filenames, function(plot_obj, file){
  ggsave(file, plot_obj$plot, dpi = 300, width = 12, height = 6)
  
  update_figure_list(figure_id = str_remove(basename(file), '.png'),
                     figure_title = plot_obj$plot$labels$title)
  
  save_figure_metadata(
    plot_obj = plot_obj,  
    figure_id = str_remove(basename(file), '.png'),
    analysis_module = '3. Exploratory',
    figure_type = 'van Krevelen diagram',
    caption = 'Van Krevelen diagram showing H/C vs O/C for all detected compounds',
    figure_file = file,
    group_aes = 'facets',
    modifiable_aesthetics = list('colour', 'facets'),
    data_source = list(data = my_data.file,
                       sample_metadata = my_metadata.file),
    r_script_path = file.path(my_outdir, 'data_exploration.R'),
    functions_used = list(plot = 'plot_van_krevelen()'),
    resolution = 300,
    height = 12,
    width = 6
  )
  
})

## Density plot ----
density_plots <- map(thermo_idx, function(idx){
  res <- plot_density(df_longer,
                      index = idx,
                      color_by = group1,
                      facet_col = group1,
                      facet_row = group2) 
  
  res$plot <- res$plot +
    scale_fill_manual(values = my_colors)
  
  return(res)
})


plot_filenames <- file.path(my_outdir, paste0('3.2.', 1:4, '_density_plot_', thermo_idx,'.png'))
walk2(density_plots, plot_filenames, function(plot_obj, file){
  ggsave(file, plot_obj$plot, dpi = 300, width = 12, height = 6)
  
  update_figure_list(figure_id = str_remove(basename(file), '.png'),
                     figure_title = plot_obj$plot$labels$title)
  
  save_figure_metadata(
    plot_obj = plot_obj,  
    figure_id = str_remove(basename(file), '.png'),
    analysis_module = '3. Exploratory',
    figure_type = 'Density plot',
    caption = 'Density plot of thermodynamic indexes in each sample group',
    figure_file = file,
    group_aes = c('facets', 'fill'),
    modifiable_aesthetics = list('fill', 'facets'),
    data_source = list(data = my_data.file,
                       sample_metadata = my_metadata.file),
    r_script_path = file.path(my_outdir, 'data_exploration.R'),
    functions_used = list(plot = 'plot_density()'),
    resolution = 300,
    height = 12,
    width = 6
  )
  
})


## Violin plot ----

violin_plots <- map(thermo_idx, function(idx){
  res <- plot_violin(df_longer,
                     index = idx,
                     color_by = group1,
                     facet_by = group2,
                     title = idx) 
  
  res$plot <- res$plot +
    scale_fill_manual(values = my_colors)
  
  return(res)
})

plot_filenames <- file.path(my_outdir, paste0('3.3.', 1:4, '_violin_plot_', thermo_idx,'.png'))
walk2(violin_plots, plot_filenames, function(plot_obj, file){
  ggsave(file, plot_obj$plot, dpi = 300, width = 12, height = 6)
  
  update_figure_list(figure_id = str_remove(basename(file), '.png'),
                     figure_title = plot_obj$plot$labels$title)
  
  save_figure_metadata(
    plot_obj = plot_obj,  
    figure_id = str_remove(basename(file), '.png'),
    analysis_module = '3. Exploratory',
    figure_type = 'Violin and box plot',
    caption = paste0('Violin and boxplots showing differences in thermodynamic indexes ',
                     'across sample groups.'),
    figure_file = file,
    group_aes = c('facets', 'fill', 'x'),
    modifiable_aesthetics = list('fill', 'facets'),
    data_source = list(data = my_data.file,
                       sample_metadata = my_metadata.file),
    r_script_path = file.path(my_outdir, 'data_exploration.R'),
    functions_used = list(plot = 'plot_violin()'),
    resolution = 300,
    height = 12,
    width = 6
  )
})

## Weighted Thermodynamical indices ----

# Magnitude-weighted values

wt_df <- map(thermo_idx, ~ calculate_weighted(df_longer, index = .x))
wt_violins <- imap(wt_df, function(idx_df, idx){
  res <- plot_violin(df = idx_df,
                     index = 'Magnitude-averaged/weigthed index',
                     color_by = group1,
                     facet_by = 'weighted_idx',
                     title = idx,
                     calculate_stat_signif = FALSE) 
  
  res$plot <- res$plot +
    scale_fill_manual(values = my_colors)
  
  return(res)
})

plot_filenames <- file.path(my_outdir, paste0('3.4.', 1:4, '_weighted_index_', thermo_idx,'.png'))
walk2(wt_violins, plot_filenames, function(plot_obj, file){
  ggsave(file, plot_obj$plot, dpi = 300, width = 12, height = 6)
  
  update_figure_list(figure_id = str_remove(basename(file), '.png'),
                     figure_title = plot_obj$plot$labels$title)
  
  save_figure_metadata(
    plot_obj = plot_obj,  
    figure_id = str_remove(basename(file), '.png'),
    analysis_module = '3. Exploratory',
    figure_type = 'Violin and box plot',
    caption = paste0('Violin and boxplots showing differences in magnitude-averaged ',
                     'thermodynamic indexes across sample groups.'),
    figure_file = file,
    group_aes = c('fill', 'x'),
    modifiable_aesthetics = list('facets'),
    data_source = list(data = my_data.file,
                       sample_metadata = my_metadata.file),
    r_script_path = file.path(my_outdir, 'data_exploration.R'),
    functions_used = list(calculate = 'calculate_weighted()',
                          plot = 'plot_violin()'),
    resolution = 300,
    height = 12,
    width = 6
  )
})

## Plot - Class comp bar ----

class_bar <- plot_comp_bar(class_comp, composition = 'Class',
                           group = c(group1, group2), title = 'Molecular Class')

filename <- file.path(my_outdir, '3.5_Composition_by_class.png')
ggsave(filename, class_bar$plot, dpi = 300, width = 8, height = 8)

update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                   figure_title = class_bar$plot$labels$title)

save_figure_metadata(
  plot_obj = class_bar,  
  figure_id = str_remove(basename(filename), '.png'),
  analysis_module = '3. Exploratory',
  figure_type = 'Composition bar',
  caption = 'Stacked bar plot showing the molecular class composition of the detected masses',
  figure_file = filename,
  group_aes = c('x'),
  modifiable_aesthetics = list('x'),
  data_source = list(data = my_classcomp.file,
                     sample_metadata = my_metadata.file),
  r_script_path = file.path(my_outdir, 'data_exploration.R'),
  functions_used = list(plot = 'plot_comp_bar()'),
  resolution = 300,
  height = 8,
  width = 8
)


## Plot - Elemental comp bar ----

el_bar <- plot_comp_bar(el_comp, composition = 'Element',
                        group = c(group1, group2), title = 'Elemental Composition')

filename <- file.path(my_outdir, '3.6_Composition_by_element.png')
ggsave(filename, el_bar$plot, dpi = 300, width = 8, height = 8)

update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                   figure_title = el_bar$plot$labels$title)

save_figure_metadata(
  plot_obj = el_bar,  
  figure_id = str_remove(basename(filename), '.png'),
  analysis_module = '3. Exploratory',
  figure_type = 'Composition bar',
  caption = 'Stacked bar plot showing the elemental composition of the detected masses',
  figure_file = filename,
  group_aes = c('x'),
  modifiable_aesthetics = list('x'),
  data_source = list(data = my_elcomp.file,
                     sample_metadata = my_metadata.file),
  r_script_path = file.path(my_outdir, 'data_exploration.R'),
  functions_used = list(plot = 'plot_comp_bar()'),
  resolution = 300,
  height = 8,
  width = 8
)

# Comparisons ----

groups <- c(group1, group2)

walk(groups, function(g){
  
  counter <- which(groups == g)
  prefix <- paste0('3.7.', counter)
  
  comparison_dir <- file.path(my_outdir, paste0(prefix, '_Comparisons-', g))
  if(!dir.exists(comparison_dir)) dir.create(comparison_dir)
  
  df_group <- df_longer %>% 
    select(Mass, HC, OC, Presence = all_of(g), Class) %>% 
    distinct() %>% 
    group_by(Mass, HC, OC, Class) %>% 
    arrange(Presence) %>% 
    summarise(Presence = paste0(Presence, collapse = ','))
  
  vk_plot <- plot_van_krevelen(df_group, color_by = 'Presence')
  
  filename <- file.path(comparison_dir, paste0(prefix, '.1_vk_', g, '_all_data.png'))
  ggsave(filename, vk_plot$plot, dpi = 300, width = 12, height = 6)
  
  update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                     figure_title = vk_plot$plot$labels$title)
  
  save_figure_metadata(
    plot_obj = vk_plot,  
    figure_id = str_remove(basename(filename), '.png'),
    analysis_module = '3. Exploratory',
    figure_type = 'van Krevelen diagram',
    caption = 'Van Krevelen diagram comparing which metabolites where found in which sample groups',
    figure_file = filename,
    group_aes = 'colour',
    modifiable_aesthetics = list('colour'),
    data_source = list(data = my_data.file,
                       sample_metadata = my_metadata.file),
    r_script_path = file.path(my_outdir, 'data_exploration.R'),
    functions_used = list(plot = 'plot_van_krevelen()'),
    resolution = 300,
    height = 12,
    width = 6
  )
  
  ## Upset plot all features group 1
  
  group_values <- unique(df_longer[[g]])
  mass_list <- map(group_values, ~select_masses(df_longer, group = g, value = .x))
  names(mass_list) <- group_values
  
  upset_p <- plot_upset(mass_list, g)
  
  filename <- file.path(comparison_dir, paste0(prefix, '.2_upset_plot_', g, '_all_data.png'))
  png(filename, width = 12, height = 6, res = 300, unit = 'in')
  print(upset_p$plot)
  dev.off()
  
  update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                     figure_title = paste0('Upset plot - ', g))
  
  save_figure_metadata(
    plot_obj = upset_p,  
    figure_id = str_remove(basename(filename), '.png'),
    analysis_module = '3. Exploratory',
    figure_type = 'Upset plot',
    caption = 'Plot shows intersections between sample sets as vertical bars',
    figure_file = filename,
    figure_title = paste0('Upset plot - ', g),
    x_axis_label = 'none',
    y_axis_label = 'Intersection size',
    group_aes = list(none = 'none'),
    modifiable_aesthetics = list('none'),
    has_legend = FALSE,
    data_source = list(data = my_data.file,
                       sample_metadata = my_metadata.file),
    r_script_path = file.path(my_outdir, 'data_exploration.R'),
    functions_used = list(plot = 'plot_upset()'),
    resolution = 300,
    height = 12,
    width = 6
  )
  
  ## Pairwise comparisons
  
  comb <- combn(unique(as.character(df_longer[[g]])), 2)
  comb_colors <- set_names(c(get_palette('Dark2', length(unique(df_longer[[g]]))), 'gray80'),
                           nm = c(unique(as.character(df_longer[[g]])), 'Shared'))
  
  ii <- 3
  walk2(comb[1,], comb[2,], function(val1, val2){
    
    prefix_i <- paste0(prefix, '.', ii)
    
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
    
    combined <- ggarrange(vk_comb$plot, 
                          ggarrange(violin_comb$plot, class_comb$plot, 
                                    nrow = 2, align = 'v', labels = c('B', 'C')),
                          ncol = 2, widths = c(1.5, 1), labels = c('A', '')) %>% 
      annotate_figure(., top = text_grob(paste0(val1, ' vs ', val2), face = 'bold', size = 20))
    
    filename <- file.path(comparison_dir, paste0(prefix_i, '.1_combined_graph_', val1, '_vs_', val2, '.png'))
    ggsave(filename, combined, dpi = 300, width = 18, height = 8)
    
    update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                       figure_title = paste0(val1, ' vs ', val2))
    
    save_figure_metadata(
      plot_obj = list(plot = combined,
                      meta = list(description = list(A = vk_comb$meta$description,
                                                     B = violin_comb$meta$description,
                                                     C = class_comb$meta$description),
                                  insight = list(vk_comb$meta$insight,
                                                 violin_comb$meta$insight,
                                                 class_comb$meta$insight))),  
      figure_id = str_remove(basename(filename), '.png'),
      analysis_module = '3. Exploratory',
      figure_type = 'Combined plot',
      caption = paste0('Plot showing the comparison between two sample groups: ',
                       'A) Van Krevelen diagram of metabolite presence; ',
                       'B) Violin plot of GFE; ',
                       'C) Molecular class composition bars'),
      figure_file = filename,
      figure_title = paste0(val1, ' vs ', val2),
      x_axis_label = list(A = 'O:C', B = 'Presence', C = 'Presence'),
      y_axis_label = list(A = 'H:C', B = 'GFE', C = 'Percentage'),
      group_aes = list(A = list(color = 'Presence'),
                       B = list(fill = 'Presence',
                                x = 'Presence'),
                       C = list(x = 'Presence')),
      modifiable_aesthetics = list('none'),
      has_legend = FALSE,
      custom_legend = list(A = list(color = 'Presence'),
                           B = list(fill = 'Presence'),
                           C = list(fill = 'Class')),
      data_source = list(data = my_data.file, 
                         data2 = my_classcomp.file,
                         sample_metadata = my_metadata.file),
      r_script_path = file.path(my_outdir, 'data_exploration.R'),
      functions_used = list(plot = list(A = 'plot_van_krevelent()',
                                        B = 'plot_violin()',
                                        C = 'plot_comp_bar()')),
      resolution = 300,
      height = 18,
      width = 8
    )
    
    # Getting masses
    mass_list <- map(c(val1, val2), ~select_masses(df_longer, group = g, value = .x))
    names(mass_list) <- c(val1, val2)
    
    # Upset
    upset_comb <- plot_upset(mass_list, g)
    
    filename <- file.path(comparison_dir, paste0(prefix_i, '.2_upset_plot_', 
                                                 val1, '_vs_', val2, '.png'))
    png(filename, width = 10, height = 8, res = 300, unit = 'in')
    print(upset_comb$plot)
    dev.off()
    
    update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                       figure_title = paste0('Upset plot - ', val1, ' vs ', val2))
    
    save_figure_metadata(
      plot_obj = upset_comb,  
      figure_id = str_remove(basename(filename), '.png'),
      analysis_module = '3. Exploratory',
      figure_type = 'Upset plot',
      caption = 'Plot shows intersections between sample sets as vertical bars',
      figure_file = filename,
      figure_title = paste0('Upset plot - ', val1, ' vs ', val2),
      x_axis_label = 'none',
      y_axis_label = 'Intersection size',
      group_aes = list(color = g),
      modifiable_aesthetics = list('none'),
      has_legend = FALSE,
      data_source = list(data = my_data.file,
                         sample_metadata = my_metadata.file),
      r_script_path = file.path(my_outdir, 'data_exploration.R'),
      functions_used = list(plot = 'plot_upset()'),
      resolution = 300,
      height = 10,
      width = 8
    )
    
    venn_plot <- plot_venn(mass_list, filt_colors, val1, val2)
    
    filename <- file.path(comparison_dir, paste0(prefix_i, '.3_venn_diagram_', 
                                                 val1, '_vs_', val2, '.png'))
    ggsave(filename, venn_plot$plot, dpi = 300, width = 7, height = 5, bg = 'white')
    
    update_figure_list(figure_id = str_remove(basename(filename), '.png'),
                       figure_title = paste0('Venn plot - ', val1, ' vs ', val2))
    
    save_figure_metadata(
      plot_obj = venn_plot,  
      figure_id = str_remove(basename(filename), '.png'),
      analysis_module = '3. Exploratory',
      figure_type = 'Venn diagram',
      caption = 'Plot shows intersections between sample sets as vertical bars',
      figure_file = filename,
      figure_title = paste0('Venn plot - ', val1, ' vs ', val2),
      x_axis_label = 'none',
      y_axis_label = 'none',
      group_aes = list(none = 'none'),
      modifiable_aesthetics = list('none'),
      has_legend = FALSE,
      data_source = list(data = my_data.file,
                         sample_metadata = my_metadata.file),
      r_script_path = file.path(my_outdir, 'data_exploration.R'),
      functions_used = list(plot = 'plot_venn()'),
      resolution = 300,
      height = 12,
      width = 6
    )
    
    ii<<- ii +1
    
  })
  
})

diag_index <- read_json(file.path(output_dir, "index.json"))

index_figures <- list(
  diag_index,
  list(
    analysis_module = '3. Exploratory',
    script_name = file.path(fix_paths(my_outdir), 'data_exploration.R'),
    last_run = Sys.time(),
    figures = figure_list
  )
) 

write_json(index_figures, 
           file.path(output_dir, 'index.json'), 
           auto_unbox = TRUE, 
           pretty = TRUE)