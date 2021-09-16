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
  library(ggpubr)
  library(rstatix)
  library(KEGGREST)
  library(vegan)
  library(UpSetR)
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

### Change all metadata columns to factors to avoid problems at plotting

for(i in 2:ncol(metadata)){
  metadata[,i] <- factor(as_vector(metadata[,i])) 
}

# Intensity data file
df_longer <- df %>%
  pivot_longer(metadata$SampleID, names_to = 'SampleID', values_to = 'NormIntensity') %>% 
  filter(NormIntensity != 0)
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

## Set samples color names

my_colors <- get_palette('Dark2', length(unique(df_longer$%group1%)))
names(my_colors) <- unique(df_longer$%group1%) 

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
    scale_color_manual(values = get_palette('d3', 8)) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
    facet_wrap(%group1% ~ %group2%)
  
  filename <- file.path(my_outdir, paste0('vk_', var, '.png'))
  ggsave(filename, vk_plot, dpi = 300, height = 8, width = 15)
}


#### Plot - Density Diagram ####

for(var in c('GFE', 'AI_mod', 'DBE')){
  var_s = syms(var)
  density_plot <- ggplot(df_longer,
                         aes(x = (!!! var_s),
                             fill = %group1%)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = my_colors) +
    theme_bw() +
    labs(title = paste(var, 'density plot'),
         x = var) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
    facet_wrap(%group1% ~ %group2%)
  
  filename <- file.path(my_outdir, paste0('density_plot_', var, '.png'))
  ggsave(filename, density_plot, dpi = 300, width = 8, height = 8)
}
#### Plot - Violin ####

for(var in c('GFE', 'AI_mod', 'DBE')){
  
  var_s = syms(var)
  formula <- formula(paste0(var, ' ~ %group1%'))
  stat.test <- df_longer %>% 
    group_by(%group2%) %>% 
    tukey_hsd(formula) %>% 
    add_y_position(step.increase = 1)
  
  violin_plot <- ggplot(df_longer,
                        aes(x = %group1%,
                            y = (!!! var_s),
                            fill = %group1%)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.2, show.legend = F) +
    scale_fill_manual(values = my_colors) +
    theme_bw() +
    stat_pvalue_manual(stat.test,
                       label = 'p.adj.signif', inherit.aes = FALSE,
                       hide.ns = TRUE) +
    labs(title = paste(var, 'violin plot'),
         y = var) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
    facet_grid(cols = vars(%group2%))
  
  filename <- file.path(my_outdir, paste0('violin_plot_', var, '.png'))
  ggsave(filename, violin_plot, dpi = 300, width = 8, height = 8)
}

#### Weighted Thermodynamical indices ####

# Magnitude-weighted values

for(var in c('GFE', 'AI_mod', 'DBE')){
  
  var_s = syms(var)
  newnames <- paste0(var, c('_weighted', '_magniture_average'))
  
  weighted_df <- df_longer %>%  
    mutate(w = (!!! var_s) * NormIntensity) %>% 
    select(Mass, SampleID, NormIntensity, !!! var_s, w) %>% 
    group_by(SampleID) %>% 
    add_tally(w, name = 'sum_w') %>% 
    add_tally(NormIntensity, name = 'sum_intensity') %>% 
    ungroup() %>% 
    mutate(m_avg = sum_w/sum_intensity) %>% 
    select(SampleID, !!! var_s, w, m_avg) %>% 
    rename(!!newnames[1] := w, !!newnames[2] := m_avg) %>% 
    pivot_longer(-SampleID, names_to = 'weighted_idx', values_to = 'value') %>% 
    left_join(metadata, by = 'SampleID')
  
  stat.test <- weighted_df %>%
    group_by(weighted_idx) %>% 
    tukey_hsd(value ~ %group1%) %>% 
    add_y_position(scales = 'free_y', step.increase = 1)
  
  weighted_plot <- ggplot(weighted_df) +
    geom_boxplot(aes(x = %group1%,
                     y = value,
                     fill = %group1%)) +
    stat_pvalue_manual(stat.test,
                       label = 'p.adj.signif', inherit.aes = FALSE,
                       hide.ns = TRUE) +
    facet_wrap(~weighted_idx, scales = 'free_y') +
    theme_bw() +
    labs(title = paste0(var, ' - Weighted abundance'),
         caption = '"_weighted" index are calculated by multiplying the index (i.e. GFE, DBE or IA) times
         the normalized intensity per each of the detected masses.
         "_magniture average" index are calculated by dividing
         the total sum of the "_weighted" index by the total sum of the normalized 
         intensities per each of the samples.') +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
  
  filename <- file.path(my_outdir, paste0('weighted_abundance_',var, '.png'))
  ggsave(filename, weighted_plot, dpi = 300, width = 8, height = 8)
  
}

#### Plot - Class comp bar####

class_colors <- brewer.pal(length(classification$Class) + 1, name = 'Set3')
names(class_colors) <- c(classification$Class, 'Other')

class_bar <- class_comp %>% 
  group_by(%group1%, %group2%, Class)  %>% 
  summarise(Count = mean(Count, na.rm = TRUE)) %>%
  mutate(Perc_count = Count/sum(Count, na.rm = TRUE)*100) %>%
  ggplot(aes(x = %group1%,
             y = Perc_count,
             fill = Class)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = class_colors) +
  labs(title = 'Class Composition',
       y = 'Percentage') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
  facet_grid(cols = vars(%group2%))

filename <- file.path(my_outdir, 'Composition_by_class.png')
ggsave(filename, class_bar, dpi = 300, width = 8, height = 8)

#### Plot - Elemental comp bar####

element_colors <- get_palette(palette = 'Accent', length(unique(el_comp$Element)))
names(element_colors) <- unique(el_comp$Element)

el_bar <- el_comp %>% 
  group_by(%group1%, %group2%, Element)  %>% 
  summarise(Count = mean(Count, na.rm = TRUE)) %>% 
  mutate(Perc_count = Count/sum(Count, na.rm = TRUE)*100) %>%
  ggplot(aes(x = %group1%,
             y = Perc_count,
             fill = Element)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = element_colors) +
  labs(title = 'Elemental Composition',
       y = 'Percentage') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5)) +
  facet_grid(cols = vars(%group2%))

filename <- file.path(my_outdir, 'Composition_by_element.png')
ggsave(filename, el_bar, dpi = 300, width = 8, height = 8)

#### Comparisons ####

### Grouping variable 1

group1 <- '%group1%'
group1_s <- syms(group1)

comparison_dir <- file.path(my_outdir, paste0('Comparisons-', group1))
dir.create(comparison_dir)

## Group comparison

df_w <- df_longer %>% 
  pivot_wider(names_from = all_of(group1), values_from = all_of(group1), names_prefix = 'GRP_') %>% 
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
  scale_color_brewer(palette = 'Paired', guide = guide_legend(order = 2)) +
  new_scale_color() +
  class_rect  +
  scale_color_brewer(palette = 'Dark2', guide = guide_legend(order = 1)) +
  labs(x = 'O:C',
       y = 'H:C',
       title = paste0('Van Krevelen Diagram by ', group1)) +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

filename <- file.path(comparison_dir, paste0('vk_', group1, '_all_data.png'))
ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)

## Upset plot all features group 1

group_list <- list()
i <- 1

for(val in unique(df_longer$'%group1%')){
  list_mass <- df_longer %>%
    filter((!!! group1_s) == val) %>% 
    pull(Mass)
  group_list[i] <- list(list_mass)
  i <- i + 1
}

names(group_list) <- unique(df_longer$'%group1%')

filename <- file.path(comparison_dir, paste0('upset_plot_', group1, '_all_data.png'))
png(filename, width = 2400, height = 2400, res = 300)
upset(fromList(group_list), order.by = "freq")
dev.off()

## Pairwise comparisons

# Get all combinations possible for group 1

comb_g1 <- combn(unique(df_longer$'%group1%'), 2)

my_colors <- c(get_palette('d3', length(unique(df_longer$'%group1%'))), 'gray80')
color_names <- c(unique(as.character(df_longer$'%group1%')), 'Shared')
names(my_colors) <- color_names

for(i in 1:ncol(comb_g1)){
  val1 <- as.character(comb_g1[1, i])
  val2 <- as.character(comb_g1[2, i])
  
  df_comb <- df_longer %>% 
    filter((!!! group1_s) == val1 | (!!! group1_s) == val2) %>% 
    pivot_wider(names_from = all_of(group1), values_from = all_of(group1), names_prefix = 'GRP_') %>% 
    select(Mass, HC, OC, Class, GFE, DBE, contains('GRP_')) %>% 
    group_by(Mass) %>% 
    fill(contains('GRP_'), .direction = 'downup')  %>% 
    distinct() %>% 
    unite(Presence, contains('GRP_'), sep = ', ', na.rm = TRUE) %>% 
    mutate(Presence = ifelse(str_detect(Presence, ', '), 'Shared', Presence))
  
  vk_plot <- ggplot(df_comb,
                    aes(x = OC,
                        y = HC,
                        color = Presence)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = my_colors, guide = guide_legend(order = 1)) +
    new_scale_color() +
    class_rect  +
    scale_color_brewer(palette = 'Dark2', guide = guide_legend(order = 1)) +
    labs(x = 'O:C',
         y = 'H:C',
         title = 'Van Krevelen Diagram') +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          legend.position = 'bottom')
  
  stat.test <- df_comb %>%
    ungroup() %>% 
    tukey_hsd(GFE ~ Presence) %>% 
    add_y_position(step.increase = 1)
  
  
  violin_plot <- ggplot(df_comb,
                        aes(x = Presence,
                            y = GFE,
                            fill = Presence)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.size = 0.2, show.legend = F) +
    scale_fill_manual(values = my_colors) +
    theme_bw() +
    labs(title = 'GFE') +
    stat_pvalue_manual(stat.test,
                       label = 'p.adj.signif', inherit.aes = FALSE,
                       hide.ns = TRUE) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          legend.position = 'None',
          axis.title.x = element_blank())
  
  class_bar <- df_comb %>% 
    group_by(Presence)  %>% 
    count(Class) %>%
    mutate(count_perc=n/sum(n)*100) %>% 
    ggplot(aes(x = Presence,
               y = count_perc,
               fill = Class)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values = class_colors) +
    labs(title = 'Class Composition') +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.title = element_blank(),
          legend.position = 'bottom')
  
  combined_plot <- ggarrange(vk_plot, ggarrange(class_bar, violin_plot, nrow =2, align = 'v'), ncol = 2, widths = c(1.5, 1))
  combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0(val1, ' vs ', val2), face = 'bold', size = 20))
  
  filename <- file.path(comparison_dir, paste0('combined_graph_', val1, '_vs_', val2, '.png'))
  ggsave(filename, combined_plot, dpi = 300, width = 18, height = 8)
  
  ## Upset plot
  
  list1 <- df_longer %>% 
    filter((!!! group1_s) == val1) %>% 
    pull(Mass)
  
  list2 <- df_longer %>% 
    filter((!!! group1_s) == val2) %>% 
    pull(Mass)
  
  group_list <- list(list1, list2)
  names(group_list) <- c(val1, val2)
  
  filename <- file.path(comparison_dir, paste0('upset_plot_', val1, '_vs_', val2, '.png'))
  png(filename, width = 2400, height = 2400, res = 300)
  print(upset(fromList(group_list), order.by = "freq"))
  dev.off()
}

## Grouping variable 2

group2 <- '%group2%'

if(group2 != 'NULL'){
group2_s <- syms(group2)

comparison_dir <- file.path(my_outdir, paste0('Comparisons-', group2))
dir.create(comparison_dir)


  df_w <- df_longer %>% 
    pivot_wider(names_from = all_of(group2), values_from = all_of(group2), names_prefix = 'GRP_') %>% 
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
         title = paste0('Van Krevelen Diagram by ', group2)) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
  
  filename <- file.path(comparison_dir, paste0('VK_', group2, '_all_data.png'))
  ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)
  
  ## Upset plot all features group 2
  
  group_list <- list()
  i <- 1
  
  for(val in unique(df_longer$'%group2%')){
    list_mass <- df_longer %>%
      filter((!!! group2_s) == val) %>% 
      pull(Mass)
    group_list[i] <- list(list_mass)
    i <- i + 1
  }
  
  names(group_list) <- unique(df_longer$'%group2%')
  
  filename <- file.path(comparison_dir, paste0('upset_plot_', group2, '_all_data.png'))
  png(filename, width = 2400, height = 2400, res = 300)
  print(upset(fromList(group_list), order.by = "freq"))
  dev.off()
  
  
  ## Pairwise comparisons
  
  # Get all combinations possible for group 2
  
  comb_g2 <- combn(unique(df_longer$'%group2%'), 2)
  
  my_colors <- c(get_palette('d3', length(unique(df_longer$'%group2%'))), 'gray80')
  color_names <- c(unique(as.character(df_longer$'%group2%')), 'Shared')
  names(my_colors) <- color_names
  
  
  for(i in 1:ncol(comb_g2)){
    val1 <- as.character(comb_g2[1, i])
    val2 <- as.character(comb_g2[2, i])
    
    df_comb <- df_longer %>% 
      filter((!!! group2_s) == val1 | (!!! group2_s) == val2) %>% 
      pivot_wider(names_from = all_of(group2), values_from = all_of(group2), names_prefix = 'GRP_') %>% 
      select(Mass, HC, OC, Class, GFE, DBE, contains('GRP_')) %>% 
      group_by(Mass) %>% 
      fill(contains('GRP_'), .direction = 'downup')  %>% 
      distinct() %>% 
      unite(Presence, contains('GRP_'), sep = ', ', na.rm = TRUE) %>% 
      mutate(Presence = ifelse(str_detect(Presence, ', '), 'Shared', Presence))
    
    vk_plot <- ggplot(df_comb,
                      aes(x = OC,
                          y = HC,
                          color = Presence)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = my_colors, guide = guide_legend(order = 2)) +
      new_scale_color() +
      class_rect  +
      scale_color_brewer(palette = 'Dark2', guide = guide_legend(order = 1)) +
      labs(x = 'O:C',
           y = 'H:C',
           title = 'Van Krevelen Diagram') +
      theme(plot.title = element_text(face = 'bold', hjust = 0.5),
            legend.position = 'bottom')
    
    stat.test <- df_comb %>%
      ungroup() %>% 
      tukey_hsd(GFE ~ Presence) %>% 
      add_y_position(step.increase = 1)
    
    
    violin_plot <- ggplot(df_comb,
                          aes(x = Presence,
                              y = GFE,
                              fill = Presence)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.1, outlier.size = 0.2, show.legend = F) +
      scale_fill_manual(values = my_colors) +
      theme_bw() +
      labs(title = 'GFE') +
      stat_pvalue_manual(stat.test,
                         label = 'p.adj.signif', inherit.aes = FALSE,
                         hide.ns = TRUE) +
      theme(plot.title = element_text(face = 'bold', hjust = 0.5),
            legend.position = 'None',
            axis.title.x = element_blank())
    
    class_bar <- df_comb %>% 
      group_by(Presence)  %>% 
      count(Class) %>%
      mutate(count_perc=n/sum(n)*100) %>% 
      ggplot(aes(x = Presence,
                 y = count_perc,
                 fill = Class)) +
      geom_col() +
      theme_bw() +
      scale_fill_manual(values = class_colors) +
      labs(title = 'Class Composition') +
      theme(plot.title = element_text(face = 'bold', hjust = 0.5),
            axis.title = element_blank(),
            legend.position = 'bottom')
    
    combined_plot <- ggarrange(vk_plot, ggarrange(class_bar, violin_plot, nrow =2, align = 'v'), ncol = 2, widths = c(1.5, 1))
    combined_plot <- annotate_figure(combined_plot, top = text_grob(paste0(val1, ' vs ', val2), face = 'bold', size = 20))
    
    filename <- file.path(comparison_dir, paste0('combined_graph_', val1, '_vs_', val2, '.png'))
    ggsave(filename, combined_plot, dpi = 300, width = 18, height = 8)
    
    ## Upset plot
    
    list1 <- df_longer %>% 
      filter((!!! group2_s) == val1) %>% 
      pull(Mass)
    
    list2 <- df_longer %>% 
      filter((!!! group2_s) == val2) %>% 
      pull(Mass)
    
    group_list <- list(list1, list2)
    names(group_list) <- c(val1, val2)
    
    filename <- file.path(comparison_dir, paste0('upset_plot_', val1, '_vs_', val2, '.png'))
    png(filename, width = 2400, height = 2400, res = 300)
    print(upset(fromList(group_list), order.by = "freq"))
    dev.off()
  }
  
}

