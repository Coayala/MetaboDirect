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
  library(pheatmap)
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

setwd('%currentdir%')

my_data.file <- file.path('%outdir%', '6_transformations', 'network_summary_statistics.csv')
my_transf_counts.file <- file.path('%outdir%', '6_transformations', 'Transformations_summary_counts.csv')
my_metadata.file <- file.path('%metadata%')
my_outdir <- file.path('%outdir%', '6_transformations')

#### Import data ####

network_data <-  read_csv(my_data.file, col_types = cols())
transf_counts <- read_csv(my_transf_counts.file, col_types = cols())
metadata <- read_csv(my_metadata.file, col_types = cols())

#### Reformat data files ####

### Change all metadata columns to factors to avoid problems at plotting

for(i in 2:ncol(metadata)){
  metadata[,i] <- factor(as_vector(metadata[,i])) 
}

network_data_longer <- network_data %>%
  pivot_longer(!SampleID, names_to = "Parameters", values_to = "Value") %>% 
  left_join(metadata, by = 'SampleID')

transf_counts_wide <- transf_counts %>% 
  select(SampleID, Transformation, Perc_Counts) %>% 
  pivot_wider(names_from = 'Transformation', values_from = 'Perc_Counts') %>%
  column_to_rownames(var = 'SampleID') %>% 
  t()

transf_counts_wide[is.na(transf_counts_wide)] <- 0

#### Heatmap of percentage of transformations ####

col_annot <- metadata %>% 
  select(SampleID, %group1%, %group2%) %>% 
  column_to_rownames(var = 'SampleID')

row_annot <- transf_counts %>% 
  select(Transformation, Group) %>% 
  distinct() %>% 
  column_to_rownames(var = 'Transformation')

mapcolor <- colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)[100:1]

filename <- filename <- file.path(my_outdir, 'percentage_transformations.pdf')
pdf(filename, height = 10, width = 8)
pheatmap(transf_counts_wide,
         annotation_col = col_annot,
         annotation_row = row_annot,
         color = mapcolor,
         main = 'Transformations (%) per sample')
dev.off()


#### Plot stats summary ####

my_colors <- get_palette('Dark2', length(unique(network_data_longer$%group1%)))
names(my_colors) <- unique(network_data_longer$%group1%) 

group1 <- '%group1%'
group1_s <- syms(group1)

group2 <- '%group2%'
group2_s <- syms(group2)

net_stats_plot <- network_data_longer %>% 
  group_by(%group1%, Parameters) %>% 
  summarise(Mean =  mean(Value), sd = sd(Value), .groups = 'drop') %>% 
  ggplot(aes(x = (!!! group1_s),
             y = Mean,
             fill = (!!! group1_s))) +
  geom_boxplot(color = 'black') +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  labs(title = 'Network statistics summary',
       x = group1,
       fill = group1) +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title.y = element_blank()) +
  facet_wrap(~Parameters, scales = "free_y", ncol = 3) 

filename <- file.path(my_outdir, 'network_statistics_summary.png')
ggsave(filename, net_stats_plot, dpi = 300)

network_parameters <- unique(network_data_longer$Parameters)

for(param in network_parameters){
  param_df <- network_data_longer %>% 
    filter(Parameters == param)
  
  kw <- param_df %>% 
    group_by(%group2%) %>% 
    kruskal_test(Value ~ %group1%)
  filename <- file.path(my_outdir, paste0(param, '_kruskal_test_by_', group1, '.csv'))
  write_csv(kw, filename)
  
  gh <- param_df %>% 
    group_by(%group2%) %>% 
    games_howell_test(Value ~ %group1%) %>% 
    add_y_position()
  
  param_plot <- ggplot(param_df,
                       aes(x = (!!! group1_s),
                           y = Value,
                           fill = (!!! group1_s))) +
    geom_boxplot() +
    scale_fill_manual(values = my_colors) +
    stat_pvalue_manual(gh,
                       label = 'p.adj.signif', inherit.aes = FALSE,
                       hide.ns = TRUE) +
    theme_bw() +
    labs(title = param,
         x = group1,
         fill = group1) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.title.y = element_blank()) +
    facet_grid(cols = vars(%group2%))
  
  filename <- file.path(my_outdir, paste0(param, '_by_', group1, '.png'))
  ggsave(filename, param_plot, dpi = 300)
}