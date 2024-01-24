# ***************************************************************
# 
# MetaboDirect
# Custom functions
# MetaboDirect version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
# 
# ***************************************************************

# Additional data ----

## Molecular class table ----
classification <- tribble(
  ~Class, ~OC_low, ~OC_high, ~HC_low, ~HC_high,
  'Lipid', 0, 0.3, 1.5, 2.5,
  'Unsat. HC', 0, 0.125, 0.8, 1.5,
  'Cond. HC', 0, 0.95, 0.2, 0.8,
  'Protein', 0.3, 0.55, 1.5, 2.3,
  'Amino sugar', 0.55, 0.7, 1.5, 2.2,
  'Carbohydrate', 0.7, 1.5, 1.5, 2.5,
  'Lignin', 0.125, 0.65, 0.8, 1.5,
  'Tannin', 0.65, 1.1, 0.8, 1.5, 
) %>% 
  mutate(label_x = (OC_low + OC_high) / 2,
         label_y = HC_high - 0.1)

## Compound class rectangles (for plotting of Van Krevelen diagrams) ----

class_rect <-  geom_rect(data = classification,
                         aes(xmin = OC_low,
                             xmax = OC_high,
                             ymin = HC_low,
                             ymax = HC_high),
                         color = 'black',
                         fill = NA,
                         linewidth = 1,
                         inherit.aes = FALSE, 
                         linetype = 'dashed')
rect_label <- geom_label(data = classification,
                         aes(x = label_x,
                             y = label_y,
                             label = Class),
                         inherit.aes = FALSE,
                         size = 3)

# Themes ----

custom_theme <- function(angle_x = 0){
  theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5),
          axis.text.x = element_text(angle = angle_x, hjust = 1))
}

# Plotting functions ----

# ******************************************************************************
plot_van_krevelen <- function(df, color_by, facet_by = NA){
  color_by <- color_by
  p <- ggplot(df,
              aes(x = OC,
                  y = HC,
                  color = !!rlang::ensym(color_by))) +
    geom_point() +
    labs(x = 'O:C',
         y = 'H:C',
         title = paste0('Van Krevelen Diagram by ', color_by)) +
    class_rect  +
    rect_label +
    custom_theme()
  
  if(!is.na(facet_by)){
    facet_formula <- as.formula(paste0('~ ', facet_by))
    p <- p +
      facet_wrap(facet_formula)
  }
  
  return(p)
}

# ******************************************************************************
plot_density <- function(df, index, color_by, facet_col, facet_row = NA){
  index <- index
  color_by <- color_by
  p <- ggplot(df_longer,
              aes(x = !!rlang::ensym(index),
                  fill = !!rlang::ensym(color_by))) +
    geom_density(alpha = 0.7) +
    theme_bw() +
    labs(title = paste(index, 'density plot'),
         x = index) +
    custom_theme()
  
  if(is.na(facet_row)){
    facet_formula <- as.formula(paste0('~ ', facet_col))
    p <- p +
      facet_wrap(facet_formula)
  } else {
    facet_formula <- as.formula(paste0(facet_row, ' ~ ', facet_col))
    p <- p +
      facet_grid(facet_formula)
  }
  
  return(p)
}

# ******************************************************************************
plot_violin <- function(df, index, color_by, facet_by = NA, title){
  index <- index
  color_by <- color_by
  
  # Testing normality
  df_test <- df %>% 
    select(Mass, all_of(index)) %>% 
    distinct()
  pval <- try(shapiro.test(df_test[[index]])$p.value)
  
  # Comparing means
  stat_formula <- as.formula(paste0(index, '~ ', color_by))
  if(!is.na(facet_by)){
    stat_df <- df %>% 
      group_by(!!rlang::ensym(facet_by))
  } else stat_df <- df
  
  if(pval < 0.05 || str_detect(pval, 'Error')){
    stat_df <- stat_df %>% 
      dunn_test(stat_formula) %>% 
      add_significance() %>% 
      add_xy_position(step.increase = 1, scales = 'free_y')
  } else {
    stat_df <- stat_df %>% 
      tukey_hsd(stat_formula) %>% 
      add_significance() %>% 
      add_xy_position(step.increase = 1, scales = 'free_y')
  }
  
  p <- ggplot(df,
              aes(x = !!rlang::ensym(color_by),
                  y = !!rlang::ensym(index),
                  fill = !!rlang::ensym(color_by))) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.2, show.legend = F) +
    theme_bw() +
    labs(title = title,
         x = color_by) +
    custom_theme(angle_x = 45)
  
  if(!is.na(facet_by)){
    facet_formula <- as.formula(paste0('~ ', facet_by))
    p <- p +
      stat_pvalue_manual(data = stat_df,
                         label = 'p.adj.signif', 
                         inherit.aes = FALSE,
                         hide.ns = TRUE) +
      facet_wrap(facet_formula)
  } else {
    p <- p +
      stat_pvalue_manual(data = stat_df,
                         label = 'p.adj.signif', 
                         inherit.aes = FALSE,
                         hide.ns = TRUE)
  }
  
  return(p)
}

# ******************************************************************************
plot_comp_bar <- function(df, composition, group, title){
  if(length(group) == 2){
    group_a <- group[1]
    group_b <- group[2]
    df <- df %>% 
      group_by(!!rlang::ensym(group_a), 
               !!rlang::ensym(group_b), 
               !!rlang::ensym(composition))
  } else {
    group_a <- group[1]
    df <- df %>% 
      group_by(!!rlang::ensym(group_a), 
               !!rlang::ensym(composition))
  }
  
  if(composition == 'Class'){
    colors <- class_colors
  } else {
    colors <- get_palette('Set3', length(unique(df$Element)))
  }
  
  plot <- df %>% 
    summarise(Count = mean(Count, na.rm = TRUE)) %>%
    mutate(Perc_count = Count/sum(Count, na.rm = TRUE)*100) %>%
    ggplot(aes(x = !!rlang::ensym(group_a),
               y = Perc_count,
               fill = !!rlang::ensym(composition))) +
    geom_col(color = 'black') +
    scale_fill_manual(values = colors) +
    labs(title = title,
         y = 'Percentage') +
    custom_theme()
  
  if(length(group) == 2){
    facet_formula <- as.formula(paste0('~ ', group_b))
    plot <- plot + 
      facet_wrap(facet_formula, scales = 'free_x')
  }
  return(plot)
}

# ******************************************************************************
plot_ordination <- function(df, x, y, color_by, col_vec, title, 
                            shape_by = NA, point_size = 2.5,
                            x_label = NA, y_label = NA){
  color_by <- color_by
  shape_by <- shape_by
  if(is.na(shape_by)){
    plot <- ggplot(df) +
      geom_point(aes(x = {{ x }},
                     y = {{ y }},
                     color = !!rlang::ensym(color_by)),
                 size = point_size)
  } else {
    plot <- ggplot(df) +
      geom_point(aes(x = {{ x }},
                     y = {{ y }},
                     color = !!rlang::ensym(color_by),
                     shape = !!rlang::ensym(shape_by)),
                 size = point_size)
  }
  
  plot <- plot + 
    scale_color_manual(values = col_vec) +
    labs(title = title) +
    custom_theme()
  
  return(plot)
}

# ******************************************************************************
plot_diversity_index <- function(df, group_by, title){
  group_by <- group_by
  
  plot <- df %>% 
    ggplot() +
    geom_boxplot(aes(x = !!rlang::ensym(group_by),
                     y = values,
                     fill = !!rlang::ensym(group_by))) +
    scale_fill_manual(values = my_colors) +
    facet_wrap(~index, scales = 'free_y') +
    labs(title = title) +
    custom_theme()
  
  return(plot)
}

# Calculate functions ----

# ******************************************************************************
select_masses <- function(df, group, value){
  group <- group
  df %>% 
    filter(!!rlang::ensym(group) == value) %>% 
    pull(Mass)
}

# ******************************************************************************
calculate_weighted <- function(df, index){
  new_names <- paste0(index, c('_weigthed', '_magnitude_average'))
  names(new_names) <- c('wt', 'm_avg')
  weighted_df <- df %>%  
    select(Mass, SampleID, NormIntensity, all_of(c(index))) %>% 
    mutate(wt = !!rlang::ensym(index) * NormIntensity) %>% 
    group_by(SampleID) %>% 
    mutate(m_avg = sum(wt, na.rm = TRUE)/sum(NormIntensity, na.rm = TRUE)) %>% 
    select(Mass,SampleID, all_of(index), wt, m_avg) %>% 
    rename_with(.cols = c(wt, m_avg), .fn = function(x) new_names[x]) %>% 
    pivot_longer(!c(SampleID, Mass), names_to = 'weighted_idx', values_to = 'value') %>% 
    left_join(metadata, by = 'SampleID') %>% 
    ungroup()
  
  return(weighted_df)
}

# ******************************************************************************
calculate_dist <- function(mat, normalized_with, dm_method = 'AUTO'){
  if(dm_method == 'AUTO'){
    if(norm_method %in% c('max', 'sum', 'none')){
      dm_method = 'bray'
    }else if(norm_method %in% c('mean', 'median', 'z_score', 'minmax')){
      dm_method = 'manhattan'
    }else {
      dm_method = 'jaccard'
    }
  }
  
  print(paste0('Using ', dm_method, ' distance'))
  
  vegdist(mat, method = dm_method)
  
}

# ******************************************************************************
calculate_nmds <- function(mat, normalized_with, 
                           dm_method = 'AUTO', seed = 123){
  dist_mat <- calculate_dist(mat, normalized_with, dm_method)
  
  set.seed(seed)
  nmds.log <- capture.output(nmds <- metaMDS(dist_mat, 
                                             k = 2,
                                             distance = dm_method,
                                             maxit = 999,
                                             trymax = 500,
                                             wascores = TRUE))
  filename <- file.path(my_outdir, '1.1_NMDS.log')
  write_lines(as.data.frame(nmds.log), filename)
  
  filename <- file.path(my_outdir, '1.2_NMDS_stressplot.png')
  
  png(filename, res = 300, height = 3600, width = 3600)
  capture.output(stressplot(nmds))
  dev.off()
  
  nmds_scores <- as.data.frame(scores(nmds, display = 'sites')) %>% 
    rownames_to_column(var = 'SampleID') %>% 
    left_join(metadata, by = 'SampleID')
  
}

# ******************************************************************************
calculate_permanova <- function(mat, normalized_with, variables,
                                dm_method = 'AUTO', seed = 123){
  
  dist_mat <- calculate_dist(mat, normalized_with, dm_method)
  set.seed(seed)
  
  formula <- as.formula(paste0('dist_mat ~', paste(variables, collapse = '+')))
  
  permanova <- adonis2(formula, 
                       data=metadata, 
                       permutations=999, 
                       method=dm.method)
  return(permanova)
}

# ******************************************************************************
calculate_pca <- function(mat){
  
  mat <- mat %>% 
    mutate(across(everything(), ~ifelse(is.na(.x), 0, .x)))
  
  pca <- prcomp(mat, center = TRUE)
  eigen <- get_eigenvalue(pca)
  
  scree_plot <- fviz_eig(pca, addlabels = TRUE) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
  
  # Extract sample coordinates for PC1 and PC2
  pca_coordinates <- as_tibble(pca$x)
  pca_coordinates$SampleID <- rownames(pca$x)
  pca_coordinates <- left_join(pca_coordinates, metadata, by ='SampleID')
  
  var_loadings <- facto_summarize(pca, element = "var", 
                                  result=c("coord","contrib","cos2"), axes=c(1,2))
  
  return(list(coordinates = pca_coordinates,
              eigenvalues = eigen,
              scree_plot = scree_plot,
              loadings = var_loadings))
}

# ******************************************************************************
calculate_richness <- function(mat){
  richness <- specaccum(mat, method = 'random', permutations = 100)
  
  richness_long <- as_tibble(richness$perm, rownames = 'Sites') %>%
    pivot_longer(!Sites, names_to = 'permutation', values_to = 'Richness') %>%
    mutate(Sites = as.numeric(Sites))
  
  return(richness_long)
}

# ******************************************************************************
calculate_rank_abundance <- function(mat){
  rank_abundance <- tibble(peak = colnames(mat),
                           intensity_sums = colSums(mat),
                           presence = colSums(mat != 0)) %>%
    arrange(intensity_sums) %>%
    mutate(position = n():1,
           presence_perc = presence/nrow(mat)*100)
  
  return(rank_abundance)
}

# ******************************************************************************
calculate_diversity_index <- function(mat, norm_mat){
  diversity <- tibble(SampleID = rownames(norm_mat),
                      Shannon = diversity(norm_mat, index = 'shannon'),
                      Gini_Simpson = diversity(norm_mat, index = 'simpson'),
                      Chao1 = estimateR(round(mat))['S.chao1',])  
  
  filename <- file.path(my_outdir, '2.1_abundance_diversity_index.csv')
  write_csv(diversity, filename)
  
  final_df <- diversity %>% 
    pivot_longer(!SampleID, names_to = 'index', values_to = 'values') %>% 
    left_join(metadata, by = 'SampleID')
  
  return(final_df)
}

# ******************************************************************************
calculate_functional_div_index <- function(fun_df, mat){
  # Elemental composition
  el_traits <- fun_df %>% 
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
  
  functional_diversity <- tibble(SampleID = rownames(mat),
                                 Elemental_composition = 
                                   rao.diversity(mat, el_traits)$FunRao,
                                 Reactivity = 
                                   rao.diversity(mat, rx_traits)$FunRao,
                                 Insaturation_and_aromaticity = 
                                   rao.diversity(mat, ia_traits)$FunRao) 
  
  filename <- file.path(my_outdir, '3.1_functional_diversity.csv')
  write_csv(functional_diversity, filename)
  
  final_df <- functional_diversity %>% 
    pivot_longer(!SampleID, names_to = 'index', values_to = 'values') %>% 
    left_join(metadata, by = 'SampleID') %>% 
  
  return(final_df)
  
}

# KEGG function ----

# ******************************************************************************
get_kegg_compound_info <- function(cpd_id, mass, formula){
  cpd_get <- keggGet(cpd_id)
  
  cpd_res <- map(cpd_get, function(entry){
    entry_df <- tibble(Mass = mass,
                       MolecularFormula = formula,
                       KEGG_id = entry$ENTRY,
                       KEGG_name = paste(entry$NAME, collapse = ';'),
                       KEGG_formula = entry$FORMULA,
                       KEGG_pathway = ifelse(!is.null(entry$PATHWAY), 
                                             paste(entry$PATHWAY, collapse = ';'), NA),
                       KEGG_module = ifelse(!is.null(entry$MODULE), 
                                            paste(entry$MODULE, collapse = ';'), NA),
                       KEGG_brite = ifelse(!is.null(entry$BRITE), 
                                           paste(entry$BRITE, collapse = ';'), NA),
                       KEGG_enzyme = ifelse(!is.null(entry$ENZYME), 
                                            paste(entry$ENZYME, collapse = ';'), NA),
                       KEGG_reaction = ifelse(!is.null(entry$REACTION), 
                                              paste(entry$REACTION, collapse = ';'), NA))
    
    return(entry_df)
  })
  
  cpd_df <- reduce(cpd_res, rbind)
  
  return(cpd_df)
}



