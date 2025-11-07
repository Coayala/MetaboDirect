# ***************************************************************
# 
# MetaboDirect
# Custom functions
# MetaboDirect for agent
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
plot_van_krevelen <- function(df, color_by, facet_col = NA, facet_row = NA){
  color_by <- color_by
  
  group_vars <- c(facet_col, facet_row)
  group_vars <- group_vars[!is.na(group_vars)]
  
  if(is.na(facet_col) && is.na(facet_row)){
    df_plot <- df %>% 
      ungroup() %>% 
      select(Mass, OC, HC, Class, all_of(color_by)) %>% 
      distinct()
  } else if(is.na(facet_row)) {
    
    df_plot <- df %>% 
      ungroup() %>% 
      select(Mass, OC, HC, Class, all_of(c(color_by, facet_col))) %>% 
      distinct()
    
  } else {
    df_plot <- df %>% 
      ungroup() %>% 
      select(Mass, OC, HC, Class, all_of(c(color_by, facet_col, facet_row))) %>% 
      distinct()
  }
  
  p <- df_plot %>% 
    ggplot(aes(x = OC,
               y = HC,
               color = !!rlang::ensym(color_by))) +
    geom_point() +
    labs(x = 'O:C',
         y = 'H:C',
         title = paste0('Van Krevelen Diagram by ', color_by)) +
    class_rect  +
    rect_label +
    custom_theme()
  
  if(!is.na(facet_col) && is.na(facet_row)){
    facet_formula <- as.formula(paste0('~ ', facet_col))
    p <- p +
      facet_wrap(facet_formula)
  } else if(!is.na(facet_row)) {
    facet_formula <- as.formula(paste0(facet_row, '~ ', facet_col))
    p <- p +
      facet_grid(facet_formula)
  }
  
  is_num <- is.numeric(df_plot[[color_by]])
  
  vk_meta <- list(description = descriptions('van_krevelen',
                                             color_by = color_by,
                                             group = group_vars),
                  insight = insights(contains = c('van_krevelen', color_by),
                                     df_plot = df_plot,
                                     color_by = color_by,
                                     facet_col = facet_col,
                                     facet_row = facet_row,
                                     color_continuous = is_num))
  
  return(list(plot = p,
              meta = vk_meta))
}

# ******************************************************************************
plot_density <- function(df, index, color_by, facet_col, facet_row = NA){
  index <- index
  color_by <- color_by
  
  group_vars <- c(facet_col, facet_row)
  group_vars <- group_vars[!is.na(group_vars)]
  
  if(is.na(facet_row)){
    df_plot <- df %>% 
      ungroup() %>% 
      select(Mass,
             all_of(c(index, color_by,
                      facet_col))) %>% 
      distinct()
  } else {
    df_plot <- df %>% 
      ungroup() %>% 
      select(Mass,
             all_of(c(index, color_by,
                      facet_col, facet_row))) %>% 
      distinct()
  }
  
  p <- ggplot(df_plot,
              aes(x = !!rlang::ensym(index),
                  fill = !!rlang::ensym(color_by))) +
    geom_density(alpha = 0.7) +
    theme_bw() +
    labs(title = paste(index, 'density plot'),
         x = index,
         y = 'Density') +
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
  
  density_meta <- list(description = descriptions('density',
                                                  thermo_idx = index,
                                                  group = group_vars),
                       insight = insights(contains = c('density', index),
                                          plot_data = ggplot_build(p),
                                          facet_col = facet_col,
                                          facet_row = facet_row,
                                          index = index))
  
  return(list(plot = p,
              meta = density_meta))
}

# ******************************************************************************
plot_violin <- function(df, index, color_by, facet_by = NA, title,
                        calculate_stat_signif = TRUE){
  index <- index
  color_by <- color_by
  
  group_vars <- c(color_by, facet_by)
  group_vars <- group_vars[!is.na(group_vars)]
  
  # Testing normality
  df_test <- df %>% 
    ungroup() %>% 
    select(Mass, all_of(index)) %>% 
    distinct()
  pval <- try(shapiro.test(df_test[[index]])$p.value)
  
  # Comparing means
  if(!is.na(facet_by)){
    df_plot <- df %>% 
      ungroup() %>%
      select(Mass, all_of(c(index, color_by, facet_by))) %>% 
      distinct() %>% 
      group_by(across(all_of(facet_by)))
  } else {
    df_plot <- df %>% 
      ungroup() %>%
      select(Mass, all_of(c(index, color_by))) %>% 
      distinct()
  }
  
  if(calculate_stat_signif){
    stat_formula <- as.formula(paste0(index, '~ ', color_by))
    if(pval < 0.05 || str_detect(pval, 'Error')){
      stat_df <- df_plot %>% 
        dunn_test(stat_formula) %>% 
        add_significance() %>% 
        add_xy_position(step.increase = 1, scales = 'free_y')
    } else {
      stat_df <- df_plot %>% 
        tukey_hsd(stat_formula) %>% 
        add_significance() %>% 
        add_xy_position(step.increase = 1, scales = 'free_y')
    }
    
    
  } else {
    stat_df <- data.frame(group1 = NA,
                          group2 = NA,
                          p.adj = 1,
                          p.adj.signif = 'ns',
                          y.position = NA)
  }
  
  p <- df_plot %>% 
    ggplot(aes(x = .data[[color_by]],
               y = .data[[index]],
               fill = .data[[color_by]])) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.2, show.legend = F) +
    theme_bw() +
    labs(title = title,
         x = color_by,
         y = index) +
    custom_theme(angle_x = 45)
  
  if(!is.na(facet_by)){
    facet_formula <- as.formula(paste0('~ ', facet_by))
    p <- p +
      stat_pvalue_manual(data = stat_df,
                         label = 'p.adj.signif', 
                         inherit.aes = FALSE,
                         hide.ns = TRUE) +
      facet_wrap(facet_formula,
                 scales = 'free_y')
  } else {
    p <- p +
      stat_pvalue_manual(data = stat_df,
                         label = 'p.adj.signif', 
                         inherit.aes = FALSE,
                         hide.ns = TRUE)
  }
  
  violin_meta <- list(description = descriptions('violin',
                                                 thermo_idx = index,
                                                 group = group_vars),
                      insight = insights(contains = c('violin', index),
                                         df_plot = df_plot, 
                                         stat_df = stat_df, 
                                         color_by = color_by, 
                                         facet_by = facet_by, 
                                         index = index))
  
  return(list(plot = p,
              meta = violin_meta))
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
  
  df_plot <- df %>% 
    summarise(Count = mean(Count, na.rm = TRUE)) %>%
    mutate(Perc_count = Count/sum(Count, na.rm = TRUE)*100)
  
  plot <- df_plot %>%
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
  
  comp_bar_meta <- list(description = descriptions('composition_bar',
                                                   group = group,
                                                   bar_group = composition),
                        insight = insights(contains = c('composition_bar', composition),
                                           df_plot = df_plot,
                                           group = group,
                                           composition = composition))
  
  
  return(list(plot = plot,
              meta = comp_bar_meta))
}


# ******************************************************************************
plot_upset <- function(mass_list, group_by){
  
  group_values <- unique(df_longer[[group_by]])
  
  plot <- upset(fromList(mass_list), order.by = "freq", 
                nsets = length(group_values))
  
  upset_meta <- list(description = descriptions('upset',
                                                group = group_by),
                     insight = insights(contains = c('upset_plot'),
                                        mass_list = mass_list))
  
  return(list(plot = plot,
              meta = upset_meta))
  
}

plot_venn <- function(mass_list, filt_colors, val1, val2){
  
  plot <- ggvenn(mass_list, fill_color = as.character(filt_colors[1:2])) +
    labs(title = 'Number of metabolites (assigned molecular formula)') +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
  
  venn_meta <- list(description = descriptions('venn',
                                               val1 = val1,
                                               val2 = val2),
                    insight = insights(contains = c('venn'),
                                       mass_list = mass_list))
  
  return(list(plot = plot,
              meta = venn_meta))
  
}

# ******************************************************************************
plot_ordination <- function(df, x, y, color_by, col_vec, title, 
                            shape_by = NA, point_size = 2.5,
                            add_info = NA, x_label = NA, y_label = NA){
  color_by <- color_by
  shape_by <- shape_by
  
  group_vars <- c(color_by, shape_by)
  group_vars <- group_vars[!is.na(group_vars)]
  
  if(is.na(shape_by)){
    plot <- ggplot(df,
                   aes(x = {{ x }},
                       y = {{ y }},
                       color = !!rlang::ensym(color_by))) +
      geom_point(size = point_size)
  } else {
    plot <- ggplot(df,
                   aes(x = {{ x }},
                       y = {{ y }},
                       color = !!rlang::ensym(color_by),
                       shape = !!rlang::ensym(shape_by))) +
      geom_point(size = point_size)
  }
  
  plot <- plot + 
    scale_color_manual(values = col_vec) +
    labs(title = title) +
    custom_theme()
  
  ord_type <- colnames(df)[2] %>% str_remove('[0-9]')
  
  ordination_meta <- list(description = descriptions(paste0(ord_type, '_ordination'),
                                                     group = group_vars),
                          insight = insights(contains = 'ordination',
                                             ord_type = ord_type,
                                             add_info = add_info,
                                             color_by = color_by,
                                             shape_by = shape_by))
  
  return(list(plot = plot,
              meta = ordination_meta))
}

# ******************************************************************************
plot_diversity_index <- function(div_table, group_by, title){
  group_by <- group_by
  
  plot <- div_table %>% 
    ggplot(aes(x = .data[[group_by]],
               y = values,
               fill = .data[[group_by]])) +
    geom_boxplot() +
    scale_fill_manual(values = my_colors) +
    facet_wrap(~index, scales = 'free_y') +
    labs(title = title,
         x = group_by,
         y = 'Index value') +
    custom_theme()
  
  diversity_meta <- list(description = descriptions('diversity', group = group_by),
                         insight = insights(contains = c('diversity',
                                                         unique(div_table$index)),
                                            df_plot = div_table,
                                            group_by = group_by))
  
  return(list(plot = plot,
              meta = diversity_meta))
}

# Calculate functions ----

# ******************************************************************************
select_masses <- function(df, group, value){
  group <- group
  df %>% 
    filter(!!rlang::ensym(group) == value) %>% 
    pull(Mass) %>% 
    unique()
}

# ******************************************************************************
calculate_weighted <- function(df, index){
  new_names <- paste0(index, c('_weighted', '_magnitude_average'))
  names(new_names) <- c('wt', 'm_avg')
  weighted_df <- df %>%  
    select(Mass, SampleID, NormIntensity, all_of(c(index))) %>%
    mutate(wt = .data[[index]] * NormIntensity) %>% 
    group_by(SampleID) %>% 
    mutate(m_avg = sum(wt, na.rm = TRUE)/sum(NormIntensity, na.rm = TRUE)) %>% 
    select(Mass,SampleID, all_of(index), wt, m_avg) %>% 
    rename_with(.cols = c(wt, m_avg), .fn = function(x) new_names[x]) %>% 
    pivot_longer(!c(SampleID, Mass), names_to = 'weighted_idx', 
                 values_to = 'Magnitude-averaged/weigthed index') %>% 
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
  filename <- file.path(my_outdir, '5.1.1_NMDS.log')
  write_lines(as.data.frame(nmds.log), filename)
  
  filename <- file.path(my_outdir, '5.1.2_NMDS_stressplot.png')
  
  png(filename, res = 300, height = 3600, width = 3600)
  capture.output(stressplot(nmds))
  dev.off()
  
  nmds_scores <- as.data.frame(scores(nmds, display = 'sites')) %>% 
    rownames_to_column(var = 'SampleID') %>% 
    left_join(metadata, by = 'SampleID')
  
  return(list(scores = nmds_scores,
              nmds = nmds))
  
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
plot_richness <- function(richness){
  
  richness_long <- as_tibble(richness$perm, rownames = 'Sites') %>%
    pivot_longer(!Sites, names_to = 'permutation', values_to = 'Richness') %>%
    mutate(Sites = as.numeric(Sites))
  
  plot <- ggplot(richness_long,
                 aes(x = Sites,
                     y = Richness,
                     group = Sites)) +
    geom_boxplot(fill = 'yellow') +
    geom_hline(yintercept = max(richness_long$Richness) * 0.85, color = 'red') +
    theme_bw() +
    labs(title = 'Richness plot',
         x = 'Number of sites',
         y = 'Richness') +
    custom_theme()
  
  richness_meta <- list(description = descriptions('richness'),
                        insight = insights(contains = ('richness'),
                                           richness = richness))
  
  return(list(plot = plot,
              meta = richness_meta))
}

# ******************************************************************************
plot_rank_abundance <- function(mat){
  df_plot <- tibble(peak = colnames(mat),
                    intensity_sums = colSums(mat),
                    presence = colSums(mat != 0)) %>%
    arrange(intensity_sums) %>%
    mutate(position = n():1,
           presence_perc = presence/nrow(mat)*100)
  
  plot <- ggplot(df_plot,
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
  
  rank_abundance_meta <- list(description = descriptions('rank_abundance'),
                              insight = insights(contains = c('rank_abundance'),
                                                 df_plot = df_plot))
  
  return(list(plot = plot,
              meta = rank_abundance_meta))
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

descriptions <- function(plot_type, thermo_idx = NA, group = NA, bar_group = NA,
                         color_by = NA, val1 = NA, val2 = NA, ord_type = NA){
  broad_des <- list(
    van_krevelen = glue::glue(
      "Van Krevelen diagrams are figures used to visualize FTICR-MS detected peaks ",
      "based on their Oxygen to Carbon (O:C) and Hydrogen to Carbon (H:C) ratios. ",
      "Plot is colored based on {thermo_idx}, and faceted by {group[1]}.",
      color_by = color_by,
      group= group
    ),
    density = glue::glue(
      "Density plot showing the distribution of values for {thermo_idx} across ",
      "all samples.",
      " Plot is colored based on {group[1]}, and faceted by {group[1]}",
      group = group
    ),
    violin = glue::glue(
      "Figure uses boxplots surrounded by violin plots to show changes in the distribution ",
      "of {thermo_idx} values across different groups. And if there is any significance differences ",
      "among the groups. Plot is colored based on {group[1]}",
      group = group
    ),
    composition_bar = glue::glue(
      "Figure is a stacked bar that shows the relative abundances of the detected metabolites based on ",
      "their {bar_group}. Each bar represent a different sample group based on {group[1]}. Colors represent ",
      "different types of metabolites.",
      bar_group = bar_group,
      group = group
    ),
    upset = glue::glue(
      "The upset plot shows how many metabolites are shared between different sets of samples. It's ",
      "particularly useful for analyzing complex relationships between many sets. Sample sets ",
      "in this plot are based on {group[1]}.",
      group = group
    ),
    venn = glue::glue(
      "This Venn diagram shows how many metabolites are shared and unique between sample sets: ",
      "{val1} and {val2}.",
      val1 = val1,
      val2 = val2
    ),
    richness = glue::glue(
      "This plot shows how the number of detected metabolites increases with sampling effort. ",
      "If sampling effort was sufficient curve is expected to gradually plateus as fewer ",
      "new metabolites are detected with additional samples."
    ),
    rank_abundance = glue::glue(
      'The rank abundance plot ranks the detected masses by their relative abundance and the number of ',
      'samples they were found in. It allows to determine wheter the most peaks samples are also those ',
      'that are found on most of the samples or not.'
    ),
    diversity = glue::glue(
      "Figure uses boxplots to visualize the differences in diversity indexes across the different ",
      "sample groups defined by the variable {group[1]}.",
      group = group
    ),
    NMDS_ordination = glue::glue(
      "NMDS ordination is a technique that arranges samples in reduced dimensional space based on ",
      "rank-order distances, preserving the relative dissimilarity between samples without ",
      "assuming linear relationships. The NMDS ordination plot show community composition patterns ",
      "where closer points represent more similar samples. Points are colored by {group[1]}",
      group = group
    ),
    PC_ordination = glue::glue(
      "Principal Component Analysis (PCA) is a linear ordination technique that creates new orthogonal axes ",
      "(principal components) that maximize the variance explained in the dataset, with PC1 explaining the most ",
      "variance and subsequent axes explaining progressively less. It assumes linear relationships between variables ",
      "PCA plots display samples along these principal component axes. Points are colored by {group[1]}",
      group = group
    )
  )
  
  if(any(is.na(group))){
    broad_des$van_krevelen <- glue::glue(
      "Van Krevelen diagrams are figures used to visualize FTICR-MS detected peaks ",
      "based on their Oxygen to Carbon (O:C) and Hydrogen to Carbon (H:C) ratios. ",
      "Plot is colored based on {color_by}.",
      color_by = color_by,
      group= group
    )
  }
  
  if(length(group) == 2){
    broad_des$van_krevelen <- glue::glue(broad_des$density,
                                         ' and {group[2]}',
                                         group = group)
    
    broad_des$density <- glue::glue(broad_des$density,
                                    ' and {group[2]}',
                                    group = group)
    
    broad_des$violin <- glue::glue(broad_des$violin,
                                   ' and faceted based on {group[2]}',
                                   group = group)
    broad_des$NMDS_ordination <- glue::glue(broad_des$NMDS_ordination,
                                            ' and shapes represent {group[2]}',
                                            group = group)
    broad_des$PC_ordination <- glue::glue(broad_des$PC_ordination,
                                          ' and shapes represent {group[2]}',
                                          group = group)
  }
  
  return(broad_des[[plot_type]])
}

vk_insight <- function(df_plot, color_by, facet_col = NA, facet_row = NA,
                       color_continuous = TRUE){
  
  group_vars <- c(facet_col, facet_row, 'Class')
  group_vars <- group_vars[!is.na(group_vars)]
  
  counts <- df_plot %>% 
    group_by(across(all_of(group_vars))) %>% 
    count()
  
  if(is.na(facet_col) && is.na(facet_row)){
    points_insight <- counts %>% 
      mutate(insight = glue::glue('{n} masses as {Class}'))
    
    insight <- paste0(points_insight$insight, collapse = ', ')
    
  } else if(is.na(facet_row)) {
    points_insight <- map(unique(df_plot[[facet_col]]), function(fc){
      temp <- counts %>%
        filter(.data[[facet_col]] == fc) %>% 
        mutate(insight = glue::glue('{n} masses as {Class}'))
      
      points_insight <- glue::glue(
        'For group {fc}: ',
        paste0(temp$insight, collapse = ', ')
      )
    }) %>% reduce(`c`) %>% paste0(collapse = '. ')
  } else {
    points_insight <- map(unique(df_plot[[facet_col]]), function(fc){
      
      temp <- map(unique(df_plot[[facet_row]]), function(fr){
        tt <- counts %>%
          filter(.data[[facet_col]] == fc,
                 .data[[facet_row]] == fr) %>% 
          mutate(insight = glue::glue('{n} masses as {Class}'))
        
        ti_1 <- glue::glue(
          'For {facet_col} {fc} and {facet_row} {fr}: ',
          paste0(tt$insight, collapse = ', ')
        )
      }) %>% reduce(`c`) %>% paste0(collapse = '. ')
      
      
    }) %>% reduce(`c`) %>% paste0(collapse = '. ')
  }
  
  if(color_continuous){
    min_val <- min(df_plot[[color_by]], na.rm = TRUE)
    max_val <- max(df_plot[[color_by]], na.rm = TRUE)
    
    color_insight <- glue::glue(
      '. Plot is colored based on the {color_by} of each peak. {color_by} ranges ',
      'from {min_val} to {max_val}.',
      color_by = color_by,
      min_val = min_val,
      max_val = max_val
    )
  } else {
    color_insight <- glue::glue(
      '. Plot is colored based on the group of samples were each of the masses were detected. ',
      'Sample groups are defined by the {color_by} variable.'
    )
  }
  
  res <- paste0(
    "The position of the masses in the van Krevelen diagram can be used to ",
    "group detected compounds into broad molecular classes. The plot shows: ",
    points_insight,
    color_insight
  )
  
  return(res)
}

density_insight <- function(plot_data, facet_col, facet_row, index){
  
  density_data <- plot_data$data[[1]]
  
  max_density <- density_data %>% 
    select(PANEL, density, x) %>% 
    group_by(PANEL) %>% 
    slice_max(order_by = density, n = 1)
  
  if(is.na(facet_row)){
    facets <- plot_data$layout$layout %>% 
      select(PANEL, all_of(facet_col)) %>% 
      left_join(max_density, by = 'PANEL') %>% 
      mutate(thermo_idx = index,
             insight = glue::glue('{thermo_idx} value with the highest density in ',
                                  '{facet_col} {.data[[facet_col]]} is {round(x, 3)}'))
  } else {
    facets <- plot_data$layout$layout %>% 
      select(PANEL, all_of(c(facet_col, facet_row))) %>% 
      left_join(max_density, by = 'PANEL') %>% 
      mutate(thermo_idx = index,
             insight = glue::glue('{thermo_idx} value with the highest density in ',
                                  '{facet_col} {.data[[facet_col]]} and {facet_row} {.data[[facet_row]]} is {round(x, 3)}'))
  }
  
  return(paste0(facets$insight, collapse = ', '))
  
}

violin_insight <- function(df_plot, stat_df, color_by, facet_by, index){
  
  if(is.na(facet_by)){
    boxplot_data <- df_plot %>% 
      group_by(across(all_of(color_by))) %>% 
      summarise(
        q1 = quantile(.data[[index]], 0.25, na.rm = TRUE),
        median = quantile(.data[[index]], 0.5, na.rm = TRUE),
        q3 = quantile(.data[[index]], 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        lower_bound = q1 - 1.5*iqr,
        upper_bound = q3 + 1.5*iqr,
        lower_whisker = min(.data[[index]][.data[[index]] >= lower_bound], na.rm = TRUE),
        upper_whisker = max(.data[[index]][.data[[index]] <= upper_bound], na.rm = TRUE)
      ) %>% 
      ungroup() %>% 
      mutate(across(c(q1, median, q3, lower_whisker, upper_whisker), ~round(.x, 3)),
             insight = glue::glue('For {color_by} {.data[[color_by]]} the boxplot has a ',
                                  'median of {median}, upper hinge of {q3}, lower hinge of {q1} ',
                                  'upper whisker of {upper_whisker}, and lower whisker of {lower_whisker}.'))
    
    boxplot_insight <- paste0(boxplot_data$insight, collapse = ', ')
    
    sig_stat <- stat_df %>% 
      filter(p.adj < 0.05)
    
    if(nrow(sig_stat) > 0){
      sig_stat <- sig_stat %>% 
        mutate(insight = glue::glue('Significant differences found between {group1} and {group2} ',
                                    '(p-value = {round(p.adj, 5)})'))
      
      sig_insight <- paste0(sig_stat$insight, collapse = '. ')
      
      boxplot_insight <- paste0(c(boxplot_insight, sig_insight), collapse = ' ')
    }
    
    
  } else {
    boxplot_data <- df_plot %>% 
      group_by(across(all_of(c(color_by, facet_by)))) %>% 
      summarise(
        q1 = quantile(.data[[index]], 0.25, na.rm = TRUE),
        median = quantile(.data[[index]], 0.5, na.rm = TRUE),
        q3 = quantile(.data[[index]], 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        lower_bound = q1 - 1.5*iqr,
        upper_bound = q3 + 1.5*iqr,
        lower_whisker = min(.data[[index]][.data[[index]] >= lower_bound], na.rm = TRUE),
        upper_whisker = max(.data[[index]][.data[[index]] <= upper_bound], na.rm = TRUE)
      )
    
    boxplot_insight <- map(unique(df_plot[[facet_by]]), function(g){
      temp <- boxplot_data %>% 
        filter(.data[[facet_by]] == g) %>% 
        mutate(across(c(q1, median, q3, lower_whisker, upper_whisker), ~round(.x, 3)),
               insight = glue::glue('For {color_by} {.data[[color_by]]} the boxplot has a ',
                                    'median of {median}, upper hinge of {q3}, lower hinge of {q1} ',
                                    'upper whisker of {upper_whisker}, and lower whisker of {lower_whisker}.'))
      
      final <- glue::glue(
        'For facet {g}: ',
        paste0(temp$insight, collapse = ', ')
      )
    }) %>% reduce(`c`) %>% paste0(collapse = '. ')
    
    sig_stat <- stat_df %>% 
      filter(p.adj < 0.05)
    
    if(nrow(sig_stat) > 0){
      sig_insight <- map(unique(sig_stat[[facet_by]]), function(g){
        
        temp <- sig_stat  %>% 
          filter(.data[[facet_by]] == g) %>%
          mutate(insight = glue::glue('Significant differences found between {group1} and {group2} ',
                                      '(p-value = {round(p.adj, 5)})'))
        
        
        final <- glue::glue(
          'For facet {g}: ',
          paste0(sig_stat$insight, collapse = '. ')
        )
        
      }) %>% reduce(`c`) %>% paste0(collapse = '. ')
      
      boxplot_insight <- paste0(c(boxplot_insight, sig_insight), collapse = ' ')
    }
    
  }
  
  return(boxplot_insight)
  
}

comp_bar_insight <- function(df_plot, group, composition){
  
  if(length(group) == 1){
    bar_insight <- map(unique(df_plot[[group]]), function(g){
      temp <- df_plot %>%
        filter(.data[[group]] == g) %>% 
        filter(!is.nan(Perc_count)) %>% 
        mutate(insight = glue::glue('{round(Perc_count, 2)}% of masses are {.data[[composition]]}'))
      
      ins <- glue::glue(
        'For group "{g}": ',
        paste0(temp$insight, collapse = ', ')
      )
    }) %>% reduce(`c`) %>% paste0(collapse = '. ')
  }  else {
    bar_insight <- map(unique(df_plot[[group[2]]]), function(f){
      
      sub_res <- map(unique(df_plot[[group[1]]]), function(g){
        temp <<- df_plot %>%
          filter(.data[[group[1]]] == g) %>% 
          filter(!is.nan(Perc_count)) %>% 
          mutate(insight = glue::glue('{round(Perc_count, 2)}% of masses are {.data[[composition]]}'))
        
        subins <- glue::glue(
          'For group {g}: ',
          paste0(temp$insight, collapse = ', ')
        )
      }) %>% reduce(`c`) %>% paste0(collapse = '. ')
      
      ins <- glue::glue(
        'For facet {f}: ',
        sub_res
      )
      
    })%>% reduce(`c`) %>% paste0(collapse = '. ')
    
  }
  
  return(bar_insight)
  
}

upset_insight <- function(mass_list){
  
  
  df_mass <- imap(mass_list, function(ml, nm){
    data.frame(Mass = ml) %>% 
      mutate(group = nm)
  }) %>% reduce(rbind) %>% 
    group_by(Mass) %>% 
    summarise(group = paste(group, collapse = ',')) %>% 
    group_by(group) %>% 
    count(name = 'count') 
  
  total_masses <- sum(df_mass$count)
  
  top <- df_mass %>% 
    ungroup() %>% 
    slice_max(order_by = count, n = 5) %>% 
    mutate(insight = glue::glue('{group} ({round(count/total_masses *100, 2)}%)'))
  
  insight <- glue::glue(
    'Most metabolites are found in sets: ',
    paste(top$insight, collapse = ', ')
  )
  
}

venn_insight <- function(mass_list){
  
  df_mass <- imap(mass_list, function(ml, nm){
    data.frame(Mass = ml) %>% 
      mutate(group = nm)
  }) %>% reduce(rbind) %>% 
    group_by(Mass) %>% 
    summarise(group = paste(group, collapse = ',')) %>% 
    group_by(group) %>% 
    count(name = 'count') %>% 
    mutate(group = ifelse(str_detect(group, ','), 'Shared', group))
  
  total_masses <- sum(df_mass$count)
  
  top <- df_mass %>% 
    mutate(insight = glue::glue('{group} ({round(count/total_masses *100, 2)}%)'))
  
  insight <- glue::glue(
    'Metabolite distribution: ',
    paste(top$insight, collapse = ', ')
  )
  
}

richness_insight <- function(richness){
  
  asymp_model <- fitspecaccum(richness, 'asymp')
  
  current_richness <- max(richness$richness)
  
  plateu_ratio <- current_richness / coef(asymp_model)[1]
  
  if(plateu_ratio > 0.9){
    is_plateu <- glue::glue("The richness curves plateus showing that sampling was sufficient to capture ",
                            "the organic matter diversity of the samples")
  } else {
    is_plateu <- glue::glue("The richness curves does not plateu showing that more sampling is needed to capture ",
                            "the organic matter diversity of the samples")
  }
  
  insight <- glue::glue(
    "The species accumulation curves show that the maximum metabolite richness achieved is ",
    "{current_richness}. {is_plateu}"
  )
  
}

rank_abundance_insight <- function(df_plot){
  
  top_5 <- df_plot %>% 
    slice_max(order_by = intensity_sums, prop = 0.05)
  
  top_5_mean_pres <- mean(top_5$presence_perc)
  
  bot_5 <- df_plot %>% 
    slice_min(order_by = intensity_sums, prop = 0.05)
  
  bot_5_mean_pres <- mean(bot_5$presence_perc)
  
  
  insight <- glue::glue(
    "The rank abundance plot shows that the top 5% more abundant metabolites appear in ",
    "around {round(top_5_mean_pres, 2)}% of the samples. The bottom 5% less abundant ",
    "metabolites appear in {round(bot_5_mean_pres, 2)}% of the samples."
  )
  
}

diversity_insight <- function(df_plot, group_by){
  boxplot_data <-  df_plot %>% 
    group_by(across(all_of(group_by)), index) %>% 
    summarise(
      q1 = quantile(values, 0.25, na.rm = TRUE),
      median = quantile(values, 0.5, na.rm = TRUE),
      q3 = quantile(values, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower_bound = q1 - 1.5*iqr,
      upper_bound = q3 + 1.5*iqr,
      lower_whisker = min(values[values >= lower_bound], na.rm = TRUE),
      upper_whisker = max(values[values <= upper_bound], na.rm = TRUE)
    ) %>% 
    ungroup() 
  
  boxplot_insight <- map(unique(df_plot$index), function(idx){
    temp <- boxplot_data %>% 
      filter(index == idx) %>% 
      mutate(across(c(q1, median, q3, lower_whisker, upper_whisker), ~round(.x, 3)),
             insight = glue::glue('For {group_by} {.data[[group_by]]} the boxplot has a ',
                                  'median of {median}, upper hinge of {q3}, lower hinge of {q1} ',
                                  'upper whisker of {upper_whisker}, and lower whisker of {lower_whisker}.'))
    
    final <- glue::glue(
      'For diversity index {idx}: ',
      paste0(temp$insight, collapse = ' ')
    )
  }) %>% reduce(`c`) %>% paste0(collapse = '. ')
}

ordination_insight <- function(ord_type, add_info, color_by, shape_by = NA){
  
  if(ord_type == 'NMDS'){
    if(is.na(shape_by)){
      if(add_info$permanova$`Pr(>F)`[1] < 0.05){
        permanova_ins <- glue::glue(
          "The permutational analysis of variance found a significant structure of the community ",
          "based on {color_by} (p_value = {round(add_info$permanova$`Pr(>F)`[1], 2)})."
        )
      } else{
        permanova_ins <- glue::glue(
          "No significant structure of the community based on {color_by}."
        )
      }
    } else {
      if(add_info$permanova$`Pr(>F)`[1] < 0.05){
        permanova_ins <- glue::glue(
          "The permutational analysis of variance found a significant structure of the community ",
          "based on {color_by} and {shape_by} (p_value = {round(add_info$permanova$`Pr(>F)`[1], 2)})."
        )
      } else{
        permanova_ins <- glue::glue(
          "No significant structure of the community based on {color_by} and {shape_by}."
        )
      }
    }
    
    stress_ins <- glue::glue(
      "NMDS stress is {round(add_info$nmds$stress,2)}",
    )
    
    final_ins <- paste(permanova_ins, stress_ins, collapse = ' ')
    
  } else {
    axis1 <- round(add_info$pca$eigenvalues$variance.percent[1], digits = 2)
    axis2 <- round(add_info$pca$eigenvalues$variance.percent[2], digits = 2)
    
    final_ins <- glue::glue(
      "PCA ordination plot shows that PC1 explains {axis1}% of the variability amongst the samples, ",
      "while PC2 explains {axis2}% of the variability."
    )
  }
  
  return(final_ins)
  
}

insights <- function(contains, ...){
  
  dots <- list(...)
  
  other_options <- list(
    GFE = paste0(
      "Gibbs Free Energy (GFE) reflects thermodynamic stability and bioavailability. ",
      "Lower GFE values indicate more reduced, energy-rich compounds (e.g., lipids), ",
      "while higher values suggest oxidized, recalcitrant materials"
    ),
    AI_mod = paste0(
      "Modified Aromaticity Index (AI_mod) quantifies aromatic character. ",
      "Higher AI_mod values (>0.5) indicate condensed aromatic structures typical of ",
      "pyrogenic or highly processed organic matter. Lower values suggest aliphatic or ",
      "less condensed structures. Variations may reflect differences in source material ",
      "or degradation history."
    ),
    DBE = paste0(
      "Double Bond Equivalents (DBE) indicate molecular unsaturation and ring structures. ",
      "Higher DBE values suggest aromatic rings and unsaturated bonds characteristic of ",
      "complex, stable organic compounds. Lower values indicate saturated, aliphatic structures"
    ),
    NOSC = paste0(
      "Nominal Oxidation State of Carbon (NOSC) reflects the average oxidation state. ",
      "Positive NOSC indicates oxidized compounds (e.g., carboxylic acids), while negative ",
      "values indicate reduced compounds (e.g., hydrocarbons, lipids). Shifts in NOSC ",
      "distributions suggest changes in redox chemistry or organic matter processing."
    ),
    Class = paste0(
      "Molecular class composition reveals the relative abundance of biochemical compound categories ",
      "(e.g., lipids, proteins, lignins, tannins, condensed hydrocarbons). Differences across sample groups",
      " indicate variations in organic matter sources, biological activity, or degradation pathways."
    ),
    Element = paste0(
      "Elemental composition (CHO, CHON, CHOS, etc.) provides insight into heteroatom content. ",
      "Nitrogen-containing compounds suggest biological/protein-like material; sulfur indicates ",
      "contributions from microbial metabolism or marine sources; phosphorus reflects biological activity. ",
      "Compositional shifts reveal changes in organic matter origin and processing."
    ),
    Shannon = paste0(
      "Measures both species (metabolite) richness and evenness by considering the abundance of each metabolite ", 
      ", with higher values indicating greater diversity. Used to compare diversity between different communities or treatments"
    ),
    Gini_Simpson = paste0(
      "Calculates the probability that two randomly selected individuals (metabolites) belong to different species, ",
      "ranging from 0-1 where 1 indicates maximum diversity. Used to compare community diversity while being less ",
      "sensitive to rare species than Shannon index, making it robust for comparing communities with different ",
      "sampling intensities."
    ),
    Chao1 = paste0(
      "Estimates the true total species richness (including unobserved species) based on the number of singletons ",
      "and doubletons in the sample. Used to assess sampling completeness and estimate how many species might ",
      "be missed, particularly valuable when comparing communities with different sampling efforts or rare species."
    ),
    Elemental_composition = paste0(
      "Using Rao's quadratic entropy to measure functional diversity based on differences in elemental composition ",
      "between molecules, where higher values indicate greater chemical composition diversity. Used to assess how ",
      "biochemically diverse a community is in terms of stoichiometry, helping understand biogeochemical cycling potential."
    ),
    Reactivity = paste0(
      "Using Rao's quadratic entropy to measure functional diversity based on the nominal oxidation state of carbon, ",
      "(NOSC) reflecting differences in thermodynamic favorability and metabolic accessibility of compounds. ",
      "Used to evaluate the range of energy-yielding potential within a community, indicating metabolic ",
      "versatility and the availability of compounds across different redox states."
    ),
    Insaturation_and_aromaticity = paste0(
      "Using Rao's quadratic entropy to measure functional diversity based on molecular structure complexity, ",
      "including double bond equivalents and aromatic character, reflecting molecular stability and bioavailability. ",
      "Used to assess structural complexity within communities, indicating the range from ",
      "labile (easily degraded) to recalcitrant (persistent) compounds, which affects decomposition rates and carbon cycling."
    )
  )
  
  res <- c()
  
  for(value in contains){
    if(value =='van_krevelen'){
      temp <- vk_insight(df_plot = dots$df_plot, 
                         color_by = dots$color_by, 
                         facet_col = dots$facet_col,
                         facet_row = dots$facet_row,
                         color_continuous = dots$color_continuous)
      
      res <- c(res, temp)
      
    } else if(value == 'density') {
      temp <- density_insight(plot_data = dots$plot_data,
                              facet_col = dots$facet_col,
                              facet_row = dots$facet_row,
                              index = dots$index)
      
      res <- c(res, temp)
      
    } else if(value == 'violin') {
      temp <- violin_insight(df_plot = dots$df_plot,
                             stat_df = dots$stat_df,
                             color_by = dots$color_by,
                             facet_by = dots$facet_by,
                             index = dots$index)
      
      res <- c(res, temp)
      
    } else if(value == 'composition_bar') {
      temp <- comp_bar_insight(df_plot = dots$df_plot,
                               group = dots$group,
                               composition = dots$composition)
      
      res <- c(res, temp)
      
    } else if(value == 'upset_plot') {
      temp <- upset_insight(mass_list = dots$mass_list)
      
      res <- c(res, temp)
      
    } else if(value == 'venn') {
      temp <- venn_insight(mass_list = dots$mass_list)
      
      res <- c(res, temp)
      
    } else if(value == 'richness') {
      temp <- richness_insight(richness = dots$richness)
      
      res <- c(res, temp)
      
    } else if(value == 'rank_abundance') {
      temp <- rank_abundance_insight(df_plot = dots$df_plot)
      
      res <- c(res, temp)
      
    } else if(value == 'diversity') {
      temp <- diversity_insight(df_plot = dots$df_plot,
                                group_by = dots$group_by)
      
      res <- c(res, temp)
      
    } else if(value == 'ordination') {
      temp <- ordination_insight(ord_type = dots$ord_type,
                                 add_info = dots$add_info,
                                 color_by = dots$color_by,
                                 shape = dots$shape_by)
      
      res <- c(res, temp)
      
    } else if(value %in% names(other_options)){
      temp <- other_options[[value]]
      
      res <- c(res, temp)
    }
    
  }
  
  res_final <- paste(res, collapse = ' ')
  
  return(res_final)
}

update_figure_list <- function(figure_id = '', 
                               figure_title = ''){
  
  list_current_idx <- length(figure_list)
  
  metadata_file <- file.path('figures_metadata', paste0(figure_id, '.json'))
  
  info <- list(figure_id = figure_id,
               figure_title = figure_title,
               metadata_file = metadata_file)
  
  figure_list[[list_current_idx + 1]] <<- info
  
}

save_figure_metadata <- function(
    plot_obj,  
    figure_id,
    analysis_module,
    figure_type,
    figure_title = NA,
    x_axis_label = NA,
    y_axis_label = NA,
    caption,
    figure_file,
    group_aes,
    modifiable_aesthetics,
    has_legend = TRUE,
    custom_legend = NULL,
    units = 'none',
    data_source,
    r_script_path,
    functions_used,
    resolution,
    height,
    width
) {
  
  # Getting information from the plot
  p <- plot_obj$plot
  
  if(is.na(figure_title)){
    figure_title <- p@labels$title
  }
  
  if(any(is.na(x_axis_label))){
    x_axis_label <- p@labels$x
  }
  
  if(any(is.na(y_axis_label))){
    y_axis_label <- p@labels$y
  }
  
  grouping_variables <- list()
  for(i in group_aes){
    if(all(i == 'facets')){
      if(!is.null(p@facet$params$facets)){
        grouping_variables[[i]] <- list(wrap = names(p@facet$params$facets))
      } else {
        grouping_variables[[i]] <- list(cols = names(p@facet$params$cols),
                                        rows = names(p@facet$params$rows))
      }
    } else if(all(i %in% c('colour', 'color', 'shape', 'fill', 'x', 'y'))) {
      grouping_variables[i] <- as.character(p@mapping[i]) %>% 
        str_remove('~') %>% 
        str_remove('.data\\[\\["') %>% 
        str_remove('"\\]\\]')
    } else {
      grouping_variables <- group_aes
    }
  }
  
  legend_variables <- list()
  if(has_legend){
    accepted_legend_vars <- c('colour', 'color', 'fill', 'shape')
    for(m in names(p@mapping)){
      if(m %in% accepted_legend_vars){
        legend_variables[m] <- as.character(p@mapping[m]) %>% 
          str_remove('~') %>% 
          str_remove('.data\\[\\["') %>% 
          str_remove('"\\]\\]')
      } 
    }
  } else {
    if(is.null(custom_legend)){
      legend_variables['none'] <- 'none'
    }
  }
  
  if(!is.null(custom_legend)){
    legend_variables <- custom_legend
  }
  
  # Saving all data as JSON
  info <- list(
    figure_id = figure_id,
    analysis_module = analysis_module,
    figure_type = figure_type,
    figure_title = figure_title,
    caption = caption,
    description = plot_obj$meta$description,
    insights = plot_obj$meta$insight,
    plot_info = list(
      x_axis_label = x_axis_label,
      y_axis_label = y_axis_label,
      grouping_variables = grouping_variables,
      modifiable_aesthetics = modifiable_aesthetics,
      units = units,
      legend = legend_variables
    ),
    data_source = data_source,
    script_path = list(script_path = r_script_path,
                       functions = file.path(dirname(dirname(r_script_path)),
                                             'custom_functions.R')),
    functions_used = functions_used,
    figure_file_info = list(
      figure_file = figure_file,
      last_modified = file.mtime(figure_file),
      resolution = resolution,
      width = width,
      height = height
    )
  )
  
  # Add to global list
  write_json(info, 
             file.path(dirname(dirname(r_script_path)),
                       'figures_metadata', 
                       paste0(figure_id, '.json')),
             auto_unbox = TRUE,
             pretty = TRUE)
}


