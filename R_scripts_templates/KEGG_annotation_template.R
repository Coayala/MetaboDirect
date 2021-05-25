# #############################################################
# 
# MetaboDirect
# Data Exploration step
# version 1.0
# by Christian Ayala
# Licensed under the MIT license. See LICENSE.md file.
# 
# #############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(KEGGREST)
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

setwd('%currentdir%')

my_data.file <- file.path('%outdir%', '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_outdir <- file.path('%outdir%', '3_exploratory')

#### Import data ####

df <-  read_csv(my_data.file, col_types = cols())

#### Molecular formula annotation with KEGG database ####

df_formula <- df %>% 
  select(Mass, MolecularFormula)

compound_df <- tibble(Mass = NA, MolecularFormula = NA, KEGG_id = NA, KEGG_name = NA, KEGG_formula = NA,
                      KEGG_pathway = NA, KEGG_module = NA, KEGG_brite = NA, .rows = 0)

for(i in 1:nrow(df_formula)){
  formula <- df_formula$MolecularFormula[i]
  cpd_id <- names(keggFind('compound', formula, 'formula'))
  if(!is.null(cpd_id)){
    for(id in cpd_id){
      cpd_info <- keggGet(id)
      for(j in length(cpd_info)){
        temp <- tibble(Mass = df_formula$Mass[i],
                       MolecularFormula = df_formula$MolecularFormula[i],
                       KEGG_id = id,
                       KEGG_name = paste(cpd_info[[j]]$NAME, collapse = ''),
                       KEGG_formula = cpd_info[[j]]$FORMULA,
                       KEGG_pathway = ifelse(!is.null(cpd_info[[j]]$PATHWAY), paste(cpd_info[[j]]$PATHWAY, collapse = ';'), NA),
                       KEGG_module = ifelse(!is.null(cpd_info[[j]]$MODULE), paste(cpd_info[[j]]$MODULE, collapse = ';'), NA),
                       KEGG_brite = ifelse(!is.null(cpd_info[[j]]$BRITE), paste(cpd_info[[j]]$BRITE, collapse = ';'), NA))
      }
    }
  } else {
    temp <- tibble(Mass = df_formula$Mass[i], MolecularFormula = df_formula$MolecularFormula[i], KEGG_id = NA, 
                   KEGG_name = NA, KEGG_formula = NA, KEGG_pathway = NA, KEGG_module = NA, KEGG_brite = NA)
  }
  compound_df <- rbind(compound_df, temp)
}

matching_compounds <- df_formula %>% 
  left_join(compound_df[compound_df$KEGG_formula %in% df_formula$MolecularFormula,], by = c('Mass', 'MolecularFormula'))

filename <- file.path(my_outdir, 'KEGG_annotation_results.csv')
write_csv(matching_compounds, filename)
