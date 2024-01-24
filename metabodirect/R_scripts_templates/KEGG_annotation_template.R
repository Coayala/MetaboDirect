# ***************************************************************
# 
# MetaboDirect
# Data Exploration step
# version 1.0
# by Christian Ayala
# Licensed under the MIT license. See LICENSE.md file.
# 
# ***************************************************************

# Loading libraries ----

library(tidyverse)
library(KEGGREST)


# Values between two '%' are to be replaced by the correct values during the python script

# Defining variables ----

current_dir <- '%currentdir%'
output_dir <- '%outdir%'

setwd('%currentdir%')

my_data.file <- file.path(output_dir, '1_preprocessing_output', 'Report_processed_MolecFormulas.csv')
my_outdir <- file.path(output_dir, '3_exploratory')

# Loading custom functions ----
source(file.path(output_dir, 'custom_functions.R'))

# Import data ----

df <-  read_csv(my_data.file, col_types = cols())

# Molecular formula annotation with KEGG database ----

kegg_searches <- map2(df$Mass, df$MolecularFormula, function(mass, formula){
  
  print(paste0('Searching formula: ', formula))
  
  cpd <- keggFind('compound', formula, 'formula')
  cpd <- cpd[which(cpd == formula)]
  cpd_id <- names(cpd)
  
  if(length(cpd_id) > 0 && length(cpd_id) <= 10){
    cpd_df <- get_kegg_compound_info(cpd_id, mass, formula)
    
  } else if(length(cpd_id) > 10){
    # Kegg get only gives back 10 compounds at the time
    id_chunks <- split(cpd_id, ceiling(seq_along(cpd_id)/10))
    
    cpd_chunks <- map(id_chunks, function(x){
      res <- get_kegg_compound_info(x, mass, formula)

      Sys.sleep(3)
    })
    
    cpd_df <- reduce(cpd_chunks, rbind)
    
  } else {
    cpd_df <- tibble(Mass = mass,
                     MolecularFormula = formula,
                     KEGG_id = NA,
                     KEGG_name = NA,
                     KEGG_formula = NA,
                     KEGG_pathway = NA,
                     KEGG_module = NA,
                     KEGG_brite = NA,
                     KEGG_enzyme = NA,
                     KEGG_reaction = NA)
  }
  
  return(cpd_df)
  
})

compound_df <- reduce(kegg_searches, rbind)

filename <- file.path(my_outdir, '7_KEGG_annotation_results.csv')
write_csv(matching_compounds, filename)
