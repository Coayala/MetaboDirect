# ##############################################################
# 
# MetaboDirect
# Data Normalization selection
# version 1.0
# by Christian Ayala
# 
# Licensed under the MIT license. See LICENSE.md file.
# 
# ##############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(pmartR)
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

setwd('%currentdir%')

my_metadata.file <- file.path('%metadata%')
my_data.file <- file.path('%data%')
source.file <- file.path('%Metabo_HOME%', 'spans.R')
log_transform <- %log%
group <- '%group%'

source(source.file)
metadata <- read.csv(my_metadata.file)
data <- read.csv(my_data.file) 
# data$Mass <- as.character(data$Mass)
  
data <- data %>% 
  select(Mass, all_of(metadata$SampleID)) %>% 
  pivot_longer(all_of(metadata$SampleID), names_to = 'SampleID', values_to = 'Intensity') %>% 
  filter(Intensity > 0) %>% 
  pivot_wider(names_from = 'SampleID', values_from = 'Intensity')

data <- as.data.frame(data)

# data[data == 0] <- 1

pepdata <- as.pepData(data, metadata, edata_cname = 'Mass', fdata_cname = 'SampleID')

if(log_transform){
  pepdata <- edata_transform(pepdata, 'log')
}

span_test <- modified_spans_procedure(pepdata, group = group, sig_thresh = 0.001)


spans <- span_test %>%
  mutate(y_axis = paste0(subset_method, ' ', parameters))

spans_hmp <- ggplot(spans) +
  geom_tile(aes(x = normalization_method,
                y = y_axis,
                fill = SPANS_score),
            color = 'white') +
  coord_fixed() +
  scale_fill_viridis_c() +
  labs(title = 'SPANS scores for different normalizations',
       x = 'Normalization methods',
       y = 'Subset method and subset parameter',
       fill = 'SPANS score')

filename <- 'SPANS_scores.png'
ggsave(filename, spans_hmp, dpi = 300)


