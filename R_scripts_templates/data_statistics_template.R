# ##############################################################
# 
# MetaboDirect
# Data Statistics step
# version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
# 
# ##############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
}) 

# Values between two '%' are to be replaced by the correct values during the python script

#### Defining paths and variables ####

workdir <- file.path(getwd(), '%outdir%', '4_statistics')
setwd(workdir)

my_data.file <- file.path('..', '1_preprocessing_output', 'matrix_features.csv')
my_metadata.file <- file.path('..', '..', '%metadata%')
my_outdir <- getwd()

#### Import data ####

df <-  read_csv(my_data.file, col_types = cols()) %>% 
  column_to_rownames(var = 'Mass')

df <- t(df)
df[is.na(df)] <- 0

metadata <- read_csv(my_metadata.file, col_types = cols())

#### NMDS ####

## Converting absolute abundance to either relative abundance or presence/absence
## Relative abundance use method "total"
## For presence/absence transformation use method "pa"

trans = "pa"
t.matrix <- decostand(df, trans)


## Calculate distance matrix
## Use "bray" for relative abundance
## Use "jaccard for presence/absence

dm.method = "euclidean"
dist.matrix <- vegdist(t.matrix , method = dm.method) 

set.seed(123)
nmds.log <- capture.output(nmds <- metaMDS(dist.matrix, 
                k = 2,
                distance = dm.method,
                maxit = 999,
                trymax = 500,
                wascores = TRUE))
write_tsv(as.data.frame(nmds.log), 'nmds.log')
stressplot(nmds)

nmds.scores <- as.data.frame(scores(nmds)) %>% 
  rownames_to_column(var = 'SampleID') %>% 
  left_join(metadata, by = 'SampleID')

## NMDS plot

nmds_plot <- ggplot(nmds.scores,
                    aes(x = NMDS1,
                        y = NMDS2,
                        color = %group1%,
                        shape = %group2%)) +
  geom_point() +
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw() +
  labs(title = 'NMDS plot') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))
nmds_plot

filename <- file.path(my_outdir, 'NMDS_plot.png')
ggsave(filename, nmds_plot, dpi = 300)

#### PERMANOVA ####

set.seed(456)

dm.method = "bray"
dist.matrix <- vegdist(df , method = dm.method) # dist matrix

permanova <- adonis(dist.matrix ~ Habitat + Sampletype, 
                    data=metadata, 
                    permutations=999, 
                    method=dm.method)

filename <- file.path(my_outdir, 'permanova.csv')
write_csv(permanova$aov.tab, filename)
