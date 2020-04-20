
# Load libraries
library(tablerDash)
library(shinyEffects)
library(echarts4r)
library(shinyWidgets)
library(shinyjs)
# library(karyoploteR)
library(DT)
library(XML)
library(gghighlight)
library(ReactomePA)
library(shinycssloaders)
library(plotly)
library(waiter)
library(DescTools)# library(hrbrthemes)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(DOSE)
library(enrichplot)
library(rentrez)
library(reactable)
library(ggridges)
library(shinymanager)
library(UpSetR)
library(randomForest) # delete in case of using an alternative model
library(chromPlot)
import::from(Gviz, "IdeogramTrack")
import::from(Gviz, 'plotTracks')
library(shiny)
library(ontologySimilarity)
library(ontologyIndex)
library(formattable)
library(valr)
library(TissueEnrich)
library(arules)
library(arulesCBA)
library(shinyalert)
library(tidyverse)


source('functions.R')

# system('gunzip -c local_data.RData.gz > local_data.RData')
load('local_data.RData')

ridges_home <- cnv_df %>%
  mutate(source =
           case_when(
             source == 'dgv' ~ 'DGV',
             source == 'decipher' ~ 'DECIPHER',
             source == 'gnomad_v2.1' ~ 'gnomAD v2.1',
             source == 'decipher_control' ~ 'DECIPHER Control'
           )) %>%
  ggplot(aes(length_cnv, y = source)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = source), alpha = 0.6, show.legend = FALSE, size = 1.25) +
  # geom_vline(aes(xintercept = size_cnv_query), linetype = 2, color = 'red', size = 1.5) +
  scale_x_log10() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis_d() +
  # scale_fill_manual(values = c('#CD5C5C','#32CD32', '#32CD32')) +
  xlab('log10(CNVs size)') +
  ylab('Database') +
  theme_ridges()



theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size=22))
  
}

human_chrom <- list('chr1' = 1, 'chr2' = 2,'chr3' = 3,'chr4' = 4,'chr5' = 5,'chr6' = 6,'chr7' = 7,'chr8' = 8,'chr9' = 9,'chr10' = 10,'chr11' = 11,'chr12' = 12,'chr13' = 13,
                    'chr14' = 14,'chr15' = 15,'chr16' = 16,'chr17' = 17,'chr18' = 18,'chr19' = 19, 'chr20' = 20, 'chr21' = 21,  'chr22' = 22,
                    'chrX' = 'X','chrY' = 'Y')