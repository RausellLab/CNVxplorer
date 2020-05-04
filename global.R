
# Load libraries
library(tablerDash)
library(shinyEffects)
library(echarts4r)
library(shinyWidgets)
library(shinyjs)
# library(karyoploteR)
library(DT)
library(XML)
library(highcharter)
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
library(ggraph)
library(igraph) 
library(widyr)
library(tidytext)
library(tidyverse)


# # Files needed to run the app
# 
# # 1. hgcn_genes - List all the genes with scores
# # 2. df_enhancers - List enhancers with phast100vert  + phast46pla + phast46pri
# # 3. lncrna_coord - coordinates lncRNA
# # 4. lncrna - Features per lncRNA
# # 5. tad - List of TADs with coordinates
# # 6. gtex - gene expression data across tissues
# # 7. hpa - human protein expression
# # 8. hpo - Human phenotype ontology
# # 9. mgi - mouse model data
# # 10. blacklist_encode - problematics regions
# # 11. mpo_dbs - Mammalian Phenotype Ontology 
# # 12. hpo_dbs - HPO ontology
# # 13. clinvar_variants - Variants from ClinVar (chrom - pos - id)
# # 14. omim - OMIM information
# # 15. orphanet_raw - ORPHANET
# # 16. dev_raw - DECIPHER (DEVELOPMENTAL DISORDER GENES)
# # 17- gnomad_sv_raw
# # 18. decipher_control_raw
# # 19. dgv_df_raw
# # 20. vector_total_terms (Dataframe of term - description - term_description)
# # 21. anato_df (Anatomies entities: description - term)
# # 22. mirtarbase - miRNA database + coordinates from miRbase
# # vector_inheritance - Terms associated with mode of inheritances
# # trrust - Transcription factor - target genes
# # drugbank
# # prot_complex - Protein complex genes
# # ohno_genes - Ohnologs genes 
# # recomb rates - recomb
# # genes_promoter - Mean phastCons46 - o/e CpG sites
# # para_genes - Paralogous genes
# # string_db - centrality measures (degree, page_rank)
# # pubmed_df - Nº hits deletion - duplication
# # ensembl_reg  - Ensembl Regulatory Build
# # region_gaps - Centromeric and telomeric regions
# # fusil_score - essential genes
# # coord_chrom_hg19 - chromosome length

# 
# file.remove('local_data.RData.gz')
# file.remove('local_data.RData')
# save(hgcn_genes, ## hg19
#     df_enhancers, ## hg19
#     tad, ## hg19
#     gtex, ## -
#     hpa, ## -
#     hpo_genes, ## -
#     cnv_df, ## hg19
#     vector_total_terms, ## -
#     gnomad_sv_raw, ## hg19
#     decipher_control_raw, ## hg19
#     dgv_df_raw, ## hg19
#     hpo_omim, ## -
#     anato_df, ## -
#     mirtarbase, ## hg19
#     vector_inheritance, ## -
#     trrust, ## hg19
#     tf_genes, ## -
#     drugbank, ## -
#     dev_raw, 
#     panel_total, 
#     omim, 
#     orphanet_raw,  
#     hpo_dbs, 
#     lncrna, ## hg19
#     lncrna_target, ## hg19
#     denovo, ## hg19
#     clinvar_variants,  ## hg19
#     plot_p100, ## -
#     plot_p46pla, ## -
#     lncrna_coord, ## hg19
#     blacklist_encode, ## hg19
#     mpo_dbs, ## -
#     gwas_variants, ## -
#     mgi, ## -
#     syndromes_total, ## hg19
#     file = "local_data.RData")

# system("gzip local_data.RData")

source('functions.R')

#  system('gunzip -c local_data.RData.gz > local_data.RData')
# load('local_data.RData')

ridges_home <- cnv_df %>%
  filter(length_cnv >= 50) %>%
  mutate(source =
           case_when(
             source == 'dgv' ~ 'DGV',
             source == 'decipher' ~ 'DECIPHER',
             source == 'gnomad_v2.1' ~ 'gnomAD v2.1',
             source == 'decipher_control' ~ 'DECIPHER Control'
           )) %>%
  ggplot(aes(length_cnv, y = source)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = source), alpha = 0.6, show.legend = FALSE, size = 1.25) +
  scale_x_log10() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis_d() +
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

coord_chrom_hg19 <- read_tsv('https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',
                             col_names = c('chrom', 'length')) %>%
  filter(nchar(chrom) < 6) %>% filter(!str_detect(chrom, 'chrM')) %>%
  mutate(chrom = str_remove(chrom, 'chr'))


# a <- c("chr1\t249250621",
#   "chr2\t243199373",
#   "chr3\t198022430",
#   "chr4\t191154276",
#   "chr5\t180915260",
#   "chr6\t171115067",
#   "chr7\t159138663",
#   "chrX\t155270560",
#   "chr8\t146364022",
#   "chr9\t141213431",
#   "chr10\t135534747",
#   "chr11\t135006516",
#   "chr12\t133851895",
#   "chr13\t115169878",
#   "chr14\t107349540",
#   "chr15\t102531392",
#   "chr16\t90354753",
#   "chr17\t81195210",
#   "chr18\t78077248",
#   "chr20\t63025520",
#   "chrY\t59373566",
#   "chr19\t59128983",
#   "chr22\t51304566",
#   "chr21\t48129895")
# 
# coord_chrom_hg19 <- read_tsv(a,
#                              col_names = c('chrom', 'length')) %>%
#   filter(nchar(chrom) < 6) %>% filter(!str_detect(chrom, 'chrM')) %>%
#   mutate(chrom = str_remove(chrom, 'chr'))
