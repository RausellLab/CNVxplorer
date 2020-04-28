

# Load packages
library(readxl)
library(biomaRt)
library(data.table)
library(tidyverse)
library(clusterProfiler)
library(data.table)
select <- dplyr::select
rename <- dplyr::rename



# TODO

# 400 genes does not have coordinates (start/end)
# OMIM
# HI%
# Percentile RVIS or original data
# ReactomePA analysis 423MB careful!
# Gene ontology based on non-coding regions
# difference between pLI and o/e
# should a put a minus + plus to the output from remot
# PIK3CD duplicated because database HI has two observations
# warning message with input outside of the maximum chromosomes
# 17q12 genes that do not have a proper chromosome (e.g. ENSG00000270806)

# QUESTIONS NECKER

# Input user, two options: 1- select chrom, g banding or select g banding where the chromosome is included

# Reference of scores

ref_scores <- tibble(score = c('gwas', # we could filter out by number of hits (no filter - 8652, >1 4498)
                               'fda',
                               'pli',
                               'deg', # 
                               'omim', # omim total - autosomal dominant - recessive - X-linked
                               'dev_disorder',
                               'haplo',
                               'triplo',
                               'vg',
                               'dev',
                               'rvis',
                               'clinvar',  # it could be extended with reference variant - omim - disease
                               'ccr',
                               'ncRVIS',
                               'ncGERP',
                               'gene panel',
                               'HI score'),
                     level = c('gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level',
                               'gene-level'),
                     description = c('-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     'Developmental disorder genes',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     'Genomics England panel genes',
                                     '-'),
                     type = c('d',
                              'd',
                              'c',
                              'd',
                              'd',
                              'd',
                              'd',
                              'd',
                              'd',
                              'c',
                              'c',
                              'd',
                              'c',
                              'c',
                              'c',
                              'd',
                              'c'),
                     date = c('-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '14-02-2020',
                                '-',
                                '-',
                                '-',
                                '-',
                              '-',
                                '-',
                                '-'),
                     any_filter = c('-',
                              '-',
                              '-',
                              '-',
                              '-',
                              '-',
                              '-',
                              '-',
                              '-',
                              'DDD category: confirmed',
                              '-',
                              '-',
                              '-',
                              '-',
                              '-',
                              'source: Expert Review Green',
                              '-')
                     

                     )


# Load datasets

## ADDED IN MADRID

# mirtarbase
# List TFs
# drugbank curl -Lfv -o filename.zip -u francisco.requena@institutimagine.org:BElerofonte93-- https://www.drugbank.ca/releases/5-1-5/downloads/target-approved-polypeptide-ids

# ------------------------------------------------------------------------------
# Dataset: Protein-coding genes with HGCN symbol 
# Source: https://www.genenames.org/download/statistics-and-files/
# Source:  ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt
# ------------------------------------------------------------------------------

# test813 <- read_tsv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt')
# hgcn_genes <- read_excel('/home/cbl02/Storage/data/gene_with_protein_product.xlsx') %>% as_tibble()

hgcn_genes <- read_tsv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt')



hgcn_genes <- hgcn_genes %>%
  filter(status == 'Approved') %>%
  select(entrez_id, ensembl_gene_id, location, symbol) %>%
  rename(gene = symbol) %>%
  na.omit() %>%
  mutate(entrez_id = as.numeric(entrez_id)) %>%
  mutate(chrom = str_remove(location, '\\q(.*)')) %>%
  mutate(chrom = str_remove(chrom, 'p(.*)')) %>%
  filter(chrom != 'mitochondria')


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

interval_genes <- getBM(attributes = c('ensembl_gene_id', 'start_position','end_position', 'chromosome_name', 
                                      'entrezgene_id', 'version'),
                        mart = human ) %>% 
                        as_tibble() %>% 
                        filter(!str_detect(chromosome_name, 'PATCH')) %>% 
                        na.omit()
                        # rename(entrez_id = entrezgene_id)

hgcn_genes <- interval_genes %>% 
  select(-entrezgene_id) %>%
  right_join(hgcn_genes, by = c('ensembl_gene_id' = 'ensembl_gene_id', 'chromosome_name' = 'chrom')) %>%
  distinct()

no_coord_with_ensembl <- hgcn_genes %>% filter(is.na(start_position)) %>% 
  select(-start_position, -end_position, -version) %>%
  left_join(interval_genes %>% select(-ensembl_gene_id), by = c('entrez_id' = 'entrezgene_id', 'chromosome_name' = 'chromosome_name')) %>%
  mutate(length = end_position - start_position + 1) %>%
  group_by(gene) %>%
  filter(version == max(version)) %>%
  filter(length == max(length)) %>%
  slice(1) %>%
  ungroup() %>% 
  distinct() %>%
  select(-length)
  
# MAKE A THIRD ROUND WITH 336 GENES WITHOUT COORDINATES AND THEREFORE REMOVED.

hgcn_genes <- hgcn_genes %>% na.omit() %>% bind_rows(no_coord_with_ensembl)


# There are 398 genes with no coordinates



# ------------------------------------------------------------------------------
# Dataset: pLI
# Source: https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
# Access:  06/06/19
# There are 46 genes with 2 canonical transcript (canonical == TRUE)
# ------------------------------------------------------------------------------


pli <- read_tsv('/home/cbl02/Storage/data/gnomad.v2.1.1.lof_metrics.by_transcript.txt') %>%
  filter(canonical == 'TRUE') %>% 
  # filter(exp_lof < 10) %>% # take into account - important for remot gw project
  group_by(gene) %>%
  slice(1) %>%
  select(gene, pLI) %>%
  mutate(pLI = round(pLI, 2)) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Dataset: Vg - Expected population variance in its dosage that is due to genetic differences among individuals
# Source: https://www.biorxiv.org/content/10.1101/632794v1.supplementary-material - Table 1
# # Access:  06/06/19
# 5 ENTREZID linked to 2 ENSEMBL GENE
# ------------------------------------------------------------------------------


vg_raw <- read_excel('/home/cbl02/Storage/data/vg.xlsx', sheet = 3)

vg_raw <- vg_raw %>% rename(vg = Avg_VG) %>% filter(vg != 'NaN') %>% mutate(vg = as.numeric(vg)) %>%
  rename(gene = GeneID)

test <- clusterProfiler::bitr(vg_raw$gene, fromType = 'ENSEMBL', toType = "ENTREZID", OrgDb="org.Hs.eg.db")

vg <- vg_raw %>% left_join(test, by = c('gene' = 'ENSEMBL')) %>% 
  select(-gene, ENTREZID, vg) %>% 
  mutate(ENTREZID = as.numeric(ENTREZID)) %>% 
  rename(entrez_id = ENTREZID) %>%
  filter(!entrez_id %in% c(3117, 5414, 64788, 148753, 100652739))



# Expression variance

# url <- 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
# 
# download.file(url, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz')
# 
# system('gunzip GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz')
# 
# sd_gtex <- read_tsv('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct', skip = 2)



# ------------------------------------------------------------------------------
# Dataset: Pubmed articles associated with G-bands and "deletion" "duplication" keywords
# ------------------------------------------------------------------------------

pubmed_bands <- chromPlot::hg_cytoBandIdeo %>% 
                  as_tibble() %>%
                  select(Name, Chrom)

get_band <- function(band_input, chrom_input) {

  
  print(glue('{chrom_input} - {band_input}'))

  # band_input <- 'q24.3'
  # chrom_input <- '5'

  band_tmp <- band_input
  band2_tmp <- paste0(chrom_input, band_input)
  chrom_tmp <- paste('chromosome', chrom_input)

  result_del <- length(entrez_search(db="pubmed", 
                                     term= paste(chrom_tmp,'AND','(', band_tmp,'OR', band2_tmp,')', 'AND deletion AND homo sapiens'), retmax = 1000 )$ids)
  result_dup <- length(entrez_search(db="pubmed", 
                                     term= paste(chrom_tmp,'AND', '(', band_tmp, 'OR', band2_tmp, ')', 'AND duplication  AND homo sapiens'), retmax = 1000 )$ids)
  
  result <- tibble(band = band_tmp, chrom = chrom_input, hits_del = result_del, hits_dup = result_dup)
  
  return(result)
}

tic()
pubmed_lists <- pmap(list(pubmed_bands$Name, 
                     pubmed_bands$Chrom),
                                  get_band)

pubmed_df <- bind_rows(lapply(pubmed_lists, as.data.frame.list)) %>% as_tibble()

toc()

pubmed_df <- pubmed_df %>% left_join( chromPlot::hg_cytoBandIdeo %>% select(Chrom, Start, End, Name), by =
                           c('band' = 'Name', 'chrom' = 'Chrom')) %>%
  rename(start = Start, end = End) %>% 
  mutate(start = start + 1)



# ------------------------------------------------------------------------------
# Annotation promoter region (-2kbs - TSS - + 2kbs)
# ------------------------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicScores)
library(Biostrings)
phast46pla <- getGScores("phastCons46wayPlacental.UCSC.hg19")
plan("multiprocess", workers = 40)


genomic_ranges <- tibble(chrom = '1', start = 1000000, end = )


  get_annot_promoter <- function(gene, chrom, tss) {
    
    # gene <- 'A'
    # chrom <- '1'
    # tss <- 58856544
    
    tmp_granges <- GRanges(seqnames= paste0('chr', chrom), 
                           IRanges(start=(tss-2000):(tss+2000), width=1))
    
    # Calculation mean PhastCons46
    phast_temp <- gscores(phast46pla, tmp_granges)
    phast_temp <- phast_temp$default %>% mean(na.rm = TRUE)
 
    
    seqs <- getSeq(Hsapiens, paste0('chr', chrom), (tss-2000), (tss+2000))
    
    cpg_temp <- length(seqs)*dinucleotideFrequency(seqs)["CG"]/
      (letterFrequency(seqs, "C")*letterFrequency(seqs, "G"))
  
    
    result_temp <- tibble(gene = gene, mean_phast = phast_temp, cpg_density = cpg_temp )
    return(result_temp)
  }
  
  
  tic()
  
  output_temp <- future_pmap(list(hgcn_genes$gene, 
                                  hgcn_genes$chrom, 
                                  hgcn_genes$start_position), 
                             get_annot_promoter)
  
  genes_promoter <- bind_rows(lapply(output_temp, as.data.frame.list)) %>% as_tibble()
  
  toc()
  

  # hgcn_genes %>% select(gene, pLI) %>%
  #   mutate(pLI = ifelse(pLI >= 90, 'yes', 'no')) %>%
  #   left_join(genes_promoter) %>%
  #   mutate(mean_phast = ntile(mean_phast, 10)) %>%
  #   mutate(cpg_density = ntile(cpg_density, 10)) %>%
  #   pivot_longer(names_to = 'score', values_to = 'value', -c(gene, pLI)) %>%
  #   mutate(value = as.factor(value)) %>%
  #   na.omit() %>%
  #   count(pLI, score, value) %>%
  #   group_by(score, value) %>%
  #   mutate(perc = 100*(n / sum(n))) %>%
  #   ggplot(aes(value, perc)) +
  #   geom_col(aes(fill = pLI)) +
  #   facet_wrap(~score)
  # 
  # hgcn_genes %>% select(gene, disease) %>%
  #   left_join(genes_promoter) %>%
  #   mutate(mean_phast = ntile(mean_phast, 10)) %>%
  #   mutate(cpg_density = ntile(cpg_density, 10)) %>%
  #   pivot_longer(names_to = 'score', values_to = 'value', -c(gene, disease)) %>%
  #   mutate(value = as.factor(value)) %>%
  #   na.omit() %>%
  #   count(disease, score, value) %>%
  #   group_by(score, value) %>%
  #   mutate(perc = 100*(n / sum(n))) %>%
  #   ggplot(aes(value, perc)) +
  #   geom_col(aes(fill = disease)) +
  #   facet_wrap(~score)
  

  
  # ------------------------------------------------------------------------------
  # Centromeric and telomeric regions
  # ------------------------------------------------------------------------------  
  
  
  session <- browserSession("UCSC")
  genome(session) <- "hg19"
  
  
  region_gaps <- getTable(ucscTableQuery(session, "Gap"))


  region_gaps <- region_gaps %>%
    as_tibble() %>%
    filter(type %in% c('telomere', 'centromere')) %>%
    mutate(chrom = str_remove(as.character(chrom), 'chr')) %>%
    rename(start = chromStart, end = chromEnd ) %>%
    select(chrom, start, end, type)
  
  
  
  
# ------------------------------------------------------------------------------
# DNAse peaks 
# ------------------------------------------------------------------------------
# 
# session <- browserSession("UCSC")
# genome(session) <- "hg19"
# 
# 
# 
# 
# query <- ucscTableQuery(session, "DNase Clusters")
# tableName(query) <- "wgEncodeRegDnaseClusteredV3"
# dhs_df <- getTable(query) %>% as_tibble()
# 
# dhs_df <- dhs_df %>%
#   rename(start = chromStart,
#          end = chromEnd) %>%
#   select(chrom, start, end, name) %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   rename(n_cells = name)


# ------------------------------------------------------------------------------
# Ensembl Regulatory Build
# Aggregation from ENCODE, Roadmap Epigenomics and Blueprint
# Version: Ensembl 95
# ------------------------------------------------------------------------------

download.file('ftp://ftp.ensembl.org/pub/grch37/release-95/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz',
              'homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz')
  
system('gunzip homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz')
  
ensembl_reg <- read_tsv('homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff', 
                        col_names = FALSE)

file.remove('homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff')  

ensembl_reg <- ensembl_reg %>% 
  select(X1, X3, X4, X5) %>%
  filter(nchar(X1) <= 2) %>%
  rename(chrom = X1,
         type = X3,
         start = X4,
         end = X5)



# ------------------------------------------------------------------------------
# Dataset: Clingen
# Source: ftp://ftp.clinicalgenome.org/
# Access: 06/06/19
# Info score: https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/help.shtml#review
# FILTER SCORE == 3
# ------------------------------------------------------------------------------

clingen_raw <- read.table('/home/cbl02/Storage/data/ClinGen_haploinsufficiency_gene_GRCh37.bed', col.names = c('chrom', 'start', 'end', 'gene', 'score'), stringsAsFactors = FALSE,
                          skip = 1)

clingen <- clingen_raw %>% 
  as_tibble() %>% 
  filter(score == 3) %>% 
  select(gene) %>% 
  pull(gene)



# ------------------------------------------------------------------------------
# Dataset: CNV intolerance score (gene-level score)
# ------------------------------------------------------------------------------

download.file('https://storage.cloud.google.com/gnomad-public/legacy/exacv1_downloads/release0.3.1/cnv/exac-final-cnv.gene.scores071316',
              'exac-final-cnv.gene.scores071316')

cnv_int <- read_tsv('exac-final-cnv.gene.scores071316')



# ------------------------------------------------------------------------------
# Dataset: CNV Syndromes - ClinGen
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_recurrent_CNV_V1.0-hg19.bed

# ------------------------------------------------------------------------------

cnv_syndromes_clingen <- read_tsv('ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_recurrent_CNV_V1.0-hg19.bed',
                                 skip = 1, col_names = FALSE)

cnv_syndromes_clingen <- cnv_syndromes_clingen %>%
  select(X1, X2, X3, X4) %>%
  rename(chrom = X1,
         start = X2,
         end = X3,
         syndrome_name = X4
         ) %>%
  mutate(start = 1 + start)  %>% # .bed format 0-based
  mutate(source = 'clingen')

# ------------------------------------------------------------------------------
# Dataset: CNV Syndromes - DECIPHER
# Source: sftp user@sftpsrv.sanger.ac.uk/decipher-agreements/pub
# ------------------------------------------------------------------------------

cnv_syndromes_decipher <- read_tsv('/home/cbl02/Storage/data/daa_decipher/decipher-syndromes-grch37-2020-01-19.txt', skip = 1)

cnv_syndromes_decipher <-cnv_syndromes_decipher %>%
  rename(syndrome_name = `# syndrome_name`,
         chrom = chr) %>%
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', ', ')) %>%
  mutate(source = 'decipher') %>%
  select(chrom, start, end, syndrome_name, variant_class, phenotypes, source)

# ------------------------------------------------------------------------------
# Dataset: DATASET - CNV syndromes from decipher and clinvar
# ------------------------------------------------------------------------------


syndromes_total <- cnv_syndromes_clingen %>%
  bind_rows(cnv_syndromes_decipher) %>%
  mutate(chrom = str_remove(chrom, 'chr'))

# ------------------------------------------------------------------------------
# Dataset: Clingen
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_triplosensitivity_gene.bed
# Access: 06/06/19
# Info scores: https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/help.shtml#review
# ------------------------------------------------------------------------------

triplo <- read.table('/home/cbl02/Storage/data/ClinGen_triplosensitivity_gene.bed', col.names = c('chrom', 'start', 'end', 'gene', 'score'), stringsAsFactors = FALSE,
                     skip = 1, sep = '\t') %>%
  as_tibble() %>%
  filter(!str_detect(score, 'Not yet evaluated')) %>%
  filter(score >= 1) %>%
  filter(score != 40) %>% 
  select(gene) %>%
  pull(gene)

# ------------------------------------------------------------------------------
# Dataset: DEVELOPMENTAL DISORDER GENES
# Source: https://decipher.sanger.ac.uk/ddd#ddgenes
# note: Filtering by: category == 'confirmed'
# Access: 06/06/19
# ------------------------------------------------------------------------------

### IMPORTANT

# THERE ARE RELEVANT FEATURES SUCH AS MODE, CONSEQUENCE, DISEASE THAT MAY BE RELEVANT
# FOR THE APP

dev_raw <- read_csv('/home/cbl02/Storage/data/DDG2P_14_2_2020.csv', 
                      col_names = TRUE) %>%
  filter(`DDD category`  == 'confirmed') %>%
  rename(gene = `gene symbol`) %>%
  select(gene, `gene mim`, `disease name`, `allelic requirement`, `organ specificity list`, pmids)
  

dev_genes <- dev_raw %>% pull(gene)

# ------------------------------------------------------------------------------
# Dataset: OMIM 
# Source: morbidmap.txt
# Access: ####
# ------------------------------------------------------------------------------


# 1 - The disorder is placed on the map based on its association with
# a gene, but the underlying defect is not known.
# 2 - The disorder has been placed on the map by linkage or other
# statistical method; no mutation has been found.
# 3 - The molecular basis for the disorder is known; a mutation has been
# found in the gene.
# 4 - A contiguous gene deletion or duplication syndrome, multiple genes
# are deleted or duplicated causing the phenotype.
#

# 
# omim <- read_excel('/home/cbl02/Storage/data/morbidmap.xlsx', skip = 4, col_names = c('phenotype', 
#                                                                                       'gene', 'gene_id', 'location'))
# omim <- omim[1:7745,] # description of the (id) 
# 
# omim$gene <-map_chr(omim$gene, function(x) str_split(x, ',')[[1]][1])
# omim$id_pheno <-map_chr(omim$phenotype, function(x) str_extract(x, '\\d{6}'))
# # omim$map_key <-map_chr(omim$phenotype, function(x) str_extract(x, '\\(([^)]+)\\)'))
# omim$map_key <-map_chr(omim$phenotype, function(x) str_extract(x, '\\(([[0-9]]+)\\)'))
# 
# omim$map_key <-map_chr(omim$map_key, function(x) str_remove(x, '\\('))
# omim$map_key <-map_chr(omim$map_key, function(x) str_remove(x, '\\)'))
# 
# omim$phenotype <-map_chr(omim$phenotype, function(x) str_replace(x, '\\d{6}', ''))
# omim$phenotype <-map_chr(omim$phenotype, function(x) str_replace(x, '\\d{6}', ''))
# omim$phenotype <-map_chr(omim$phenotype, function(x) gsub( '*\\(.*?\\) *', '', x)) 
# omim$phenotype <-map_chr(omim$phenotype, function(x) str_remove( x, ',  ')) 
# 
# 
# omim  <- omim %>% 
#   filter(map_key %in% c(1:4)) %>%
#   pull(gene)

# omim %>% count(map_key)
# A tibble: 4 x 2
# map_key       n
#   1          68
#   2        1083
#   3        6444
#   4         148


# ------------------------------------------------------------------------------
# Dataset: OMIM 
# Source: Barthelemy script (data/OMIM/2019_06_10)
# Most of the symbols are alias
# ------------------------------------------------------------------------------

omim_filtered <- read_tsv('/home/cbl02/Storage/data/matched_genes_all_pheno.tsv')

omim <- omim_filtered %>%
  filter(quality == 3,
         somatic == 0,
         cancer == 0,
         complexity == 0) %>%
  rename(gene = Gene_symbols) %>%
  separate_rows(gene, sep = ', ')


omim_genes <- omim_filtered %>%
  filter(quality == 3,
         somatic == 0,
         cancer == 0,
         complexity == 0) %>%
  rename(gene = Gene_symbols) %>%
  select(gene) %>% 
  separate_rows(gene, sep = ', ') %>%
  distinct() %>%
  pull()
# ------------------------------------------------------------------------------
# Dataset: FDA 
# Source: MacArthur lab github
# Access: ####
# ------------------------------------------------------------------------------

fda <- read.table('https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/fda_approved_drug_targets.tsv') %>%
  as_tibble() %>%
  select(V1) %>%
  mutate(V1 = as.character(V1)) %>%
  rename(gene = V1) %>%
  pull(gene) 

# ------------------------------------------------------------------------------
# Dataset: RVIS score 
# Note: The newest RVIS release (v4) is based on the ExAC v2 standing variation 
# Source: http://genic-intolerance.org/about.jsp
# ------------------------------------------------------------------------------

rvis <- read.table('http://genic-intolerance.org/data/GenicIntolerance_v3_12Mar16.txt', sep = '\t', header = TRUE)
rvis <- rvis %>% select(GENE, X.ExAC_0.05.popn) %>%
  as_tibble() %>% 
  na.omit() %>%
  mutate(GENE = as.character(GENE)) %>%
  rename(rvis =X.ExAC_0.05.popn, gene = GENE)


# ------------------------------------------------------------------------------
# Dataset: ClinVar
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2019/
# Access: 
# ------------------------------------------------------------------------------




# MAYBE INCLUDE RS_REFERENCE (dbSNP)
clinvar_raw <- read_tsv('/home/cbl02/Storage/data/clinvar_20191219.vcf', skip = 28, 
         col_types = list(`#CHROM` = col_character()))

clinvar_genes <- clinvar_raw %>% 
  filter(str_detect(V8, 'CLNSIG=Pathogenic') | str_detect(V8, 'CLNSIG=Likely_pathogenic')) %>%
  select(V8) %>%
  mutate(gene = str_extract(V8, pattern =  'GENEINFO[^;]*')) %>%
  select(-V8) %>%
  mutate(gene = str_remove(gene, 'GENEINFO=')) %>%
  mutate(gene = str_remove(gene, '\\:.*')) %>%
  distinct() %>%
  na.omit() %>%
  pull()


clinvar_variants <- clinvar_raw %>%
  filter(str_detect(V8, 'CLNSIG=Pathogenic') | str_detect(V8, 'CLNSIG=Likely_pathogenic')) %>%
  mutate(V8 = as.character(V8)) %>%
  mutate(disease_identifier = str_extract(V8, 'CLNDISDB([^;]+)')) %>%
  mutate(disease_identifier = str_remove(disease_identifier, 'CLNDISDB=')) %>%
  mutate(disease_name = str_extract(V8, 'CLNDN([^;]+)')) %>%
  mutate(disease_name = str_remove(disease_name, 'CLNDN=')) %>%
  mutate(clinical_sign = str_extract(V8, 'CLNSIG([^;]+)')) %>%
  mutate(clinical_sign = str_remove(clinical_sign, 'CLNSIG=')) %>%
  mutate(gene = str_extract(V8, 'GENEINFO([^:]+)')) %>% 
  mutate(gene = str_remove(gene, 'GENEINFO=')) %>%
  rename(chrom = V1, pos = V2, id = V3, reference = V4, alternative = V5) %>%
  select(chrom, pos, reference, alternative, gene, clinical_sign, disease_identifier, disease_name, id) %>%
  mutate(chrom = as.character(chrom))




# ------------------------------------------------------------------------------
# Dataset: Haploinsufficiency prediction
# Source: https://decipher.sanger.ac.uk/about#downloads/data
# Access: 6/6/19
# ------------------------------------------------------------------------------


hi <- read.table('/home/cbl02/Storage/data/HI_Predictions_Version3.bed', sep = '\t', skip = 1) %>%
  as_tibble() %>%
  select(V4) %>%
  mutate(V4 = as.character(V4)) %>%
  separate(V4, into = LETTERS[1:10]) %>%
  mutate(A = paste(A, B, sep = '-')) %>%
  mutate(A = str_replace(A, '-0', '')) %>%
  select(A, D, E) %>%
  mutate(E = if_else(E == '', '00', E )) %>%
  mutate(hi = paste(D, E, sep = '.')) %>%
  mutate(hi = as.numeric(hi)) %>%
  select(-D, -E) %>%
  rename(gene = A)

# ------------------------------------------------------------------------------
# Dataset: GWAS genes
# Source:https://www.ebi.ac.uk/gwas/docs/file-downloads
# Note: file All associations v1.0
# Access: 6/6/19
# ------------------------------------------------------------------------------


gwas_raw <- read.table('/home/cbl02/Storage/data/gwas_catalog_v1.0-associations_e96_r2019-05-03.tsv', header = TRUE, sep = '\t',
                       fill = TRUE)

gwas <- gwas_raw %>% as_tibble() %>% 
  select(CHR_ID, CHR_POS, INTERGENIC, REPORTED.GENE.S.) %>% 
  mutate(CHR_ID = as.character(CHR_ID),
         CHR_POS = as.character(CHR_POS),
         INTERGENIC = as.numeric(as.character(INTERGENIC)),
         REPORTED.GENE.S. = as.character(REPORTED.GENE.S.)) %>%
  filter(CHR_ID %in% c(1:22,'X')) %>%
  separate(REPORTED.GENE.S., into = as.character(1:150), sep = ',') %>% 
  select(CHR_ID, CHR_POS, INTERGENIC, '1') %>% 
  rename('gene' = '1') %>%
  filter(INTERGENIC == 0) %>%
  select(gene) %>%
  distinct() %>%
  na.omit() %>%
  pull()


gwas_variants <- gwas_raw %>% 
  as_tibble() %>% 
  select(CHR_ID, CHR_POS, INTERGENIC, REPORTED.GENE.S., DISEASE.TRAIT, LINK) %>% 
  mutate(CHR_ID = as.character(CHR_ID),
         CHR_POS = as.character(CHR_POS),
         INTERGENIC = as.numeric(as.character(INTERGENIC)),
         REPORTED.GENE.S. = as.character(REPORTED.GENE.S.)) %>%
  filter(CHR_ID %in% c(1:22,'X')) %>%
  separate(REPORTED.GENE.S., into = as.character(1:150), sep = ',') %>% 
  select(CHR_ID, CHR_POS, INTERGENIC, DISEASE.TRAIT,'1', LINK) %>% 
  rename('gene' = '1') %>%
  mutate(CHR_POS = as.numeric(CHR_POS)) %>%
  mutate(gene = if_else(gene == 'NR', '-', gene)) %>%
  mutate(INTERGENIC = if_else(INTERGENIC == 1, 'Yes', 'No'))
  

# ------------------------------------------------------------------------------
# Dataset: CCR score - Nº of regions located in a gene that are above of P99
# Source: https://www.nature.com/articles/s41588-018-0294-6#Sec28
# Source2: https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM3_ESM.txt
# Access: 6/6/19
# ------------------------------------------------------------------------------

ccr <- read.table('https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0294-6/MediaObjects/41588_2018_294_MOESM3_ESM.txt', 
                  header = TRUE) %>%
  as_tibble() %>%
  rename(ccr = number_of_99th_percentile_CCRs) %>%
  mutate(gene = as.character(gene))

# ------------------------------------------------------------------------------
# Dataset: ncRVIS and ncGERP
# Source: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005492#sec025
# Table: S1 Data. Collection of RVIS and GERP scores and their corresponding percentile
# ------------------------------------------------------------------------------

nc_raw <- read_excel('/home/cbl02/Storage/data/journal.pgen.1005492.s011.XLSX.xlsx', sheet = 1)

nc <- nc_raw %>%
  select( `CCDS release 9` ,`CCDS release 15`, ncRVIS, `ncGERP      [Average GERP++]` ) %>%
  rename(gene1 =  `CCDS release 9`, gene2 = `CCDS release 15`, ncrvis = ncRVIS, ncgerp = `ncGERP      [Average GERP++]` ) %>%
  mutate(ncgerp = as.numeric(ncgerp)) %>%
  select(-gene1) %>%
  rename(gene = gene2) %>%
  mutate(ncrvis = round(ncrvis, 3),
         ncgerp = round(ncgerp, 3))



# ------------------------------------------------------------------------------
# Dataset: GTEx
# Source: https://gtexportal.org/home/datasets
# File name: 	GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
# ------------------------------------------------------------------------------

gtex <- read.table('/home/cbl02/Storage/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct',
                   sep = '\t', header = TRUE, skip = 2)

gtex <- gtex %>% 
  as_tibble() %>%
  rename(gene = Description) %>%
  select(-gene_id) %>%
  gather('tissue', 'value', - gene) %>%
  mutate(tissue = str_replace_all(tissue, '\\.\\.\\.', '-')) %>%
  mutate(tissue = str_replace_all(tissue, '\\.\\.', '-')) %>%
 mutate(tissue = str_replace_all(tissue, '\\.', '-'))



# ------------------------------------------------------------------------------
# Dataset: HPA
# Source: https://gtexportal.org/home/datasets
# File name: 	GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
# ------------------------------------------------------------------------------

hpa <- read_tsv('/home/cbl02/Storage/data/normal_tissue.tsv')

hpa <- hpa %>%
  rename(gene = `Gene name`, tissue = Tissue, cell_type = `Cell type`) %>%
  select(-Gene) %>%
  filter(Level != 'Not detected')

# ------------------------------------------------------------------------------
# Dataset: Mouse phenotype
# Source: http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
# ------------------------------------------------------------------------------

mgi <- read_tsv('http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt', col_names = FALSE)

mgi <- mgi %>%
  as_tibble() %>%
  select(-X8, -X3) %>%
  rename(gene = X1, entrez_id = X2, gene_mouse = X5, pheno = X7, mgi = X6) %>%
  mutate(mgi = str_remove(mgi, pattern = '  ')) %>%
  select(-X4) %>%
  filter(pheno != '')

# ------------------------------------------------------------------------------
# Dataset: OMIM
# Source: https://data.omim.org/downloads/cpJFEdlrQ5qqPk0TOzvVBA/morbidmap.txt
# ------------------------------------------------------------------------------

# 1 - The disorder is placed on the map based on its association with
# a gene, but the underlying defect is not known.
# 2 - The disorder has been placed on the map by linkage or other
# statistical method; no mutation has been found.
# 3 - The molecular basis for the disorder is known; a mutation has been
# found in the gene.
# 4 - A contiguous gene deletion or duplication syndrome, multiple genes
# are deleted or duplicated causing the phenotype.
#

# morbidmap <- read_xlsx('/home/cbl02/Storage/data/morbidmap.xlsx', skip = 3)
# 
# morbidmap <- morbidmap %>% 
#   as_tibble() %>%
#   rename(pheno = `# Phenotype`, gene = `Gene Symbols`, mim_gene = `MIM Number`) %>%
#   select(- `Cyto Location`) %>%
#   mutate(mapping = 
#            case_when(
#              str_detect(pheno, '\\([1]\\)') ~ "1",
#              str_detect(pheno, '\\([2]\\)') ~ "2",
#              str_detect(pheno, '\\([3]\\)') ~ "3",
#              str_detect(pheno, '\\([4]\\)') ~ "4"
#            )) %>%
#   filter(!isNA(mapping)) %>%  # eliminate description at the bottom
#   mutate(pheno = str_remove(pheno, '\\([1-4]\\)' )) %>%
#   mutate(mim_disease = str_extract(pheno, '[0-9]{6}')) %>%
#   mutate(pheno = str_remove(pheno, '[0-9]{6}')) %>%
#   mutate(pheno = str_remove(pheno, ',  ')) %>%
#   select(gene, pheno, mim_gene, mim_disease, mapping)


# ------------------------------------------------------------------------------
# Dataset: Paralogous genes
# Source: biomart
# ------------------------------------------------------------------------------


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")


para_genes <- getBM(attributes = c('external_gene_name', 'hsapiens_paralog_ensembl_gene'),
                        mart = human )

para_genes <- para_genes %>% 
  as_tibble() %>%
  rename(gene = external_gene_name, para = hsapiens_paralog_ensembl_gene) %>%
  filter(para != '') %>%
  count(gene)


# ------------------------------------------------------------------------------
# LiftOver
# Instructions: https://www.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
# ------------------------------------------------------------------------------


library(rtracklayer)
library(liftOver)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
from_hg38_to_hg19 = import.chain(path)

# ------------------------------------------------------------------------------
# Dataset: Recombination rates
# Source: https://science.sciencemag.org/content/suppl/2019/01/23/363.6425.eaau1043.DC1
# Genome assembly: hg38
# ------------------------------------------------------------------------------



url <- "https://science.sciencemag.org/highwire/filestream/721792/field_highwire_adjunct_files/4/aau1043_DataS3.gz"
download.file(url, 'aau1043_DataS3.gz')
system("gunzip  aau1043_DataS3.gz")  
 recomb <- read_tsv('aau1043_DataS3', skip = 8, col_names = c('chrom', 'start', 'end', 'cm_mb', 'cm'))

file.remove('aau1043_DataS3')

tmp_granges_recomb <- recomb %>% GRanges()

seqlevelsStyle(tmp_granges_recomb) = "UCSC"  # necessary
tmp_granges_recomb = liftOver(tmp_granges_recomb, from_hg38_to_hg19)

recomb <- tmp_granges_recomb %>% as_tibble() %>% select(seqnames, start, end, cm_mb, cm) %>%
  rename(chrom = seqnames) %>% mutate(chrom = str_remove(as.character(chrom), 'chr'))


# ------------------------------------------------------------------------------
# Dataset: Structural Variants (SV)
# Source: https://gnomad.broadinstitute.org/downloads
# Name file: https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2_sv.sites.bed.gz
# More info: https://macarthurlab.org/2019/03/20/structural-variants-in-gnomad/
# SVLEN a threshold?
# A tibble: 8 x 2
# SVTYPE      n
# BND     72411 # Breakends
# CPX      5249 # complex SV
# CTX         9 # Translocation
# DEL    199498 # deleted
# DUP     51428 # duplicated
# INS    115407 # insertion
# INV       707 # inversion
# MCNV     1148 # Multiallelic CNV
# Populations:
# AMR (Latino)
# EUR (European)
# OTH (Other)
# AFR (African)
# EAS (East Asian)

# ------------------------------------------------------------------------------





gnomad_sv_raw <-  read_tsv('/home/cbl02/Storage/data/gnomad_v2.1_sv.sites.bed', col_types = list(`#chrom` = col_character())) %>%
  filter(SVTYPE %in% c('DEL', 'DUP')) %>%
  filter(FILTER == 'PASS') %>%
  rename(id = name) %>%
  mutate(id  = str_remove(id, 'gnomAD-SV_v2.1_')) %>%
  mutate(source = 'gnomad_v2.1') %>%
  rename(chrom = `#chrom`) %>%
  mutate(chrom = as.character(chrom)) %>%
  mutate(start = start + 1) %>%
  select(id, chrom, start, end, svtype, AF)



gnomad_sv <- gnomad_sv_raw %>%
  select(chrom, start, end, id, source) 
  
    


# ------------------------------------------------------------------------------
# Dataset: Structural Variants (SV) - Population CNVs
# Source: https://decipher.sanger.ac.uk/about#downloads/data
# Population Copy-Number Variation Frequencies
# Variables: duplicated - deletion - general // observations - frequency - se
# ------------------------------------------------------------------------------
# 
decipher_control_raw <- read_tsv('/home/cbl02/Storage/data/population_cnv.txt', col_names = TRUE) %>%
  mutate(source = 'decipher_control') %>%
  rename(id = population_cnv_id, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  mutate(chrom = as.character(chrom)) %>%
  mutate(chrom = if_else(chrom == 23, 'X', chrom))

decipher_control <- decipher_control_raw %>%
  select(id, chrom, start, end, source)

# ------------------------------------------------------------------------------
# Dataset: CNVs PATHOGENIC Decipher
# Source: sftp user@sftpsrv.sanger.ac.uk/decipher-agreements/pub
# ------------------------------------------------------------------------------

decipher_sv_raw <- read_tsv('/home/cbl02/Storage/data/daa_decipher/decipher-cnvs-grch37-2020-01-19.txt', skip = 1) %>%
  as_tibble() %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  mutate(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  # DE NOVO CONSTITUTIVE!!! 
  filter(pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>%
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>%
  # add
  select(id, chrom, start, end, source, pathogenicity, genotype, variant_class, phenotypes)



# ------------------------------------------------------------------------------
# Dataset: DGV
# Source: http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt
# ------------------------------------------------------------------------------
## CHECK POSSIBLE MISTAKE READING DATA

dgv_df_raw <- read_tsv('/home/cbl02/Storage/data/GRCh37_hg19_variants_2016-05-15.txt',
                   col_types = list(chr = col_character())) %>%
  filter(varianttype == 'CNV') %>%
  filter(!variantsubtype %in% c('deletion', 'duplication')) %>%
  rename(id = variantaccession, chrom = chr) %>%
  mutate(source = 'dgv',
         chrom = as.character(chrom))  %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  select( -source, -varianttype, -mergedorsample,
         -mergedvariants, -supportingvariants, -samples, -platform, -cohortdescription, 
         -frequency)
  

dgv_df <- dgv_df_raw %>%
  select(id, chrom, start, end, source)

  


# ------------------------------------------------------------------------------
# Dataset: Aggregation of data from DECIPHER, gnomAD and DGV
# ------------------------------------------------------------------------------

cnv_df <- decipher_sv_raw %>% 
  bind_rows(gnomad_sv) %>% 
  bind_rows(dgv_df) %>%
  bind_rows(decipher_control) %>%
  mutate(length_cnv = end - start + 1) %>% 
  mutate(phenotypes = replace_na(phenotypes, '-'))



# ------------------------------------------------------------------------------
# Dataset: Protein-Protein interaction network
# Source: https://www.intomics.com/inbio/map.html#downloads
# Name file: inBio_Map_core_2016_09_12.tar.gz
# ------------------------------------------------------------------------------

# inbio_network_raw <- read.table('/home/cbl02/Storage/data/core.psimitab', sep = '\t', header = FALSE)
# 
# inbio_network <- inbio_network_raw %>%
#   as_tibble() %>%
#   slice(1:100)















# ------------------------------------------------------------------------------
# Dataset: TADs
# Name file: http://promoter.bx.psu.edu/hi-c/publications.html
# ------------------------------------------------------------------------------

tad <- read_tsv('/home/cbl02/Storage/data/H1-ESC_Dixon2015-raw_TADs.txt', col_names = FALSE)

tad <- tad %>%
  rename(chrom = X1, start = X2, end = X3) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  mutate(id = row_number()) %>% select(id, chrom, start, end)


# ------------------------------------------------------------------------------
# Dataset: Genehancer
# Version: 4.11
# Name file: Email from Marilyn (genenhancer team)
# Genome reference: hg38
# ------------------------------------------------------------------------------

genehancer <- read.table('/home/cbl02/Storage/data/genehancer_V4_11.gff', header = FALSE, sep = '\t',
                  stringsAsFactors = FALSE)

df <- genehancer %>%
  as_tibble() %>%
  filter(V3 == 'Enhancer') %>%
  select(-V2, -V3, -V7, -V8) %>%
  rename(chrom = V1, start = V4, end = V5, score_enh = V6, attributes = V9)



a <- str_split(df$attributes, ';')
a <- lapply(a, function(x) x[[1]][1])
a <- unlist(a)
a <- str_remove(a, 'genehancer_id=')
df$id <- a
df_prev <- df %>% select(attributes)
df <- df %>% select(-attributes)

df_assoc <- df_prev

a <- str_split(df_assoc$attributes, ';')
df_assoc$id <- unlist(lapply(a, function(x) x[[1]]))
df_assoc <- df_assoc %>% select(id, attributes)
df_assoc$id <- str_remove(df_assoc$id, 'genehancer_id=')

df_assoc <- separate_rows(df_assoc, attributes, sep = ';')
df_assoc <- df_assoc %>% filter(!str_detect(attributes, "genehancer_id="))

df_assoc_gene <- df_assoc %>% filter(str_detect(attributes, 'connected_gene'))
df_assoc_score <- df_assoc %>% filter(str_detect(attributes, 'score'))

df_assoc_gene$score <- df_assoc_score$attributes
colnames(df_assoc_gene) <- c('id', 'gene', 'score')

df_assoc_gene$gene <- str_remove(df_assoc_gene$gene, 'connected_gene=')
df_assoc_gene$score <- str_remove(df_assoc_gene$score, 'score=')
df_assoc_gene$score <- as.numeric(df_assoc_gene$score)

#
glimpse(df)
glimpse(df_assoc_gene)

#
df_assoc_enh <- df %>% filter(score_enh > 1)
df_assoc_gene <- df_assoc_gene %>% filter(score > 1)

df_ge <- df_assoc_enh %>% left_join(df_assoc_gene, by = 'id')
df_ge <- df_ge %>% select(gene, everything())
df_ge <- df_ge %>% na.omit()
df_ge <- df_ge %>% filter(chrom != 'chrY')
# 
# write.table(df_ge %>% select(chrom, start, end, id) %>% distinct(), 'enhancer_cnvxplore', quote = FALSE, row.names = FALSE,
#             col.names = FALSE, sep = '\t', append = TRUE)

# Used liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver). Default parameters
# from Grch38 to Grch37
# Succesfully converted 27514 records
# Conversion failed on 8 records.

enhancer_raw <- read.table('/home/cbl02/Storage/data/hglft_genome_4c90f_a0b3d0.bed', header = FALSE, sep = '\t')

df_enhancers <- enhancer_raw %>% as_tibble() %>% rename(chrom = V1, start = V2, end = V3, id = V4) %>% mutate(chrom = str_remove(chrom, 'chr')) %>%
  mutate(id = as.character(id))

df_enhancers <- df_ge %>% select(gene, id, score_enh, score) %>% left_join(df_enhancers, by = 'id')

# write.table(enhancer_raw %>% filter(!str_detect(V1, 'PATCH')) %>% distinct(),
#             '/home/cbl02/Storage/data/enhancer_cnvxplorer', quote = FALSE, row.names = FALSE,
#             col.names = FALSE, sep = '\t', append = TRUE)


from_python <- read.table('/home/cbl02/Storage/data/enhancer_cnvxplorer_cons', header = TRUE, 
                     sep = '\t', stringsAsFactors = FALSE) %>% as_tibble()

from_remot <- read.table('/home/cbl02/Storage/data/enhancer_cnvxplorer_3_result.txt', header = TRUE,
                         sep = '\t', stringsAsFactors = FALSE) %>% 
                as_tibble()  %>%
                mutate(meanCov =  as.double(meanCov)) %>%
                mutate(oe = obs / expProba) %>%
                mutate(oe = ntile(oe, 100)) %>%
                select(-meanCov, -expProba, -obs)

 

from_python <- from_python %>% 
  mutate(phast100 = round(phast100, 3)) %>% 
  mutate(phast46pla = round(phast46pla, 3)) %>% 
  mutate(phast46pri = round(phast46pri, 3))


df_enhancers <- df_enhancers %>% 
  left_join(from_python %>% select(id, phast100, phast46pla, phast46pri), by = 'id') %>%
  left_join(from_remot, by = c('chrom', 'start', 'end'))

df_enhancers <- df_enhancers %>% mutate(phast100 = round(phast100, 2)) %>%
  mutate(phast46pla = round(phast46pla, 2)) %>%
  mutate(phast46pri = round(phast46pri, 2))
  

plot_p100 <- df_enhancers %>% select(id, phast100) %>% distinct() %>% ggplot(aes(phast100)) + geom_density() +
  xlab('Phast100way score')
plot_p46pla <-  df_enhancers %>% 
  select(id, phast46pla) %>% 
  distinct() %>% 
  ggplot(aes(phast46pla)) + 
  geom_density() +
  xlab('Phast46way score')
  
# 
# enh_post <- mod_remot('from_remot/enhancer_apolo_crossmap_cleaned_7_result.txt', 'EUR', 7, TRUE)
# 
# enh_post <- remove_duplicated_regions(enh_post, seg_dup, self_chain)
# 
# crossmap_id <- read.table('to_remot/enhancer_apolo_crossmap_cleaned')
# 
# enh_def <- enh_post %>% left_join(crossmap_id, by = c('chrom' = 'V1', 'start' = 'V2', 'end' = 'V3')) %>%
#   left_join(df_ge %>% select(gene, id), by = c('V4' = 'id'))
# 
# 
# test <- bitr(enh_def$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% 
#   as_tibble()
# 
# enh_def <- enh_def %>% 
#   left_join(test, by = c('gene' = 'SYMBOL')) %>% 
#   select(-gene) %>% 
#   na.omit() %>% 
#   rename(gene = ENTREZID) %>%
#   mutate(gene = as.numeric(gene)) %>%
#   mutate(length = end - start) %>%
#   group_by(gene) %>%
#   sample_n(1) %>%
#   # filter(length == max(length)) %>%
#   ungroup() %>%
#   mutate(pos = round((end + start) /2, 0))


# ------------------------------------------------------------------------------
# Dataset: STRING
# Source: https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz
# Version: 11
# ------------------------------------------------------------------------------


library(tidygraph)

download.file('https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz', 
              '9606.protein.links.v11.0.txt.gz')
system('gunzip 9606.protein.links.v11.0.txt.gz')
string_db <- read_delim('9606.protein.links.v11.0.txt', delim = ' ', skip = 1, col_names = c('p1', 'p2', 'score'))
file.remove('9606.protein.links.v11.0.txt')


string_db <- string_db %>%
  filter(score >= 700) %>%
  mutate(p1 = str_remove(p1, '9606.'),
         p2 = str_remove(p2, '9606.')) %>%
          select(-score)

human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

interval_genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_peptide_id'),
                        mart = human ) %>% 
  as_tibble()


string_db <- as_tbl_graph(string_db, directed = FALSE) %>%
  mutate(page_rank = centrality_pagerank(),
         degree = centrality_degree()) %>%
  activate(nodes) %>%
  as_tibble()

string_db <- string_db %>% left_join(interval_genes, by = c('name' = 'ensembl_peptide_id')) %>%
  select(hgnc_symbol, page_rank, degree, -name) %>%
  rename(gene = hgnc_symbol)

# ------------------------------------------------------------------------------
# Dataset: lncRNA
# Source: Filtered data (species == 'Human') - https://apps.kaessmannlab.org/lncRNA_app/
# Genome reference: hg19
# ------------------------------------------------------------------------------

lncrna <- read.table('/home/cbl02/Storage/data/lncRNA/filtered_data.csv', header = TRUE, sep = ',',
                     stringsAsFactors = FALSE)

lncrna <- lncrna %>%
  as_tibble() %>%
  select(-X, -Species ) %>%
  rename(id = XLOC.id, chrom = chromosome, dynamic = Dynamic, ensembl75_id = ENSEMBL75.id, name = LncRNA.name, conservation = Minimum.age, genomic_class = Genomic.class,
         nearest_coding_gene = Nearest.coding.gene)


lncrna_coord <- read.table('/home/cbl02/Storage/data/lncRNA/human.lncRNA.gtf', sep = '\t',
                           stringsAsFactors = FALSE)

lncrna_coord <- lncrna_coord %>%
  as_tibble() %>%
  separate(V9, into = LETTERS[8:14], sep = ';') %>%
  select(-V3, -M,  -K, -V6, -V8, -N, -L, -V2, -I) %>%
  mutate(H = str_remove(H, 'gene_id ')) %>%
  rename(chrom = V1, start = V4, end = V5, sense = V7, n_exon = J, id = H ) %>%
  mutate(n_exon = as.numeric(str_remove(n_exon, 'exon_number '))) %>%
  select(id, chrom, start, end, sense, n_exon)



## 7 organs: Brain, Cerebellum, Heart, Kidney, Liver, Ovary, Testis
## 26 developmental step: "10wpc, 11wpc, 12wpc, 13wpc, 16wpc, 18wpc, 19wpc, 20wpc, 4wpc, 5wpc, 6wpc, 7wpc, 
# 8wpc, 9wpc, infant, newborn, olderMidAge, oldTeenager, school, senior, Senior, teenager, toddler, youngAdult, 
# youngMidAge, youngTeenager"

# lncrna_expression <-  read.table('/home/cbl02/Storage/data/lncRNA/HumanRPKMs.txt', header = TRUE, sep = '\t',
#                                  stringsAsFactors = FALSE) %>%
#   rownames_to_column(var = 'id') %>%
#   as_tibble() %>%
#   gather('type', 'rpkm', -id) %>%
#   separate(col = type, into = c('organ', 'dev_step', 'dup'), sep = '\\.') %>%
#   mutate(dev_step = if_else(dev_step == 'Senior', 'senior', dev_step)) %>% # there was a typo in 85037 rows (Senior) and the rest, senior (510222)
#   filter(!str_detect(string = id, pattern = 'ENSG')) %>%
#   filter(str_detect(id, 'XLOC'))



# ------------------------------------------------------------------------------
# Dataset: seg_dup and self_chain
# Description: Bed files with segmental duplications and self-chain of the genome(hg19)
# Source: http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
# Source: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chainSelf.txt.gz.
# Source with name of each field: http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgsid=1432484_CeWhnfzDUkBhJCkZ00iaiTQcyb3E&hgta_doSchemaDb=hg19&hgta_doSchemaTable=chainSelf
# Comment: the score 90 is obtained from Quinlab (A map of constraint regions in coding part)
# ------------------------------------------------------------------------------

# seg_dup <- read.table('http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab', header = TRUE) %>% 
#   as_tibble()
# 
# seg_dup1 <- seg_dup %>% select(chrom, chromStart, chromEnd) %>% 
#   rename(start = chromStart, end = chromEnd)
# seg_dup2 <- seg_dup %>% select(otherChrom, otherStart, otherEnd) %>%
#   rename(chrom = otherChrom, start = otherStart, end = otherEnd)
# seg_dup <- seg_dup1 %>% bind_rows(seg_dup2) %>% distinct()  %>% mutate(chrom = str_remove(chrom, 'chr')) %>%
#   makeGRangesFromDataFrame()
# 
# # Self-chain coordinates
# 
# self_chain <- read.table('data/chainSelf.txt', header = FALSE) %>% as_tibble()
# colnames(self_chain) <- c('bin', 'score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStrand',
#                           'qStart', 'qEnd', 'id', 'normScore')
# 
# # NormScore = score / length bases. This score threshold is based on the paper: 
# self_chain <- self_chain %>% filter(normScore >= 90)
# 
# # Clean dataframe and convert it to a GRanges object
# 
# self_chain1 <- self_chain %>% select(tName, tStart, tEnd) %>% 
#   rename(chrom = tName, start = tStart, end = tEnd)
# self_chain2 <- self_chain %>% select(qName, qStart, qEnd) %>%
#   rename(chrom = qName, start = qStart, end = qEnd)
# self_chain <- self_chain1 %>% bind_rows(self_chain2) %>% distinct()  %>% mutate(chrom = str_remove(chrom, 'chr')) %>%
#   makeGRangesFromDataFrame()
# 

# ------------------------------------------------------------------------------
# Dataset: blacklist.v2.bed
# Source: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz
# ------------------------------------------------------------------------------

blacklist_encode <- read_tsv('/home/cbl02/Storage/data/hg19-blacklist.v2.bed', col_names = FALSE) %>%
  as_tibble() %>%
  rename(chrom = X1, start = X2, end = X3, class = X4) %>%
mutate(start = start + 1) %>%
    mutate(chrom = str_remove(chrom, 'chr'))
  


# ------------------------------------------------------------------------------
# Dataset: Essential genes (mouse model + cell line)
# Source: https://www.nature.com/articles/s41467-020-14284-2
# ------------------------------------------------------------------------------

# # 3262 genes found it / 3326 total nº genes
# e_mgi <- read.table('/home/cbl02/Storage/data/essential_gene_lists/mgi', header = FALSE,
#                     stringsAsFactors = FALSE) %>% as_tibble()
# 
# # 1931 genes found it / 2010 total nº genes
# 
# e_invitro <- read.table('/home/cbl02/Storage/data/essential_gene_lists/invitro', header = FALSE,
#                     stringsAsFactors = FALSE) %>% as_tibble()
# 
# e_intersect <- e_mgi %>% filter(V1 %in% e_invitro$V1)

# Cellular lethal (CL)
# Developmental lethal (DL)
# Subviable (SV)
# Viable with phenotype (VP)
# Viable with no phenotype (VN)

human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

interval_genes <- getBM(attributes = c('hgnc_symbol', 'hgnc_id'),
                        mart = human ) %>% 
  as_tibble() %>% 
  na.omit()

listAttributes(human)

fusil_score <- read_tsv('https://www.ebi.ac.uk/biostudies/files/S-BSST293/u/SourceDataFile1_FUSIL_bins.txt')


fusil_score <- fusil_score %>% select(HGNC_ID, FUSIL) %>%
  mutate(HGNC_ID = as.integer(str_remove(HGNC_ID, 'HGNC:'))) %>%
  left_join(interval_genes, by = c('HGNC_ID' = 'hgnc_id')) %>%
  filter(FUSIL != '-')


# ------------------------------------------------------------------------------
# Dataset: GDI score
# Source: Paper: The human gene damage index as a gene-level approach to prioritizing exome variants
# Explanation: High GDI values reflect highly damaged genes. 
# ------------------------------------------------------------------------------

gdi <- read.table('/home/cbl02/Storage/data/GDI_full_10282015.txt', skip = 1, stringsAsFactors = FALSE) %>% as_tibble()
colnames(gdi) <- c('gene', 'gdi', 'gdi_phred', 'all_diseases', 'all_mendelian', 'mendelian_ad', 
                   'mendelian_ar', 'all_pid', 'pid_ad', 'pid_ar', 'all_cancer', 'cancer_dominant',
                   'cancer_recessive')

gdi <- gdi %>% 
  select(gene, gdi) %>% 
  mutate(gdi = round(gdi, 2))

# ------------------------------------------------------------------------------
# Dataset: SNIPre
# Source: antonio's email
# two genes are duplicated and have wrong symbol (01-Mar, 02-Mar)
# ------------------------------------------------------------------------------

snipre <- read_tsv('/home/cbl02/Storage/data/Sniprefallgenes.txt', col_names = TRUE)

snipre <- snipre %>% select(Gene, SnIPRE.f) %>% 

    rename(gene = Gene, snipre = SnIPRE.f) %>% 
  mutate(snipre = str_replace(snipre, ',', '.')) %>%
  mutate(snipre = as.numeric(snipre)) %>%
  mutate(snipre = round(snipre, 2)) %>%
  filter(!gene %in% c('01-Mar', '02-Mar'))


# ------------------------------------------------------------------------------
# Dataset: Phenotype terms associated with OMIM diseases
# Source: https://hpo.jax.org/app/download/annotation
# Explanation: phenotype_annotation_hpoteam.tab: contains annotations made explicitly
# and manually by the HPO-team (mostly referring to OMIM entries)
# https://hpo-annotation-qc.readthedocs.io/en/latest/annotationFormat.html#phenotype-hpoa-format
# Annotation: https://hpo.jax.org/app/help/annotations
# ------------------------------------------------------------------------------


# what happens with the evience?
# what happens with the all frequent - frequent?
# clinical modifier?
# onset?
# make the intersect with the list of omim filtered genes
# update the list with the account new version from OMIM
# what is a qualifier?

url <- 'http://compbio.charite.de/jenkins/job/hpo.annotations.current/lastSuccessfulBuild/artifact/current/phenotype.hpoa'

hpo_omim <- read_tsv(url, col_names = TRUE, skip = 4) %>%
  rename(mim_disease = DatabaseID,
         desc = DiseaseName,
         hp = HPO_ID) %>%
  mutate(disease_source = str_extract(mim_disease, '^[^:]+:\\s*')) %>%
  mutate(disease_source = str_remove(disease_source, ':')) %>%
  mutate(mim_disease = str_remove(mim_disease, 'OMIM:')) %>%
  mutate(mim_disease = str_remove(mim_disease, 'DECIPHER:')) %>%
  mutate(mim_disease = str_remove(mim_disease, 'ORPHA:')) %>%
    mutate(desc = str_remove(desc, '[0-9]{6}'),
           desc = str_remove(desc, '#')) %>%
  mutate(mim_disease = as.numeric(mim_disease)) %>%
  rename(identifier = mim_disease) 

  
# ------------------------------------------------------------------------------
# Dataset: Human Phenotype Ontology
# Source: http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
# Name file: phenotype_annotation.tab
# ------------------------------------------------------------------------------

# ERROR - 10 ROWS MISSING!!

url <- 'http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/genes_to_phenotype.txt'
hpo_genes <- read_tsv(url, skip = 1, col_names = FALSE)

hpo_genes <- hpo_genes %>% select(X2, X3, X4, X9) %>%
  mutate(disease_source = str_extract(X9, '^[^:]+:\\s*')) %>%
  mutate(X9 = str_remove(X9, '^[^:]+:\\s*')) %>%
  mutate(disease_source = str_remove(disease_source, ':')) %>%
  rename(gene = X2, hp = X3, desc = X4, identifier = X9) %>%
  mutate(identifier = as.double(identifier))

# hpo_omim <- hpo_omim %>%
#   filter(X1 == 'OMIM') %>%
#   # select(-X2) %>%
#   select(-X1, -X10, -X12, -X13, -X14) %>% # remove x1 because all values are omim - change if we include orphanet/omim - the other file
#   rename(mim_disease = X2,
#          desc = X3,
#          pace = X4,
#          hp = X5,
#          mim_gene = X6,
#          evidence = X7,
#          onset = X8,
#          frequency = X9,
#          clinical_modifier = X11) %>%
#   mutate(desc = str_remove(desc, '[0-9]{6}'),
#          desc = str_remove(desc, '#'))
# ------------------------------------------------------------------------------
# Dataset: Gene panel
# Source: hhttps://panelapp.genomicsengland.co.uk
# header:  web - subhead: Downloading Gene Panels
# We got 306 (5 less than expected) because these panels were empty:
# "gene_panel_218" "gene_panel_520" "gene_panel_530" "gene_panel_561" "gene_panel_82" 
# ------------------------------------------------------------------------------
# 
# library(rvest)
# website <- "https://panelapp.genomicsengland.co.uk/panels/"
# page <- read_html(website)
# 
# c_ref <- page %>%
#   html_nodes("a") %>%       # find all links
#   html_attr("href")
# 
# df_ref <- tibble(ref = c_ref, id = NA) %>%
#   filter(str_detect(ref, 'download')) %>%
#   mutate(ref = str_remove(ref, '/panels/')) %>%
#   mutate(id = ref) %>%
#   mutate(id = str_remove(id, '/download/01234/'))
# 
setwd('/home/cbl02/Storage/data/gene_panel')
# 
# walk2(df_ref$ref, df_ref$id, function(a, b)
#   download.file(url = paste0(website, a), destfile = paste0('gene_panel_', b))
#   )
files_panel <- list.files()
panel_total <- tibble()

for (i in 1:length(files_panel)) {
  print(i)
  df_tmp <- read_tsv(paste0('/home/cbl02/Storage/data/gene_panel/', files_panel[i]))
  df_tmp <- df_tmp %>% mutate(source = files_panel[i])
  panel_total <- rbind(panel_total, df_tmp)
}

# we filtered out those genes not containing a "review ranking" ( 3,233 out 45488)
# Filtering out genes with a evidence level (red - amber)
panel_total <- panel_total %>%
  rename(entity_name = `Entity Name`,
         entity_type = `Entity type`,
         gene = `Gene Symbol`,
         sources = `Sources(; separated)`) %>%
  filter(entity_type == 'gene') %>%  # optional - we can include regions in our analysis
  filter(str_detect(sources, 'Expert Review')) %>%
  separate_rows(sources, sep = ';') %>%
  filter(str_detect(sources, 'Expert Review Green')) %>%
  select(gene, Level4, -sources, source, Phenotypes)

panel_total_genes <- panel_total %>% pull(gene) %>% unique()

setwd('/home/cbl02/Storage/cnvxplore')

# ------------------------------------------------------------------------------
# Dataset: Orphadata - RARE DISEASES WITH THEIR ASSOCIATED GENES
# Source: http://www.orphadata.org/data/xml/en_product6.xml
# ------------------------------------------------------------------------------
library(xml2)

orphanet_raw <- read_tsv("/home/cbl02/Storage/data/orphanet_test.tsv" )

orphanet_raw <- orphanet_raw %>%
  select(Symbol, OrphaNumber5, Name6, OrphaNumber, Name, SourceOfValidation, Name15) %>%
  distinct() %>%
  rename(gene = Symbol)

orphanet_genes <- orphanet_raw %>% pull(gene) %>% unique()
  

# ------------------------------------------------------------------------------
# Source: /home/cbl02/Storage/data/curated_gene_disease_associations.tsv
# http://www.disgenet.org/static/disgenet_ap1/files/downloads/readme.txt
# Readme file: https://www.disgenet.org/static/disgenet_ap1/files/downloads/readme.txt
# ------------------------------------------------------------------------------

# geneId 		-> NCBI Entrez Gene Identifier
# geneSymbol	-> Official Gene Symbol
# DSI		-> The Disease Specificity Index for the gene
# DPI		-> The Disease Pleiotropy Index for the gene
# diseaseId 	-> UMLS concept unique identifier
# diseaseName 	-> Name of the disease	
# diseaseType  	-> The DisGeNET disease type: disease, phenotype and group
# diseaseClass	-> The MeSH disease class(es)
# diseaseSemanticType	-> The UMLS Semantic Type(s) of the disease
# score		-> DisGENET score for the Gene-Disease association
# EI		-> The Evidence Index for the Gene-Disease association
# YearInitial	-> First time that the Gene-Disease association was reported
# YearFinal	-> Last time that the Gene-Disease association was reported
# NofPmids	-> Total number of publications reporting the Gene-Disease association
# NofSnps		-> Total number of SNPs associated to the Gene-Disease association
# source		-> Original source reporting the Gene-Disease association


disgenet <- read_tsv('/home/cbl02/Storage/data/curated_gene_disease_associations.tsv')

disgenet %>%
  filter(source %in% c('CLINGEN', 'GENOMICS_ENGLAND', 'ORPHANET'))

genes_disgenet <- disgenet %>% select(geneSymbol) %>% distinct() %>% pull()


# ------------------------------------------------------------------------------
# Dataset: Imprinting genes
# Source: http://www.geneimprint.com/site/genes-by-species
# ------------------------------------------------------------------------------
# 
# imp_genes <- read_tsv('/home/cbl02/Storage/data/imprinted_genes')
# 
# imp_genes <- imp_genes %>% 
#   separate(Aliases, into = LETTERS[1:25], sep = ',') %>%
#   gather('remove', 'gene', -Location, -Status, -`Expressed Allele`) %>%
#   select(-remove) %>%
#   rename(expressed_allele = `Expressed Allele` ) %>%
#   select(gene, everything()) %>%
#   distinct()

# ------------------------------------------------------------------------------
# Dataset: De novo variants (denovo-db)
# Source: http://denovo-db.gs.washington.edu/denovo-db/Download.jsp
# Version: 1.6.1
# non-SSC Samples	
# ------------------------------------------------------------------------------

denovo_raw <- read_tsv('/home/cbl02/Storage/data/denovo-db.non-ssc-samples.variants.v.1.6.1.tsv', skip = 1, col_names = TRUE,
                       col_types = list(Chr = col_character()) )

denovo <- denovo_raw %>% 
  as_tibble() %>%
  filter(FunctionClass %in% c('frameshift', 'splice-donor', 'missense', 'stop-gained', 'start-lost', 'splice-acceptor')) %>%
  select(Chr, Position, Gene, PrimaryPhenotype, StudyName, PubmedID, FunctionClass, CaddScore, LofScore) %>%
  rename(chrom = Chr) %>%
  mutate(LofScore = na_if(LofScore, -1)) %>%
  mutate(CaddScore = na_if(CaddScore, -1))

# ------------------------------------------------------------------------------
# Dataset: .vcf file
# Source: https://github.com/HudsonAlpha/UDN_SV_export/blob/master/data/UDN_dump_20190320.vcf.gz
# Genome_assembly: hg19
# ------------------------------------------------------------------------------

# https://github.com/HudsonAlpha/UDN_SV_export/blob/master/data/UDN_dump_20190320.vcf.gz




# ------------------------------------------------------------------------------
# Dataset: Structural Variants (dbVar)
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/
# Guide: https://github.com/ncbi/dbvar/blob/master/Structural_Variant_Sets/Nonredundant_Structural_Variants/ToolGuide.md
# Genome reference: hg38
# Sources: 1000 Genomes Consortium + gnomAD + Sudmant2015 + Genome_in_a_Bottle
# ------------------------------------------------------------------------------

dbvar_deletion_patho <- fread('ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh38.nr_deletions.pathogenic.tsv.gz',
                                 skip = 1) %>% as_tibble()
dbvar_deletion_common <- fread('ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh38.nr_deletions.common.tsv.gz', 
                                   skip = 1) %>% as_tibble()
dbvar_duplication_patho <- fread('ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/duplications/GRCh38.nr_duplications.pathogenic.tsv.gz',
                                 skip = 1) %>% as_tibble()
dbvar_duplication_common <- fread('ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/duplications/GRCh38.nr_duplications.common.tsv.gz', 
                                  skip = 1) %>% as_tibble()

dbvar <- dbvar_deletion_patho %>%
  rbind(dbvar_deletion_common, dbvar_duplication_patho, dbvar_duplication_common)

dbvar <- dbvar %>%
  rename(chrom = `#chr`, start = outermost_start, end = outermost_stop) %>%
  filter(clinical_assertion %in% c('Pathogenic', NA)) %>%
  filter(!str_detect(study, ';')) %>%
  filter(variant_type != 'copy_number_variation') %>%
  mutate(clinical_assertion = if_else(is.na(clinical_assertion), 'control', 'pathogenic')) %>%
  mutate(category_variant = case_when(
    str_detect(variant_type, 'deletion') ~ 'deletion',
    str_detect(variant_type, 'gain') ~ 'duplication',
    str_detect(variant_type, 'duplication') ~ 'duplication',
    str_detect(variant_type, 'loss') ~ 'deletion'
  )) %>%
  select(chrom, start, end, clinical_assertion, category_variant)


# ------------------------------------------------------------------------------
# Dataset: ACMG 59 GENES - I guess 100% overlapping with OMIM genes
# Source: https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/
# ------------------------------------------------------------------------------

# acmg_genes <- read_tsv('/home/cbl02/Storage/data/acmg_genes.tsv', col_names = FALSE) %>%
#   rename(disease_name = X1, gene = X3) %>%
#   select(gene, disease_name) %>%
#   mutate(MIM_gene = str_extract(gene, '\\(([^\\)]+)\\)')) %>%
#   mutate(gene = str_remove(gene, '\\(([^\\)]+)\\)')) %>%
#   filter(disease_name != 'LDLR(MIM 606945)') %>% # field with two genes
#   bind_rows(tibble(gene = 'LDLR', disease_name = 'Familial hypercholesterolemia (MIM 143890)', MIM_gene = '(MIM 606945)' ))



# ------------------------------------------------------------------------------
# Get ontologies (Gene ontology (GO), The Mammalian Phenotype Ontology (mp))
# https://github.com/obophenotype/upheno
# ------------------------------------------------------------------------------

mpo_dbs <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/mp.obo')
mpo_dbs <- mpo_dbs$name %>% as_tibble(rownames = 'term') %>% rename(description = value)
hpo_dbs <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/hp.obo')
uberon_dbs <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/uberon.obo')

test <- ontologyIndex::get_OBO('https://raw.githubusercontent.com/obophenotype/upheno-dev/master/modules/hp-mp/upheno.obo')

# ------------------------------------------------------------------------------
# Match HP - MP
# https://github.com/obophenotype/upheno
# ------------------------------------------------------------------------------

match_hp_mp <- read_tsv('https://raw.githubusercontent.com/obophenotype/upheno/master/mappings/hp-to-mp-bestmatches.tsv',
                        col_names = c('hp_term', 'hp_label','mp_term','mp_label', 'equivalence_score', 'subclass_score'))

# ------------------------------------------------------------------------------
# Dataset: Human Phenotype Ontology
# Description: Input vector chosen by the user (tab - Phenotypic analysis)
# ------------------------------------------------------------------------------

vector_total_terms <- hpo_dbs$name %>% 
  enframe(name = 'term') %>%  
  mutate(term_desc = paste(value, '-', term)) %>%
  filter(!term %in% 'HP:0000005')

mode_inheritance <- c("Autosomal dominant inheritance", "Gonosomal inheritance", "Multifactorial inheritance", "Uniparental disomy", "Contiguous gene syndrome", "Genetic anticipation", "Somatic mutation", "Autosomal recessive inheritance", "Heterogeneous", "Sporadic", "Mitochondrial inheritance", "Semidominant mode of inheritance")


vector_inheritance <- vector_total_terms %>% filter(value %in% mode_inheritance) %>% select(-term_desc)

vector_inheritance <- split(c('Any', vector_inheritance$term), c('Any', vector_inheritance$value))

vector_total_terms <- vector_total_terms %>% filter(!value %in% c(mode_inheritance, 'All'))

# ------------------------------------------------------------------------------
# Dataset: Anatomy entities from HPO database
# ------------------------------------------------------------------------------

vector_main_hpo <- c("Abnormality of the skeletal system", "Abnormality of limbs", "Abnormality of the nervous system", "Abnormality of metabolism/homeostasis", "Abnormality of head or neck", "Abnormality of the cardiovascular system", "Abnormality of the eye", "Abnormality of the genitourinary system", "Abnormality of the integument", "Abnormality of the immune system", "Abnormality of the digestive system", "Abnormality of blood and blood-forming tissues", "Neoplasm", "Abnormality of the musculature", "Abnormality of the respiratory system", "Abnormality of the endocrine system", "Abnormality of the ear", "Abnormal cellular phenotype", "Abnormality of connective tissue", "Abnormality of prenatal development or birth", "Growth abnormality", "Constitutional symptom", "Abnormality of the breast", "Abnormality of the voice", "Abnormality of the thoracic cavity")


anato_df <- tibble(value = vector_main_hpo) %>% left_join(enframe(hpo_dbs$name), by = 'value') %>%
  rowwise() %>%
  mutate(value = str_remove(value, 'Abnormality of the '),
         value = str_remove(value, 'Abnormal of '),
         value = str_remove(value, 'Abnormality of '),
         value = paste(toupper(substring(value, 1,1)), substring(value, 2), 
                       sep="", collapse=" ")) %>% ungroup()

# ------------------------------------------------------------------------------
# Dataset: CiViC
# Source : https://civicdb.org/releases
# ------------------------------------------------------------------------------

civic_raw <- read_tsv('https://civicdb.org/downloads/01-Mar-2020/01-Mar-2020-ClinicalEvidenceSummaries.tsv')

civic_raw %>% count(variant_types) %>% arrange(desc(n)) %>% View()


variant_summary <- read_tsv('https://civicdb.org/downloads/01-Mar-2020/01-Mar-2020-VariantSummaries.tsv')

# ------------------------------------------------------------------------------
# mirtarbase
# version 8
# ------------------------------------------------------------------------------


vector_colnames <- c('id', 'name','specie_mirna', 'gene_symbol', 'gene_entrezid', 'specie_target', 'experiment', 'support_type',
                     'references')

mirna_raw <- read_excel('C:/Users/Requena/Desktop/miRTarBase_MTI.xlsx', 
                        col_names = vector_colnames, skip = 1)


mirna <- mirna_raw %>%
  filter(specie_mirna == 'Homo sapiens' & specie_target == 'Homo sapiens') %>%
  filter(support_type == 'Functional MTI') %>%
  select(-support_type, -specie_mirna , -specie_target)

mirna_coord <- read_tsv('C:/Users/Requena/Desktop/hsa.gff3', skip = 13, col_names = LETTERS[1:9]) %>% 
  mutate(name = str_extract(I, 'Name=[^;]*')) %>%
  mutate(name = str_remove(name, 'Name=')) %>%
  # mutate(name = tolower(name)) %>%
  # filter(C == 'miRNA_primary_transcript') %>%
  select(A, D, E, name) %>%
  rename(chrom = A, start = D, end = E) %>%
  mutate(start = start - 1) # Convert to 0-based input liftover

# Failed 7 coordinates
# Final result: 4794 coordinates

from_liftover <- read_tsv('C:/Users/Requena/Desktop/hglft_genome_3633f_8e5a50.bed', 
                          col_names = c('chrom', 'start', 'end', 'name')) %>%
  mutate(start = start + 1) %>% distinct()

mirtarbase <- mirna %>% 
  left_join(from_liftover, by = 'name') %>%
  # select(chrom, start, end) %>%
  select(id, name, chrom, start, end, gene_symbol, experiment, references) %>%
  na.omit() %>%
  mutate(chrom = str_remove(chrom, 'chr'))



# ------------------------------------------------------------------------------
# Dataset: COSMIC CNVs
# ------------------------------------------------------------------------------

non_coding_rna <- read_tsv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/RNA_micro.txt')


# ------------------------------------------------------------------------------
# Dataset: List TFs human
# Source: 
# ------------------------------------------------------------------------------
library(httr)

url1 <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fnature13182/MediaObjects/41586_2014_BFnature13182_MOESM86_ESM.xlsx'

GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))

tf_genes <- read_excel(tf, sheet = 8) %>% select(SYMBOL) %>% distinct() %>% pull()



# ------------------------------------------------------------------------------
# Dataset: DrugBank
# Source:https://www.drugbank.ca/releases/latest#protein-identifiers - Approved - Pharma. active
# Version: 2020-01-03	
# ------------------------------------------------------------------------------

drugbank <- read_csv('pharmacologically_active.csv')

drugbank <- drugbank %>% 
  filter(Species == 'Humans') %>%
  select(`Gene Name`, Name, `UniProt ID`, `Drug IDs`) %>%
  rename(gene = `Gene Name`, 
         uniprot_id = `UniProt ID`,
         # pdb_id = `PDB ID`,
         drug_id = `Drug IDs`) %>%
  separate_rows(drug_id, sep = ';') %>%
  mutate(drug_id = str_replace(drug_id, ' ', ''))

download.file('https://www.drugbank.ca/releases/5-1-5/downloads/all-drugbank-vocabulary', 'drugbank_all_drugbank_vocabulary.csv.zip')
file.remove('drugbank_all_drugbank_vocabulary.csv.zip')
names_drugs <- read_csv('drugbank_all_drugbank_vocabulary.csv.zip')
names_drugs <- names_drugs %>%
  select(`DrugBank ID`, `Common name`, Synonyms) %>%
  rename(drug_id = `DrugBank ID`,
         name = `Common name`,
         syn = Synonyms)


drugbank <- drugbank %>% left_join(names_drugs, by = 'drug_id')

# ------------------------------------------------------------------------------
# geno2MP database
# Version  January 06, 2020
# hg19
# ------------------------------------------------------------------------------

geno2mp <- read_tsv("http://geno2mp.gs.washington.edu/download/Geno2MP.variants.vcf.gz", 
                    skip = 14,  col_types = list(`#CHROM` = col_character()))


# ------------------------------------------------------------------------------
# Protein complex genes
# CORUM - Mammals Protein complexes 
# ------------------------------------------------------------------------------


download.file('http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip', 'allComplexes.txt.zip')
prot_complex <- read_tsv('allComplexes.txt.zip')

prot_complex <- prot_complex %>% 
  filter(Organism == 'Human') %>%
  select(ComplexID, ComplexName, `subunits(Gene name)`) %>%
  rename(id = ComplexID, name = ComplexName, gene =  `subunits(Gene name)`) %>%
  separate_rows(gene, sep = ';')
  

file.remove('allComplexes.txt.zip')


# ------------------------------------------------------------------------------
# Ohnologs genes
# OHNOLOGS database
# Filter: Strict
# ------------------------------------------------------------------------------


ohno <- read_tsv('http://ohnologs.curie.fr/cgi-bin/DownloadBrowse.cgi?crit=[0]&org=hsapiens&opt=pairs&wgd=2R')


ohno_genes <- ohno %>% select(Symbol1, Symbol2) %>% 
  pivot_longer(values_to = 'gene', names_to = 'delete', cols = c(Symbol1, Symbol2)) %>% 
  select(gene) %>%
  distinct()

# ------------------------------------------------------------------------------
# ncRNAs associated with diseases
# OHNOLOGS database
# Filter: Strict
# ------------------------------------------------------------------------------
# 
# url <- 'https://www.rna-society.org/mndr/php/download.php?file=All%20ncRNA-disease.zip'
# download.file(url, 'all_ncrna.zip')
# system('unzip all_ncrna.zip')
# 
# ncrna_disease <- read_tsv('All ncRNA-disease.txt', col_names = TRUE)
# 
# ncrna_disease %>%
#   filter(Species == 'Homo sapiens') %>%
#   filter(`ncRNA Category` == 'lncRNA') %>% 
#   filter(`Related Gene` != '') %>%
#   count(Methods) %>% arrange(desc(n))
# 
# ncrna_disease %>% glimpse()
# 
# file.remove('All ncRNA-disease.txt')
# file.remove('all_ncrna.zip')


# ------------------------------------------------------------------------------
# lncRNA2target database
# Only 153 unique lncRNAs 
# ------------------------------------------------------------------------------
library(readxl)

url <- 'http://123.59.132.21/lncrna2target/data/lncRNA_target_from_low_throughput_experiments.xlsx'

download.file(url, 'lncRNA_target_from_low_throughput_experiments.xlsx')
file.remove('lncRNA_target_from_low_throughput_experiments.xlsx')
lncrna_raw <- read_xlsx('lncRNA_target_from_low_throughput_experiments.xlsx', sheet = 1)

lncrna_target <- lncrna_raw %>% 
  filter(Species == 'Homo sapiens') %>%
  # filter(Ensembl_ID != 'NA') %>%
  select(LncRNA_official_symbol,Ensembl_ID, GENCODE_gene_name,Entrez_ID,
         LncRNA_experiment, Target_official_symbol,Tissue_Origin, Disease_state, PMID)

human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

interval_genes <- getBM(attributes = c('ensembl_gene_id', 'start_position','end_position', 'chromosome_name', 
                                       'entrezgene_id','hgnc_symbol', 'version'),
                        mart = human ) %>% 
  as_tibble() %>% 
  filter(!str_detect(chromosome_name, 'PATCH')) %>% 
  na.omit()

coord_lncrna <- lncrna_target %>% 
  select(Ensembl_ID) %>% 
  distinct() %>%
  left_join(interval_genes, by = c('Ensembl_ID' = 'ensembl_gene_id')) %>% 
  select(-entrezgene_id, -version, -hgnc_symbol) %>%
  na.omit() %>%
  distinct()


lncrna_target <- lncrna_target %>%
  select(Ensembl_ID, LncRNA_official_symbol, Target_official_symbol, Tissue_Origin, Disease_state, PMID) %>%
  left_join(coord_lncrna, by = 'Ensembl_ID') %>%
  na.omit() %>%
  rename(start = start_position, end = end_position, chrom = chromosome_name,
         official_symbol = LncRNA_official_symbol, target_symbol = Target_official_symbol) %>%
  select(official_symbol, chrom, start, end, everything())


# ------------------------------------------------------------------------------
# AGGREGATE ALL THE INFORMATION
# ------------------------------------------------------------------------------


hgcn_genes <- hgcn_genes %>%
  left_join(pli) %>% # pli score
  left_join(vg) %>% # variance gene expression
  mutate(haplo = as.factor(if_else(gene %in% clingen, 'Yes', 'No'))) %>% # haploinsufficiency genes
  mutate(triplo = as.factor(if_else(gene %in% triplo, 'Yes', 'No'))) %>% # triploinsufficiency genes 
  mutate(dev = as.factor(if_else(gene %in% dev_genes, 'Yes', 'No'))) %>% # developmental disorder genes - it can be extended with mode, consecuence and disease
  mutate(fda = as.factor(if_else(gene %in% fda, 'Yes', 'No'))) %>% #  Mechanistic targets of FDA-approved drugs 
  mutate(clinvar = as.factor(if_else(gene %in% clinvar_genes, 'Yes', 'No'))) %>% # List of genes with likely pathogenic and pathogenic variants
  mutate(gwas = as.factor(if_else(gene %in% gwas, 'Yes', 'No'))) %>% # GWAS genes
  mutate(omim = as.factor(if_else(gene %in% omim_genes, 'Yes', 'No'))) %>% # OMIM genes
  mutate(orphanet = as.factor(if_else(gene %in% orphanet_genes, 'Yes', 'No'))) %>% # OMIM genes
  mutate(genomics_england = as.factor(if_else(gene %in% panel_total_genes, 'Yes', 'No'))) %>% # Genomics England panel 
  mutate(essent = as.factor(if_else(ensembl_gene_id %in% e_intersect$V1 , 'Yes', 'No'))) %>% # essential genes (intersection mgi_invitro)
  # mutate(disease = as.factor(if_else(gene %in% genes_disgenet , 'Yes', 'No'))) %>%
  left_join(ccr) %>% # Genes with CCRs in the 99th percentile or higher 
  left_join(nc) %>% # non-coding scores RVIS and ncGERP - 5UTR + 3UTR + 250bp upstream
  left_join(rvis) %>% # RVIS score based
  left_join(hi) %>% # Haploinsufficiency Score (HI index)
  left_join(snipre) %>% # SNIPre score
  left_join(gdi) # GDI score (High GDI values reflect highly damaged genes. )

hgcn_genes <- hgcn_genes %>% 
  rename(chrom = chromosome_name) %>%
  mutate(clingen = if_else(haplo == 'Yes' | triplo == 'Yes', 'Yes', 'No')) %>%
  mutate(ccr = ifelse(is.na(ccr), 0, ccr)) %>%
  rename(band = location) %>%
  mutate(disease = if_else(clingen == 'Yes' | genomics_england == 'Yes' | dev == 'Yes' | orphanet == 'Yes' | omim == 'Yes', 'Yes', 'No' ) )


hgcn_genes <- hgcn_genes %>% 
  mutate(pLI = ntile(pLI, 100)) %>% # high pLI = 1 Likely Pathogenic
  mutate(rvis = ntile(-(rvis), 100)) %>% # low rvis = Likely Pathogenic
  mutate(ncrvis = ntile(-(ncrvis), 100)) %>% # low ncrvis = Likely Pathogenic
  mutate(ncgerp = ntile(ncgerp, 100)) %>% # high ncgerp = Likely Pathogenic
  mutate(gdi = ntile(-(gdi), 100)) %>% # low gdi = Likely Pathogenic
  mutate(hi = ntile(-(hi), 100)) %>% # low hi = Likely Pathogenic
  mutate(snipre = ntile(-(snipre), 100)) # low snipre = Likely Pathogenic

hgcn_genes <- hgcn_genes %>% rename(start = start_position, end = end_position)

# ------------------------------------------------------------------------------
# Dataset: TRRUST v.2
# Version 2
# ------------------------------------------------------------------------------

trrust <- read_tsv('https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv', 
                   col_names = c('tf', 'target', 'mechanism', 'reference'))

trrust <- trrust %>% separate_rows(reference, sep = ';')

trrust <- trrust %>% 
  left_join(hgcn_genes %>% select(gene, chrom, start_position, end_position) %>% 
              rename(start = start_position, end = end_position)
  , by = c('tf' = 'gene')) %>%
  left_join(hgcn_genes %>% select(gene, chrom, start_position, end_position) %>% 
              rename(target_chrom = chrom, 
                     target_start = start_position, 
                     target_end = end_position)
            , by = c('target' = 'gene')) %>%
  na.omit() %>%
  select(tf, chrom, start, end, target, mechanism, everything())

# ------------------------------------------------------------------------------
# Length chromosomes (hg19)
# ------------------------------------------------------------------------------  

coord_chrom_hg19 <- read_tsv('https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',
                             col_names = c('chrom', 'length')) %>%
  filter(nchar(chrom) < 6) %>% filter(!str_detect(chrom, 'chrM')) %>%
  mutate(chrom = str_remove(chrom, 'chr'))

result_tbl <- tibble()

for (i in 1:nrow(coord_chrom_hg19)) {
  
  last_nt <- coord_chrom_hg19 %>% slice(i) %>% pull(length)
  
  tmp_tbl <- tibble('chrom' = coord_chrom_hg19 %>% slice(i) %>% pull(chrom), 
         'start' = seq(1, (last_nt-10**6), 10**6 ), 
         'end' = seq(1000000, last_nt, 10**6,  ))
  
  tmp2_tbl <- tibble('chrom' = coord_chrom_hg19 %>% slice(i) %>% pull(chrom),
         'start' = tmp_tbl %>% tail(1) %>% pull(end) + 1,
         'end' = last_nt)
  
  
  result_tbl <- result_tbl %>% bind_rows(tmp_tbl, tmp2_tbl)
  
  
}

count_genes <- function(chrom, start, end) {

  result_tmp <- hgcn_genes %>% 
    rename(start = start_position, end = end_position) %>%
    bed_intersect(tibble('chrom' = chrom, 'start' = start, 'end' = end)) %>%
    nrow()
  
  return(result_tmp)
}


gene_density_tbl <- result_tbl %>% 
  mutate(is_end = if_else((end - start + 1) < 10**6, 'yes', 'no')) %>%
  mutate(gene_density = pmap_dbl(list(chrom, start, end), count_genes))

