

# Load packages
library(tidyverse)
library(readxl)
library(biomaRt)
library(data.table)

# TODO

# 398 genes does not have coordinates (start/end)
# OMIM
# HI%
# Percentile RVIS or original data
# Functional analysis - Gene Ontology - Pathway analysis (which db?)
# ReactomePA analysis 423MB careful!
# Gene ontology based of non-coding regions
# difference between pLI and o/e
# should a put a minus + plus to the output from remot
# PIK3CD duplicated because database HI have two observations
# warning message with input outside of the maximum chromosomes

# QUESTIONS NECKER

# Input user, two options: 1- select chrom, g banding or select g banding where the chromosome is included

# Reference of scores

ref_scores <- tibble(score = c('gwas', # we could filter out by number of hits (no filter - 8652, >1 4498)
                               'fda',
                               'pli',
                               'deg', # CRISPR - mice K.O - aggregated dataset (DEG) (last source from 2013)
                               'omim', # omim total - autosomal dominant - recessive - X-linked
                               'dev_disorder',
                               'haplo',
                               'triplo',
                               'vg',
                               'rvis',
                               'clinvar',  # it could be extended with reference variant - omim - disease
                               'ccr',
                               'ncRVIS',
                               'ncGERP',
                               'HI score'),
                     level = c('gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene',
                               'gene'),
                     description = c('-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-',
                                     '-'),
                     type = c('d',
                              'd',
                              'c',
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
                              'c'),
                     source = c('-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-',
                                '-'))


# Load datasets

# ------------------------------------------------------------------------------
# Dataset: Genes with HGCN symbol and 
# Source:  ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt
# ------------------------------------------------------------------------------

hgcn_genes <- read_excel('/home/cbl02/Storage/data/gene_with_protein_product.xlsx') %>% as_tibble()

hgcn_genes <- hgcn_genes %>%
  select(entrez_id, ensembl_gene_id, location, symbol) %>%
  rename(gene = symbol) %>%
  na.omit()


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

# 
# 'external_gene_name' - SYMBOL
interval_genes <- getBM(attributes = c('ensembl_gene_id_version', 'start_position','end_position', 'chromosome_name'), 
                        mart = human ) %>% as_tibble() %>% filter(!str_detect(chromosome_name, 'PATCH'))

hgcn_genes <- interval_genes %>% 
  rename(ensembl_gene_id = ensembl_gene_id_version) %>%
  mutate(ensembl_gene_id = str_remove(ensembl_gene_id, '\\..*')) %>%
  right_join(hgcn_genes, by = c('ensembl_gene_id'))


hgcn_genes %>% count(start_position) %>% arrange(desc(n))
hgcn_genes %>% filter(is.na(start_position)) %>% select(ensembl_gene_id) %>% pull()

# There are 398 genes with no coordinates

hgcn_genes <- hgcn_genes %>%
  left_join(pli) %>% # pli score
  left_join(vg) %>% # variance gene expression
  mutate(haplo = if_else(gene %in% clingen, 1, 0)) %>% # haploinsufficiency genes
  mutate(triplo = if_else(gene %in% triplo, 1, 0)) %>% # triploinsufficiency genes 
  mutate(dev = if_else(gene %in% dev_genes, 1, 0)) %>% # developmental disorder genes - it can be extended with mode, consecuence and disease
  mutate(fda = if_else(gene %in% fda, 1, 0)) %>% #  Mechanistic targets of FDA-approved drugs 
  left_join(rvis) %>% # RVIS score based
  mutate(clinvar = as.factor(if_else(gene %in% clinvar_raw, 1, 0))) %>% # List of genes with likely pathogenic and pathogenic variants
  mutate(gwas = if_else(gene %in% gwas, 1, 0)) %>% # GWAS genes
  left_join(ccr) %>% # Genes with CCRs in the 99th percentile or higher 
  left_join(nc) %>% # non-coding scores RVIS and ncGERP - 5UTR + 3UTR + 250bp upstream
  left_join(hi) # Haploinsufficiency Score (HI index)

hgcn_genes <- hgcn_genes %>% rename(chrom = chromosome_name)
# ------------------------------------------------------------------------------
# Dataset: pLI
# Source: https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
# Access:  06/06/19
# ------------------------------------------------------------------------------

pli <- read.table('data/gnomad.v2.1.1.lof_metrics.by_gene.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  # filter(canonical == 'true') %>% 
  select(gene, pLI, transcript, oe_lof, oe_lof_lower, oe_lof_upper ) %>%
  mutate(pLI = round(pLI, 2))

# ------------------------------------------------------------------------------
# Dataset: Vg - Expected population variance in its dosage that is due to genetic differences among individuals
# Source: https://www.biorxiv.org/content/10.1101/632794v1.supplementary-material - Table 1
# # Access:  06/06/19
# ------------------------------------------------------------------------------


vg_raw <- read_excel('data/vg.xlsx', sheet = 3)

vg_raw <- vg_raw %>% rename(vg = Avg_VG) %>% filter(vg != 'NaN') %>% mutate(vg = as.numeric(vg)) %>%
  rename(gene = GeneID)

test <- bitr(vg_raw$gene, fromType = 'ENSEMBL', toType = "ENTREZID", OrgDb="org.Hs.eg.db")

vg <- vg_raw %>% left_join(test, by = c('gene' = 'ENSEMBL')) %>% select(-gene) %>% select(ENTREZID, vg) %>%
  mutate(ENTREZID = as.numeric(ENTREZID)) %>% rename(entrez_id = ENTREZID)

# ------------------------------------------------------------------------------
# Dataset: Clingen
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_haploinsufficiency_gene_GRCh37.bed
# Access: 06/06/19
# ------------------------------------------------------------------------------

clingen_raw <- read.table('data/ClinGen_haploinsufficiency_gene_GRCh37.bed', col.names = c('chrom', 'start', 'end', 'gene', 'score'), stringsAsFactors = FALSE,
                          skip = 1)

clingen <- clingen_raw %>% as_tibble() %>% filter(score == 3) %>% select(gene) %>% pull(gene)

# ------------------------------------------------------------------------------
# Dataset: Clingen
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_triplosensitivity_gene.bed
# Access: 06/06/19
# ------------------------------------------------------------------------------

triplo <- read.table('data/ClinGen_triplosensitivity_gene.bed', col.names = c('chrom', 'start', 'end', 'gene', 'score'), stringsAsFactors = FALSE,
                     skip = 1, sep = '\t') %>%
  as_tibble() %>%
  filter(score == 1) %>%
  select(gene) %>%
  pull(gene)

# ------------------------------------------------------------------------------
# Dataset: DEVELOPMENTAL DISORDER GENES
# Source: https://decipher.sanger.ac.uk/ddd#ddgenes
# note: Filtering by: category == 'confirmed'
# Access: 06/06/19
# ------------------------------------------------------------------------------

dev_genes <- read.table('data/DDG2P_6_6_2019.csv', 
                        header = TRUE, stringsAsFactors = FALSE,
                        sep = ',') %>%
  as_tibble() %>%
  filter(DDD.category == 'confirmed') %>%
  select(gene.symbol) %>%
  rename(gene = gene.symbol) %>%
  pull(gene)

# ------------------------------------------------------------------------------
# Dataset: OMIM 
# Source: https://decipher.sanger.ac.uk/ddd#ddgenes
# Access: ####
# ------------------------------------------------------------------------------

## there are some genes with duplicated symbol on the same row

omim <- read.table('data/', sep = '\t', header = TRUE,
                   stringsAsFactors = F)
omim <- omim[-3087,]
colnames(omim) <- c('pheno', 'gene','id_gene', 'location')

omim$gene <-map_chr(omim$gene, function(x) str_split(x, ',')[[1]][1])
omim$id_pheno <-map_chr(omim$pheno, function(x) str_extract(x, '\\d{6}'))
omim$pheno <-map_chr(omim$pheno, function(x) str_replace(x, '\\d{6}', ''))
omim$pheno <-map_chr(omim$pheno, function(x) str_replace(x, '\\d{6}', ''))
omim$pheno <-map_chr(omim$pheno, function(x) gsub( '*\\(.*?\\) *', '', x)) 



# ------------------------------------------------------------------------------
# Dataset: FDA 
# Source: https://decipher.sanger.ac.uk/ddd#ddgenes
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
# Access: 
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

clinvar_raw <- read.table('data/clinvar_20190527.vcf', skip = 28)

clinvar_raw <- clinvar_raw %>% 
  as_tibble() %>%
  filter(!str_detect(V8, 'CLNSIG=Pathogenic') | str_detect(V8, 'CLNSIG=likely_pathogenic')) %>%
  select(V8) %>%
  mutate(gene = str_extract(V8, pattern =  'GENEINFO[^;]*')) %>%
  select(-V8) %>%
  mutate(gene = str_remove(gene, 'GENEINFO=')) %>%
  mutate(gene = str_remove(gene, '\\:.*')) %>%
  distinct() %>%
  na.omit() %>%
  pull()

# ------------------------------------------------------------------------------
# Dataset: Haploinsufficiency prediction
# Source: https://decipher.sanger.ac.uk/about#downloads/data
# Access: 6/6/19
# ASK ANTONIO IF IT'S WORTHY
# ------------------------------------------------------------------------------


double_genes  <- read.table('data/HI_Predictions_Version3.bed', sep = '\t', skip = 1) %>%
  as_tibble() %>%
  select(V4) %>%
  mutate(V4 = as.character(V4)) %>%
  filter(str_detect(V4, '-')) %>%
  separate(V4, into = as.character(1:10)) %>%
  select(`1`, `5`, `6`) %>%
  rename(A = `1`, D = `5`, E = `6`)






hi <- read.table('data/HI_Predictions_Version3.bed', sep = '\t', skip = 1) %>%
  as_tibble() %>%
  select(V4) %>%
  mutate(V4 = as.character(V4)) %>%
  filter(!str_detect(V4, '-')) %>%
  separate(V4, into = LETTERS[1:10]) %>% 
  select(A, D, E) %>%
  rbind(double_genes) %>%
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


gwas_raw <- read.table('data/gwas_catalog_v1.0-associations_e96_r2019-05-03.tsv', header = TRUE, sep = '\t',
                       fill = TRUE)

gwas <- gwas_raw %>% as.tibble() %>% 
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


nc_raw <- read_excel('data/journal.pgen.1005492.s011.XLSX.xlsx', sheet = 1)

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

# gtex <- read.table('https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt',
#                    sep = '\t', header = TRUE)

gtex <- read.table('data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct',
                   sep = '\t', header = TRUE, skip = 2)

gtex <- gtex %>% 
  as_tibble() %>%
  rename(gene = Description) %>%
  select(-gene_id) %>%
  gather('tissue', 'value', - gene)

# ------------------------------------------------------------------------------
# Dataset: HPA
# Source: https://gtexportal.org/home/datasets
# File name: 	GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
# ------------------------------------------------------------------------------

hpa <- read.table('data/normal_tissue.tsv', header = TRUE, sep = '\t')

hpa <- hpa %>%
  as_tibble() %>%
  rename(gene = Gene.name, tissue = Tissue, cell_type = Cell.type) %>%
  # mutate(Level = as.numeric(case_when(
  #   Level == 'High' ~ "3",
  #   Level == 'Medium' ~ "2",
  #   Level == 'Low' ~ "1",
  #   Level == 'Not detected' ~ "0"
  # ))) %>%
  select(-Gene)

# ------------------------------------------------------------------------------
# Dataset: Mouse phenotype
# Source: http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
# ------------------------------------------------------------------------------

mgi <- read.table('http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt', sep = '\t')

mgi <- mouse_p %>%
  as_tibble() %>%
  select(-V8, -V3) %>%
  rename(gene = V1, entrez_id = V2, gene_mouse = V5, pheno = V7, mgi = V6) %>%
  mutate(mgi = str_remove(mgi, pattern = '  ')) %>%
  select(-V4) %>%
  filter(pheno != '')

# ------------------------------------------------------------------------------
# Dataset: OMIM
# Source: https://data.omim.org/downloads/cpJFEdlrQ5qqPk0TOzvVBA/morbidmap.txt
# ------------------------------------------------------------------------------


morbidmap <- read_xlsx('data/morbidmap.xlsx', skip = 3)

morbidmap <- morbidmap %>% 
  as_tibble() %>%
  rename(pheno = `# Phenotype`, gene = `Gene Symbols`, mim_gene = `MIM Number`) %>%
  select(- `Cyto Location`) %>%
  mutate(mapping = 
           case_when(
             str_detect(pheno, '\\([1]\\)') ~ "1",
             str_detect(pheno, '\\([2]\\)') ~ "2",
             str_detect(pheno, '\\([3]\\)') ~ "3",
             str_detect(pheno, '\\([4]\\)') ~ "4"
           )) %>%
  filter(!isNA(mapping)) %>%  # eliminate description at the bottom
  mutate(pheno = str_remove(pheno, '\\([1-4]\\)' )) %>%
  mutate(mim_disease = str_extract(pheno, '[0-9]{6}')) %>%
  mutate(pheno = str_remove(pheno, '[0-9]{6}')) %>%
  mutate(pheno = str_remove(pheno, ',  ')) %>%
  select(gene, pheno, mim_gene, mim_disease, mapping)


# ------------------------------------------------------------------------------
# Dataset: Paralogous genes
# Source: biomart
# ------------------------------------------------------------------------------


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")


para_genes <- getBM(attributes = c('external_gene_name', 'hsapiens_paralog_ensembl_gene'), 
                        mart = human )

para_genes %>% as_tibble() %>%
  rename(gene = external_gene_name, para = hsapiens_paralog_ensembl_gene) %>%
  filter(para != '') %>%
  count(gene) %>%
  filter(gene == 'DMD')


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
# PCRPLUS_DEPLETED I don't know what it is....
# PESR_GT_OVERDISPERSION this neither...


gnomad_sv_raw <- read.table('/home/cbl02/Storage/data/gnomad_v2_sv.sites.bed', sep = '\t', header = TRUE)


gnomad_sv_raw <- gnomad_sv_raw %>%
  as_tibble() %>%
  filter(SVTYPE %in% c('DEL', 'DUP')) %>%
  filter(FILTER == 'PASS') %>%
  select(CHROM, START, END, NAME) %>%
  rename(id = NAME) %>%
  mutate(id  = str_remove(id, 'gnomAD_v2_')) %>%
  mutate(source = 'gnomad_v2') %>%
  rename(chrom = CHROM, start = START, end = END)
    


# ------------------------------------------------------------------------------
# Dataset: Structural Variants (SV)
# Source: https://decipher.sanger.ac.uk/about#downloads/data
# Population Copy-Number Variation Frequencies
# Variables: duplicated - deletion - general // observations - frequency - se
# ------------------------------------------------------------------------------

decipher_sv_raw <- read.table('/home/cbl02/Storage/data/population_cnv.txt', sep = '\t', header = TRUE) %>%
  as_tibble() %>%
  mutate(source = 'decipher') %>%
  rename(id = population_cnv_id, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  select(id, chrom, start, end, source)
  
  
# ------------------------------------------------------------------------------
# Dataset: Aggregation of data from DECIPHER, gnomAD and DGV
# ------------------------------------------------------------------------------

  cnv_df <- decipher_sv_raw %>% rbind(gnomad_sv_raw) %>% rbind(dgv_df)

# ------------------------------------------------------------------------------
# Dataset: Protein-Protein interaction network
# Source: https://www.intomics.com/inbio/map.html#downloads
# Name file: inBio_Map_core_2016_09_12.tar.gz
# ------------------------------------------------------------------------------

inbio_network_raw <- read.table('/home/cbl02/Storage/data/core.psimitab', sep = '\t', header = FALSE)

inbio_network <- inbio_network_raw %>%
  as_tibble() %>%
  slice(1:100)


# ------------------------------------------------------------------------------
# Dataset: Human Phenotype Ontology
# Source: http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
# Name file: phenotype_annotation.tab
# ------------------------------------------------------------------------------

# ERROR - LESS ROWS THAN  FILE!!

url <- 'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt'
hpo_raw <- read.table(url, sep = '\t', skip = 1, stringsAsFactors = FALSE)

hpo_genes <- hpo_raw %>% as_tibble() %>% rename(entrez_id = V1, gene = V2, term = V3, hp = V4)


# ------------------------------------------------------------------------------
# Dataset: TADs
# Source: 30765865 - 25693564
# Name file: https://raw.githubusercontent.com/JacobSpectorMD/ClinTAD/master/home/files/boundary.txt
# ------------------------------------------------------------------------------

tad <- read.table('https://raw.githubusercontent.com/JacobSpectorMD/ClinTAD/master/home/files/boundary.txt', header = FALSE, sep = '\t', 
                  stringsAsFactors = FALSE)

tad <- tad %>%
  as_tibble() %>%
  rename(chrom = V1, start = V2, end = V3) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  mutate(id = row_number()) %>% select(id, chrom, start, end)


# ------------------------------------------------------------------------------
# Dataset: Genehancer
# Version: 4.11
# Name file: Email from Marilyn (genenhancer team)
# Genome reference: hg19
# ------------------------------------------------------------------------------

genehancer <- read.table('/home/cbl02/Storage/data/genehancer_V4_11.gff', header = FALSE, sep = '\t', 
                  stringsAsFactors = FALSE)

df <- genehancer %>% 
  as.tibble() %>%
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
df_assoc_gene <- df_assoc_gene %>% filter(score_gene > 1)

df_ge <- df_assoc_enh %>% left_join(df_assoc_gene, by = 'id')
df_ge <- df_ge %>% select(gene, everything())
df_ge <- df_ge %>% na.omit()
df_ge <- df_ge %>% filter(chrom != 'chrY')

write.table(df_ge %>% select(chrom, start, end, id) %>% distinct(), 'enhancer_cnvxplore', quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = '\t', append = TRUE)

# Used liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver). Default parameters
# from Grch38 to Grch37
# Succesfully converted 27514 records
# Conversion failed on 8 records.

enhancer_raw <- read.table('/home/cbl02/Storage/data/hglft_genome_4c90f_a0b3d0.bed', header = FALSE, sep = '\t')

df_enhancers <- enhancer_raw %>% as_tibble() %>% rename(chrom = V1, start = V2, end = V3, id = V4) %>% mutate(chrom = str_remove(chrom, 'chr')) %>%
  mutate(id = as.character(id))

df_enhancers <- df_ge %>% select(gene, id, score_enh, score_gene) %>% left_join(df_enhancers, by = 'id')
# 
# write.table(enhancer_raw %>% filter(!str_detect(V1, 'PATCH')) %>% distinct(),
#             'enhancer_cnvxplore_crossmap_cleaned', quote = FALSE, row.names = FALSE,
#             col.names = FALSE, sep = '\t', append = TRUE)
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


lncrna_expression_raw <- read.table('/home/cbl02/Storage/data/lncRNA/HumanRPKMs.txt', header = TRUE, sep = '\t',
                     stringsAsFactors = FALSE)

## 7 organs: Brain, Cerebellum, Heart, Kidney, Liver, Ovary, Testis
## 26 developmental step: "10wpc, 11wpc, 12wpc, 13wpc, 16wpc, 18wpc, 19wpc, 20wpc, 4wpc, 5wpc, 6wpc, 7wpc, 
# 8wpc, 9wpc, infant, newborn, olderMidAge, oldTeenager, school, senior, Senior, teenager, toddler, youngAdult, 
# youngMidAge, youngTeenager"

lncrna_expression <- lncrna_expression_raw %>%
  rownames_to_column(var = 'id') %>%
  as_tibble() %>%
  gather('type', 'rpkm', -id) %>%
  separate(col = type, into = c('organ', 'dev_step', 'dup'), sep = '\\.') %>%
  mutate(dev_step = if_else(dev_step == 'Senior', 'senior', dev_step)) %>% # there was a typo in 85037 rows (Senior) and the rest, senior (510222)
  filter(!str_detect(string = id, pattern = 'ENSG')) %>%
  filter(str_detect(id, 'XLOC'))



# ------------------------------------------------------------------------------
# Dataset: seg_dup and self_chain
# Description: Bed files with segmental duplications and self-chain of the genome(hg19)
# Source: http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
# Source: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chainSelf.txt.gz.
# Source with name of each field: http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgsid=1432484_CeWhnfzDUkBhJCkZ00iaiTQcyb3E&hgta_doSchemaDb=hg19&hgta_doSchemaTable=chainSelf
# Comment: the score 90 is obtained from Quinlab (A map of constraint regions in coding part)
# ------------------------------------------------------------------------------

seg_dup <- read.table('http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab', header = TRUE) %>% 
  as_tibble()

seg_dup1 <- seg_dup %>% select(chrom, chromStart, chromEnd) %>% 
  rename(start = chromStart, end = chromEnd)
seg_dup2 <- seg_dup %>% select(otherChrom, otherStart, otherEnd) %>%
  rename(chrom = otherChrom, start = otherStart, end = otherEnd)
seg_dup <- seg_dup1 %>% bind_rows(seg_dup2) %>% distinct()  %>% mutate(chrom = str_remove(chrom, 'chr')) %>%
  makeGRangesFromDataFrame()

# Self-chain coordinates

self_chain <- read.table('data/chainSelf.txt', header = FALSE) %>% as_tibble()
colnames(self_chain) <- c('bin', 'score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStrand',
                          'qStart', 'qEnd', 'id', 'normScore')

# NormScore = score / length bases. This score threshold is based on the paper: 
self_chain <- self_chain %>% filter(normScore >= 90)

# Clean dataframe and convert it to a GRanges object

self_chain1 <- self_chain %>% select(tName, tStart, tEnd) %>% 
  rename(chrom = tName, start = tStart, end = tEnd)
self_chain2 <- self_chain %>% select(qName, qStart, qEnd) %>%
  rename(chrom = qName, start = qStart, end = qEnd)
self_chain <- self_chain1 %>% bind_rows(self_chain2) %>% distinct()  %>% mutate(chrom = str_remove(chrom, 'chr')) %>%
  makeGRangesFromDataFrame()


# ------------------------------------------------------------------------------
# Dataset: blacklist.v2.bed
# Source: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz
# ------------------------------------------------------------------------------

blacklist_encode <- read.table('/home/cbl02/Storage/data/hg19-blacklist.bed', sep = '\t', stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(chrom = V1, start = V2, end = V3, class = V4) %>%
  select(-V5, -V6) %>%
  mutate(chrom = str_remove(chrom, 'chr'))
  
# ------------------------------------------------------------------------------
# Dataset: DGV
# Source: http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt
# ------------------------------------------------------------------------------

dgv_df <- read.table('/home/cbl02/Storage/data/GRCh37_hg19_variants_2016-05-15.txt', sep = '\t', stringsAsFactors = FALSE,
                     header = TRUE) %>%
  as_tibble() %>%
  filter(varianttype == 'CNV') %>%
  filter(variantsubtype %in% c('deletion', 'duplication')) %>%
  rename(id = variantaccession, chrom = chr) %>%
  select(id, chrom, start, end) %>%
  mutate(source = 'dgv')


# ------------------------------------------------------------------------------
# Dataset: Imprinting genes
# Source: http://www.geneimprint.com/site/genes-by-species
# ------------------------------------------------------------------------------



  