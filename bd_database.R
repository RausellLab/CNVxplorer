

# Load packages
library(readxl)
library(biomaRt)
library(data.table)
library(tidyverse)
library(clusterProfiler)
library(data.table)
select <- dplyr::select


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
                               'rvis',
                               'clinvar',  # it could be extended with reference variant - omim - disease
                               'ccr',
                               'ncRVIS',
                               'ncGERP',
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
# Dataset: Protein-coding genes with HGCN symbol 
# Source:  ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt
# ------------------------------------------------------------------------------

hgcn_genes <- read_excel('/home/cbl02/Storage/data/gene_with_protein_product.xlsx') %>% as_tibble()

hgcn_genes <- hgcn_genes %>%
  select(entrez_id, ensembl_gene_id, location, symbol) %>%
  rename(gene = symbol) %>%
  na.omit() %>%
  mutate(entrez_id = as.numeric(entrez_id))


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

# 
# 'external_gene_name' - SYMBOL
# interval_genes <- getBM(attributes = c('ensembl_gene_id_version', 'start_position','end_position', 'chromosome_name'),
#                         mart = human ) %>% as_tibble() %>% filter(!str_detect(chromosome_name, 'PATCH'))
interval_genes <- getBM(attributes = c('entrezgene_id', 'start_position','end_position', 'chromosome_name'),
                        mart = human ) %>% 
                        as_tibble() %>% 
                        filter(!str_detect(chromosome_name, 'PATCH')) %>% 
                        na.omit() %>%
                        rename(entrez_id = entrezgene_id) %>%
                        mutate(entrez_id = as.numeric(entrez_id))

hgcn_genes <- interval_genes %>% 
  # rename(ensembl_gene_id = ensembl_gene_id_version) %>%
  # mutate(entrez_id = str_remove(entrez_id, '\\..*')) %>%
  right_join(hgcn_genes, by = c('entrez_id'))


hgcn_genes %>% count(start_position) %>% arrange(desc(n))
hgcn_genes %>% filter(is.na(start_position)) %>% select(ensembl_gene_id) %>% pull()

# There are 398 genes with no coordinates



# ------------------------------------------------------------------------------
# Dataset: pLI
# Source: https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
# Access:  06/06/19
# ------------------------------------------------------------------------------

pli <- read.table('/home/cbl02/Storage/data/gnomad.v2.1.1.lof_metrics.by_gene.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  # filter(canonical == 'true') %>% 
  select(gene, pLI, transcript, oe_lof, oe_lof_lower, oe_lof_upper ) %>%
  mutate(pLI = round(pLI, 2))

# ------------------------------------------------------------------------------
# Dataset: Vg - Expected population variance in its dosage that is due to genetic differences among individuals
# Source: https://www.biorxiv.org/content/10.1101/632794v1.supplementary-material - Table 1
# # Access:  06/06/19
# ------------------------------------------------------------------------------


vg_raw <- read_excel('/home/cbl02/Storage/data/vg.xlsx', sheet = 3)

vg_raw <- vg_raw %>% rename(vg = Avg_VG) %>% filter(vg != 'NaN') %>% mutate(vg = as.numeric(vg)) %>%
  rename(gene = GeneID)

test <- clusterProfiler::bitr(vg_raw$gene, fromType = 'ENSEMBL', toType = "ENTREZID", OrgDb="org.Hs.eg.db")

vg <- vg_raw %>% left_join(test, by = c('gene' = 'ENSEMBL')) %>% select(-gene) %>% select(ENTREZID, vg) %>%
  mutate(ENTREZID = as.numeric(ENTREZID)) %>% rename(entrez_id = ENTREZID)

# ------------------------------------------------------------------------------
# Dataset: Clingen
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_haploinsufficiency_gene_GRCh37.bed
# Access: 06/06/19
# Info score: https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/help.shtml#review
# ------------------------------------------------------------------------------

clingen_raw <- read.table('/home/cbl02/Storage/data/ClinGen_haploinsufficiency_gene_GRCh37.bed', col.names = c('chrom', 'start', 'end', 'gene', 'score'), stringsAsFactors = FALSE,
                          skip = 1)

clingen <- clingen_raw %>% 
  as_tibble() %>% 
  filter(score == 3) %>% 
  select(gene) %>% 
  pull(gene)

# ------------------------------------------------------------------------------
# Dataset: Clingen
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_recurrent_CNV_V1.0-hg19.bed

# ------------------------------------------------------------------------------

cnv_syndrome_clingen <- read_tsv('ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_recurrent_CNV_V1.0-hg19.bed',
                                 skip = 1, col_names = FALSE)

cnv_syndrome_clingen <- cnv_syndrome_clingen %>%
  select(X1, X2, X3) %>%
  rename(chrom = X1,
         start = X2,
         end = X3) %>%
  mutate(start = 1 + start) # .bed format 0-based



# ------------------------------------------------------------------------------
# Dataset: Clingen
# Source: ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen/ClinGen_triplosensitivity_gene.bed
# Access: 06/06/19
# Info scores: https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/help.shtml#review
# ------------------------------------------------------------------------------

triplo <- read.table('/home/cbl02/Storage/data/ClinGen_triplosensitivity_gene.bed', col.names = c('chrom', 'start', 'end', 'gene', 'score'), stringsAsFactors = FALSE,
                     skip = 1, sep = '\t') %>%
  as_tibble() %>%
  filter(score >= 1) %>%
  filter(score != 40) %>% 
  filter(!str_detect(score, 'Not yet evaluated')) %>%
  select(gene) %>%
  pull(gene)

# ------------------------------------------------------------------------------
# Dataset: CNV syndromes
# Source: Decipher
# ------------------------------------------------------------------------------

cnv_syndrome_decipher <- read_tsv('/home/cbl02/Storage/data/decipher_cnv_syndromes_13_12_19', col_names = TRUE)


# ------------------------------------------------------------------------------
# Dataset: DEVELOPMENTAL DISORDER GENES
# Source: https://decipher.sanger.ac.uk/ddd#ddgenes
# note: Filtering by: category == 'confirmed'
# Access: 06/06/19
# ------------------------------------------------------------------------------

### IMPORTANT

# THERE ARE RELEVANT FEATURES SUCH AS MODE, CONSEQUENCE, DISEASE THAT MAY BE RELEVANT
# FOR THE APP

dev_genes <- read.table('/home/cbl02/Storage/data/DDG2P_6_6_2019.csv', 
                        header = TRUE, stringsAsFactors = FALSE,
                        sep = ',') %>%
  as_tibble() %>%
  filter(DDD.category == 'confirmed') %>%
  select(gene.symbol) %>%
  rename(gene = gene.symbol) %>%
  pull(gene)

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


omim <- read_excel('/home/cbl02/Storage/data/morbidmap.xlsx', skip = 4, col_names = c('phenotype', 
                                                                                      'gene', 'gene_id', 'location'))
omim <- omim[1:7745,] # description of the (id) 

omim$gene <-map_chr(omim$gene, function(x) str_split(x, ',')[[1]][1])
omim$id_pheno <-map_chr(omim$phenotype, function(x) str_extract(x, '\\d{6}'))
# omim$map_key <-map_chr(omim$phenotype, function(x) str_extract(x, '\\(([^)]+)\\)'))
omim$map_key <-map_chr(omim$phenotype, function(x) str_extract(x, '\\(([[0-9]]+)\\)'))

omim$map_key <-map_chr(omim$map_key, function(x) str_remove(x, '\\('))
omim$map_key <-map_chr(omim$map_key, function(x) str_remove(x, '\\)'))

omim$phenotype <-map_chr(omim$phenotype, function(x) str_replace(x, '\\d{6}', ''))
omim$phenotype <-map_chr(omim$phenotype, function(x) str_replace(x, '\\d{6}', ''))
omim$phenotype <-map_chr(omim$phenotype, function(x) gsub( '*\\(.*?\\) *', '', x)) 
omim$phenotype <-map_chr(omim$phenotype, function(x) str_remove( x, ',  ')) 


omim  <- omim %>% 
  filter(map_key %in% c(1:4)) %>%
  pull(gene)

# omim %>% count(map_key)
# A tibble: 4 x 2
# map_key       n
#   1          68
#   2        1083
#   3        6444
#   4         148

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

clinvar_raw <- read.table('/home/cbl02/Storage/data/clinvar_20190527.vcf', skip = 28)

clinvar <- clinvar_raw %>% 
  as_tibble() %>%
  filter(!str_detect(V8, 'CLNSIG=Pathogenic')) %>%
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

# gtex <- read.table('https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt',
#                    sep = '\t', header = TRUE)

gtex <- read.table('/home/cbl02/Storage/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct',
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

hpa <- read.table('/home/cbl02/Storage/data/normal_tissue.tsv', header = TRUE, sep = '\t')

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

mgi <- mgi %>%
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

# 1 - The disorder is placed on the map based on its association with
# a gene, but the underlying defect is not known.
# 2 - The disorder has been placed on the map by linkage or other
# statistical method; no mutation has been found.
# 3 - The molecular basis for the disorder is known; a mutation has been
# found in the gene.
# 4 - A contiguous gene deletion or duplication syndrome, multiple genes
# are deleted or duplicated causing the phenotype.
#

morbidmap <- read_xlsx('/home/cbl02/Storage/data/morbidmap.xlsx', skip = 3)

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
# Dataset: DGV
# Source: http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt
# ------------------------------------------------------------------------------

dgv_df <- read_tsv('/home/cbl02/Storage/data/GRCh37_hg19_variants_2016-05-15.txt') %>%
  as_tibble() %>%
  filter(varianttype == 'CNV') %>%
  filter(variantsubtype %in% c('deletion', 'duplication')) %>%
  rename(id = variantaccession, chrom = chr) %>%
  #select(id, chrom, start, end) %>%
  mutate(source = 'dgv')  


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

vector_term <- hpo_genes %>% select(term, hp) %>% distinct() %>% pull(term)
vector_hp <- hpo_genes %>% select(term, hp) %>% distinct() %>% pull(hp)


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
# Genome reference: hg38
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
  # select(-X5) %>%
  mutate(chrom = str_remove(chrom, 'chr'))
  


# ------------------------------------------------------------------------------
# Dataset: Essential genes (intersection mgi + invitro (paper: ))
# Source: Email from Bartha
# ------------------------------------------------------------------------------

# 3262 genes found it / 3326 total nº genes
e_mgi <- read.table('/home/cbl02/Storage/data/essential_gene_lists/mgi', header = FALSE,
                    stringsAsFactors = FALSE) %>% as_tibble()

# 1931 genes found it / 2010 total nº genes

e_invitro <- read.table('/home/cbl02/Storage/data/essential_gene_lists/invitro', header = FALSE,
                    stringsAsFactors = FALSE) %>% as_tibble()

e_intersect <- e_mgi %>% filter(V1 %in% e_invitro$V1)



# ------------------------------------------------------------------------------
# Dataset: GDI score
# Source: Paper: The human gene damage index as a gene-level approach to prioritizing exome variants
# Explanation: High GDI values reflect highly damaged genes. 
# ------------------------------------------------------------------------------

gdi <- read.table('/home/cbl02/Storage/data/GDI_full_10282015.txt', skip = 1, stringsAsFactors = FALSE) %>% as_tibble()
colnames(gdi) <- c('gene', 'gdi', 'gdi_phred', 'all_diseases', 'all_mendelian', 'mendelian_ad', 
                   'mendelian_ar', 'all_pid', 'pid_ad', 'pid_ar', 'all_cancer', 'cancer_dominant',
                   'cancer_recessive')

gdi <- gdi %>% select(gene, gdi) %>% mutate(gdi = round(gdi, 2))

# ------------------------------------------------------------------------------
# Dataset: SNIPre
# Source: antonio's email
# ------------------------------------------------------------------------------

snipre <- read.table('/home/cbl02/Storage/data/Sniprefallgenes.txt', header = TRUE, stringsAsFactors = FALSE) %>%
  as_tibble()

snipre <- snipre %>% select(Gene, SnIPRE.f) %>% 
  rename(gene = Gene, snipre = SnIPRE.f) %>% 
  mutate(snipre = str_replace(snipre, ',', '.')) %>%
  mutate(snipre = as.numeric(snipre)) %>%
  mutate(snipre = round(snipre, 2))


# ------------------------------------------------------------------------------
# Dataset: Paralogous genes
# Source: http://ogee.medgenius.info/downloads/
# ------------------------------------------------------------------------------

para_genes <- read_tsv('/home/cbl02/Storage/data/dup_genes/dupgenes_ogee.txt')

# ------------------------------------------------------------------------------
# Dataset: Phenotype terms associated with OMIM diseases
# Source: https://hpo.jax.org/app/download/annotation
# Explanation: phenotype_annotation_hpoteam.tab: contains annotations made explicitly 
# and manually by the HPO-team (mostly referring to OMIM entries)
# Annotation: https://hpo.jax.org/app/help/annotations
# ------------------------------------------------------------------------------

hpo_omim <- read_tsv('/home/cbl02/Storage/data/phenotype_annotation_hpoteam.tab', col_names = FALSE)

hpo_omim <- hpo_omim %>% 
  filter(X1 == 'OMIM') %>%
  # select(-X2) %>%
  select(-X1, -X2, -X10, -X12, -X13, -X14) %>% # remove x1 because all values are omim - change if we include orphanet/omim - the other file
  rename(desc = X3, 
         pace = X4, 
         hpo = X5, 
         term = X6, 
         evidence = X7, 
         onset = X8, 
         frequency = X9, 
         clinical_modifier = X11) %>%
  mutate(desc = str_remove(desc, '[0-9]{6}'),
         desc = str_remove(desc, '#'))


# ------------------------------------------------------------------------------
# Dataset: Gene panel
# Source: hhttps://panelapp.genomicsengland.co.uk
# header:  web - subhead: Downloading Gene Panels
# We got 306 (5 less than expected) because these panels were empty:
# "gene_panel_218" "gene_panel_520" "gene_panel_530" "gene_panel_561" "gene_panel_82" 
# ------------------------------------------------------------------------------

library(rvest)
website <- "https://panelapp.genomicsengland.co.uk/panels/"
page <- read_html(website)

c_ref <- page %>%
  html_nodes("a") %>%       # find all links
  html_attr("href")

df_ref <- tibble(ref = c_ref, id = NA) %>%
  filter(str_detect(ref, 'download')) %>%
  mutate(ref = str_remove(ref, '/panels/')) %>% 
  mutate(id = ref) %>%
  mutate(id = str_remove(id, '/download/01234/'))
  
setwd('/home/cbl02/Storage/data/gene_panel')

walk2(df_ref$ref, df_ref$id, function(a, b) 
  download.file(url = paste0(website, a), destfile = paste0('gene_panel_', b))
  )
files_panel <- list.files()
panel_total <- tibble()

for (i in 1:length(files_panel)) {
  print(i)
  df_tmp <- read_tsv(paste0('/home/cbl02/Storage/data/gene_panel/', files_panel[i]))
  df_tmp <- df_tmp %>% mutate(source = files_panel[i])
  panel_total <- rbind(panel_total, df_tmp)
}

# we filtered out those genes not containing a "review ranking" ( 3,259 out 45488)
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
  select(gene, Level4, -sources, source) 
  

setwd('/home/cbl02/Storage/cnvxplore')




# ------------------------------------------------------------------------------
# Source: /home/cbl02/Storage/data/curated_gene_disease_associations.tsv
# http://www.disgenet.org/static/disgenet_ap1/files/downloads/readme.txt
# ------------------------------------------------------------------------------

disgenet <- read_tsv('/home/cbl02/Storage/data/curated_gene_disease_associations.tsv')

# ------------------------------------------------------------------------------
# Source: https://atlas.ctglab.nl/
# Version: Release 2 v20190117
# ------------------------------------------------------------------------------

gwas_atlas <- read_tsv('/home/cbl02/Storage/data/gwasATLAS_v20190117.txt')

# ------------------------------------------------------------------------------
# Dataset: Imprinting genes
# Source: http://www.geneimprint.com/site/genes-by-species
# ------------------------------------------------------------------------------

imp_genes <- read_tsv('/home/cbl02/Storage/data/imprinted_genes')

imp_genes <- imp_genes %>% 
  separate(Aliases, into = LETTERS[1:25], sep = ',') %>%
  gather('remove', 'gene', -Location, -Status, -`Expressed Allele`) %>%
  select(-remove) %>%
  rename(expressed_allele = `Expressed Allele` ) %>%
  select(gene, everything()) %>%
  distinct()

# ------------------------------------------------------------------------------
# Dataset: De novo variants (denovo-db)
# Source: http://denovo-db.gs.washington.edu/denovo-db/Download.jsp
# Version: 1.6.1
# non-SSC Samples	
# ------------------------------------------------------------------------------

denovo <- read_tsv('/home/cbl02/Storage/data/denovo-db.non-ssc-samples.variants.v.1.6.1.tsv', skip = 1)

denovo <- denovo %>% 
  select(Chr, Position, Gene, PrimaryPhenotype, StudyName, PubmedID, FunctionClass) %>%
  rename(chrom = Chr)


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
# AGGREGATE ALL THE INFORMATION
# ------------------------------------------------------------------------------

hgcn_genes <- hgcn_genes %>%
  left_join(pli) %>% # pli score
  left_join(vg) %>% # variance gene expression
  mutate(haplo = as.factor(if_else(gene %in% clingen, 'Yes', 'No'))) %>% # haploinsufficiency genes
  mutate(triplo = as.factor(if_else(gene %in% triplo, 'Yes', 'No'))) %>% # triploinsufficiency genes 
  mutate(dev = as.factor(if_else(gene %in% dev_genes, 'Yes', 'No'))) %>% # developmental disorder genes - it can be extended with mode, consecuence and disease
  mutate(fda = as.factor(if_else(gene %in% fda, 'Yes', 'No'))) %>% #  Mechanistic targets of FDA-approved drugs 
  mutate(clinvar = as.factor(if_else(gene %in% clinvar, 'Yes', 'No'))) %>% # List of genes with likely pathogenic and pathogenic variants
  mutate(gwas = as.factor(if_else(gene %in% gwas, 'Yes', 'No'))) %>% # GWAS genes
  mutate(omim = as.factor(if_else(gene %in% omim, 'Yes', 'No'))) %>% # GWAS genes
  mutate(essent = as.factor(if_else(ensembl_gene_id %in% e_intersect$V1 , 'Yes', 'No'))) %>% # essential genes (intersection mgi_invitro)
  left_join(ccr) %>% # Genes with CCRs in the 99th percentile or higher 
  left_join(nc) %>% # non-coding scores RVIS and ncGERP - 5UTR + 3UTR + 250bp upstream
  left_join(rvis) %>% # RVIS score based
  left_join(hi) %>% # Haploinsufficiency Score (HI index)
  left_join(snipre) %>% # SNIPre score
  left_join(gdi) # GDI score (High GDI values reflect highly damaged genes. )

hgcn_genes <- hgcn_genes %>% 
  rename(chrom = chromosome_name) %>%
  mutate(ccr = ifelse(is.na(ccr), 0, ccr)) %>%
  rename(band = location)