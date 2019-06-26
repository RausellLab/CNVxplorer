

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

hgcn_genes <- read_excel('data/gene_with_protein_product.xlsx') %>% as_tibble()

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
# ------------------------------------------------------------------------------


gnomad_sv_raw <- read.table('/home/cbl02/Storage/data/gnomad_v2_sv.sites.bed', sep = '\t', header = TRUE)


gnomad_sv_raw %>%
  as_tibble()


# ------------------------------------------------------------------------------
# Dataset: Structural Variants (SV)
# Source: https://decipher.sanger.ac.uk/about#downloads/data
# ------------------------------------------------------------------------------

decipher_sv_raw <- read.table('/home/cbl02/Storage/data/population_cnv.txt', sep = '\t', header = TRUE)

decipher_sv_raw %>%
  as_tibble()



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

# ERROR!!
hpo_raw <- read.table('http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab', 
                      sep = '\t', header = FALSE, fill = TRUE)

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
