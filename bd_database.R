

# Load packages
library(tidyverse)
library(readxl)
library(biomaRt)

# Reference of scores

tibble(gene_symbol = NA, gene_entrezid = NA,
  score = c('gwas', # we could filter out by number of hits (no filter - 8652, >1 4498)
                 'fda',
                 'pli',
                 'deg', # CRISPR - mice K.O - aggregated dataset (DEG) (last source from 2013)
                 'omim', # omim total - autosomal dominant - recessive - X-linked
                 'dev_disorder',
                 'haplo',
                 'triplo',
                 'vg',
            'rvis',
            'clinvar',
            'ccr'), # it could be extended with reference variant - omim - disease
       level = c('gene',
                 'gene',
                 'gene',
                 'gene',
                 'gene',
                 'gene',
                 'gene',
                 'gene',
                 'gene',
                 'gene'),
       source = c('http:',
                  'http:',
                  'http:',
                  'http:',
                  'http:',
                  'http:',
                  'http:',
                  'http:',
                  'http:',
                  'http:'))


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

hgcn_genes %>%
  left_join(pli) %>% # pli score
  left_join(vg) %>% # variance gene expression
  mutate(haplo = if_else(gene %in% clingen, 1, 0)) %>% # haploinsufficiency genes
  mutate(triplo = if_else(gene %in% triplo, 1, 0)) %>% # triploinsufficiency genes 
  mutate(dev = if_else(gene %in% dev_genes, 1, 0)) %>% # developmental disorder genes - it can be extended with mode, consecuence and disease
  mutate(fda = if_else(gene %in% fda, 1, 0)) %>% #  Mechanistic targets of FDA-approved drugs 
  left_join(rvis) %>% # RVIS score based
  mutate(clinvar = if_else(gene %in% clinvar_raw, 1, 0)) %>% # List of genes with likely pathogenic and pathogenic variants
  mutate(gwas = if_else(gene %in% gwas, 1, 0)) %>% # GWAS genes
  left_join(ccr) %>%
  count(ccr)# Genes with CCRs in the 99th percentile or higher

# ------------------------------------------------------------------------------
# Dataset: pLI
# Source: https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
# Access:  06/06/19
# ------------------------------------------------------------------------------

pli <- read.table('data/gnomad.v2.1.1.lof_metrics.by_gene.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  # filter(canonical == 'true') %>% 
  select(gene, pLI, transcript, oe_lof, oe_lof_lower, oe_lof_upper )

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
  na.omit()
  
# ------------------------------------------------------------------------------
# Dataset: Haploinsufficiency prediction
# Source: https://decipher.sanger.ac.uk/about#downloads/data
# Access: 6/6/19
# ASK ANTONIO IF IT'S WORTHY
# ------------------------------------------------------------------------------

a <- read.table('data/HI_Predictions_Version3.bed', sep = '\t', skip = 1) %>%
  as_tibble() %>%
  select(V4) %>%
  mutate(V4 = as.character(V4))



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
  na.omit()

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
