
library(chromPlot)

data(hg19_cytoBandIdeo)

chromPlot::hg_cytoBandIdeo

##########

hgcn_genes <- read_excel('data/gene_with_protein_product.xlsx') %>% as_tibble()

hgcn_genes <- hgcn_genes %>%
  select(entrez_id, ensembl_gene_id, location, symbol) %>%
  rename(gene = symbol)


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

test <- useEnsembl(biomart="ensembl",
                   GRCh=37,
                   dataset = "hsapiens_gene_ensembl",
                   version=96)


interval_genes <- getBM(attributes = c('ensembl_gene_id_version', 'start_position','end_position', 'chromosome_name'), 
                        mart = human ) %>% as_tibble() %>% filter(!str_detect(chromosome_name, 'PATCH'))

hgcn_genes <- interval_genes %>% 
  rename(ensembl_gene_id = ensembl_gene_id_version) %>%
  mutate(ensembl_gene_id = str_remove(ensembl_gene_id, '\\..*')) %>%
  right_join(hgcn_genes, by = c('ensembl_gene_id'))

hgcn_genes %>% count(start_position) %>% arrange(desc(n))

interval_genes <- getBM(attributes = c('external_gene_name', 'start_position','end_position', 'chromosome_name'), 
                        mart = human ) %>% as_tibble() %>% filter(!str_detect(chromosome_name, 'PATCH'))

hgcn_genes <- interval_genes %>% 
  rename(gene = external_gene_name) %>%
  right_join(hgcn_genes, by = c('gene'))

hgcn_genes %>% count(start_position) %>% arrange(desc(n))




interval_genes <- getBM(attributes = c('entrezgene', 'start_position','end_position', 'chromosome_name'), 
                        mart = human ) %>% as_tibble() %>% filter(!str_detect(chromosome_name, 'PATCH')) %>% na.omit()

hgcn_genes <- interval_genes %>% 
  rename(entrez_id = entrezgene) %>%
  right_join(hgcn_genes)

hgcn_genes %>% count(start_position) %>% arrange(desc(n))
