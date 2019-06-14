
library(chromPlot)

data(hg19_cytoBandIdeo)

hg_cytoBandIdeo <- chromPlot::hg_cytoBandIdeo %>% filter(Name == 'p14')

##########

library(rentrez)
test <- entrez_search(db="omim", term="DMD") %>% as_tibble()

###

test <- GET('https://ghr.nlm.nih.gov/search?query=DMD&show=xml')
xmlParse(test)
XML::xmlToDataFrame(getNodeSet(test, path='//row'))
###
library('heatmaply')

data <- data.frame(Value = round(rnorm(50, 10, 2), 0))
ggplot(data) + 
  geom_histogram(aes(x = Value, fill = Value == 13))

hgcn_genes %>%
  select(gene, pLI, rvis) %>%
  na.omit() %>%
  slice(1:50) %>%
  column_to_rownames('gene') %>%
  as.matrix() %>%
  heatmaply()


test <- hgcn_genes %>% slice(1:100) %>% select(entrez_id) %>% pull()  %>% as.character()
univ <- hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character

go_analysis <- enrichGO(gene  = test,
                        universe      = univ,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05) %>% as.data.frame()

go_analysis %>% 
  # filter(pvalue <= 0.05) %>%
  ggplot(aes(reorder(Description, p.adjust), p.adjust)) +
  geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
  scale_fill_viridis_d() +
  # scale_y_log10() +
  coord_flip() +
  xlab('') +
  ylab('p-adjusted') +
  geom_vline(xintercept = 0.05) +
  theme_minimal()

tablerCard(
  title = "Plots",
  zoomable = FALSE,
  closable = FALSE,
  plotOutput('func_analysis'),
  options = tagList(
    switchInput(
      inputId = "enable_distPlot",
      label = "Plot?",
      value = TRUE,
      onStatus = "success",
      offStatus = "danger"
    )
  )),


test <- hgcn_genes %>% slice(1:100) %>% select(entrez_id) %>% pull()  %>% as.character()
univ <- hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character

  

ego <- enrichGO(gene          = test,
                universe      = univ,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

barplot(ego, showCategory=30)


df <- mtcars %>% mutate(name = row.names(.))
df %>% 
  ggplot(aes(mpg, disp)) +
  geom_point(col = "darkred") +
  gghighlight(name == 'Camaro Z28',
              unhighlighted_colour = alpha("steelblue", 0.4),
              use_direct_label = TRUE,
              label_key = name,
              label_params = list(size = 5)) +
  geom_point(col = "darkred", size = 2.5) 

a <- matrix(NA, nrow = 2, ncol = 500)
a[1,] <- hgcn_genes$pLI[1:500]
a[2,] <- hgcn_genes$rvis[1:500]
a[1,][as.numeric(which(is.na(a[1,])))] <- 0
a[2,][as.numeric(which(is.na(a[2,])))] <- 0
colnames(a) <- hgcn_genes$gene[1:500]
pheatmap(a, cluster_rows = FALSE, cluster_cols = FALSE)
  
  df <- data.frame(x = c("a", "b"), y = c(3, 4), z = c(5, 6))
  df %>% spread(x, y) 
  

test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

library(pheatmap)


c(hgcn_genes$start_position[28], data_raw$end_position[28]) %overlaps% c(1000, 60000000)


de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)



data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edo, foldChange=geneList)
barplot(edo, showCategory=20)


for (i in 1:19000) {
  
 b <- c(hgcn_genes$start_position[i], hgcn_genes$end_position[i]) %overlaps% c(10, 45)
 print(b)
}

hgcn_genes %>% mutate(testeo =  c(start_position, end_position) %overlaps% c(10, 45)) %>% count(testeo)

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


plotKaryotype(chromosomes = 'chr1', plot.type = 2) %>%
  kpDataBackground(data.panel = 1) %>%
  kpAddBaseNumbers() %>%
  kpRect(chr="chr1", x0=c(100000000, 110300000), x1=c(120000000, 150050000), y0=0.2, y1=0.3)
