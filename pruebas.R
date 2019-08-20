

library(rjson)

json_file <- 'http://myvariant.info/v1/query?q=chr1:69000-70000'

json_data <- fromJSON(file=json_file)





snv_df <- read.table('/home/cbl02/Desktop/input_snv.bed', sep = '\t', header = FALSE) %>%
  rename(chrom = V1, start = V2, end = V3)


for (i in 1:nrow(snv_df)) {
  
  df_add <- hgcn_genes %>% 
    filter(chrom == snv_df$chrom[i]) %>%
    mutate(keep = NA) %>%
    rowwise() %>%
    mutate(keep = c(start_position, end_position) %overlaps% c(snv_df$start[i], snv_df$end[i])) %>%
    filter(keep == TRUE) %>% 
    select(-keep) %>%
    ungroup()
  # check with snvs that overlap in more than 1 gene
  if (!df_add$gene %in% data_raw$gene) {
    data_raw <-  bind_rows(hgcn_genes, df_add)
  }
}







test99 %>%
  count(source) %>%
  mutate(total = sum(n)) %>%
  rowwise() %>%
  mutate((percentage = n / total) * 100)
  

get_perc_overlap(test99 %>% rename(start_position = start, end_position = end), 1000, 10000000) %>% count(p_overlap)


cnv_df %>%
  ggplot(aes(length_cnv)) +
  geom_density(aes(fill = source), alpha = 0.5, color = 'black') +
  coord_cartesian(xlim = c(0,250000)) +
  geom_vline(aes(xintercept = 50000), linetype = 2) +
  theme_minimal()

cnv_df %>%
  ggplot(aes(length_cnv, y = source)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = source), alpha = 0.6, show.legend = FALSE) +
  geom_vline(aes(xintercept = 50000), linetype = 2, color = 'red', size = 1.5) +
  scale_x_log10() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis_d() +
  xlab('CNVs size') +
  ylab('Database') +
  theme_ridges()

cnv_df %>%
  ggplot(aes(x = source, length_cnv)) +
  geom_boxplot(aes(fill = source), alpha = 0.8) +
  theme_minimal() +
  scale_y_log10()



1000000

cnv_df %>%
  filter(length_cnv < 100000)
  


cnv_df %>%
  ggplot(aes(length_cnv)) +
  geom_density(aes(fill = source), alpha = 0.4) +
  xlim(0, 100000)





cnv_region <- GRanges(
  seqnames = "1",
  ranges = IRanges(1000, 10000000))

genes_region <- test1 %>% 
  select(chrom, start_position, end_position) %>%
  rename(start = start_position, end = end_position) %>%
  makeGRangesFromDataFrame()

overlaps_hits <- pintersect(genes_region, cnv_region)

percentOverlap <- width(gr2) /  width(overlaps_hits)


library(tidyverse)
library(randomForest)



url <- 'https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data'

df <- read_csv(url, col_names = c('age', 'sex', 'cp', 'testbps', 'chol', 'fbs', 'restecg',
                                  'thalach', 'exang', 'oldpeak', 'slope', 'ca',
                                  'thal', 'hd'))

df[df == '?'] <- NA

df <- df %>%
  mutate(sex = as.factor(if_else(sex == 1, 'M', 'F'))) %>%
  mutate(cp = as.factor(cp),
         ca = as.factor(ca),
         fbs = as.factor(fbs),
         restecg = as.factor(restecg),
         exang = as.factor(exang),
         slope = as.factor(slope),
         exang = as.factor(as.integer(ca)),
         thal = as.factor(thal)) %>%
  mutate(hd = as.factor(if_else(hd == 0, 'Healthy', 'Unhealthy')))


test66
go <- Ontology("hp")


chrom_tmp <- '1'
start_tmp <- 1000
end_tmp <- 10000


chromPlot::hg_cytoBandIdeo %>%
  filter(Chrom %in% chrom_tmp) %>%
  mutate(keep = map2_chr(Start, End, function(x,y) c(1000, 243700000) %overlaps% c(x,y)))



test66 %>%
  mutate(description = map_chr(hp, function(x) termDesc(term(go, x))))
  


set.seed(42)
  
df <- rfImpute(hd ~ ., data = df, iter = 6)
model  <- randomForest(hd ~ ., data = df, proximity = TRUE, ntree = 10000)


tmp_plot <- tibble(
  trees = rep(1:nrow(model$err.rate), 3),
  type = rep(c('OOB', 'healthy', 'unhealthy'), each = nrow(model$err.rate)),
  error = c(model$err.rate[,'OOB'],
            model$err.rate[,'Healthy'],
            model$err.rate[,'Unhealthy'])
)

ggplot(tmp_plot, aes(x = trees, y = error)) + geom_point(aes(color = type)) + geom_line(aes(group = type, color = type)) +
  theme_minimal()


glimpse(df)

tmp_df <- read.table('/home/cbl02/Storage/remot/7mer_mutation_rate_nonCodon.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, dec = ',')

cor(tmp_df$European.autosomal.substitution.probability.from.reference.to.alternative, tmp_df$Asian.autosomal.substitution.probability.from.reference.to.alternative)

summary(tmp_df)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org.Hs.eg.db,
                         gr=gr)


data_raw <- cnv_df %>% mutate(keep = NA) %>% slice(1:10)


cnv_df

map(data_raw, function(x) x[,3] - x[,4])

data_raw %>%
  mutate(keep = c(end, start) %overlaps% c(10529, 10000000))



df_output <- cnv_df %>%
  filter(chrom == '1') %>%
mutate(keep = map2_lgl(start, end, function(x, y) c(x, y) %overlaps% c(1000, 1000000))) %>%
  filter(keep == TRUE) %>%
  select(-keep)






test1



data_tmp <- lncrna_coord %>% filter(chrom == '1') %>%
  mutate(keep = 0)

for (i in 1:nrow(data_tmp)) {
  data_tmp$keep[i] <- c(data_tmp$start[i], data_tmp$end[i]) %overlaps% c(1000, 10000000)
}

data_tmp %>% filter(keep == 1) %>% select(id) %>% distinct() %>% pull(id)


ideoTrack <- IdeogramTrack(genome="hg19", chromosome= 'chr1')
plotTracks(ideoTrack, from= 1000 , to= 100000, showBandId=TRUE,
           cex.bands=0.5)
ideoTrack <- IdeogramTrack(genome="hg19", chromosome="chrX")
a <- plotTracks(ideoTrack, from=85e6, to=129e6)
  
  filtered_genes <- test1 %>% select(entrez_id) %>% pull()  %>% as.character()

  ggo <- groupGO(gene     = filtered_genes,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 13,
                 readable = TRUE),
  
  
  
test20 %>%
    as_tibble() %>%
    filter(Count != 0) %>%
    arrange(desc(Count)) %>%
    separate(geneID, sep = '/', into = as.character(1:1000)) %>%
    gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio) %>%
    select(-delete) %>%
    na.omit() %>%
    filter(Description == 'axonogenesis') %>%
    distinct()
  
  
  library(Gviz)
  
  ideoTrack <- IdeogramTrack(genome="hg19", chromosome="chr4")
  plotTracks(ideoTrack, from=85e6, to=129e6)
  
  
  tad 
  
  
  lncrna_coord
  
  
  a <- check_regions('1', 1000, 1000000)
  
  lncrna %>% filter(id %in%)
  

check_tads <- function(chrom = NULL, start = NULL, end = NULL) {
  

        chrom_tmp <-  chrom
        start_tmp <-  start
        end_tmp <-  end
    

  tmp_df <- tad %>%
      filter(chrom == !!chrom_tmp) %>%
      mutate(check_tad = if_else(start_tmp <= start & end_tmp >= start, 1, 0 )) %>%
      mutate(check_tad2 = if_else(start_tmp <= end & end_tmp >= end, 1, 0 )) %>%
      mutate(check_tad3 = if_else(start_tmp <= start & end_tmp >= end, 1, 0 )) %>%
      mutate(check_final_tad = if_else(check_tad  + check_tad2 + check_tad3 > 0, 1, 0 )) %>%
      select(id, check_final_tad) %>%
      filter(check_final_tad > 0)
  
  if (nrow(tmp_df) > 0) {
    
    list_tads <- as.vector(tmp_df %>% pull(id))
    
    
  } else {
    
    list_tads <- 0
  
  }
  
  return(list_tads)
  
  
}


check_regions('3', 1020000, 3900000 )
  
start = 1720000
end = 3400000

  
if_else(start_tmp >= start & end_tmp >= end, 1, 0)
    
  
       
  tad %>% slice(1027)
  
    
  test_df %>% filter(tad_check > 0)
  
  
  
  
  barplot(ggo)
  
  query_tmp <- entrez_summary(db="pubmed", id= test21[['ids']])
  
  title <- unname(map_chr(query_tmp, function(x) x[["title"]]))
  n_cites <- unname(map_chr(query_tmp, function(x) x[["pmcrefcount"]]))
  
  df_output <- tibble(title = title, n_cites = n_cites)
  
  
  datatable(v, options = list(
    pageLength = 5, autoWidth = TRUE, style = 'bootstrap'))
  
  ggo %>% 
    as_tibble() %>%
    filter(Count != 0) %>%
    arrange(desc(Count)) %>%
    separate(geneID, sep = '/', into = as.character(1:1000)) %>%
    gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio) %>%
    select(-delete) %>%
    na.omit()
    
  
  ggo %>%
    as_tibble() %>%
    arrange(desc(Count)) %>%
    slice(1:10) %>%
    as_tibble() %>% 
    # mutate(p.adjust = -log10(p.adjust)) %>%
    ggplot(aes(reorder(Description, Count), Count)) +
    geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
    scale_fill_viridis_d() +
    coord_flip() +
    xlab('') +
    ylab('-log10(p-adjusted)') +
    ggtitle('Pathway analysis') +
    geom_text(
      aes(label = GeneRatio, y = Count + 0.05),
      position = position_stack(vjust = 0.5),
      vjust = 0
    ) +
    theme_minimal() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  

a <- entrez_search(db="pubmed", term="22q11.2", retmax = 200 )


a <- entrez_summary(db="pubmed", id=a$ids)



b <- unname(map_chr(a, function(x) x[["title"]]))
b <- unname(map_chr(a, function(x) x[["pmcrefcount"]]))


##


##
library(httr)

a <- GET('https://ghr.nlm.nih.gov/condition/alzheimer-disease?report=json')

hello <- test1

go <- Ontology("mp")

test_tmp <- hello %>% select(gene) %>% pull()
mgi_tmp <- mgi %>%  filter(gene %in% test_tmp)
mgi_tmp <- mgi_tmp %>% separate(pheno, into = LETTERS[1:230], sep = ' ') %>%
  gather('delete', 'mpo_id', -gene, -entrez_id, -gene_mouse, -mgi) %>%
  filter(mpo_id != '') %>%
  select(-delete)

vector_mpo <- mgi_tmp %>% select(mpo_id) %>% unique() %>%
  mutate(description = map_chr(mpo_id, function(x) termLabel(term(go, x)))) %>%
  mutate(description = str_remove(description, ' phenotype'))

mgi_tmp <- mgi_tmp %>% left_join(vector_mpo)

###

go <- Ontology("mp")

test_tmp <- test15 %>% select(gene) %>% pull()
mgi_test <- mgi %>% filter(gene %in% test_tmp)
mgi_test <- mgi_test %>% separate(pheno, into = LETTERS[1:230], sep = ' ') %>%
  gather('delete', 'mpo_id', -gene, -entrez_id, -gene_mouse, -mgi) %>%
  filter(mpo_id != '') %>%
  select(-delete)

vector_mpo <- mgi_test %>% select(mpo_id) %>% unique() %>%
  mutate(description = map_chr(vector_mpo$mpo_id, function(x) termLabel(term(go, x)))) %>%
  mutate(description = str_remove(description, ' phenotype'))

mgi_test <- mgi_test %>% left_join(vector_mpo)

a <- tibble(a = 25, b = '<img src="/download.jpg" height="52"></img>')

datatable(a, escape = FALSE)

  e_charts(description) %>% 
  e_bar(n, name = "Phenotype") %>%
  e_flip_coords() # flip axis


  mgi_test %>% count(description) %>% arrange(n) %>%
  e_charts() %>% 
  e_treemap(description, description, n) %>% 
  e_title("Description")

library(rols)

entrez_db_summary()

entrez_db_searchable("nuccore")

  manolito <- entrez_search(db="gene", term="(DMD[Gene Name]) AND Homo sapiens[Organism]")

  all_recs <- entrez_fetch(db="gene", id=1756, rettype="gene_table")
  
a <- entrez_summary(db = 'gene', id = 1752)
### 

  hgcn_genes %>%
  slice(1:100) %>%
  select(gene, pLI) %>%
  mutate(haplo = if_else(haplo == 1, 'Yes', 'No'))
  datatable() %>%
  formatStyle(
      'pLI',
      background = styleColorBar(c(0,1), '#ca7171'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    )



###


hgMale <- gganatogram(data=hgMale_key, fillOutline='#a6bddb', organism='human', sex='male', fill="value") + theme_void()





gtex %>%
  filter(gene == 'DDX11L1') %>%
  ggplot(aes(reorder(tissue, -value), value)) +
  geom_col(aes(fill = tissue), color = 'black', show.legend = FALSE) +
  theme_minimal() +
  xlab('Tissue') +
  ylab('Median TPM')


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
