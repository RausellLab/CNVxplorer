
#3

hgcn_genes


test3111 %>%
  rowwise() %>%
  mutate(is_overlap = c(hgcn_genes$start_position[1:100], hgcn_genes$end_position[1:100]) %overlaps% c(start, end))


hgcn_genes %>%
  rowwise() %>%
  mutate(is_overlap = c(start_position, end_position) %overlaps%  c(test3111$start[1], test3111$end[1])) %>%
  filter(isTRUE(is_overlap)) %>%
  nrow()


for (i in 1:nrow(test3111)){
  
  
  tmp_start <- test3111$start[1]
  tmp_end <- test3111$end[1]
  tmp_genes <- hgcn_genes %>%
    rowwise() %>%
    mutate(is_overlap = c(start_position, end_position) %overlaps%  c(tmp_start, tmp_end)) %>%
    filter(isTRUE(is_overlap)) %>%
    pull(gene)
  
  genes_not_cnv <- tmp_genes[!tmp_genes %in% (test2019 %>% pull(gene))]
  
  test3111$n_genes[i] <- tmp_n
  test3111$n_genes_not_cnv[i] <- length(genes_not_cnv)
  
}


##


test9211 <<- filtered_genes_cnv
filtered_tissue <- 'Adipose...Subcutaneous'

gtex %>%
  filter(tissue == !!filtered_tissue) %>%
  filter(gene %in% !!test9211) %>%
  ggplot(aes(reorder(gene, -value), value)) +
  geom_col(aes(fill = tissue), color = 'black', show.legend = FALSE) +
  # theme_fancy() +
  xlab('Tissue') +
  ylab(paste('Median TPM')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position='none') +
  ggtitle(paste('Gene expression:', filtered_tissue))





df <- read_csv('/home/cbl02/Storage/remot/test_to_r')

df2 <- hgcn_genes %>% filter(chrom == 20)

df %>%
  select(id_chunk, start, end, rank) %>%
  na.omit() %>%
  rowwise() %>%
  mutate(is_gene =  c(start, end) %overlaps% c(df2$start_position, df2$end_position)) %>%
  ungroup() %>%
  # count(is_gene)
  ggplot(aes(id_chunk, rank)) + 
  geom_point()
  ggplot(aes(is_gene, rank)) + 
  geom_boxplot()
  


test231312321$ID

bp_go <- godata('org.Hs.eg.db', ont="BP")
mf_go <- godata('org.Hs.eg.db', ont="MF")
cc_go <- godata('org.Hs.eg.db', ont="CC")

goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Jiang")

a <- GOSemSim::termSim(test231312321$ID, test231312321$ID, semData = hsGO, method = 'Resnik')
###




b <- GOSemSim::mgeneSim(tes912$entrez_id, mf_go, measure = 'Resnik')
c <- GOSemSim::mgeneSim(tes912$entrez_id, mf_go, measure = 'Resnik', combine = 'avg')

b == c


go1 = c("GO:0004022","GO:0004024","GO:0004174")
go2 = c("GO:0009055","GO:0005515")
mgoSim(go1, go2, semData=hsGO, measure="Wang", combine= 'BMA')



dgv_df


######
# Given a genomic interval chrom:start-end, identify genes that are overlapping


final_result <- tibble()


for (i in 1:nrow(dgv_df)) {

print(i)
tmp_cnv_id <-  dgv_df$id[i]
tmp_cnv_start <- dgv_df$start[i]
tmp_cnv_end <- dgv_df$end[i]
tmp_cnv_chrom <- dgv_df$chrom[i]

hpo_down <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/hp.obo')

tmp_result <- hgcn_genes %>%
  filter(chrom == tmp_cnv_chrom) %>%
  rowwise() %>%
  mutate(overlap = c(start_position, end_position) %overlaps% c(tmp_cnv_start, tmp_cnv_end)) %>%
  ungroup() %>%
  filter(isTRUE(overlap)) %>%
  select(entrez_id) %>%
  mutate(id_cnv = tmp_cnv_id)

final_result <- final_result %>% bind_rows(tmp_result)

}

library(tidyverse)
library(readxl)

df_paper <- read_xlsx('/home/cbl02/Desktop/test_paper_antonio.xlsx', sheet = 1, skip = 2)

df_paper <- df_paper %>% 
  filter(Consequence_severity == 'High probability LoF') %>%
  filter(Manual.annotation == 'LoF') %>%
  filter(!str_detect(GeneName, '^OR')) %>%
  pull(GeneName)
  
test_df <- hgcn_genes %>%
  filter(gene %in% df_paper)





###






df_windows <- read_tsv('/home/cbl02/Storage/remot/test_chunks/chr20_window')
  #na.omit() %>%
  #mutate(id = row_number()) %>%
  #arrange(start) %>%
  #mutate(is_gene = NA)
df_windows %>% filter(obs == 2) %>% slice(1) %>% pull(exp)
df_windows %>% filter(start == 65901)  %>% pull(exp)

refGR <- makeGRangesFromDataFrame(df_windows)
testGR <- makeGRangesFromDataFrame(hgcn_genes %>% filter(chrom == '20') %>%
                                     rename(start = start_position,
                                            end = end_position))


hits <- findOverlaps(refGR, testGR) %>% as_tibble() %>% rename(id = queryHits)

hits <- hits %>% count(id) %>% arrange(desc(n)) %>% rename(n_hits = n) %>% right_join(df_windows)

hits  %>% mutate(n_hits = ifelse(isNA(n_hits), FALSE, TRUE)) %>%
  mutate(n_hits = as.factor(n_hits)) %>%
  ggplot(aes(n_hits, ratio)) +
  geom_boxplot(aes(fill = n_hits)) +
  coord_cartesian(ylim = c(0, 10^7))

vector_start_genes <- hgcn_genes %>% filter(chrom == '20') %>% pull(start_position)
vector_end_genes <- hgcn_genes %>% filter(chrom == '20') %>% pull(end_position)

df_windows %>%
  time_decompose(ratio, method = 'stl') %>%
  #anomalize(remainder, method = 'ged', alpha = 0.05, max_anoms = 0.2) %>%
  plot_anomaly_decomposition()


for (i in 1:nrow(df_windows)) {
  print(i)
  i <- 3000
  tmp_start = df_windows$start[i]
  tmp_end =  df_windows$end[i]
  df_windows$is_gene[i] <- c(tmp_start, tmp_end) %overlaps% c(vector_start_genes, vector_end_genes)

}
  print(i)

  df_windows %>% ggplot(aes(is_gene, ratio)) +
    geom_boxplot() +
    coord_cartesian(ylim = c(0,10^7))
  
  df_windows %>% ggplot(aes(ratio)) +
    geom_density(aes(fill = is_gene), alpha = 0.3) + 
    coord_cartesian(xlim = c(0,10^8))
  
  
t.test(df_windows %>% filter(is_gene == TRUE) %>% pull(ratio),
       df_windows %>% filter(is_gene == FALSE) %>% pull(ratio))

test <- df_windows %>% 

  rowwise() %>%
  mutate(is_gene = c(start, end) %overlaps% c(hgcn_genes %>% filter(chrom == '20') %>% pull(start_position), 
                                              hgcn_genes %>% filter(chrom == '20') %>% pull(end_position))) %>%
  ungroup()

test %>% count(is_gene)

    tmp_gr <- GRanges(
      seqnames = paste0('chr',1),
      ranges = IRanges(10272:10275, width=1))
    a <- gscores(cadd_1_3, tmp_gr)


test1931 <- test312
    
cnv_patho <- test1931 %>% filter(source == 'decipher')
cnv_nonpatho <- test1931 %>% filter(source == 'dgv')


    
start_pos <- 34813719
end_pos <- 36278623  


interval_values <- seq(from = start_pos, to = end_pos, length.out = 200)

df_interval <- matrix(interval_values, ncol = 2, byrow = TRUE)
colnames(df_interval) <- c('start', 'end')
df_interval <- as_tibble(df_interval)

query <- IRanges(df_interval$start, df_interval$end)
subject_patho <- IRanges(cnv_patho$start, cnv_patho$end)
hits_intervals_patho <- countOverlaps(query, subject_patho)

subject_nonpatho<- IRanges(cnv_nonpatho$start, cnv_nonpatho$end)
hits_intervals_nonpatho <- countOverlaps(query, subject_nonpatho)


df_interval <- df_interval %>% mutate(n_patho = hits_intervals_patho, 
                                      n_nonpatho = hits_intervals_nonpatho) %>%
  mutate(id = row_number()) %>%
  gather('category', 'n_overlap', -start, -end, -id)


df_interval %>%
  mutate(category = if_else(category == 'n_patho', 'Pathogenic CNVs', 'Non-pathogenic CNVs')) %>%
  ggplot(aes(id, n_overlap)) +
  geom_col(aes(fill = category), color = 'black', show.legend = FALSE) + 
  #geom_point() +
 # geom_smooth(aes(color = category, group = category))
  #geom_line(aes(color = category, group = category)) + 
  theme_minimal() + 
  facet_wrap(~category, nrow = 2)

df_interval %>%
  mutate(category = if_else(category == 'n_patho', 'Pathogenic CNVs', 'Non-pathogenic CNVs')) %>%
  ggplot(aes(id, n_overlap)) +
  geom_density(aes(fill = category), color = 'black', show.legend = FALSE) + 
  #geom_point() +
  # geom_smooth(aes(color = category, group = category))
  #geom_line(aes(color = category, group = category)) + 
  theme_minimal() + 
  facet_wrap(~category, nrow = 2)



get_model_score <- function(chrom, start_coord, end_coord, df_genes, model1) {
  
  # chrom_coord <- '1'
  # start_coord <- 34813719
  # end_coord <- 36278623
  chrom_coord <- str_remove(chrom, 'chr')
  start_coordinates <- start_coord
  end_coordinates <- end_coord
  length_input_cnv <- end_coordinates - start_coordinates + 1
  tmp_df <- df_genes %>%
            filter(chrom == chrom_coord) %>%
            mutate(keep = NA) %>%
            rowwise() %>%
            mutate(keep = c(start_position, end_position) %overlaps% c(start_coordinates, end_coordinates)) %>%
            filter(keep == TRUE) %>% 
            select(-keep) %>%
            ungroup()
  
  input_model <- tibble(n_genes = NA, pli = NA, rvis = NA, omim = NA, 
                        length_cnv =  length_input_cnv)
  
  if (nrow(tmp_df) == 0) {
    
    input_model$n_genes[1]  <- 0
    input_model$pli[1]  <- 0
    input_model$rvis[1]  <- 0
    input_model$omim[1] <- 0
    
  } else {
    
    # Variable: number of genes disrupted
    input_model$n_genes[1]  <- nrow(tmp_df)
    # Variable: sum pLI score
    input_model$pli[1]  <- tmp_df %>% pull(pLI) %>% sum(na.rm = TRUE)
    # Variable: sum rvis score
    input_model$rvis[1]  <- tmp_df %>% pull(rvis) %>% sum(na.rm = TRUE)
    # Variable: nº genes included in OMIM
    tmp_omim <- tmp_df %>% count(omim) %>% filter(omim == 'Yes') %>% pull(n)
    if (length(tmp_omim) == 0) {
      input_model$omim[1] <- 0
    } else {
      input_model$omim[1] <- tmp_omim
      
    }
  }
  

  score_predicted <- predict(model1, input_model, type = "prob") %>% 
    as_tibble() %>% pull(decipher) %>% as.numeric() * 100
  list_result <- list(score_predicted, input_model$n_genes)
  return(list_result)
}



a <- c("DOID:14095", "DOID:5844", "DOID:2044", "DOID:8432", "DOID:9146",
       "DOID:10588", "DOID:3209", "DOID:848", "DOID:3341", "DOID:252")
b <- c("DOID:9409", "DOID:2491", "DOID:4467", "DOID:3498", "DOID:11256")
doSim(a, b, measure="Resnik")


g1 <- c("84842", "2524", "10590", "3070", "91746")
g2 <- c("84289", "6045", "56999", "9869")

geneSim(g1, g2, measure="Resnik", combine="BMA")

hpo_omim 

library(HPOSim) 

list_genes <- test2019 %>% pull(gene)
list_hpo <- hpo_genes %>% filter(gene %in% list_genes) %>% pull(hp)
# check omim terms with pmid
list_omim_diseases <- hpo_omim %>% select(term) %>% filter(str_detect(term, 'OMIM')) %>% distinct() %>% pull() 

result_sim <- tibble()

for (i in 1:length(list_omim_diseases)) {
  
  print(i)
  tmp_list <- hpo_omim %>% filter(term == !!list_omim_diseases[i]) %>% pull(hpo)
  tmp_score <- getTermListSim(list_hpo, tmp_list, method = 'Resnik', combinemethod = 'funSimMax',  IC = IC)
  
  tmp_df <- tibble(omim_disease = list_omim_diseases[i], 
                   score = tmp_score)
  
  result_sim <- rbind(result_sim, tmp_df)
}

a <- c('HP:0100807', 'HP:0001166')

b <- c('HP:0011793', 'HP:0001166')

d <- c('HP:0011792', 'HP:0002861')

anno1<-c("HP:0011793", "HP:0011795")
anno2<-c("HP:0011792", "HP:0011794")


# HPOSim::calcTermSim('HP:0011793', 'HP:0001909', IC = IC, method = 'Resnik')

  
library(xml)


your.ids <- c("26386083","26273372","26066373","25837167","25466451","25013473")
# rentrez function to get the data from pubmed db
fetch.pubmed <- entrez_fetch(db = "pubmed", id = your.ids,
                             rettype = "xml", parsed = T)
# Extract the Abstracts for the respective IDS.  
abstracts = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
  xmlValue(xmlChildren(x)$Abstract))


entrez_db_links()


tmp <- entrez_link(dbfrom='pubmed', id= 25581431, db='gene')

tmp$links
#### CHI-SQUARED TEST

list_panel_names <- panel_total %>% select(Level4) %>% distinct() %>% pull()
input_genes <- test2019 %>% select(gene) %>% pull()

fisher_result <- tibble()
for (i in 1:length(list_panel_names)) {
print(i)
input_test <- matrix(rep(NA, 4), ncol = 2, dimnames = list(c('in_panel', 'no_panel'), c('col1', 'col2') ))

input_test[1,][2] <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% pull(gene) %in% input_genes %>% sum()
input_test[1,][1] <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% nrow() - input_test[1,][2] 
input_test[2,][2] <- length(input_genes) - input_test[1,][2] 
input_test[2,][1] <- 19146 - input_test[2,][2] - input_test[1,][1] - input_test[1,][2]

tmp_tibble <- tibble(name_panel = list_panel_names[i], 
                     p_value = fisher.test(input_test, 'greater')$p.value,
                     gene_ratio = paste(as.character(input_test[1,][2]), '/',
                                        as.character(input_test[1,][2] + input_test[1,][1]))
)

fisher_result <- rbind(fisher_result, tmp_tibble)
}

tmp <- matrix(rep(NA, 4), ncol = 2)

tmp2 <- matrix(rep(NA, 4), ncol = 2)


tmp[1,][2] <-  2
tmp[2,][2] <-  19 - tmp[1,][2]
tmp[1,][1] <- 10 - tmp[1,][2]
tmp[2,][1] <- 19146 - tmp[1,][1] - tmp[2,][2] - tmp[1,][2] 

tmp[1,][2] <- 3 
tmp[2,][2] <-  9
tmp[1,][1] <- 13
tmp[2,][1] <-  4


fisher.test(tmp, alternative = 'greater')$p.value


tmp2[1,][2] <- 9 
tmp2[2,][2] <-  4
tmp2[1,][1] <- 3
tmp2[2,][1] <-  13

fisher.test(tmp2, alternative = 'greater')$p.value



#### POISSON DISTRIBUTION

tmp1 <- test1492 %>% filter(source == 'decipher') %>% as_tibble()



median_sliding_w <- tmp1 %>% filter(source == 'decipher') %>%
  pull(length_cnv) %>%
  median()

total_length_query = 36278623 - 34813719 + 1
n_windows = round(total_length_query / median_sliding_w, 0)
34813719


36278623


df_result <- tibble(id = 1, from = 34813719, to = 34813719 + median_sliding_w)

for (i in 2:n_windows) {
  print(i)
  
  from_tmp <- df_result[,3][i-1,] %>% pull()
  to_tmp <- df_result[,3][i-1,] %>% pull() + median_sliding_w
  
  tmp_tibble <- tibble(id = i, 
                       from =  from_tmp,
                       to = to_tmp)
  
  df_result <- rbind(df_result, tmp_tibble)
  
}
  

gr_input <- GRanges(
  seqnames= test2019$chrom,
  ranges= IRanges(start= test2019$start_position, end = test2019$end_position)
)

gr_cnvs <- GRanges(
  seqnames= test821321313$chrom,
  ranges= IRanges(start= test821321313$start_position, end = test821321313$end_position)
)

gr_intersect <- intersect(gr_input, gr_cnvs)


plot(dpois(1:19, nrow(tmp1)/n_windows))


test711

data_raw <- hgcn_genes %>% filter(chrom == 1)

jajaa <- data_raw  %>% mutate(keep = NA) %>%
  rowwise() %>%
  mutate(keep = c(start_position, end_position) %overlaps% c(test711$start, test711$end)) %>%
  filter(keep == TRUE) %>% 
  select(-keep) %>%
  ungroup()

jajaa <- jajaa %>% mutate(pos = round((end_position + start_position)/2), 0) %>%
  mutate(id = rep('gene', n())) %>%
  select(id, pos)

test711 <- test711 %>%
  mutate(id = as.factor(seq(1:n()))) %>%
  select(id, start, end) %>%
  gather('rm', 'pos', -id) %>%
  select(-rm) %>%
  rbind(jajaa) %>%
  mutate(id_color = if_else(id == 'gene(s)', 'steelblue', 'red'))


  ggplot() +
  geom_point(data = test711, aes(pos, id), color = 'black', shape = 21) + 
  geom_path(data =test711 %>% filter(id != 'gene'), aes(pos, id, group = id)) + 
  coord_cartesian(xlim = c(34813719, 36278623 )) + 
  theme_fancy()

######

df <- test771





df_tmp <- df %>% 
  select(term) %>%
  distinct() %>%
  mutate(id_row = row_number())
  

vector_hpo <- df %>% select(term) %>% distinct() %>% pull()


validate(
  need(length(vector_hpo) > 1, "Please, select more than one phenotype term.")
)

vector_genes <- df %>% select(gene) %>% distinct() %>% pull()

list_result <- replicate(length(vector_hpo), NA, simplify = FALSE)


for (i in 1:length(vector_hpo)) {
  
  list_result[[i]] <- df %>% filter(term == !!vector_hpo[i]) %>% pull(gene) %in% vector_genes %>% which()
  names(list_result)[i] <- vector_hpo[i]
  
}

upset(fromList(list_result), empty.intersections = "on", order.by = "freq",
      point.size = 3.5, line.size = 2, number.angles = 0,
      mainbar.y.label = "Phenotype Terms Intersections", sets.x.label = "Genes Associated Per Phenotype Term",
      text.scale = c(1.3, 1.3, 1, 1, 2, 2))















start_p <- 34813719
end_p <- 36278623

points_check <- seq(start_p, end_p, by = 100)

df_test <- tibble(id = 1:length(points_check), pos = points_check, n_times = NA)
n_rows <- nrow(df_test)

ole_df <- test1492 %>% filter(source == 'decipher')

for (i in 1:n_rows) {
  score = 0
  for (j in 1:nrow(ole_df)) {
  keep <- df_test$pos[i] %overlaps% c(ole_df$start_position[j], ole_df$end_position[j])
  if (isTRUE(keep)) {
    score <- score +  1
  }
  }
  df_test$n_times[i] <- score
  print(paste0(i, '/', n_rows))
}


theo_value <- ole_df %>% mutate(diff = end_position - start_position + 1) %>% pull(diff) %>% sum()    / (end_p - start_p + 1)

df_test %>% ggplot(aes(pos, n_times)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = theo_value, color = 'red')








hp_chosen <- test322
genes_chosen <- test2019 %>% select(entrez_id) %>% pull()


df_tmp <- hpo_genes %>%
  filter(hp %in% hp_chosen) %>%
  filter(entrez_id %in% genes_chosen)

df_tmp2 <- df_tmp %>%
  select(hp) %>%
  distinct() %>%
  mutate(description = ifelse(is_null(map_chr(hp, function(x) termDesc(term(go, x)))), '', 
                              map_chr(hp, function(x) termDesc(term(go, x)))
                              ))

df_tmp <- df_tmp %>% 
  left_join(df_tmp2, by = 'hp')

df_tmp







test1936 %>% ggplot(aes(phast100)) + geom_histogram()

save_object <- save_object[! save_object %in% test912$gene]


genes_outside <- hgcn_genes %>% filter(ensembl_gene_id %in% save_object)


test1 <- test766 %>% select(gene, term) %>% mutate(id_row = row_number())

vector_hpo <- test1 %>% select(term) %>% distinct() %>% pull()
vector_genes <- test1 %>% select(gene) %>% distinct() %>% pull()

list_result <- replicate(length(vector_hpo), NA, simplify = FALSE)


for (i in 1:length(vector_hpo)) {
  
  list_result[[i]] <- test1 %>% filter(term== !!vector_hpo[i]) %>% pull(id_row)
  names(list_result)[i] <- vector_hpo[i]
  
}


upset(fromList(list_result), empty.intersections = "on", order.by = "freq")


library(biomaRt)


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

interval_genes <- getBM(attributes = c('ensembl_gene_id_version', 'start_position','end_position', 'chromosome_name'), 
                        mart = human ) %>% as_tibble()



interval_genes %>% filter(str_detect(ensembl_gene_id_version, 'ENSG00000270806'))


m <- matrix(1:100, ncol = 5, nrow = 10) %>% as.data.frame()

plus_1 <- function(x, y = 1) {
  x <- x + y + 1 
  return(x)
}


m %>%
map_df(function(x) plus_1(x, 2))

enrichDO(gene  = test4577,
         universe      = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character(),
         # pAdjustMethod = "BH",
         pvalueCutoff  = 0.05,
         readable = TRUE)

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
