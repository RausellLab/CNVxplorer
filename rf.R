
library(randomForest)

hgcn_genes2 <- hgcn_genes %>% drop_na(start_position, end_position) 



cnv_df2 <- cnv_df %>% 
  mutate(n_genes = NA) %>% 
  mutate(source = as.factor(source)) %>%
  filter(source == 'decipher' | source == 'gnomad_v2') %>%
  sample_n(10000)


start_time <- Sys.time()

for ( i in 1:nrow(cnv_df2)) {
  
  print(i)
  
  tmp_df <- hgcn_genes2 %>% 
    filter(chrom == cnv_df2$chrom[i]) %>%
    rowwise() %>%
    mutate(keep = c(start_position, end_position) %overlaps% c(cnv_df2$start[i], cnv_df2$end[i])) %>%
    filter(keep == TRUE)
  
  if (nrow(tmp_df) == 0) {
    
    cnv_df2$n_genes[i]  <- 0
    cnv_df2$pli[i]  <- 0
    cnv_df2$rvis[i]  <- 0
    
  } else {
    
    # Variable: number of genes disrupted
    cnv_df2$n_genes[i]  <- nrow(tmp_df)
    # Variable: sum pLI score
    cnv_df2$pli[i]  <- tmp_df %>% pull(pLI) %>% sum(na.rm = TRUE)
    # Variable: sum rvis score
    cnv_df2$rvis[i]  <- tmp_df %>% pull(rvis) %>% sum(na.rm = TRUE)

    
  }
 
  
  tmp_omim <- tmp_df %>% count(omim) %>% filter(omim == 'Yes') %>% pull(n)
  
  # Variable: nº genes included in OMIM
  if (length(tmp_omim) == 0) {
    cnv_df2$omim[i] <- 0
    
  } else {
    cnv_df2$omim[i] <- tmp_omim
  }
  
}
end_time <- Sys.time()
end_time - start_time

cnv_df2$source <- droplevels(cnv_df2$source)

set.seed(100)
train <- sample(nrow(cnv_df2), 0.7*nrow(cnv_df2), replace = FALSE)
TrainSet <- cnv_df2[train,]
ValidSet <- cnv_df2[-train,]


model1 <- randomForest(source ~ length_cnv + n_genes + pli + rvis + omim, data = TrainSet, importance = TRUE)

predict(model1, cnv_df)

importance(model1, type = 1)
