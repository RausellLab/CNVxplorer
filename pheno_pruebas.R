

library(ontologyIndex)
library(ontologySimilarity)
library(GOSemSim)




##




  
ddd <- read_tsv('C:/Users/Requena/Desktop/ddd169_phenotypes.tsv') %>%
  mutate(child_hpo = child_hpo %>% str_replace_all('\\|', ', ') %>% str_split(', ')) %>%
  rename(hpo_list = child_hpo)





total_genes <- base::split(hpo_genes$hp, hpo_genes$gene)

counter = 0

df_result <- tibble(id_patient = 'NA', rank_resnik_avg = NA, rank_lin_avg = NA)

for (i in 1:nrow(ddd)) {
  print(i)
  causal_gene <- ddd$causal_gene[i]
  patient_id <- ddd$patient_ID[i]
  hpo_list <- ddd$hpo_list[i]
  
  not_included <- hpo_list[[1]][!hpo_list[[1]] %in% hpo_dbs$id]
  any_dup <- hpo_list[duplicated(hpo_list)]
  
  if (length(not_included) > 0) {
    if (length(not_included) == length(hpo_list[[1]])) {
      counter = counter + 1
      print('obsolete')
      next
    }
    next
  }
  


  tmp_result <- get_profile_sims(ontology = hpo_dbs, profile = ddd$hpo_list[i], term_sets = total_genes,
                 term_sim_method = 'resnik')

  output_resnik <- tmp_result %>% 
    enframe() %>%
    rename(gene = name) %>%
    arrange(desc(value)) %>%
    mutate(rank = row_number()) %>%
    filter(gene == causal_gene) %>%
    pull(rank)
  
  tmp_result <- get_profile_sims(ontology = hpo_dbs, profile = ddd$hpo_list[i], term_sets = total_genes,
                                 term_sim_method = 'lin')
  
  output_lin <- tmp_result %>% 
    enframe() %>%
    rename(gene = name) %>%
    arrange(desc(value)) %>%
    mutate(rank = row_number()) %>%
    filter(gene == causal_gene) %>%
    pull(rank)
  
  tmp_df <- tibble(id_patient =  patient_id,
                   rank_resnik_avg = output_resnik,
                   rank_lin_avg = output_lin)
  
  df_result <- df_result %>% bind_rows(tmp_df)
}

tmp2 <- df_result %>% 
  na.omit() %>% 
  arrange(rank_resnik_avg) %>%
  mutate(top_resnik = case_when(
    rank_resnik_avg <= 50 ~ 'top_50',
    rank_resnik_avg <= 100 ~ 'top_100',
    rank_resnik_avg <= 150 ~ 'top_150',
    rank_resnik_avg <= 200 ~ 'top_200',
    rank_resnik_avg <= 250 ~ 'top_250',
    rank_resnik_avg <= 300 ~ 'top_300',
    rank_resnik_avg <= 350 ~ 'top_350',
    rank_resnik_avg > 350 ~ 'total'
    
),
top_lin = case_when(
  rank_lin_avg <= 50 ~ 'top_50',
  rank_lin_avg <= 100 ~ 'top_100',
  rank_lin_avg <= 150 ~ 'top_150',
  rank_lin_avg <= 200 ~ 'top_200',
  rank_lin_avg <= 250 ~ 'top_250',
  rank_lin_avg <= 300 ~ 'top_300',
  rank_lin_avg <= 350 ~ 'top_350',
  rank_lin_avg > 350 ~ 'total'

))

my_order <- c('top_50', 'top_100', 'top_150', 'top_200', 'top_250', 'top_300', 'top_350', 'total')

tmp_resnik <- tmp2 %>%
  count(top_resnik) %>%
  right_join(tibble(top_resnik = my_order, n = rep(0, 8)), by = 'top_resnik') %>%
  select(-n.y) %>%
  rename(n = n.x) %>%
  replace_na(list(n = 0)) %>%
  arrange(match(top_resnik, my_order)) %>%
  mutate(perc_resnik = n / sum(n) * 100) %>%
  select(top_resnik, perc_resnik) %>%
  mutate(cum_resnik = cumsum(perc_resnik))
  


tmp_lin <- tmp2 %>%
  count(top_lin) %>%
  right_join( tibble(top_lin = my_order, n = rep(0, 8)), by = 'top_lin') %>%
  select(-n.y) %>%
  rename(n = n.x) %>%
  replace_na(list(n = 0)) %>%
  arrange(match(top_lin, my_order)) %>%
  mutate(perc_lin = n / sum(n) * 100) %>%
  select(top_lin, perc_lin) %>%
  arrange(match(top_lin, my_order)) %>%
  mutate(cum_lin = cumsum(perc_lin))
  


tmp_resnik %>%
  left_join(tmp_lin, by = c('top_resnik' = 'top_lin')) %>%
  pivot_longer(cols = starts_with('cum')) %>%
  ggplot(aes(reorder(top_resnik, value), value)) +
  geom_line(aes(group = name, color = name)) +
  geom_point(aes(fill =name, group = name), shape = 21) +
  xlab('Rank of causative gene') +
  ylab('% of patients') +
  ggtitle('Performance with different similarity measures') +
  theme_bw()
  



jaccard_mp <- read_csv('D:/mp-hp-mappings/upheno_equivalence_model_jaccard.csv')

jaccard_mp <- jaccard_mp %>%
  mutate(c1 = str_remove(c1, '<http://purl.obolibrary.org/obo/'),
         c2 = str_remove(c2, '<http://purl.obolibrary.org/obo/')) %>%
  mutate(c1 = str_remove(c1, '>')) %>%
  mutate(c2 = str_remove(c2, '>')) %>%
  mutate(c1 = str_replace(c1,'_', ':')) %>%
  mutate(c2 = str_replace(c2,'_', ':')) %>%
  rename(mp = c1, hp = c2) %>%
  select(-contains('label'))


mgi_tmp <-  mgi %>% separate(pheno, into = LETTERS[1:230], sep = ' ') %>%
  gather('delete', 'mpo_id', -gene, -entrez_id, -gene_mouse, -mgi) %>%
  filter(mpo_id != '') %>%
  select(gene, mpo_id) %>%
  filter(gene %in% (hpo_genes$gene %>% unique()))


ddd_mg <- ddd %>% filter(causal_gene %in% (mgi_tmp %>% pull(gene) %>% unique()))

total_genes_mp <- base::split(mgi_tmp$mpo_id, mgi_tmp$gene)



hpo_list <- hpo_list[[1]]

empty_v <- c()

for (i in 1:length(total_genes_mp)){

tmp_gene <- total_genes_mp[[i]]
print(i)

jeje <- jaccard_mp %>%
  filter(hp %in% hpo_list,
         mp %in% tmp_gene) %>%
  nrow()
empty_v <- c(empty_v, jeje)

}















get_sim_score <- function(genes_vector, patient_terms) {
  
  genes_vector <- c('KIAA0319L', 'GJB3', 'GJB4')
  patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
  
  tmp_df <- hpo_genes %>% 
    filter(gene %in% genes_vector) %>%
    select(gene, hp)
  
  to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
  
  genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
  
  hpo_patient <- patient_terms
  
  names(hpo_patient) <- c('patient')  
  total_set <- c(genes_sample, hpo_patient)
  
  
  
  mat_test <- get_sim_grid(ontology = hpo_down, 
                           term_sets = total_set,
                           # term_sets2 = hpo_patient,
                           # combine = 'average',
                           term_sim_method = 'resnik',
                           combine = 'average') %>%
    as_tibble(rownames = 'gene') %>%
    filter(gene != 'patient') %>%
    select(gene, patient) %>%
    mutate(p_value = NA) %>%
    left_join(to_p_value)
    # rowwise() %>%
     # mutate(por_jajas = pmap_int(get_p))
    # mutate(p_value = pmap_int(jaja, ~get_p(patient_terms, patient, n_freq)))
    
    mat_te
    
    library(ontologyIndex)
    library(ontologySimilarity)
    library(GOSemSim)
    
    
    library(ontologyIndex)
    
    go_test <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/go.obo')
    
    get_descendants(go_test, 'GO:0045737')
    
    library(GOSemSim)
    
    hsGO <- godata('org.Hs.eg.db', ont="MF")
    
    
    
    # Read .obo file
    # Source: http://purl.obolibrary.org/obo/hp.obo
    
    hpo_down <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/hp.obo')
    
    genes_sample <- replicate(simplify=FALSE, n= 1, expr=minimal_set(hpo, sample(hpo$id, size= 1)))
    names(genes_sample) <- paste0('gene', seq(1,1))
    
    hpo_patient <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo, sample(hpo$id, size= 1)))
    names(hpo_patient) <- c('patient1')
    
    total_set <- c(genes_sample, hpo_patient)
    
    
    c('HP:0025013', 'HP:0003791')
    
    # Calculate Resnik
    
    anc_1 <- hpo_down$ancestors['HP:0025013'][[1]]
    anc_2 <- hpo_down$ancestors['HP:0003791'][[1]]
    anc_int <- intersect(anc_1, anc_2)
    
    desc_term <- hpo_down$children['HP:0000118'][[1]]
    
    resnik_value <- -log(length(get_descendants(hpo_down, 'HP:0000118')) / length(hpo_down$id))
    
    # Calculate Lin 
    
    ic_1 <- -log(length(get_descendants(hpo_down, 'HP:0025013')) / length(hpo_down$id))
    ic_2 <- -log(length(get_descendants(hpo_down, 'HP:0003791')) / length(hpo_down$id))
    
    lin_value <- 2 * (resnik_value / (ic_1 + ic_2))
    
    get_sim_grid(ontology = hpo_down, 
                 term_sets = total_set,
                 term_sim_method = 'lin',
                 combine = 'average')
    
    
    # resnik method = 1.234995
    # HP:0005289
    # HP:0004438
    
    
    get_sim_score <- function(genes_vector, patient_terms) {
      
      genes_vector <- c('KIAA0319L', 'GJB3', 'GJB4')
      patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
      
      tmp_df <- hpo_genes %>% 
        filter(gene %in% genes_vector) %>%
        select(gene, hp)
      
      to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
      
      genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
      
      hpo_patient <- patient_terms
      
      names(hpo_patient) <- c('patient')  
      total_set <- c(genes_sample, hpo_patient)
      
      
      
      mat_test <- get_sim_grid(ontology = hpo_down, 
                               term_sets = total_set,
                               # term_sets2 = hpo_patient,
                               # combine = 'average',
                               term_sim_method = 'resnik',
                               combine = 'average') %>%
        as_tibble(rownames = 'gene') %>%
        filter(gene != 'patient') %>%
        select(gene, patient) %>%
        mutate(p_value = NA) %>%
        left_join(to_p_value)
      # rowwise() %>%
      # mutate(por_jajas = pmap_int(get_p))
      # mutate(p_value = pmap_int(jaja, ~get_p(patient_terms, patient, n_freq)))
      
      mat_test$p_value <-  c(apply(mat_test, 1, get_p))
      
      return(mat_test)
      
      # get_sim_p(mat_test, group = ncol(mat_test))
    }
    
    
    get_p <- function(.patient_terms, .real_value, .n_hpo_gene) {
      
      # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
      # names(patient_terms) <- 'patient1'
      # n_hpo_gene <- 4
      fake_gene <- replicate(simplify=FALSE, n= 100, expr=minimal_set(hpo_down, sample(hpo_down$id, size= n_hpo_gene)))
      
      
      total_set <- c(fake_gene, hpo_patient)
      
      mat_test <- get_sim_grid(ontology = hpo_down, 
                               term_sets = total_set,
                               term_sim_method = 'resnik',
                               combine = 'average')
      
      p_value <- mat_test %>% 
        as_tibble() %>%
        select(patient) %>%
        filter(patient >= real_value) %>%
        nrow()
      
      p_value <- p_value / 100
      
      return(p_value)
      
    }
    
    
    
    
    
    
    
    hsGO <- godata('org.Hs.eg.db', ont="MF")
    goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
    
    GOSemSim:::termSim('GO:0004022', 'GO:0005515', semData = hsGO, method = 'Wang')
    GOSemSim:::infoContentMethod('GO:0004022', 'GO:0005515', method = 'Wang', hsGO)
    GOSemSim:::infoContentMethod_cpp('GO:0004022', 'GO:0005515', '' )
    
    
    ontologySimilarity::get_term_sim_mat(ontology = hpo_down, 
                                         information_content = c(0)
                                         method = 'resnik')
    
    
    library(ontologyIndex)
    library(ontologySimilarity)
    library(GOSemSim)
    
    
    library(ontologyIndex)
    
    go_test <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/go.obo')
    
    get_descendants(go_test, 'GO:0045737')
    
    library(GOSemSim)
    
    hsGO <- godata('org.Hs.eg.db', ont="MF")
    
    
    
    # Read .obo file
    # Source: http://purl.obolibrary.org/obo/hp.obo
    
    hpo_down <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/hp.obo')
    
    genes_sample <- replicate(simplify=FALSE, n= 1, expr=minimal_set(hpo, sample(hpo$id, size= 1)))
    names(genes_sample) <- paste0('gene', seq(1,1))
    
    hpo_patient <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo, sample(hpo$id, size= 1)))
    names(hpo_patient) <- c('patient1')
    
    total_set <- c(genes_sample, hpo_patient)
    
    
    c('HP:0025013', 'HP:0003791')
    
    # Calculate Resnik
    
    anc_1 <- hpo_down$ancestors['HP:0025013'][[1]]
    anc_2 <- hpo_down$ancestors['HP:0003791'][[1]]
    anc_int <- intersect(anc_1, anc_2)
    
    desc_term <- hpo_down$children['HP:0000118'][[1]]
    
    resnik_value <- -log(length(get_descendants(hpo_down, 'HP:0000118')) / length(hpo_down$id))
    
    # Calculate Lin 
    
    ic_1 <- -log(length(get_descendants(hpo_down, 'HP:0025013')) / length(hpo_down$id))
    ic_2 <- -log(length(get_descendants(hpo_down, 'HP:0003791')) / length(hpo_down$id))
    
    lin_value <- 2 * (resnik_value / (ic_1 + ic_2))
    
    get_sim_grid(ontology = hpo_down, 
                 term_sets = total_set,
                 term_sim_method = 'lin',
                 combine = 'average')
    
    
    # resnik method = 1.234995
    # HP:0005289
    # HP:0004438
    
    
    get_sim_score <- function(genes_vector, patient_terms) {
      
      genes_vector <- c('KIAA0319L', 'GJB3', 'GJB4')
      patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
      
      tmp_df <- hpo_genes %>% 
        filter(gene %in% genes_vector) %>%
        select(gene, hp)
      
      to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
      
      genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
      
      hpo_patient <- patient_terms
      
      names(hpo_patient) <- c('patient')  
      total_set <- c(genes_sample, hpo_patient)
      
      
      
      mat_test <- get_sim_grid(ontology = hpo_down, 
                               term_sets = total_set,
                               # term_sets2 = hpo_patient,
                               # combine = 'average',
                               term_sim_method = 'resnik',
                               combine = 'average') %>%
        as_tibble(rownames = 'gene') %>%
        filter(gene != 'patient') %>%
        select(gene, patient) %>%
        mutate(p_value = NA) %>%
        left_join(to_p_value)
      # rowwise() %>%
      # mutate(por_jajas = pmap_int(get_p))
      # mutate(p_value = pmap_int(jaja, ~get_p(patient_terms, patient, n_freq)))
      
      mat_test$p_value <-  c(apply(mat_test, 1, get_p))
      
      return(mat_test)
      
      # get_sim_p(mat_test, group = ncol(mat_test))
    }
    
    
    get_p <- function(.patient_terms, .real_value, .n_hpo_gene) {
      
      # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
      # names(patient_terms) <- 'patient1'
      # n_hpo_gene <- 4
      fake_gene <- replicate(simplify=FALSE, n= 100, expr=minimal_set(hpo_down, sample(hpo_down$id, size= n_hpo_gene)))
      
      
      total_set <- c(fake_gene, hpo_patient)
      
      mat_test <- get_sim_grid(ontology = hpo_down, 
                               term_sets = total_set,
                               term_sim_method = 'resnik',
                               combine = 'average')
      
      p_value <- mat_test %>% 
        as_tibble() %>%
        select(patient) %>%
        filter(patient >= real_value) %>%
        nrow()
      
      p_value <- p_value / 100
      
      return(p_value)
      
    }
    
    
    
    
    
    
    
    hsGO <- godata('org.Hs.eg.db', ont="MF")
    goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
    
    GOSemSim:::termSim('GO:0004022', 'GO:0005515', semData = hsGO, method = 'Wang')
    GOSemSim:::infoContentMethod('GO:0004022', 'GO:0005515', method = 'Wang', hsGO)
    GOSemSim:::infoContentMethod_cpp('GO:0004022', 'GO:0005515', '' )
    
    
    ontologySimilarity::get_term_sim_mat(ontology = hpo_down, 
                                         information_content = c(0)
                                         method = 'resnik')
    
    
    library(ontologyIndex)
    library(ontologySimilarity)
    library(GOSemSim)
    
    
    library(ontologyIndex)
    
    go_test <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/go.obo')
    
    get_descendants(go_test, 'GO:0045737')
    
    library(GOSemSim)
    
    hsGO <- godata('org.Hs.eg.db', ont="MF")
    
    
    
    # Read .obo file
    # Source: http://purl.obolibrary.org/obo/hp.obo
    
    hpo_down <- ontologyIndex::get_OBO('http://purl.obolibrary.org/obo/hp.obo')
    
    genes_sample <- replicate(simplify=FALSE, n= 1, expr=minimal_set(hpo, sample(hpo$id, size= 1)))
    names(genes_sample) <- paste0('gene', seq(1,1))
    
    hpo_patient <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo, sample(hpo$id, size= 1)))
    names(hpo_patient) <- c('patient1')
    
    total_set <- c(genes_sample, hpo_patient)
    jaja <- ontologySimilarity::descendants_IC(hpo_down)
    
    c('HP:0025013', 'HP:0003791')
    
    # Calculate Resnik
    
    anc_1 <- hpo_down$ancestors['HP:0025013'][[1]]
    anc_2 <- hpo_down$ancestors['HP:0003791'][[1]]
    anc_int <- intersect(anc_1, anc_2)
    
    desc_term <- hpo_down$children['HP:0000118'][[1]]
    
    resnik_value <- -log(length(get_descendants(hpo_down, 'HP:0000118')) / length(hpo_down$id))
    
    # Calculate Lin 
    
    ic_1 <- -log(length(get_descendants(hpo_down, 'HP:0025013')) / length(hpo_down$id))
    ic_2 <- -log(length(get_descendants(hpo_down, 'HP:0003791')) / length(hpo_down$id))
    
    lin_value <- 2 * (resnik_value / (ic_1 + ic_2))
    
    get_sim_grid(ontology = hpo_down, 
                 term_sets = total_set,
                 term_sim_method = 'lin',
                 combine = 'average')
    
    
    # resnik method = 1.234995
    # HP:0005289
    # HP:0004438
    
    
get_sim_score <- function(genes_vector, patient_terms) {
      
      # genes_vector <- c('KIAA0319L', 'GJB3', 'GJB4')
      # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
      
      tmp_df <- hpo_genes %>% 
        filter(gene %in% genes_vector) %>%
        select(gene, hp)
      
      to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
      
      genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
      
      hpo_patient <- patient_terms
      
      names(hpo_patient) <- c('patient')  
      total_set <- c(genes_sample, hpo_patient)
      
      
      
      mat_test <- get_sim_grid(ontology = hpo_down, 
                               term_sets = total_set,
                               # term_sets2 = hpo_patient,
                               # combine = 'average',
                               term_sim_method = 'resnik',
                               combine = 'average') %>%
        as_tibble(rownames = 'gene') %>%
        filter(gene != 'patient') %>%
        select(gene, patient) %>%
        mutate(p_value = NA) %>%
        left_join(to_p_value)
      
        mat_test$p_value <-  c(apply(mat_test, 1, get_p))
      
      return(mat_test)
      
      # get_sim_p(mat_test, group = ncol(mat_test))
    }
    
    
    get_p <- function(.patient_terms, .real_value, .n_hpo_gene) {
      
      # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
      # names(patient_terms) <- 'patient1'
      # n_hpo_gene <- 4
      fake_gene <- replicate(simplify=FALSE, n= 100, expr=minimal_set(hpo_down, sample(hpo_down$id, size= n_hpo_gene)))
      
      
      total_set <- c(fake_gene, hpo_patient)
      
      mat_test <- get_sim_grid(ontology = hpo_down, 
                               term_sets = total_set,
                               term_sim_method = 'resnik',
                               combine = 'average')
      
      p_value <- mat_test %>% 
        as_tibble() %>%
        select(patient) %>%
        filter(patient >= real_value) %>%
        nrow()
      
      p_value <- p_value / 100
      
      return(p_value)
      
    }
    
    
    
    
    
    
    
    hsGO <- godata('org.Hs.eg.db', ont="MF")
    goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")
    
    GOSemSim:::termSim('GO:0004022', 'GO:0005515', semData = hsGO, method = 'Wang')
    GOSemSim:::infoContentMethod('GO:0004022', 'GO:0005515', method = 'Wang', hsGO)
    GOSemSim:::infoContentMethod_cpp('GO:0004022', 'GO:0005515', '' )
    
    
    ontologySimilarity::get_term_sim_mat(ontology = hpo_down, 
                                         information_content = c(0)
                                         method = 'resnik')
    st$p_value <-  c(apply(mat_test, 1, get_p))

  return(mat_test)
  
  # get_sim_p(mat_test, group = ncol(mat_test))
}


get_p <- function(.patient_terms, .real_value, .n_hpo_gene) {
  
  # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
  # names(patient_terms) <- 'patient1'
  # n_hpo_gene <- 4
  fake_gene <- replicate(simplify=FALSE, n= 100, expr=minimal_set(hpo_down, sample(hpo_down$id, size= n_hpo_gene)))
  
  
  total_set <- c(fake_gene, hpo_patient)
  
  mat_test <- get_sim_grid(ontology = hpo_down, 
                           term_sets = total_set,
                           term_sim_method = 'resnik',
                           combine = 'average')
  
  p_value <- mat_test %>% 
    as_tibble() %>%
    select(patient) %>%
    filter(patient >= real_value) %>%
    nrow()
  
  p_value <- p_value / 100
  
  return(p_value)

}







hsGO <- godata('org.Hs.eg.db', ont="MF")
goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")

GOSemSim:::termSim('GO:0004022', 'GO:0005515', semData = hsGO, method = 'Wang')
GOSemSim:::infoContentMethod('GO:0004022', 'GO:0005515', method = 'Wang', hsGO)
GOSemSim:::infoContentMethod_cpp('GO:0004022', 'GO:0005515', '' )


ontologySimilarity::get_term_sim_mat(ontology = hpo_down, 
                                     information_content = c(0)
                                     method = 'resnik')
