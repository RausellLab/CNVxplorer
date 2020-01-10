

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
