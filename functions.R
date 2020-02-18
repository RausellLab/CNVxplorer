
# library(GenomicScores)
select <- dplyr::select

# ------------------------------------------------------------------------------
# Name function: check_tads
# Description: Check if a specific interval (chrom, start, end) overlap with any TADs limit
# in that case, the function retrieves a vector with the TADs disrupted.
# ------------------------------------------------------------------------------


check_tads <- function(chrom = NULL, start = NULL, end = NULL, tad_object = NULL) {
  # 
  # chrom_tmp <- test655
  # start_tmp <- test555
  # end_tmp <- test666
  # 
  chrom_tmp <-  chrom
  start_tmp <-  start
  end_tmp <-  end
  
  
  
  
  tmp_df <- tad_object %>%
    filter(chrom == !!chrom_tmp) %>%
    mutate(check_tad = if_else(start_tmp <= start & end_tmp >= start, 1, 0 )) %>%
    mutate(check_tad2 = if_else(start_tmp <= end & end_tmp >= end, 1, 0 )) %>%
    mutate(check_tad3 = if_else(start_tmp <= start & end_tmp >= end, 1, 0 )) %>%
    mutate(check_final_tad = if_else(check_tad  + check_tad2 + check_tad3 > 0, 1, 0 )) %>%
    filter(check_final_tad > 0) %>%
    select(-contains('check'))

  
  
  if (nrow(tmp_df) > 0) {
    
    df_tads <- tmp_df
    
    
  } else {
    
    df_tads <- 0
    
  }
  
  return(df_tads)
  
  
}


# ------------------------------------------------------------------------------
# Name function: check_regions
# Description: Check if a specific interval (chrom, start, end) overlap with any lncRNA limit
# in that case, the function retrieves a vector with the id lncRNA that overlap the region.
# ------------------------------------------------------------------------------


check_regions <- function(chrom = NULL, start = NULL, end = NULL) {
  
  
  chrom_tmp <-  chrom
  start_tmp <-  start
  end_tmp <-  end
  
  
  tmp_df <- lncrna_coord %>%
    filter(chrom == !!chrom_tmp) %>%
    mutate(check_tad = if_else(start_tmp <= start & end_tmp >= start, 1, 0 )) %>%
    mutate(check_tad2 = if_else(start_tmp <= end & end_tmp >= end, 1, 0 )) %>%
    mutate(check_tad3 = if_else(start_tmp <= start & end_tmp >= end, 1, 0 )) %>%
    mutate(check_final_tad = if_else(check_tad  + check_tad2 + check_tad3 > 0, 1, 0 )) %>%
    select(id, check_final_tad) %>%
    filter(check_final_tad > 0) %>%
    distinct()
  
  if (nrow(tmp_df) > 0) {
    
    list_tads <- as.vector(tmp_df %>% pull(id))
    
    
  } else {
    
    list_tads <- 0
    
  }
  
  return(list_tads)
  
  
}


# ------------------------------------------------------------------------------
# Name: get_perc_overlap
# Description: Calculate the percentage of overlapping between regions(chrom:start-end) and the CNV interval
# Input: dataframe
# Output: dataframe - p_overlap  (column)
# ------------------------------------------------------------------------------



get_perc_overlap <- function(df, chrom_cnv, start_cnv, end_cnv) {
  
  # start_cnv <- 34813719
  # end_cnv <- 36278623
  # 

  
  tmp_df <- df %>% dplyr::select(chrom, start, end)
  
  df <- bed_intersect(tibble(chrom = chrom_cnv, start = start_cnv, end = end_cnv ), tmp_df ) %>%
    mutate(p_overlap = (.overlap /(end.x - start.x + 1))*100) %>% arrange(p_overlap) %>%
    mutate(p_overlap = round(p_overlap, 2)) %>%
    select(start.y, end.y, p_overlap) %>%
    right_join(df, by = c('start.y' = 'start', 'end.y' = 'end')) %>%
    rename(start = start.y, end = end.y) %>%
    select(-p_overlap, p_overlap)
  
  # df <- df %>%
  #   rename(start_gene = start_position, end_gene = end_position) %>%
  #   mutate(type_overlap = case_when(
  #     start_cnv <  start_gene & end_cnv > start_gene ~ "left",
  #     start_cnv > start_gene & end_cnv < end_gene  ~ "center",
  #     start_cnv <  end_gene & end_cnv > end_gene ~ "right",
  #     start_cnv < start_gene & end_cnv > end_gene ~ 'all'
  #   )) %>%
  #   mutate(p_overlap = case_when(
  #     type_overlap ==  'center' ~ 1,
  #     type_overlap ==  'all' ~ (end_gene - start_gene + 1) / (end_cnv - start_cnv + 1) ,
  #     type_overlap ==  'right'  ~ (end_gene - start_cnv + 1) / (end_cnv - start_cnv + 1) ,
  #     type_overlap ==  'left'  ~ (end_cnv - start_gene + 1) / (end_cnv - start_cnv + 1)
  #   )) %>%
  #   mutate(p_overlap = round(p_overlap * 100, 0)) %>%
  #   rename(start_position = start_gene, end_position = end_gene)
  # select(-type_overlap)
  
  

  
  return(df)
}



# ------------------------------------------------------------------------------
# Name: plot_upset
# Description: Plot a upset graphic which depicts overlapping of phenotype terms in genes
# Input: dataframe
# Output: plot
# ------------------------------------------------------------------------------

get_upset <- function(df, gene = FALSE) {
  

  # df <- test41114
  
  
  df_tmp <- df %>% 
    dplyr::select(term) %>%
    distinct() %>%
    mutate(id_row = row_number())
  
  
  vector_hpo <- df %>% select(term) %>% distinct() %>% pull()
  
  
  validate(
    need(length(vector_hpo) > 1, "No intersection found")
  )
  
  vector_genes <- df %>% select(gene) %>% distinct() %>% pull()
  
  list_result <- replicate(length(vector_hpo), NA, simplify = FALSE)
  
  
  for (i in 1:length(vector_hpo)) {
    
    list_result[[i]] <- df %>% filter(term == !!vector_hpo[i]) %>% pull(gene) %in% vector_genes %>% which()
    names(list_result)[i] <- vector_hpo[i]
    
  }
  
 name_axis <-  if_else(gene,'Disease database', 'Phenotype term')
  
        upset(fromList(list_result), order.by = "freq" ,
        point.size = 3.5, line.size = 2, number.angles = 0,
        mainbar.y.label = paste(name_axis, "Intersections"), sets.x.label = paste('Genes Associated Per', name_axis),
        text.scale = c(1.3, 1.3, 1, 1, 2, 2))
  
  
}


# ------------------------------------------------------------------------------
# Calculate distance between set of genes and patient's HPO terms
# ------------------------------------------------------------------------------

get_sim_score <- function(genes_vector, patient_terms, hpo_list_genes, hpo_dbs) {
  
  
  # hpo_list_genes <- hpo_genes
  # hpo_dbs <- hpo_down
  # genes_vector <- c('KIAA0319L', 'GJB3', 'GJB4')
  # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
  
  # genes_vector <- test001
  # patient_terms <- test002
  # hpo_list_genes <- test003
  # hpo_dbs <- test004

  # patient_terms[[1]] <- c(patient_terms[[1]] , 'HP:0001166' )
  
  
  hpo_patient <- patient_terms[[1]]
  # names(hpo_patient) <- c('value_patient')  
  
  
  
  tmp_df <- hpo_list_genes %>% 
    filter(gene %in% genes_vector) %>%
    select(gene, hp)
  
  to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
  
  genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
  
  

  mat_test <- get_profile_sims(ontology = hpo_dbs, profile = hpo_patient, term_sets = genes_sample,
                   term_sim_method = 'resnik') %>%
    as_tibble(rownames = 'gene') %>%
    mutate(p_value = NA) %>%
    # mutate(patient_terms = hpo_patient) %>%
    left_join(to_p_value, by = "gene")
  

  # for (i in 1:nrow(mat_test)) {
  #   
  #   mat_test$value_patient[i] 
  #   mat_test$patient_terms[i] 
  #   mat_test$n_freq[i] 
  #   
  #   mat_test$p_value[i] <- get_p(mat_test$patient_terms[i], 
  #                                mat_test$value_patient[i],
  #                                mat_test$n_freq[i] 
  #                                )
  #   
  # }

  return(mat_test)
  
  # get_sim_p(mat_test, group = ncol(mat_test))
}



# ------------------------------------------------------------------------------
# Calculate p.value
# ------------------------------------------------------------------------------



get_p <- function(patient_terms, value_patient, n_freq) {
  
  
  # patient_terms <- mat_test$patient_terms[i]
  # value_patient <- mat_test$value_patient[i]
  # n_freq <- mat_test$n_freq[i] 
  # 
  
  # print(n_freq)

  fake_gene <- replicate(simplify=FALSE, n= 100, expr=minimal_set(hpo_down, sample(hpo_down$id, size= n_freq)))
  hpo_patient <-  patient_terms
  names(hpo_patient) <- 'patient'

  total_set <- c(fake_gene, hpo_patient)

  mat_test <- get_sim_grid(ontology = hpo_down,
                           term_sets = total_set,
                           term_sim_method = 'resnik',
                           combine = 'average')

  p_value <- mat_test %>%
    as_tibble() %>%
    select(patient) %>%
    filter(patient >= value_patient) %>%
    nrow()

  p_value <- p_value / 100

  return(p_value)

}


# ------------------------------------------------------------------------------
# GENE SIMILARITY BASED ON GENE ONTOLOGY
# THINGS TO PAY ATTENTION:
# GENES > 1 ONTOLOGY TERMS
# GENES == 0 ONTOLOGY TERMS
# ------------------------------------------------------------------------------


mgeneSim_mod <- function (genes, semData, measure = "Wang", drop = "IEA", combine = "BMA", 
                          verbose = TRUE) 
{
  
  
  # genes <- tmp_genes
  # measure = 'Resnik'
  # drop = NULL
  # semData = mf_go
  
  genes <- unique(as.character(genes))
  n <- length(genes)
  scores <- matrix(NA, nrow = n, ncol = n)
  rownames(scores) <- genes
  colnames(scores) <- genes
  gos <- lapply(genes, gene2GO, godata = semData, dropCodes = NULL)
  uniqueGO <- unique(unlist(gos))
  
  
  if (length(uniqueGO) == 0) {
    return(NA)
  }
  
  
  go_matrix <- mgoSim(uniqueGO, uniqueGO, semData, measure = measure, 
                      combine = NULL)
  if (verbose) {
    cnt <- 1
    pb <- txtProgressBar(min = 0, max = sum(1:n), style = 3)
  }
  for (i in seq_along(genes)) {
    for (j in seq_len(i)) {
      if (verbose) {
        setTxtProgressBar(pb, cnt)
        cnt <- cnt + 1
      }
      scores[i, j] <- combineScores(go_matrix[gos[[i]], 
                                              gos[[j]]], combine = combine)
      scores[j, i] <- scores[i, j]
    }
  }
  if (verbose) 
    close(pb)
  removeRowNA <- apply(!is.na(scores), 1, sum) > 0
  removeColNA <- apply(!is.na(scores), 2, sum) > 0
  return(scores[removeRowNA, removeColNA, drop = FALSE])
}

# ------------------------------------------------------------------------------
# Name: get_score
# Description: Conservation score given a model
# Input: dataframe (bed format)
# Output: input + output_name
# ------------------------------------------------------------------------------
# 
# phastcons100 <- getGScores("phastCons100way.UCSC.hg19")
# phylop100 <- getGScores("phyloP100way.UCSC.hg19")
# phastcons46pla <- getGScores("phastCons46wayPlacental.UCSC.hg19")
# phastcons46primates <- getGScores("phastCons46wayPrimates.UCSC.hg19")
# cadd_1_3 <- getGScores("cadd.v1.3.hg19")
# 
# get_score <- function(df, model_conserv, output_name) {
#   
#   df <- df_enhancers[1:10,]
#   model_conserv <- phastcons46pla
#   output_name <- 'caca'
#   
#   df <- df %>% mutate(remove_later = NA)
#   n_df <- nrow(df)
#   for (i in 1:n_df){
#     print(paste0(i, '/', n_df))
#     tmp_gr <- GRanges(
#       seqnames = paste0('chr',1),
#       ranges = IRanges(1100000:1100010, width=1))
#     a <- gscores(phastcons100, tmp_gr)
# #     
#     df$remove_later[i] <- mean(a$default)
#     
#   }
#   
#   df <- df %>%
#     rename(!!quo_name(output_name) := remove_later) 
#   
#   
#   return(df)
# }
# 


