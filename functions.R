
select <- dplyr::select
count <- dplyr::count
rename <- dplyr::rename
slice <- dplyr::slice

# ------------------------------------------------------------------------------
# Name function: check_tads
# Description: Check if a specific interval (chrom, start, end) overlap with any TADs limit
# in that case, the function retrieves a vector with the TADs disrupted.
# ------------------------------------------------------------------------------


check_tads <- function(coord_tbl, tad_object = NULL) {

 
  
  select_tads_disrupted <-  tad_object %>%
    pivot_longer(-c(id,chrom), names_to = 'coord', values_to = 'start') %>%
    mutate(end = start) %>%
    select(-coord) %>%
    bed_intersect(coord_tbl) %>%
    pull(id.x)
  
 df_tads <- tad_object %>% filter(id %in% select_tads_disrupted)
  
  if (nrow(df_tads) == 0) df_tads <- 0
  
  return(df_tads)
  
  
}


# ------------------------------------------------------------------------------
# Name function: check_regions
# Description: Check if a specific interval (chrom, start, end) overlap with any lncRNA limit
# in that case, the function retrieves a vector with the id lncRNA that overlap the region.
# ------------------------------------------------------------------------------
# 
# 
# check_regions <- function(chrom = NULL, start = NULL, end = NULL) {
#   
#   
#   chrom_tmp <-  chrom
#   start_tmp <-  start
#   end_tmp <-  end
#   
#   
#   tmp_df <- lncrna_coord %>%
#     filter(chrom == !!chrom_tmp) %>%
#     mutate(check_tad = if_else(start_tmp <= start & end_tmp >= start, 1, 0 )) %>%
#     mutate(check_tad2 = if_else(start_tmp <= end & end_tmp >= end, 1, 0 )) %>%
#     mutate(check_tad3 = if_else(start_tmp <= start & end_tmp >= end, 1, 0 )) %>%
#     mutate(check_final_tad = if_else(check_tad  + check_tad2 + check_tad3 > 0, 1, 0 )) %>%
#     select(id, check_final_tad) %>%
#     filter(check_final_tad > 0) %>%
#     distinct()
#   
#   if (nrow(tmp_df) > 0) {
#     
#     list_tads <- as.vector(tmp_df %>% pull(id))
#     
#     
#   } else {
#     
#     list_tads <- 0
#     
#   }
#   
#   return(list_tads)
#   
#   
# }


# ------------------------------------------------------------------------------
# Name: get_perc_overlap
# Description: Calculate the percentage of overlapping between regions(chrom:start-end) and the CNV interval
# Input: dataframe
# Output: dataframe - p_overlap  (column)
# ------------------------------------------------------------------------------



get_perc_overlap <- function(df, input_tbl, 
                             is_a_gene = FALSE,
                             is_patho = FALSE) {

  
  df <- df %>%
    mutate(id_tmp = row_number())

  
  if (is_a_gene == TRUE) {
    
    tmp_df <- df %>%
      dplyr::select(chrom, start, end, id_tmp) 
    
  df <- bed_intersect(tmp_df, input_tbl) %>%
    group_by(chrom, start.x, end.x) %>%
    filter(.overlap == max(.overlap)) %>%
    distinct() %>%
    ungroup() %>%
    mutate(p_overlap = ((.overlap + 1)  /(end.x - start.x + 1))*100) %>% 
    mutate(p_overlap = round(p_overlap, 2)) %>%
    select(start.x, end.x, p_overlap, id_tmp.x) %>%
    right_join(df, by = c('start.x' = 'start', 'end.x' = 'end',
                           'id_tmp.x' = 'id_tmp')) %>%
    rename(start = start.x, end = end.x) %>%
    select(-p_overlap, p_overlap, -id_tmp.x) %>%
    arrange(desc(p_overlap)) %>%
    distinct()

  } else {
    
    if (isTRUE(is_patho)) {
      
      df_over_tmp <- df %>% 
        dplyr::select(chrom, start, end, id_tmp) %>%
        bed_intersect(input_tbl)
      
      
      df <- df_over_tmp %>%
        group_by(chrom, start.x, end.x) %>%
        filter(.overlap == max(.overlap)) %>%
        ungroup() %>%
        mutate(p_overlap = ((.overlap + 1) /(end.x - start.x + 1))*100) %>% 
        arrange(p_overlap) %>%
        mutate(p_overlap = round(p_overlap, 2)) %>%
        select(start.x, end.x, id_tmp.x, p_overlap) %>%
        right_join(df, by = c('start.x' = 'start', 'end.x' = 'end', 'id_tmp.x' = 'id_tmp')) %>%
        rename(start = start.x, end = end.x) %>%
        select(-p_overlap, p_overlap, -id_tmp.x) %>%
        arrange(desc(p_overlap)) %>%
        distinct()
      
      
      
    } else {
      
      tmp_df <- df %>% 
        dplyr::select(chrom, start, end, id_tmp)
      
      df_over_tmp <-  bed_intersect(input_tbl, tmp_df)
      
      
      df <- df_over_tmp %>%
        group_by(chrom, start.y, end.y) %>%
        filter(.overlap == max(.overlap)) %>%
        ungroup() %>%
        mutate(p_overlap = ((.overlap + 1) /(end.x - start.x + 1))*100) %>% 
        arrange(p_overlap) %>%
        mutate(p_overlap = round(p_overlap, 2)) %>%
        select(start.y, end.y, id_tmp.y, p_overlap) %>%
        right_join(df, by = c('start.y' = 'start', 'end.y' = 'end', 'id_tmp.y' = 'id_tmp')) %>%
        rename(start = start.y, end = end.y) %>%
        select(-p_overlap, p_overlap, -id_tmp.y) %>%
        arrange(desc(p_overlap)) %>%
        distinct()
      
      
    }


  }

  return(df)
}



# ------------------------------------------------------------------------------
# Name: plot_upset
# Description: Plot a upset graphic which depicts overlapping of phenotype terms in genes
# Input: dataframe
# Output: plot
# ------------------------------------------------------------------------------

get_upset <- function(df) {
  

  df_tmp <- df %>% 
    dplyr::select(term) %>%
    distinct() %>%
    mutate(id_row = row_number())
  
  
  vector_hpo <- df %>% dplyr::select(term) %>% distinct() %>% pull()
  
  
  validate(
    need(length(vector_hpo) > 1, "No intersection found.")
  )
  
  vector_genes <- df %>% select(gene) %>% distinct() %>% pull()
  
  list_result <- replicate(length(vector_hpo), NA, simplify = FALSE)
  
  
  for (i in 1:length(vector_hpo)) {
    
    order_genes <- df %>% filter(term == !!vector_hpo[i]) %>% pull(gene)
    list_result[[i]] <-vector_genes %in%  order_genes %>% which()
    names(list_result)[i] <- vector_hpo[i]
    
  }
  

        upset(fromList(list_result), order.by = "freq" ,
        point.size = 3.5, line.size = 2, number.angles = 0,
        mainbar.y.label = 'Gene intersections', sets.x.label = 'Number of genes',
        text.scale = c(1.3, 1.3, 2, 2, 2, 2))
  
  
}


# ------------------------------------------------------------------------------
# Calculate distance between set of genes and patient's HPO terms
# ------------------------------------------------------------------------------

# get_sim_score <- function(genes_vector, patient_terms, hpo_list_genes, hpo_dbs) {
#   
#   
#   # hpo_list_genes <- hpo_genes
#   # hpo_dbs <- hpo_down
#   # genes_vector <- c('KIAA0319L', 'GJB3', 'GJB4')
#   # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
#   
#   genes_vector <- test001
#   patient_terms <- test002
#   hpo_list_genes <- test003
# 
# 
#   
#   hpo_patient <- patient_terms[[1]]
# 
#   tmp_df <- hpo_list_genes %>% 
#     filter(gene %in% genes_vector) %>%
#     select(gene, hp)
#   
#   to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
#   
#   genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
#   
#   test00021 <<- genes_sample
#   test00022 <<- hpo_patient
# 
#   # mat_test <- get_profile_sims(ontology = hpo_dbs, profile = hpo_patient, term_sets = genes_sample,
#   #                  term_sim_method = 'resnik') %>%
#   #   as_tibble(rownames = 'gene') %>%
#   #   mutate(p_value = NA) %>%
#   #   # mutate(patient_terms = hpo_patient) %>%
#   #   left_join(to_p_value, by = "gene")
#   
#   mat_test <- get_sim_grid(ontology=hpo_dbs, 
#                                 term_sets= list('patient' = hpo_patient),
#                                 term_sets2 = genes_sample,
#                                 term_sim_method = 'resnik') %>%
#     as_tibble(rownames = 'gene') %>%
#     mutate(p_value = NA) %>%
#     left_join(to_p_value, by = "gene")
#   
# 
#   # for (i in 1:nrow(mat_test)) {
#   #   
#   #   mat_test$value_patient[i] 
#   #   mat_test$patient_terms[i] 
#   #   mat_test$n_freq[i] 
#   #   
#   #   mat_test$p_value[i] <- get_p(mat_test$patient_terms[i], 
#   #                                mat_test$value_patient[i],
#   #                                mat_test$n_freq[i] 
#   #                                )
#   #   
#   # }
# 
#   return(mat_test)
#   
#   # get_sim_p(mat_test, group = ncol(mat_test))
# }



# ------------------------------------------------------------------------------
# Calculate p.value
# ------------------------------------------------------------------------------


# 
# get_p <- function(patient_terms, value_patient, n_freq) {
#   
#   
#   # patient_terms <- mat_test$patient_terms[i]
#   # value_patient <- mat_test$value_patient[i]
#   # n_freq <- mat_test$n_freq[i] 
#   # 
#   
#   # print(n_freq)
# 
#   fake_gene <- replicate(simplify=FALSE, n= 100, expr=minimal_set(hpo_down, sample(hpo_down$id, size= n_freq)))
#   hpo_patient <-  patient_terms
#   names(hpo_patient) <- 'patient'
# 
#   total_set <- c(fake_gene, hpo_patient)
# 
#   mat_test <- get_sim_grid(ontology = hpo_down,
#                            term_sets = total_set,
#                            term_sim_method = 'resnik',
#                            combine = 'average')
# 
#   p_value <- mat_test %>%
#     as_tibble() %>%
#     select(patient) %>%
#     filter(patient >= value_patient) %>%
#     nrow()
# 
#   p_value <- p_value / 100
# 
#   return(p_value)
# 
# }
