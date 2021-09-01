
select <- dplyr::select
count <- dplyr::count
rename <- dplyr::rename
slice <- dplyr::slice


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
        mainbar.y.label = 'Number of genes', sets.x.label = 'Number of genes',
        text.scale = c(1.3, 1.3, 2, 2, 2, 2))
  
  
}

