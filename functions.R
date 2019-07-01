
# ------------------------------------------------------------------------------
# Name function: check_tads
# Description: Check if a specific interval (chrom, start, end) overlap with any TADs limit
# in that case, the function retrieves a vector with the TADs disrupted.
# ------------------------------------------------------------------------------


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


check_regions('1', 1000, 1000000)
