
library(GenomicScores)

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


# ------------------------------------------------------------------------------
# Name: get_perc_genes
# Description: Calculate the percentage of overlapping between genes and the CNV interval
# Input: dataframe
# Output: dataframe - p_overlap  (column)
# ------------------------------------------------------------------------------



get_perc_genes <- function(df, start_cnv, end_cnv) {
  
  # start_cnv <- 1000
  # end_cnv <- 100000
  # 
  # df <- test1
  
  df <- df %>%
    rename(start_gene = start_position, end_gene = end_position) %>%
    mutate(type_overlap = case_when(
      start_gene <  start_cnv & end_gene < end_cnv ~ "left",
      start_gene > start_cnv & end_gene < end_cnv  ~ "center",
      start_gene >  start_cnv & end_gene > end_cnv ~ "right"
    )) %>%
    mutate(p_overlap = case_when(
      type_overlap ==  'center' ~ 1,
      type_overlap ==  'right'  ~ (end_cnv - start_gene + 1) / (end_gene - start_gene + 1) ,
      type_overlap ==  'left'  ~ (end_gene - start_gene + 1) / (end_gene - start_gene + 1)
    )) %>%
    mutate(p_overlap = round(p_overlap * 100, 0)) %>%
    rename(start_position = start_gene, end_position = end_gene) %>%
    select(-type_overlap)
  
  return(df)
}



# ------------------------------------------------------------------------------
# Name: get_score
# Description: Conservation score given a model
# Input: dataframe (bed format)
# Output: input + output_name
# ------------------------------------------------------------------------------

phastcons100 <- getGScores("phastCons100way.UCSC.hg19")
phylop100 <- getGScores("phyloP100way.UCSC.hg19")
phastcons46pla <- getGScores("phastCons46wayPlacental.UCSC.hg19")
phastcons46primates <- getGScores("phastCons46wayPrimates.UCSC.hg19")
cadd_1_3 <- getGScores("cadd.v1.3.hg19")

get_score <- function(df, model_conserv, output_name) {
  
  df <- df_enhancers[1:10,]
  model_conserv <- phastcons46pla
  output_name <- 'caca'
  
  df <- df %>% mutate(remove_later = NA)
  n_df <- nrow(df)
  for (i in 1:n_df){
    print(paste0(i, '/', n_df))
    tmp_gr <- GRanges(
      seqnames = paste0('chr',df$chrom[i]),
      ranges = IRanges(df$start[i]:df$end[i], width=1))
    a <- gscores(model_conserv, tmp_gr)
    
    df$remove_later[i] <- mean(a$default)
    
  }
  
  df <- df %>%
    rename(!!quo_name(output_name) := remove_later) 
  
  
  return(df)
}


test1 <- get_score(df_enhancers[1:10,], phastcons46primates, 'phast46pri' )

