
# library(GenomicScores)
select <- dplyr::select

# ------------------------------------------------------------------------------
# Name function: check_tads
# Description: Check if a specific interval (chrom, start, end) overlap with any TADs limit
# in that case, the function retrieves a vector with the TADs disrupted.
# ------------------------------------------------------------------------------


check_tads <- function(chrom = NULL, start = NULL, end = NULL, tad_object = NULL) {
  
  
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



get_perc_overlap <- function(df, start_cnv, end_cnv) {
  
  # start_cnv <- 1000
  # end_cnv <- 100000
  # 
  # df <- test1
  
  start_cnv <- as.numeric(start_cnv)
  end_cnv <- as.numeric(end_cnv)
  
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
# Name: plot_upset
# Description: Plot a upset graphic which depicts overlapping of phenotype terms in genes
# Input: dataframe
# Output: plot
# ------------------------------------------------------------------------------

get_upset <- function(df) {
  

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


