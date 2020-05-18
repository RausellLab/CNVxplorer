
# library(GenomicScores)
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



get_perc_overlap <- function(df, input_tbl, 
                             is_a_gene = FALSE,
                             is_patho = FALSE) {


  
  input_tbl <- test2020
  
  
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

get_sim_score <- function(genes_vector, patient_terms, hpo_list_genes, hpo_dbs) {
  
  
  # hpo_list_genes <- hpo_genes
  # hpo_dbs <- hpo_down
  # genes_vector <- c('KIAA0319L', 'GJB3', 'GJB4')
  # patient_terms <- replicate(simplify=FALSE, n=1, expr=minimal_set(hpo_down, sample(hpo_down$id, size=10)))
  
  genes_vector <- test001
  patient_terms <- test002
  hpo_list_genes <- test003


  
  hpo_patient <- patient_terms[[1]]

  tmp_df <- hpo_list_genes %>% 
    filter(gene %in% genes_vector) %>%
    select(gene, hp)
  
  to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
  
  genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
  
  test00021 <<- genes_sample
  test00022 <<- hpo_patient

  # mat_test <- get_profile_sims(ontology = hpo_dbs, profile = hpo_patient, term_sets = genes_sample,
  #                  term_sim_method = 'resnik') %>%
  #   as_tibble(rownames = 'gene') %>%
  #   mutate(p_value = NA) %>%
  #   # mutate(patient_terms = hpo_patient) %>%
  #   left_join(to_p_value, by = "gene")
  
  mat_test <- get_sim_grid(ontology=hpo_dbs, 
                                term_sets= list('patient' = hpo_patient),
                                term_sets2 = genes_sample,
                                term_sim_method = 'resnik') %>%
    as_tibble(rownames = 'gene') %>%
    mutate(p_value = NA) %>%
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



# ------------------------------------------------------------------------------
# Name: check_app_cnv
# Description: CNV annotation -> arules model
# ------------------------------------------------------------------------------


check_app_cnv <- function(input_id, input_clinical, input_variant, input_inheritance,
                      input_chrom, input_start, input_end) {
  
  
  
  id_tmp <- input_id
  clinical_tmp <- input_clinical
  type_variant_tmp <- input_variant
  type_inheritance_tmp <- input_inheritance
  chrom_tmp <- input_chrom
  start_tmp <- input_start
  end_tmp <- input_end
  length_tmp <- end_tmp - start_tmp + 1
  threshold_30_tmp <- round((length_tmp / 100)*30,0)
  
  # just_test <- input_check_cnv %>% filter(id == '1621')
  # id_tmp <- just_test$id
  # clinical_tmp <-  just_test$pathogenicity
  # type_variant_tmp <-  just_test$variant_class
  # type_inheritance_tmp <-  just_test$inheritance
  # chrom_tmp <-  just_test$chrom
  # start_tmp <-  just_test$start
  # end_tmp <-  just_test$end
  # length_tmp <- end_tmp - start_tmp + 1
  # threshold_30_tmp <- round((length_tmp / 100)*30,0)
  
  tmp_cnv <- tibble(chrom = chrom_tmp, start = start_tmp, end = end_tmp)
  
  # 0. Basic information
  # # N.N Number of genes
  # # N.N Length
  # 1. CNV databases
  # # 1.1 OVERLAP CNV SYNDROMES
  # # 1.2 OVERLAP ANNOTATED PATHOGENIC
  # # 1.3 OVERLAP NON-PATHOGENIC CNVs
  # # 1.4 BLACKLIST REGION
  # 2. Disease annotation
  # # 2.1 OVERLAP DISEASE GENES
  # # 2.2 OVERLAP DISEASE VARIANTS
  # 3. Model mouse information
  # # 3.1 RULE - OVERLAP WITH ORTHOLOGS GENES ASSOCIATED WITH LETHALITY
  # 4. Biological categories
  ## 4.1 RULE - SIGNIFICATIVE FUNCTIONAL SIMILARITY
  ## 4.2 RULE - SIGNIFICATIVE HITS - GENE ONTOLOGY
  ## 4.3 RULE - SIGNIFICATIVE HITS - PATHWAYS
  # 
  
  vector_genes <- hgcn_genes %>%
    rename(start = start_position, end = end_position) %>%
    bed_intersect(tmp_cnv ) %>%
    mutate(length_gene = end.x - start.x + 1) %>%
    mutate(perc_overlap = .overlap / length_gene) %>%
    filter(perc_overlap >= 0.3) %>% pull(gene.x)
  
  vector_entrez <- hgcn_genes %>%
    rename(start = start_position, end = end_position) %>%
    bed_intersect(tmp_cnv ) %>%
    mutate(length_gene = end.x - start.x + 1) %>%
    mutate(perc_overlap = .overlap / length_gene) %>%
    filter(perc_overlap >= 0.3) %>% 
    pull(entrez_id.x)
  
  
  
  
  # # NUMBER OF GENES
  
  result_n_genes <- length(vector_genes)
  
  # # 1º RULE - OVERLAP CNV SYNDROMES
  
  tmp_df <- syndromes_total %>%
    filter(chrom == chrom_tmp)
  
  
  result_n_overlap_cnv_syndrome <- bed_intersect(tmp_df, tmp_cnv) %>%
    filter(.overlap > threshold_30_tmp) %>%
    nrow()
  
  
  
  # 2º RULE - OVERLAP ANNOTATED PATHOGENIC
  tmp_df <- cnv_df %>%
    filter(! id %in% id_tmp) %>%
    filter(pathogenicity == 'Pathogenic') %>%
    filter(chrom == chrom_tmp)
  
  result_n_overlap_patho <- bed_intersect(tmp_df, tmp_cnv) %>%
    filter(.overlap > threshold_30_tmp) %>%
    nrow()
  
  
  # 3º RULE - OVERLAP NON-PATHOGENIC CNVs
  tmp_df <- cnv_df %>% filter(! id %in% id_tmp) %>%
    filter(source != 'decipher') %>%
    filter(chrom == chrom_tmp)
  
  result_n_overlap_nonpatho <- bed_intersect(tmp_df, tmp_cnv) %>%
    filter(.overlap > threshold_30_tmp) %>%
    nrow()
  
  
  # Nº RULE - OVERLAP BLACKLIST REGIONS
  
  
  result_n_blacklist <- bed_intersect(tmp_cnv, blacklist_encode) %>%
    filter(.overlap > threshold_30_tmp) %>%
    nrow()
  
  
  
  # 4º RULE - OVERLAP DISEASE GENES
  
  result_n_disease_genes <- hgcn_genes %>%
    filter(gene %in% vector_genes) %>%
    filter(disease == 'Yes') %>%
    filter(chrom == chrom_tmp) %>%
    rename(start = start_position, end = end_position) %>%
    nrow()
  
  
  # 5º RULE - OVERLAP DISEASE VARIANTS
  
  tmp_df <- clinvar_variants %>%
    filter(clinical_sign == "Pathogenic") %>%
    filter(chrom == chrom_tmp) %>%
    mutate(start = pos, end = pos) %>%
    select(chrom, start, end) %>%
    mutate(chrom = as.character(chrom))
  
  
  
  n_clinvar <- bed_intersect(tmp_cnv, tmp_df) %>% nrow()
  
  
  tmp_df <- gwas_variants %>%
    filter(INTERGENIC == "No") %>%
    filter(CHR_ID == chrom_tmp) %>%
    mutate(start = CHR_POS, end = CHR_POS) %>%
    rename(chrom = CHR_ID) %>%
    select(chrom, start, end)
  
  
  n_gwas <- bed_intersect(tmp_cnv, tmp_df) %>% nrow()
  
  
  result_n_variants_genes <- n_clinvar + n_gwas
  
  
  # # 3.1 RULE - OVERLAP WITH ORTHOLOGS GENES ASSOCIATED WITH EMBRYONIC 
  # PHENOTYPE IN MOUSE MODEL
  
  # MP:0010768 - mortality/aging
  # MP:0005380 - embryo
  
  result_n_mouse_embryo <- mgi %>% filter(gene %in% vector_genes, str_detect(pheno, 'MP:0010768')) %>% nrow()
  
  # # MAXIMUM pLI score
  temporal_df <- hgcn_genes %>% filter(gene %in% vector_genes) %>%  
    select(pLI) %>%
    na.omit() %>%
    arrange(desc(pLI)) %>%
    slice(1) %>% 
    pull(pLI)
  
  maximum_pli <- ifelse(length(temporal_df) == 0, 0, temporal_df)

  

  
  
  # Number of genes associated with at least 1 HPO term
  
  
  n_genes_one_hpo <- vector_genes %in% hpo_genes$gene %>% sum()

  
  # Number ohnologs genes
  
  
  result_n_ohno <- vector_genes %in% ohno_genes$gene %>% sum()

  
  # TFs
  
  result_n_tf <- length(vector_genes[vector_genes %in% tf_genes])
  
  
  # Number of genes associated with target genes of approved drugs
  
  
  result_n_target_drugs <- vector_genes %in% drugbank$gene %>% sum()
  
  
  # Nº RULE - OVERLAP BLACKLIST REGIONS
  
  
  result_n_blacklist <- bed_intersect(tmp_cnv, blacklist_encode) %>%
    filter(.overlap > threshold_30_tmp) %>%
    nrow()

  # AGGREGATION
  
  result_tmp <- tibble(
    'id' = id_tmp,
    'clinical' = clinical_tmp,
    'type_variant' = type_variant_tmp, #
    'type_inheritance' =   type_inheritance_tmp, #
    'n_cnv_syndromes' = result_n_overlap_cnv_syndrome, #
    'patho_cnv' = result_n_overlap_patho, #
    'nonpatho_cnv' = result_n_overlap_nonpatho, #
    'disease_genes' = result_n_disease_genes, #
    'disease_variants' = result_n_variants_genes, #
    'embryo_mouse' = result_n_mouse_embryo, #
    'max_pli' = maximum_pli, #
    'n_genes_hpo' = n_genes_one_hpo, #
    'n_genes' = result_n_genes, #
    # 'n_ohno' = result_n_ohno,
    'n_tf' = result_n_tf, #
    'n_target_drugs' = result_n_target_drugs, #
    'n_blacklist' = result_n_blacklist, #

  )
  
  return(result_tmp)
  
}

