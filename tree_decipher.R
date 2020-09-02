# library(enrichR)
# library(DOSE)
# library(ReactomePA)
# library(clusterProfiler)
# library(UpSetR)

library(valr)
library(future)
library(grid)
library(tidyverse)
library(tictoc)
library(ontologyIndex)
library(corrr)
library(furrr)
library(tune)
library(parsnip)
library(yardstick)
library(glue)
library(patchwork)
library(rsample)
library(xrf)
library(rules)
library(dials)
library(workflows)
library(rstanarm)
library(bayestestR)


set.seed(123)


rename <- dplyr::rename
slice <- dplyr::slice


# save(hgcn_genes, df_enhancers, lncrna_coord, lncrna, tad, gtex, hpa, hpo_genes, cnv_df, vector_total_terms,
#      gnomad_sv_raw, decipher_control_raw, dgv_df_raw, hpo_omim, anato_df, mirtarbase, pubmed_df,
#       vector_inheritance, trrust, tf_genes, drugbank,prot_complex,ohno_genes,recomb,
#       genes_promoter,para_genes,string_db,region_gaps,ensembl_reg, gene_density_tbl,
#     select, dev_raw, panel_total, omim, orphanet_raw,  hpo_dbs, model1, denovo, clinvar_variants, ridges_home, plot_p100, plot_p46pla, blacklist_encode, mpo_dbs, gwas_variants,mgi, syndromes_total,
# file = "env_annot_cnvs.RData")
# 
load('env_annot_cnvs.RData')
load('local_data.RData')


theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size=22))
  
}

input_check_cnv <- read_tsv('/data-cbl/frequena_data/from_workstation/daa_decipher/decipher-cnvs-grch37-2020-01-19.txt', skip = 1) %>%
  as_tibble() %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  select(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  filter(pathogenicity %in% c('Unknown')) %>%
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>% 
  # filter(inheritance == 'De novo constitutive') %>% 
  filter(variant_class %in% c('Deletion', 'Duplication')) %>%
  mutate(id_tmp = row_number())


coord_cytobands <- coord_cytobands %>% rename(chrom = Chrom, start = Start, end = End)

gen_random <- function(input_check_cnv, n_rep = 9) {
  
  result_df <- tibble()
  
  for (i in 1:nrow(input_check_cnv)) {
    
    print(i)
    
    tmp_chrom <- input_check_cnv %>% slice(i) %>% pull(chrom)
    tmp_start <- input_check_cnv %>% slice(i) %>% pull(start)
    tmp_end <- input_check_cnv %>% slice(i) %>% pull(end)
    tmp_id <- input_check_cnv %>% slice(i) %>% pull(id_tmp)
    tmp_cyto <- bed_intersect(coord_cytobands,
                            tibble('chrom' = tmp_chrom,'start' = tmp_start, 'end' = tmp_end)) %>% 
      pull(Name.x)
    
    # tmp_cyto %in% (coord_cytobands %>% filter(chrom == tmp_chrom) %>% pull(Name))
    
    tmp_n_genes <- hgcn_genes %>% 
      bed_intersect(tibble('chrom' = tmp_chrom,'start' = tmp_start, 'end' = tmp_end)) %>% nrow()
  
    
    result_tmp_cnv <- tibble()
    
    n_try <- 0

    while (nrow(result_tmp_cnv) < n_rep) {
      
      n_try <- n_try + 1
      
      if (n_try < 500) {
        
      selected_chrom <- coord_chrom_hg19 %>% filter(chrom == tmp_chrom)
      random_chrom <- selected_chrom  %>% pull(chrom)
      random_chrom_end <- selected_chrom  %>% pull(length)
      random_chrom_end <- random_chrom_end - (tmp_end - tmp_start + 1)
      random_start <- sample(1:random_chrom_end, 1)
      random_end <- random_start + (tmp_end - tmp_start + 1)
      random_cyto <- bed_intersect(coord_cytobands,
                                   tibble('chrom' = random_chrom,
                                          'start' = random_start, 
                                          'end' = random_end)) %>%
      pull(Name.x)
      
      # Filter 1: eliminate cnvs mapping the original one
      tmp_check <- random_cyto %in% tmp_cyto %>% unique()
      
      if (length(tmp_check) == 2 | (length(tmp_check) == 1 & isTRUE(tmp_check))) next
      
      random_genes <- hgcn_genes %>% bed_intersect(tibble('chrom' = random_chrom,
                            'start' = random_start, 
                            'end' = random_end)) %>% nrow()
      
      # Filter 2: eliminate cnvs with no genes if the original cnv has genes
      if (tmp_n_genes > 0) {
        if (random_genes == 0 ) next
      }
      
      
      tmp_df <- tibble('chrom' = random_chrom,
                       'start' = random_start,
                       'end' = random_end)
      
      result_tmp_cnv <- result_tmp_cnv %>% bind_rows(tmp_df)
      } else {
        
        result_tmp_cnv <- tibble('chrom' = 'error', 'start' = 0, 'end' = 0)
        break
      }
      
      }
    
    result_tmp_cnv <- result_tmp_cnv %>% 
      mutate(id_tmp = row_number()) %>%
      mutate(id_tmp = paste0(tmp_id,'_', id_tmp))
    
    result_df <- result_df %>% bind_rows(result_tmp_cnv)
  }
  
  return(result_df)
  
 }

random_cnvs <- gen_random(input_check_cnv, n_rep = 9)
remove_cnvs_error <- random_cnvs %>% 
  filter(chrom == 'error') %>% 
  pull(id_tmp) %>%
  as.character()

input_check_cnv <- input_check_cnv %>% filter(! id_tmp %in% remove_cnvs_error)
random_cnvs <- random_cnvs %>% 
  filter(chrom != 'error') %>%
  mutate(inheritance = NA, variant_class = NA )

input_check_cnv <- input_check_cnv %>% mutate(id_tmp = as.character(id_tmp)) %>% bind_rows(random_cnvs)

check_cnv <- function(input_id, input_clinical, input_variant, input_inheritance,
                      input_chrom, input_start, input_end) {
  
  
  # input_id <- input_check_cnv %>% slice(1) %>% pull(id_tmp)
  # input_clinical <- input_check_cnv %>% slice(1) %>% pull(pathogenicity)
  # input_variant <- input_check_cnv %>% slice(1) %>% pull(variant_class)
  # input_inheritance <- input_check_cnv %>% slice(1) %>% pull(inheritance)
  # input_chrom <- input_check_cnv %>% slice(1) %>% pull(chrom)
  # input_start <- input_check_cnv %>% slice(1) %>% pull(start)
  # input_end <- input_check_cnv %>% slice(1) %>% pull(end)
  
  
  
  
  

  id_tmp <- input_id
  clinical_tmp <- input_clinical
  type_variant_tmp <- input_variant
  type_inheritance_tmp <- input_inheritance
  chrom_tmp <- input_chrom
  start_tmp <- input_start
  end_tmp <- input_end
  length_tmp <- end_tmp - start_tmp + 1
  # threshold_30_tmp <- round((length_tmp / 100)*30,0)
  threshold_30_tmp <- 0

  # just_test <- input_check_cnv %>% filter(id == '131') %>% slice(1)
  # id_tmp <- just_test$id
  # clinical_tmp <-  just_test$pathogenicity
  # type_variant_tmp <-  just_test$variant_class
  # type_inheritance_tmp <-  just_test$inheritance
  # chrom_tmp <-  just_test$chrom
  # start_tmp <-  just_test$start
  # end_tmp <-  just_test$end
  # length_tmp <- end_tmp - start_tmp + 1
  # threshold_30_tmp  <- 0

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
    bed_intersect(tmp_cnv ) %>%
    mutate(length_gene = end.x - start.x + 1) %>% 
    # mutate(perc_overlap = .overlap / length_gene)
    # filter(perc_overlap >= 0.3) %>%
    pull(gene.x)
  
  vector_entrez <- hgcn_genes %>%
    bed_intersect(tmp_cnv) %>%
    mutate(length_gene = end.x - start.x + 1) %>%
    # mutate(perc_overlap = .overlap / length_gene) %>%
    # filter(perc_overlap >= 0.3) %>% 
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


# # MAXIMUM CCR score

temporal_df <- hgcn_genes %>% filter(gene %in% vector_genes) %>%  
  select(ccr) %>%
  na.omit() %>%
  arrange(desc(ccr)) %>%
  slice(1) %>% 
  pull(ccr)

maximum_ccr <- ifelse(length(temporal_df) == 0, 0, temporal_df)



# # MAXIMUM HI score

temporal_df <- hgcn_genes %>% filter(gene %in% vector_genes) %>%  
  select(hi) %>%
  na.omit() %>%
  arrange(desc(hi)) %>%
  slice(1) %>% 
  pull(hi)

maximum_hi <- ifelse(length(temporal_df) == 0, 0, temporal_df)


# # Min Vg score

temporal_df <- hgcn_genes %>% 
  mutate(p_vg = ntile(vg, 100)) %>% 
  filter(gene %in% vector_genes) %>%  
  select(p_vg) %>%
  na.omit() %>%
  arrange(p_vg) %>%
  slice(1) %>% 
  pull(p_vg)

minimum_vg <- ifelse(length(temporal_df) == 0, 0, temporal_df)


# Number of genes associated with at least 1 HPO term


n_genes_one_hpo <- vector_genes %in% hpo_genes$gene %>% sum()


# Protein complex gene


result_n_prot_complex <- vector_genes %in% prot_complex$gene %>% sum()


# Number ohnologs genes


result_n_ohno <- vector_genes %in% ohno_genes$gene %>% sum()


# Number of genes associated with target genes of approved drugs


result_n_target_drugs <- vector_genes %in% drugbank$gene %>% sum()


# Number of physiological systems affected

hpo_from_gene <- hpo_genes %>% filter(gene %in% vector_genes) %>% pull(hp)

result_n_systems <- unlist(map(hpo_from_gene, function(x) get_ancestors(hpo_dbs, x))) %>% 
  enframe() %>%
  filter(value %in% anato_df$name) %>%
  count(value) %>% 
  nrow()

# Maximum recombination rate

result_recomb = valr::bed_closest(tmp_cnv, recomb) %>% pull(cm_mb.y) %>% max()


# Distance (Mb) closest centromeric and telomeric region



dist_cent <- bed_closest(tmp_cnv, region_gaps %>% 
              filter(type == 'centromere')) %>%
  pull(.dist) %>% abs()


dist_tel <- bed_closest(tmp_cnv, region_gaps %>% 
              filter(type == 'telomere')) %>%
  pull(.dist) %>% abs() %>% min()

dist_cent <- dist_cent / 10**6
dist_tel <- dist_tel / 10**6

if (dist_tel == Inf) dist_tel <- NA

# Gene density


result_gene_density <- gene_density_tbl %>% bed_intersect(tmp_cnv) %>% filter(is_end.x != 'yes') %>%
  pull(gene_density.x) %>% max()

if (result_gene_density == -Inf) result_gene_density <- NA


# Pubmed 

result_max_del <- bed_intersect(tmp_cnv, pubmed_df) %>% pull(hits_del.y)
if (length(result_max_del) == 0) {
  result_max_del <- NA
} else {
  
  result_max_del <- result_max_del %>% max()
  
}

result_max_dup <- bed_intersect(tmp_cnv, pubmed_df) %>% pull(hits_dup.y)
if (length(result_max_dup) == 0) {
  result_max_dup <- NA
} else {
  
  result_max_dup <- result_max_dup %>% max()
  
}

# if (websiteLive) {
# 
# 
# if (length(vector_genes) == 0) {
# 
#   result_n_hits_go_bp <- 0
# 
# } else  {
# 
#   Sys.sleep(sample(1:2,1))
# 
#   enriched <- enrichr(vector_genes, dbs)
# 
#   Sys.sleep(sample(1:2,1))
#   
# 
#   enriched_n_row <- enriched$GO_Biological_Process_2018 %>% as_tibble() %>% nrow()
# 
#     if (enriched_n_row == 0) {
# 
#     result_n_hits_go_bp <- 0
#   } else {
# 
#   result_n_hits_go_bp <- enriched$GO_Biological_Process_2018 %>%
#     as_tibble() %>%
#     filter(P.value <= 0.05) %>%
#     rowwise() %>%
#     mutate(n_genes = as.integer(str_split(Overlap, pattern = '/')[[1]][1])) %>%
#     ungroup() %>%
#     filter(n_genes > 1) %>%
#     nrow()
# 
# }
# }
#   
# } else {
#   
#   result_n_hits_go_bp <- 9999
#   
# }

# result_n_hits_go_bp <- 0

# if (length(vector_entrez) == 0) {
#   
#   result_n_hits_go_bp <- 0
# } else {
#   Sys.sleep(sample(1:3,1))
# result_n_hits_go_bp <- enrichGO(gene = vector_entrez,
#                 # universe      = names(geneList),
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "BP",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE) %>%
#       
#       as_tibble() %>% 
#       filter(Count > 1) %>%
#       nrow()
# 
# }


# Maximum PhastCons46 score promoter region

result_max_phast <- genes_promoter %>% filter(gene %in% vector_genes) %>% arrange(desc(mean_phast)) %>%
  slice(1) %>% pull(mean_phast)


if (length(result_max_phast) == 0) result_max_phast <- 0

# Maximum CpG density

result_max_cpg <- genes_promoter %>% filter(gene %in% vector_genes) %>% arrange(desc(cpg_density)) %>%
  slice(1) %>% pull(cpg_density)


if (length(result_max_cpg) == 0) result_max_cpg <- 0


# Maximum degree
result_max_degree <- string_db %>% filter(gene %in% vector_genes) %>% arrange(desc(degree)) %>%
  slice(1) %>% pull(degree)

if (length(result_max_degree) == 0) result_max_degree <- 0


# Maximum page_rank

result_max_page_rank <- string_db %>% filter(gene %in% vector_genes) %>% arrange(desc(page_rank)) %>%
  slice(1) %>% pull(page_rank)

if (length(result_max_page_rank) == 0) result_max_page_rank <- 0


# Paralogous genes

result_max_par <- para_genes %>% filter(gene %in% vector_genes) %>% arrange(desc(n)) %>%
  slice(1) %>% pull(n)

if (length(result_max_par) == 0) result_max_par <- 0

# TFs

result_n_tf <- length(vector_genes[vector_genes %in% tf_genes])


# TFs - disease genes

result_tfs_gene_disease <-  bed_intersect(trrust, tmp_cnv) %>%
  filter(.overlap > threshold_30_tmp) %>%
  select(target.x) %>% rename(gene = target.x) %>% distinct() %>% filter(!gene %in% vector_genes) %>%
  filter(gene %in% (hgcn_genes %>% filter(disease == 'Yes') %>% pull(gene))) %>% nrow()


# Number open chromatin regions overlapping

result_n_open <- bed_intersect(tmp_cnv, ensembl_reg %>% filter(type == 'open_chromatin_region')) %>% nrow()

# Number TFBs overlapping

result_n_tfbs <- bed_intersect(tmp_cnv, ensembl_reg %>% filter(type == 'TF_binding_site')) %>% nrow()

# Number CTCF regions overlapping

result_n_ctcf <- bed_intersect(tmp_cnv, ensembl_reg %>% filter(type == 'CTCF_binding_site')) %>% nrow()

# Number essential genes


result_n_essent_genes_cl <- hgcn_genes %>% 
  filter(gene %in% vector_genes) %>% 
  filter(str_detect(fusil, 'CL')) %>% 
  nrow()

result_n_essent_genes_dl <- hgcn_genes %>% 
  filter(gene %in% vector_genes) %>% 
  filter(str_detect(fusil, 'DL')) %>% 
  nrow()


# Number enhancers


result_n_enhancers <- bed_intersect(df_enhancers, tmp_cnv) %>%
  filter(.overlap > threshold_30_tmp) %>%
  select(id.x) %>% distinct() %>% nrow()

genes_enhancers_no_cnv <- bed_intersect(df_enhancers, tmp_cnv) %>%
  filter(.overlap > threshold_30_tmp) %>%
  select(gene.x) %>% rename(gene = gene.x) %>% distinct() %>% filter(!gene %in% vector_genes)

result_enhancer_gene_disease <-  genes_enhancers_no_cnv %>%
  filter(gene %in% (hgcn_genes %>% filter(disease == 'Yes') %>% 
                      pull(gene))) %>% 
  nrow()
  

n_genes_one_hpo <- genes_enhancers_no_cnv %>% filter(gene  %in% hpo_genes$gene) %>% nrow()

if(length(n_genes_one_hpo) == 0) {
  vector_hpo_genes <- NA
} else {
  
  vector_hpo_genes <- hpo_genes %>% filter(gene %in% (genes_enhancers_no_cnv %>% pull(gene))) %>% 
    pull(gene) %>% unique() %>% paste(collapse = ', ')
  
}


# Number miRNAs

result_n_mirna <-bed_intersect(mirtarbase, tmp_cnv) %>%
  filter(.overlap > threshold_30_tmp) %>%
  select(id.x) %>% distinct() %>% nrow()

result_mirnas_gene_disease <-  bed_intersect(mirtarbase, tmp_cnv) %>%
  filter(.overlap > threshold_30_tmp) %>%
  select(gene_symbol.x) %>% rename(gene = gene_symbol.x) %>% distinct() %>% filter(!gene %in% vector_genes) %>%
  filter(gene %in% (hgcn_genes %>% filter(disease == 'Yes') %>% pull(gene))) %>% nrow()







# # FUNCTIONAL SIMILARITY OF GENES
# 
# g1 <- hgcn_genes %>% filter(gene %in% vector_genes) %>% pull(entrez_id) %>% as.character()
# 
# result_similarity_wang <- ifelse(length(g1) == 0, NA, 
#                                  geneSim(g1, g1, measure="Wang", combine="BMA") %>% mean(na.rm = TRUE))

# AGGREGATION

result_tmp <- tibble(
  'id' = id_tmp,
  'chrom' = chrom_tmp,
  'start' = start_tmp,
  'end' = end_tmp,
  'clinical' = clinical_tmp,
  'type_variant' = type_variant_tmp,
  'type_inheritance' =   type_inheritance_tmp,
  'n_cnv_syndromes' = result_n_overlap_cnv_syndrome,
  'patho_cnv' = result_n_overlap_patho,
  'nonpatho_cnv' = result_n_overlap_nonpatho,
  'disease_genes' = result_n_disease_genes,
  'disease_variants' = result_n_variants_genes,
  'embryo_mouse' = result_n_mouse_embryo,
  'max_pli' = maximum_pli,
  'max_hi' = maximum_hi,
  'max_ccr' = maximum_ccr,
  'length_cnv' = length_tmp,
  'n_genes_hpo' = n_genes_one_hpo,
  'n_genes' = result_n_genes,
  'n_ohno' = result_n_ohno,
  'n_tf' = result_n_tf,
  'n_target_drugs' = result_n_target_drugs,
  'n_blacklist' = result_n_blacklist,
  'n_prot_complex' = result_n_prot_complex,
  'min_vg' = minimum_vg,
  'max_par' = result_max_par,
  'max_degree' = result_max_degree,
  'max_page_rank' = result_max_page_rank,
  'max_phast' = result_max_phast,
  'max_cpg' = result_max_cpg,
  'n_systems' = result_n_systems,
  'n_open' = result_n_open,
  'n_tfbs' = result_n_tfbs,
  'n_ctcf' = result_n_ctcf,
  'recomb_rate' = result_recomb,
  'dist_cent' = dist_cent,
  'dist_tel' = dist_tel,
  'pubmed_del' = result_max_del,
  'pubmed_dup' = result_max_dup,
  'essent_cl' = result_n_essent_genes_cl,
  'essent_dl' = result_n_essent_genes_dl,
  'n_enhancers' = result_n_enhancers,
  'n_mirnas' = result_n_mirna,
  'enh_gene_disease' = result_enhancer_gene_disease,
  'enh_gene_n_hpo' = n_genes_one_hpo,
  'enh_gene_hpo_vector' = vector_hpo_genes,
  'mirna_gene_disease' = result_mirnas_gene_disease,
  'tf_gene_disease' = result_tfs_gene_disease,
  'gene_density' = result_gene_density
)

return(result_tmp)

}

# ------------------------------------------------------------------------------
# PARALLEL ANNOTATION VISUALIZATION
# ------------------------------------------------------------------------------

plan("multiprocess", workers = 40)

a <- check_cnv(input_check_cnv$id[1],
          input_check_cnv$pathogenicity[1],
          input_check_cnv$variant_class[1],
          input_check_cnv$inheritance[1],
          input_check_cnv$chrom[1],
          input_check_cnv$start[1],
          input_check_cnv$end[1]
          )

tic()

output_list <- future_pmap(list(input_check_cnv$id_tmp, 
                         input_check_cnv$pathogenicity, 
                         input_check_cnv$variant_class,
                         input_check_cnv$inheritance,
                         input_check_cnv$chrom, 
                         input_check_cnv$start, 
                         input_check_cnv$end), 
                       check_cnv)

output_df <- bind_rows(lapply(output_list, as.data.frame.list)) %>% as_tibble()

toc()

# ------------------------------------------------------------------------------
# RMARKDOWN VISUALIZATION
# ------------------------------------------------------------------------------


output_df %>% write_tsv('output_df.tsv')


library(prettydoc)
rmarkdown::render("output_annotation_cnvs.Rmd", 'html_document')


# ------------------------------------------------------------------------------
# LOGISTIC REGRESSION
# ------------------------------------------------------------------------------

# temporal


model_df <- output_df %>% 
  mutate(clinical = if_else(is.na(clinical), 'unlabeled', clinical)) %>%
  separate(id, sep = '_', c('id', 'random_n')) %>%
  mutate(clinical = factor(clinical))
  
remove_ids <- model_df %>% group_by(id) %>% count() %>% filter(n != 10) %>% pull(id)
model_df <- model_df %>% filter(! id %in% remove_ids)


total_ids <- model_df %>%
  select(id) %>% 
  distinct() %>% 
  pull()

training_ids <- total_ids[sample(length(total_ids)*0.7)]

training_tbl <- model_df %>% filter(id %in% training_ids)
test_tbl <- model_df %>% filter(!id %in% training_ids)


training_cv <- vfold_cv(training_tbl, prop = 0.9, strata = clinical)


logistic_rs <- fit_resamples(logistic_reg(mode = 'classification') %>% 
                               set_engine(engine = 'glm'),
                             clinical ~ n_genes + disease_genes + n_genes_hpo ,
                             control = control_resamples(save_pred = TRUE),
                             training_cv)

xgboost_rs <- fit_resamples(boost_tree(mode = 'classification') %>% 
                              set_engine(engine = 'xgboost'),
                            clinical ~ n_genes + disease_genes + n_genes_hpo ,
                            control = control_resamples(save_pred = TRUE),
                            training_cv)

# ------------------------------------------------------------------------------
# 10CV - ROC AND PR CURVES
# ------------------------------------------------------------------------------

from_cv <- function(model, name_model) {
  
  # ROC AUC
  tmp_mean <- model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    roc_auc(clinical, .pred_unlabeled) %>% 
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <-model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    roc_auc(clinical, .pred_unlabeled) %>% 
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p1 <- model %>% 
    unnest(.predictions) %>%
    group_by(id) %>%
    roc_curve(clinical, .pred_unlabeled) %>%
    ggplot(aes(1-specificity, sensitivity)) +
    geom_path(aes(group = id, color = id),  show.legend = FALSE) +
    theme_bw() +
    labs(title = glue('{name_model}  (Mean AUC: {tmp_mean} ± {tmp_err})'))
  
  # PR AUC
  tmp_mean <- model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    pr_auc(clinical, .pred_unlabeled) %>% 
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <-model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    pr_auc(clinical, .pred_unlabeled) %>% 
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p2 <- model %>% 
    unnest(.predictions) %>%
    group_by(id) %>%
    pr_curve(clinical, .pred_unlabeled) %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = id, color = id), show.legend = FALSE) +
    theme_bw() +
    labs(title = glue('{name_model} (Mean AUCpr: {tmp_mean} ± {tmp_err})'))
  
  list(p1, p2)
  
}


a <- from_cv(xgboost_rs, name_model = 'Logistic regression')

a1 <- a[[1]]
a2 <- a[[2]]

a1 + a2

# ------------------------------------------------------------------------------
# TEST DATASET
# ------------------------------------------------------------------------------
# xrf_model <- xrf(formule_models, data = training_tbl_to_model, family = "binomial")


training_tbl_to_model <- training_tbl %>% select(-id, -random_n, -chrom, -start, -end, -type_variant,
                        -type_inheritance, -enh_gene_hpo_vector)
# n_cnv_syndromes + n_genes_hpo + max_ccr + max_hi + max_pli + embryo_mouse
# essent_cl + max_cpg + max_phast  + n_prot_complex
formule_models <- clinical ~   n_cnv_syndromes + n_genes_hpo + max_ccr + max_hi + max_pli + embryo_mouse +
  essent_cl + max_cpg + max_phast  + n_prot_complex + patho_cnv + nonpatho_cnv + n_systems
  



xgboost_model <- boost_tree() %>%
  set_mode('classification') %>%
  set_engine('xgboost') %>%
  fit(formule_models, data = training_tbl_to_model)

forest_model <- rand_forest() %>%
  set_mode('classification') %>%
  set_engine('ranger') %>%
  fit(formule_models, data = training_tbl_to_model)


logistic_model <- logistic_reg() %>%
  set_mode('classification') %>%
  set_engine('glm') %>%
  fit(formule_models, data = training_tbl_to_model)

rulefit_model1 <- rule_fit(trees = 100, tree_depth = 3, penalty = 0.3) %>%
  set_mode('classification') %>%
  set_engine('xrf') %>%
  fit(formule_models, data = training_tbl_to_model)

rulefit_model2 <- xrf(formule_models,
                     data = training_tbl_to_model,
                     # xgb_control = list(nrounds = 100, max_depth = 5, min_child_weight = 3),
                     family = "binomial")



perc_xgboost <- predict(xgboost_model, test_tbl, type = 'prob') %>%
  rename(pred_patho = .pred_Pathogenic) %>%
  select(pred_patho) %>%
  bind_cols(test_tbl) %>%
  select(id, pred_patho, random_n) %>%
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  # filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))

perc_forest <- predict(forest_model, test_tbl, type = 'prob') %>%
    rename(pred_patho = .pred_Pathogenic) %>%
  select(pred_patho) %>% 
  bind_cols(test_tbl) %>%
  select(id, pred_patho, random_n) %>% 
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  # filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))


perc_xrf <- tibble('pred_patho' = 1 - as.vector(predict(rulefit_model1, test_tbl, type = 'response'))) %>%
  predict(rulefit_model1, test_tbl, type = 'prob') %>%
  rename(pred_patho = .pred_Pathogenic) %>%
  select(pred_patho) %>%
  bind_cols(test_tbl) %>%
  select(id, pred_patho, random_n) %>%
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))


perc_xrf <- predict(rulefit_model, test_tbl, type = 'prob') %>%
  rename(pred_patho = .pred_Pathogenic) %>%
  select(pred_patho) %>% 
  bind_cols(test_tbl) %>%
  select(id, pred_patho, random_n) %>% 
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  # filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))


perc_logistic <- predict(logistic_model, test_tbl, type = 'prob') %>%
  rename(pred_patho = .pred_Pathogenic) %>%
  select(pred_patho) %>% 
  bind_cols(test_tbl) %>%
  select(id, pred_patho, random_n) %>% 
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  # filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))


##

# ------------------------------------------------------------------------------
# PLOT COMPARISON
# ------------------------------------------------------------------------------

  ggplot() +
    geom_path(data = tibble('rank' = seq(1,10), 'accum_perc' = seq(10,100,10)), aes(rank, accum_perc, color = 'red')) +
    geom_point(data = tibble('rank' = seq(1,10), 'accum_perc' = seq(10,100,10)), aes(rank, accum_perc)) +
    geom_path(data = perc_xgboost, aes(rank, accum_perc, color = 'blue'), ) +
    geom_point(data = perc_xgboost, aes(rank, accum_perc)) +
    geom_path(data = perc_xrf, aes(rank, accum_perc, color = 'green') ) +
    geom_point(data = perc_xrf, aes(rank, accum_perc)) +
    geom_path(data = perc_forest, aes(rank, accum_perc, color = 'yellow') ) +
    geom_point(data = perc_forest, aes(rank, accum_perc)) +
    geom_path(data = perc_logistic, aes(rank, accum_perc, color = 'orange') ) +
    geom_point(data = perc_logistic, aes(rank, accum_perc)) +
    scale_colour_manual(name = 'Model',
                      values =c('red'='red','blue'='blue', 'green' = 'green', 'yellow' = 'yellow', 'orange' = 'orange'), 
                      labels = c('xgboost','rulefit', 'logistic regression', 'random','random forest')) +
    scale_x_continuous(breaks = scales::pretty_breaks(10)) +
    theme_bw() +
    labs(x = 'Rank', y = 'Cumulative Frequency (%)')

# ------------------------------------------------------------------------------
# HYPERPARAMETER TUNING
# ------------------------------------------------------------------------------

# penalty = 0.003
# tree_depth = 8


xgboost_model <- rule_fit(trees = 100, tree_depth = tune(), penalty = tune()) %>%
  set_mode('classification') %>%
  set_engine('xrf')

tree_grid <- grid_regular(tree_depth(), penalty(), levels = 5)

training_folds <- vfold_cv(training_tbl_to_model[1:5000,], v =  2)


tree_wf <- workflow() %>%
  add_model(xgboost_model) %>%
  add_formula(formule_models)

tree_res <- 
  tree_wf %>% 
  tune_grid(
    resamples = training_folds,
    grid = tree_grid
  )



tree_res %>%
  collect_metrics() %>%
  # mutate(tree_depth = factor(tree_depth)) %>%
  ggplot(aes(tree_depth, mean)) +
  geom_line(size = 1.5, alpha = 0.6,  aes(color = penalty )) +
  geom_point(size = 2) +
  facet_wrap(~ .metric, scales = "free", nrow = 2) +
  scale_x_log10() +
  scale_color_viridis_d(option = "plasma", begin = .9, end = 0)


# trees = 500, tree_depth = 4
best_tree <- tree_res %>% select_best(metric = 'roc_auc')

final_wf <- 
  tree_wf %>% 
  finalize_workflow(best_tree)

# ------------------------------------------------------------------------------
# FINAL MODEL 
# ------------------------------------------------------------------------------

final_tree <- 
  final_wf %>%
  fit(data = training_tbl_to_model) 


final_xgboost <- predict(final_tree, test_tbl, type = 'prob') %>%
  rename(pred_patho = .pred_Pathogenic) %>%
  select(pred_patho) %>% 
  bind_cols(test_tbl) %>%
  select(id, pred_patho, random_n) %>% 
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  # filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))

ggplot() +
  geom_path(data = tibble('rank' = seq(1,10), 'accum_perc' = seq(10,100,10)), aes(rank, accum_perc, color = 'red')) +
  geom_point(data = tibble('rank' = seq(1,10), 'accum_perc' = seq(10,100,10)), aes(rank, accum_perc)) +
  geom_path(data = final_xgboost, aes(rank, accum_perc, color = 'orange') ) +
  geom_point(data = final_xgboost, aes(rank, accum_perc)) +
  # scale_colour_manual(name = 'Model',
  #                     values =c('red'='red','blue'='blue', 'green' = 'green', 'yellow' = 'yellow', 'orange' = 'orange'), 
  #                     labels = c('xgboost','rulefit', 'logistic regression', 'random','random forest')) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  theme_bw() +
  labs(x = 'Rank', y = 'Cumulative Frequency (%)')

# ------------------------------------------------------------------------------
# END
# ------------------------------------------------------------------------------

output_xrf <- xrf(formule_models,
                  training_tbl_to_model, family = 'binomial',
                  xgb_control = list(nrounds = 100, max_depth = 2))

data_pre_lasso <- output_xrf$full_data

rules_tbl <- output_xrf %>%
  coef(s = "lambda.min") %>%
  as_tibble() %>%
  select(term, rule) %>%
  na.omit()
# separate_rows(rule, sep = '&') %>%
# slice(1) %>%
# pull(rule)

 



# plan("multiprocess", workers = 40)
# 
# tic()
# 
# output_list <- future_pmap(list(input_check_cnv$id_tmp, 
#                                 input_check_cnv$pathogenicity, 
#                                 input_check_cnv$variant_class,
#                                 input_check_cnv$inheritance,
#                                 input_check_cnv$chrom, 
#                                 input_check_cnv$start, 
#                                 input_check_cnv$end), 
#                            check_cnv)
# 
# output_df <- bind_rows(lapply(output_list, as.data.frame.list)) %>% as_tibble()
# 
# toc()

bayesian_model <- rstanarm::stan_glm(formule_models,
                                     family = 'binomial',
                                     data = data_pre_lasso,
                                     cores = 4,
                                     iter = 4000,
                                     chains = 4)
                                     # QR = TRUE,
                                     # prior = laplace())


ready_test <- tibble()

for (i in 1:nrow(test_tbl)) {
  print(i)
  tmp_test_row <- test_tbl %>% slice(i)
  para_test <- rules_tbl %>%
    rowwise() %>%
    mutate(ok = tmp_test_row %>% filter_(rule) %>% nrow()) %>%
    pivot_wider(-rule,values_from = ok, names_from = term)
  ready_test_tmp <- para_test %>% bind_cols(tmp_test_row)
  ready_test <- ready_test %>% bind_rows(ready_test_tmp)
}

# library(ggridges)
# posterior_linpred(bayesian_model, newdata = ready_test, transform = TRUE) %>%
#   as_tibble() %>%
#   pivot_longer(everything(), names_to = 'obs', values_to = 'prob') %>%
#   ggplot(aes(prob, obs)) +
#     geom_density_ridges()

pred_bayesian  <- posterior_linpred(bayesian_model, newdata = ready_test, transform = TRUE) %>%
  as_tibble() %>%
  map_dbl(~ mean(.x))

pred_lasso <- predict(output_xrf, ready_test, type = 'response')


check_pred <- ready_test %>%
  select(clinical) %>%
  mutate(pred_lasso = pred_lasso %>% as.vector(),
         pred_bayesian = pred_bayesian)


roc_lasso <- check_pred %>% roc_curve(clinical, pred_lasso) %>% mutate(model = 'lasso')
roc_bayesian <- check_pred %>% roc_curve(clinical, pred_bayesian) %>% mutate(model = 'bayesian')

auc_lasso <- check_pred %>% roc_auc(clinical, pred_lasso) %>% pull(.estimate) %>% round(3)
auc_bayesian <- check_pred %>% roc_auc(clinical, pred_bayesian) %>% pull(.estimate) %>% round(3)

roc_both <- roc_lasso %>% bind_rows(roc_bayesian)


roc_both %>%
  ggplot(aes((1-specificity), sensitivity)) +
  geom_line(aes(color = model), size = 1) +
  theme_bw() +
  geom_abline(linetype = 3) +
  labs(title = 'Comparison performance logistic and Bayesian model',
       subtitle = glue('AUC (lasso) = {auc_lasso} - AUC (Bayesian) = {auc_bayesian}'))



perc_xrf <- tibble('pred_patho' = 1 - as.vector(predict(output_xrf, ready_test, type = 'response'))) %>%
  # predict(rulefit_model1, test_tbl, type = 'prob') %>%
  # rename(pred_patho = .pred_Pathogenic) %>%
  # select(pred_patho) %>%
  bind_cols(ready_test) %>%
  select(id, pred_patho, random_n) %>%
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  # filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))

perc_xrf_bayesian <- tibble('pred_patho' = 1 - as.vector(pred_bayesian)) %>%
  # predict(rulefit_model1, test_tbl, type = 'prob') %>%
  # rename(pred_patho = .pred_Pathogenic) %>%
  # select(pred_patho) %>%
  bind_cols(ready_test) %>%
  select(id, pred_patho, random_n) %>%
  mutate(random_n = if_else(is.na(random_n), 'original', random_n )) %>%
  group_by(id) %>%
  arrange(desc(pred_patho)) %>%
  mutate(rank = row_number()) %>%
  # filter(id == '2169') %>%
  ungroup() %>%
  group_by(random_n, rank) %>%
  count() %>%
  filter(random_n == 'original') %>%
  ungroup() %>%
  mutate(perc = 100*(n / sum(n))) %>%
  mutate(accum_perc = cumsum(perc)) %>%
  mutate(same_group = 'yes') %>%
  mutate(rank = as.integer(rank))

ggplot() +
  geom_path(data = tibble('rank' = seq(1,10), 'accum_perc' = seq(10,100,10)), aes(rank, accum_perc, color = 'red')) +
  geom_point(data = tibble('rank' = seq(1,10), 'accum_perc' = seq(10,100,10)), aes(rank, accum_perc)) +
  geom_path(data = perc_xrf, aes(rank, accum_perc, color = 'orange') ) +
  geom_point(data = perc_xrf, aes(rank, accum_perc)) +
  geom_path(data = perc_xrf_bayesian, aes(rank, accum_perc, color = 'red') ) +
  geom_point(data = perc_xrf_bayesian, aes(rank, accum_perc)) +
  # scale_colour_manual(name = 'Model',
  #                     values =c('red'='red','blue'='blue', 'green' = 'green', 'yellow' = 'yellow', 'orange' = 'orange'), 
  #                     labels = c('xgboost','rulefit', 'logistic regression', 'random','random forest')) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  theme_bw() +
  labs(x = 'Rank', y = 'Cumulative Frequency (%)')



# ------------------------------------------------------------------------------
# END
# ------------------------------------------------------------------------------


# plan("multiprocess", workers = 2)
# 
# tic()
# 
# 
# xd <- output_df %>%
#   select(id, chrom, start, end) %>%
#   mutate(n_hits = future_pmap_chr(list(id, chrom, start, end),  get_enrich))
# 
# toc()














# 
# output_df %>%
# 
# to_cba <- output_df %>% 
#   # mutate(clinical = if_else(clinical == 'Likely benign', 'Benign', clinical)) %>%
#   filter(clinical %in% c('Pathogenic', 'Benign', 'Likely benign', 'Likely pathogenic')) %>%
#   mutate(clinical = as.factor(clinical)) %>%
#   mutate(type_variant = if_else(type_variant == 'Deletion', 1, 0)) %>%
#   mutate(type_inheritance = if_else(type_inheritance == 'De novo constitutive', 1, 0))
#   # mutate(max_pli = if_else(maximum_pli >= 90, 1, 0)) %>%
#   # mutate(embryo_mouse = if_else(embryo_mouse >= 1, 1, 0))
# 
# 
# # to_cba[,-c(1:2,12)] <- map_df(to_cba[,-c(1:2,12 )], function(x) if_else(x > 0, 1, 0))
# 
# 
# 
# 
# model_cba <- arulesCBA::CBA(clinical ~ n_cnv_syndromes + type_variant + embryo_mouse + 
#                               maximum_pli + n_genes_hpo + n_genes + n_tf + n_blacklist + n_target_drugs +
#                               type_inheritance + patho_cnv + nonpatho_cnv + disease_genes + disease_variants, 
#                             data = to_cba %>% select(-length_cnv) ,
#                             # method = 'first',
#                             support = 0.005, 
#                             confidence = 0.7)
# 
# 


# 
# output_df %>%
#   # mutate(n_systems = if_else(n_systems > 1, 'Yes', 'No')) %>%
#   mutate(max_pli = if_else(max_pli >= 90, 'Yes', 'No')) %>%
#   mutate(max_hi = if_else(max_hi >= 90, 'Yes', 'No')) %>%
#   mutate_if(is.integer, ~ if_else(. > 0, 'Yes', 'No')) %>%
#   mutate_if(is.double, ~ if_else(. > 0, 'Yes', 'No')) %>%
#   # mutate(clinical = if_else(clinical == 'Likely benign', 'Benign', clinical)) %>%
#   # mutate(clinical = if_else(clinical == 'Likely pathogenic', 'Pathogenic', clinical)) %>%
# 
#   count(clinical,
#         n_cnv_syndromes,
#         disease_genes,
#         patho_cnv,
#         n_ohno,
#         n_tf,
#         n_target_drugs,
#         n_prot_complex,
#         n_tfbs,
#         n_ctcf,
#         pubmed_del,
#         pubmed_dup,
#         n_open,
#         max_pli,
#         max_ccr,
#         max_hi,
#         disease_variants,
#         nonpatho_cnv,
#         essent_cl,
#         essent_dl,
#         enh_gene_disease,
#         mirna_gene_disease,
#         tf_gene_disease,
#         n_genes_hpo,
#         n_blacklist) %>%
#   mutate(clinical = factor(clinical,levels = c("Pathogenic",'Likely pathogenic', "Unknown", "Uncertain",
#                                                'Likely benign', "Benign"))) %>%
#   group_by(clinical) %>%
#   mutate(perc = n / sum(n)*100) %>%
#   select(-n) %>%
#   pivot_longer(-c(clinical,perc),  names_to = 'rule', values_to = 'yes_no') %>%
#   ggplot(aes(clinical, perc)) +
#   geom_col(aes(fill = yes_no)) +
#   theme_fancy() +
#   facet_wrap(~ rule) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   xlab('Clinical significance') +
#   ylab('Percentage (%)')
# 
# Plotting distances (telomeric and centromeric regions)
# output_df %>% ggplot(aes(dist_cent)) + 
#   geom_histogram(binwidth = 5, fill = 'steelblue', color = 'black') + facet_grid(~ clinical) +
#   theme_bw()
# 
# output_df %>% ggplot(aes(dist_tel)) + geom_histogram(binwidth = 5, fill = 'steelblue', color = 'black') + 
#   facet_grid(~ clinical) +
#   theme_bw()
# 
# output_df %>% ggplot(aes(dist_cent)) + 
#   geom_density(aes(fill = clinical), color = 'black', show.legend = FALSE, alpha = 0.7) + facet_grid(~ clinical) +
#   theme_bw()
# 
# output_df %>% ggplot(aes(dist_tel)) + 
#   geom_density(aes(fill = clinical), color = 'black', show.legend = FALSE) + facet_grid(~ clinical) +
#   theme_bw()
# 
# # Gene density ~ clinical
# output_df %>% ggplot(aes(gene_density)) + 
#   geom_density(aes(fill = clinical), color = 'black', alpha = 0.6) + facet_grid(~ clinical) +
#   theme_bw()
# 
# output_df %>% ggplot(aes(clinical, gene_density)) + 
#   geom_boxplot(aes(fill = clinical), color = 'black') +
#   theme_bw()
# 
# # Plotting distance
# output_df %>% 
#   ggplot(aes(length_cnv, clinical)) + 
#   stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = clinical), 
#                       alpha = 0.6, show.legend = FALSE, size = 1.25) +
#   theme_ridges() +
#   scale_x_log10() +
#   xlab('log10(Length CNV)')





# Plotting correlation with length_cnv
# 
# output_df %>% select(-clinical, -type_variant, -type_inheritance, -id) %>% correlate(method = 'pearson') %>% select(rowname, length_cnv) %>% na.omit()  %>% ggplot(aes(reorder(rowname, length_cnv), length_cnv)) + 
#   geom_col(fill = 'steelblue', color = 'black') + coord_flip() + theme_bw()
# 
# # PCA
# pca_output <- output_df %>% select(-clinical, -type_variant, -type_inheritance, -id, -dist_tel, -gene_density) %>% 
#   prcomp(scale = TRUE, center = TRUE)
#   
# a$x %>% 
#   as_tibble() %>% 
#   select(PC1, PC2) %>%  
#   bind_cols(output_df %>% select(clinical)) %>% 
#   ggplot(aes(PC1, PC2)) + 
#   geom_point(aes(fill = factor(clinical)), shape = 21, color = 'black', show.legend = FALSE) +
#   facet_grid(~ clinical) +
#   theme_bw()
# 
# a %>%  fviz_contrib(choice = "var", axes = 1, top = 15)
# a %>%  fviz_contrib(choice = "var", axes = 2, top = 15)
# a %>%  fviz_eig(addlabels = TRUE, ylim = c(0, 50))





# 
# 
# 
# 
# output_df %>% 
#   mutate(max_pli = if_else(max_pli >= 90, 'Yes', 'No')) %>%
#   mutate(max_hi = if_else(max_hi >= 90, 'Yes', 'No')) %>%
#   mutate_if(is.integer, ~ if_else(. > 0, 'Yes', 'No')) %>%
#   mutate_if(is.double, ~ if_else(. > 0, 'Yes', 'No')) %>%
#   mutate(clinical = if_else(clinical == 'Likely benign', 'Benign', clinical)) %>%
#   mutate(clinical = if_else(clinical == 'Likely pathogenic', 'Pathogenic', clinical)) %>%
#   filter(clinical == 'Benign') %>%
#   count(n_cnv_syndromes, patho_cnv, disease_genes, disease_variants)
# 
#   

# result_df <- first_df %>%
#   mutate(patho_cnv = patho_cnv + n_cnv_syndromes) %>%
#   mutate(disease_genes = disease_genes + disease_variants) %>%
#   select(-type_inheritance, -type_variant, -disease_variants, -n_cnv_syndromes)
# 
# result_df[,-c(1:2)] <- map_df(result_df[,-c(1:2)], function(x) if_else(x > 0, 1, 0))
# 

# 
# 
# tmp_plot <- result_df %>% filter(clinical %in% c('Pathogenic', 'Likely pathogenic')) %>% select(-clinical )
# n_no_intersect <- tmp_plot %>% filter(patho_cnv == 0, nonpatho_cnv == 0, disease_genes == 0) %>% nrow()
# title_plot <- paste('Pathogenic - Likely pathogenic CNVs','\n', 'No intersection (', n_no_intersect, 'CNVs)')
# 
# upset(tmp_plot %>% as.data.frame() , order.by = "freq" ,
#       point.size = 3.5, line.size = 2, number.angles = 0, sets.x.label = 'Number of CNVs',
#       text.scale = c(1.3, 1.3, 2, 2, 2, 2))
# grid.text(title_plot,x = 0.75, y=0.95, gp=gpar(fontsize=16))
# 
# 
# 
# 
# tmp_plot <- result_df %>% filter(!clinical %in% c('Pathogenic', 'Likely pathogenic')) %>% select(-clinical )
# n_no_intersect <- tmp_plot %>% filter(patho_cnv == 0, nonpatho_cnv == 0, disease_genes == 0) %>% nrow()
# title_plot <- paste('Benign/Likely benign/Uncertain/Unknown CNVs','\n', 'No intersection (', n_no_intersect, 'CNVs)')
# 
# upset(tmp_plot %>% as.data.frame() , order.by = "freq" ,
#       point.size = 3.5, line.size = 2, number.angles = 0, sets.x.label = 'Number of CNVs',
#       text.scale = c(1.3, 1.3, 2, 2, 2, 2))
# grid.text(title_plot,x = 0.75, y=0.95, gp=gpar(fontsize=16))
