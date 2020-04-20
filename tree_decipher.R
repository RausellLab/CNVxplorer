# library(enrichR)
# library(DOSE)
# library(ReactomePA)
# library(clusterProfiler)
# library(UpSetR)

library(valr)
library(future)
library(tictoc)
library(furrr)
library(grid)
library(tidyverse)
library(tictoc)
library(ontologyIndex)
library(corrr)



rename <- dplyr::rename
slice <- dplyr::slice

# save(hgcn_genes, df_enhancers, lncrna_coord, lncrna, tad, gtex, hpa, hpo_genes, cnv_df, vector_total_terms,
#      gnomad_sv_raw, decipher_control_raw, dgv_df_raw, hpo_omim, anato_df, mirtarbase, pubmed_df,
#       vector_inheritance, trrust, tf_genes, drugbank,prot_complex,ohno_genes,recomb,
#       genes_promoter,para_genes,string_db,region_gaps, fusil_score,ensembl_reg,
#     select, dev_raw, panel_total, omim, orphanet_raw,  hpo_dbs, model1, denovo, clinvar_variants, ridges_home, plot_p100, plot_p46pla, blacklist_encode, mpo_dbs, gwas_variants,mgi, syndromes_total,
# file = "env_annot_cnvs.RData")
# 
# load('env_annot_cnvs.RData')

theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size=22))
  
}

input_check_cnv <-  read_tsv('/data-cbl/frequena_data/from_workstation/daa_decipher/decipher-cnvs-grch37-2020-01-19.txt', skip = 1) %>%
  as_tibble() %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  select(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  # filter(pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>% 
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>% 
  # filter(inheritance == 'De novo constitutive') %>% 
  filter(variant_class %in% c('Deletion', 'Duplication'))
  # mutate(id_tmp = row_number())


input_check_overlap <- input_check_cnv %>% select(id,pathogenicity,variant_class, 
                                                  inheritance,chrom,start, end) %>%
  bind_rows(gnomad_sv_raw %>% mutate(pathogenicity = 'benign_gnomad',
                                     inheritance = NA) %>%
              rename(variant_class = svtype) %>%
              select(id,pathogenicity,
                     chrom, start, end))


check_cnv <- function(input_id, input_clinical, input_variant, input_inheritance,
                      input_chrom, input_start, input_end) {

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

  just_test <- input_check_cnv %>% filter(id == '131') %>% slice(1)
  id_tmp <- just_test$id
  clinical_tmp <-  just_test$pathogenicity
  type_variant_tmp <-  just_test$variant_class
  type_inheritance_tmp <-  just_test$inheritance
  chrom_tmp <-  just_test$chrom
  start_tmp <-  just_test$start
  end_tmp <-  just_test$end
  length_tmp <- end_tmp - start_tmp + 1
  threshold_30_tmp  <- 0

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
    # mutate(perc_overlap = .overlap / length_gene)
    # filter(perc_overlap >= 0.3) %>%
    pull(gene.x)
  
  vector_entrez <- hgcn_genes %>%
    rename(start = start_position, end = end_position) %>%
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

# N Number of genes identified as Cellular letal (FUSIL score)


# result_n_genes_cl <- vector_genes %in% fusil_scoregene %>% sum()


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


result_n_essent_genes_cl <- length(vector_genes[vector_genes %in% 
                                               (fusil_score %>% filter(FUSIL == 'CL') %>% pull(hgnc_symbol))]) 
  
  
result_n_essent_genes_dl <- length(vector_genes[vector_genes %in% 
                                                  (fusil_score %>% filter(FUSIL == 'DL') %>% pull(hgnc_symbol))]) 



# Number enhancers


result_n_enhancers <- bed_intersect(df_enhancers, tmp_cnv) %>%
  filter(.overlap > threshold_30_tmp) %>%
  select(id.x) %>% distinct() %>% nrow()


result_enhancer_gene_disease <-  bed_intersect(df_enhancers, tmp_cnv) %>%
  filter(.overlap > threshold_30_tmp) %>%
  select(gene.x) %>% rename(gene = gene.x) %>% distinct() %>% filter(!gene %in% vector_genes) %>%
  filter(gene %in% (hgcn_genes %>% filter(disease == 'Yes') %>% pull(gene))) %>% nrow()
  
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
  'mirna_gene_disease' = result_mirnas_gene_disease,
  'tf_gene_disease' = result_tfs_gene_disease
)

return(result_tmp)

}

plan("multiprocess", workers = 40)


# test_before <- input_check_cnv %>% slice(sample(1:nrow(.), 100)) 




tic()

output_list <- future_pmap(list(input_check_cnv$id, 
                         input_check_cnv$pathogenicity, 
                         input_check_cnv$variant_class,
                         input_check_cnv$inheritance,
                         input_check_cnv$chrom, 
                         input_check_cnv$start, 
                         input_check_cnv$end), 
                       check_cnv)

output_df <- bind_rows(lapply(output_list, as.data.frame.list)) %>% as_tibble()

toc()

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
output_df %>%
  # mutate(n_systems = if_else(n_systems > 1, 'Yes', 'No')) %>%
  mutate(max_pli = if_else(max_pli >= 90, 'Yes', 'No')) %>%
  mutate(max_hi = if_else(max_hi >= 90, 'Yes', 'No')) %>%
  mutate_if(is.integer, ~ if_else(. > 0, 'Yes', 'No')) %>%
  mutate_if(is.double, ~ if_else(. > 0, 'Yes', 'No')) %>%
  # mutate(clinical = if_else(clinical == 'Likely benign', 'Benign', clinical)) %>%
  # mutate(clinical = if_else(clinical == 'Likely pathogenic', 'Pathogenic', clinical)) %>%

  count(clinical,
        n_cnv_syndromes,
        disease_genes,
        patho_cnv,
        n_ohno,
        n_tf,
        n_target_drugs,
        n_prot_complex,
        n_tfbs,
        n_ctcf,
        pubmed_del,
        pubmed_dup,
        n_open,
        max_pli,
        max_ccr,
        max_hi,
        disease_variants,
        nonpatho_cnv,
        essent_cl,
        essent_dl,
        enh_gene_disease,
        mirna_gene_disease,
        tf_gene_disease,
        n_genes_hpo,
        n_blacklist) %>%
  mutate(clinical = factor(clinical,levels = c("Pathogenic",'Likely pathogenic', "Unknown", "Uncertain",
                                               'Likely benign', "Benign"))) %>%
  group_by(clinical) %>%
  mutate(perc = n / sum(n)*100) %>%
  select(-n) %>%
  pivot_longer(-c(clinical,perc),  names_to = 'rule', values_to = 'yes_no') %>%
  ggplot(aes(clinical, perc)) +
  geom_col(aes(fill = yes_no)) +
  theme_fancy() +
  facet_wrap(~ rule) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Clinical significance') +
  ylab('Percentage (%)')
# 
# Plotting distances (telomeric and centromeric regions)
output_df %>% ggplot(aes(dist_cent)) + 
  geom_histogram(binwidth = 5, fill = 'steelblue', color = 'black') + facet_grid(~ clinical) +
  theme_bw()

output_df %>% ggplot(aes(dist_tel)) + geom_histogram(binwidth = 5, fill = 'steelblue', color = 'black') + 
  facet_grid(~ clinical) +
  theme_bw()

output_df %>% ggplot(aes(dist_cent)) + 
  geom_density(aes(fill = clinical), color = 'black', show.legend = FALSE) + facet_grid(~ clinical) +
  theme_bw()

output_df %>% ggplot(aes(dist_tel)) + 
  geom_density(aes(fill = clinical), color = 'black') + facet_grid(~ clinical) +
  theme_bw()

# Plotting distance
output_df %>% 
  ggplot(aes(length_cnv, clinical)) + 
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = clinical), 
                      alpha = 0.6, show.legend = FALSE, size = 1.25) +
  theme_ridges() +
  scale_x_log10() +
  xlab('log10(Length CNV)')


# Plotting correlation with length_cnv

output_df %>% select(-clinical, -type_variant, -type_inheritance, -id) %>% correlate(method = 'pearson') %>% select(rowname, length_cnv) %>% na.omit()  %>% ggplot(aes(reorder(rowname, length_cnv), length_cnv)) + 
  geom_col(fill = 'steelblue', color = 'black') + coord_flip() + theme_bw()



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
