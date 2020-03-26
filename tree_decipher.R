

library(valr)
library(future)
library(tictoc)
library(furrr)
library(UpSetR)
library(grid)
library(tidyverse)
library(tictoc)


load('local_data.RData')

result_df <- tibble()


input_check_cnv <-  read_tsv('C:/Users/Requena/Desktop/decipher-cnvs-grch37-2020-01-19.txt', skip = 1) %>%
  as_tibble() %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  mutate(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  # filter(pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>% 
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>% 
  # filter(inheritance == 'De novo constitutive') %>% 
  filter(variant_class %in% c('Deletion', 'Duplication'))


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
  threshold_30_tmp <- round((length_tmp / 100)*30,0)


  # id_tmp <- only_patho$id[1]
  # clinical_tmp <-  only_patho$pathogenicity[1]
  # type_variant_tmp <-  only_patho$variant_class[1]
  # type_inheritance_tmp <-  only_patho$inheritance[1]
  # chrom_tmp <-  only_patho$chrom[1]
  # start_tmp <-  only_patho$start[1]
  # end_tmp <-  only_patho$end[1]
  # length_tmp <- end_tmp - start_tmp + 1

  tmp_cnv <- tibble(chrom = chrom_tmp, start = start_tmp, end = end_tmp)

  # 1. CNV databases
  # # 1.1 RULE - OVERLAP CNV SYNDROMES
  # # 1.2 RULE - OVERLAP ANNOTATED PATHOGENIC
  # # 1.3 RULE - OVERLAP NON-PATHOGENIC CNVs
  # 2. Disease annotation
  # # 2.1 RULE - OVERLAP DISEASE GENES
  # # 2.2 RULE - OVERLAP DISEASE VARIANTS
  # 3. Model mouse information
  # # 3.1 RULE - OVERLAP WITH ORTHOLOGS GENES ASSOCIATED WITH LETHALITY
  # 4. Biological categories

  vector_genes <- hgcn_genes %>%
    rename(start = start_position, end = end_position) %>%
    bed_intersect(tmp_cnv ) %>%
    mutate(length_gene = end.x - start.x + 1) %>%
    mutate(perc_overlap = .overlap / length_gene) %>%
    filter(perc_overlap >= 0.3) %>% pull(gene.x)


# # 1º RULE - OVERLAP CNV SYNDROMES

  tmp_df <- syndromes_total %>%
    filter(chrom == chrom_tmp)


  result_n_overlap_cnv_syndrome <- bed_intersect(tmp_cnv, tmp_df) %>%
    filter(.overlap > threshold_30_tmp) %>%
    nrow()



# 2º RULE - OVERLAP ANNOTATED PATHOGENIC
tmp_df <- cnv_df %>%
  filter(! id %in% id_tmp) %>%
  filter(pathogenicity == 'Pathogenic') %>%
  filter(chrom == chrom_tmp)

result_n_overlap_patho <- bed_intersect(tmp_cnv, tmp_df) %>%
  filter(.overlap > threshold_30_tmp) %>%
  nrow()


# 3º RULE - OVERLAP NON-PATHOGENIC CNVs
tmp_df <- cnv_df %>% filter(! id %in% id_tmp) %>%
  filter(source != 'decipher') %>%
  filter(chrom == chrom_tmp)

result_n_overlap_nonpatho <- bed_intersect(tmp_cnv, tmp_df) %>%
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

result_n_mouse_embryo <- mgi %>% filter(gene %in% vector_genes, str_detect(pheno, 'MP:0005380')) %>% nrow()

# # MAXIMUM pLI score

temporal_df <- hgcn_genes %>% filter(gene %in% vector_genes) %>% select(pLI) %>% na.omit() 
n_rows <- temporal_df %>% nrow()

maximum_pli <- if_else(n_rows == 0, 0, as.double(max(temporal_df$pLI)))




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
  'maximum_pli' = maximum_pli,
  'length_cnv' = length_tmp
)

return(result_tmp)

}

plan("multiprocess", workers = 2)



only_patho <- input_check_cnv %>% slice(1:100)

tic()


output_check <- pmap(list(only_patho$id, only_patho$pathogenicity, only_patho$variant_class,
                            only_patho$inheritance,
                            only_patho$chrom, only_patho$start, only_patho$end), 
                       check_cnv)
toc()

output_check <- bind_rows(lapply(output_check, as.data.frame.list)) %>% as_tibble()

write_tsv(output_check, 'output_check.tsv')













# result_df <- first_df %>%
#   mutate(patho_cnv = patho_cnv + n_cnv_syndromes) %>%
#   mutate(disease_genes = disease_genes + disease_variants) %>%
#   select(-type_inheritance, -type_variant, -disease_variants, -n_cnv_syndromes)
# 
# result_df[,-c(1:2)] <- map_df(result_df[,-c(1:2)], function(x) if_else(x > 0, 1, 0))
# 
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
