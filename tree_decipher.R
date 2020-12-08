
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
library(xrftest)
library(xrf)
library(rules)
library(dials)
library(workflows)
library(rstanarm)
library(bayestestR)
library(prettydoc)
library(applicable)
library(reticulate)


set.seed(123)


rename <- dplyr::rename
slice <- dplyr::slice


# save(hgcn_genes, df_enhancers, lncrna_coord, lncrna, tad, gtex, hpa, hpo_genes, cnv_df, vector_total_terms,
#      gnomad_sv_raw, decipher_control_raw, dgv_df_raw, hpo_omim, anato_df, mirtarbase, pubmed_df,
#       vector_inheritance, trrust, tf_genes, drugbank,prot_complex,ohno_genes,recomb, interactions_db,
#       genes_promoter,para_genes,string_db,region_gaps,ensembl_reg, gene_density_tbl,
#     select, dev_raw, panel_total, omim, orphanet_raw,  hpo_dbs, model1, denovo, clinvar_variants, ridges_home, plot_p100, plot_p46pla, blacklist_encode, mpo_dbs, gwas_variants,mgi, syndromes_total,
# file = "env_annot_cnvs.RData")
# 
load('env_annot_cnvs.RData')
# load('local_data.RData')


theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size=22))
  
}

# ------------------------------------------------------------------------------
# READ CNVs DATASET
# ------------------------------------------------------------------------------

# List of filters decipher
# 1. >= 50b.p
# 2. sv type = deletion or duplication
# 3. pathogenicity != 'benign' or 'likely benign'
# 4. de novo
# 5. heterozygous
# 6. contribution != None

# decipher pathogenic labels (pathogenic - likely pathogenics - unknown)
tmp_decipher <- read_tsv('/data-cbl/frequena_data/from_workstation/daa_decipher/decipher-cnvs-grch37-2020-01-19.txt', skip = 1) %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  select(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>% 
  filter(pathogenicity %in% c('Pathogenic', 
                              'Likely pathogenic',
                              # 'Benign', 
                              # 'Likely benign', 
                              'Unknown') & (! contribution %in% c('None'))) %>%  
  filter(inheritance == 'De novo constitutive') %>%
  filter(genotype == 'Heterozygous') %>%
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>% 
  # filter(inheritance == 'De novo constitutive') %>% 
  filter(variant_class %in% c('Deletion', 'Duplication')) %>%
  mutate(variant_class = tolower(variant_class))
  # distinct(chrom, start, end, .keep_all = TRUE)  # remove all the equal CNVs (exact chrom:start-end)

# gnomAD
tmp_gnomad <-  read_tsv('data/gnomad_v2.1_sv.sites.bed', col_types = list(`#chrom` = col_character())) %>%
  filter(SVTYPE %in% c('DEL', 'DUP')) %>%
  filter(FILTER == 'PASS') %>%
  rename(id = name) %>%
  mutate(id  = str_remove(id, 'gnomAD-SV_v2.1_')) %>%
  mutate(source = 'gnomad_v2.1') %>%
  rename(chrom = `#chrom`) %>%
  mutate(chrom = as.character(chrom)) %>%
  mutate(start = start + 1) %>%
  select(id, chrom, start, end, svtype, -AF) %>% # we can filter by allele frequency
  rename(variant_class = svtype) %>%
  mutate(variant_class = if_else(str_detect(variant_class, 'DEL'), 'Deletion', 'Duplication')) %>%
  mutate(source = 'gnomad_v2.1') %>%
  mutate(variant_class = tolower(variant_class))

# decipher Control
tmp_decipher_control <- read_tsv('https://decipher.sanger.ac.uk/files/downloads/population_cnv.txt.gz', col_names = TRUE) %>%
  rename(id = `#population_cnv_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  mutate(chrom = as.character(chrom)) %>%
  mutate(chrom = if_else(chrom == 23, 'X', chrom)) %>%
  mutate(source = 'decipher_control')

# dgv check gold standard
tmp_dgv <- read_tsv('http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt',
                                  col_types = list(chr = col_character())) %>%
  filter(varianttype == 'CNV') %>%
  filter(!variantsubtype %in% c('deletion', 'duplication', 'loss', 'gain')) %>%
  rename(id = variantaccession, chrom = chr, variant_type = variantsubtype) %>%
  filter(!str_detect(reference, 'gnomAD')) %>%
  mutate(source = 'dgv',
         chrom = as.character(chrom))  %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  mutate(variant_type = if_else(str_detect(variant_type, 'loss'), 'deletion', 'duplication')) %>%
  select(id, chrom, start, end, source, variant_type)

  


input_check_cnv <- tmp_decipher %>%
  bind_rows(tmp_gnomad,
            tmp_decipher_control,
            tmp_dgv)  %>%
  mutate(length_cnv = end - start  + 1) %>%
  mutate(id_tmp = row_number()) %>%
  mutate(pathogenicity = if_else(is.na(pathogenicity), 'benign', 'pathogenic')) %>%
  rename(clinical = pathogenicity)
  



# ------------------------------------------------------------------------------
# FILTER 1 - REMOVE IDENTICAL CNVs (exact chrom:start-end)
# ------------------------------------------------------------------------------

print(glue('Number of rows before identical filter: {nrow(input_check_cnv)}'))
input_check_cnv <- input_check_cnv %>%
  distinct(chrom, start, end, .keep_all = TRUE)

print(glue('Number of rows after identical filter: {nrow(input_check_cnv)}'))



# ------------------------------------------------------------------------------
# FILTER 2 - MATCH BY LENGTH
# ------------------------------------------------------------------------------

bin_length <- 100


tbl_bins <- tibble('start' = seq(1, 1e7, bin_length), 'end' = seq(bin_length, 1e7, bin_length)) %>%
  mutate(chrom = 'chr1') 

tmp_input_check_cnv <- input_check_cnv %>% mutate(chrom = 'chr1', 
                                                  start = length_cnv,
                                                  end = length_cnv) %>%
  select(chrom, start, end, source, id_tmp) %>%
  mutate(source = if_else(source == 'decipher', 'decipher', 'control'))

tmp_input_check_cnv <- tmp_input_check_cnv %>% 
  bed_intersect(tbl_bins %>% mutate(id = paste(start, end, sep = '-')))

bins_with_decipher_and_control <- tmp_input_check_cnv %>% 
  count(id.y, source.x) %>%
  count(id.y) %>%
  filter(n > 1) %>%
  pull(id.y)

selected_decipher_cnvs <- tmp_input_check_cnv %>% 
  filter(id.y %in% bins_with_decipher_and_control) %>%
  filter(source.x == 'decipher') %>%
  pull(id_tmp.x)

# only 35 intervals with > 1 decipher cnv
number_decipher_cnvs <- tmp_input_check_cnv %>% 
  filter(id.y %in% bins_with_decipher_and_control) %>%
  filter(source.x == 'decipher') %>%
  count(id.y) %>%
  rename(n_control = n)

selected_control_cnvs <- c()
for (i in 1:nrow(number_decipher_cnvs)) {
  print(glue('{i} / {nrow(number_decipher_cnvs)}'))
  tmp_control_id_tmp <- tmp_input_check_cnv %>% 
    filter(id.y %in% bins_with_decipher_and_control) %>%
    filter(source.x != 'decipher') %>%
    filter(id.y == number_decipher_cnvs$id.y[i]) %>%
    slice_sample(n = number_decipher_cnvs$n_control[i]) %>%
    pull(id_tmp.x)
  
  selected_control_cnvs <- c(selected_control_cnvs, tmp_control_id_tmp)
  
}

print(glue('Number decipher CNVs: {length(selected_decipher_cnvs)}'))
print(glue('Number control CNVs: {length(selected_control_cnvs)}'))


p_distribution0_1 <- plot_length_distribution(input_check_cnv)

ggsave(glue("cnvscore_results/plot_model/length_distribution_pre.png"), p_distribution0_1, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

  
input_check_cnv <- input_check_cnv %>% 
  filter(id_tmp %in% c(selected_control_cnvs, selected_decipher_cnvs))


# ------------------------------------------------------------------------------
# FILTER 3 - REMOVE RECIPROCAL CNVS (>= 90 %)
# ------------------------------------------------------------------------------

keep_ids <- reciprocal_overlap(input_check_cnv)

input_check_cnv <- input_check_cnv %>% filter(id_tmp %in% keep_ids)
  
# ------------------------------------------------------------------------------
# LENGTH DISTRIBUTION
# ------------------------------------------------------------------------------

p_distribution0_2 <- plot_length_distribution(input_check_cnv)

ggsave(glue("cnvscore_results/plot_model/length_distribution_post.png"), p_distribution0_2, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

# ------------------------------------------------------------------------------
# DISTRIBUTION OF CNVs ACROSS CHROMOSOMES
# ------------------------------------------------------------------------------

# coord_cytobands <- coord_cytobands %>% rename(chrom = Chrom, start = Start, end = End)
  
tmp_distribution <- coord_cytobands %>%
  rowwise() %>%
  mutate(n_cnvs = bed_intersect(tibble('chrom' = chrom, 'start' = start, 'end' = end), input_check_cnv) %>% nrow()) %>%
  ungroup()

tmp_distribution %>% filter(n_cnvs != 0) # 115 CNVs


p_distribution1 <- tmp_distribution %>% 
  group_by(chrom) %>%
  summarise(total_cnvs = sum(n_cnvs), .groups = 'keep') %>%
  ggplot(aes(reorder(chrom, -total_cnvs), total_cnvs)) +
  geom_col() +
  theme_minimal() +
  xlab('Chromosome') +
  ylab('Count CNVs') +
  ggtitle('Number of CNVs (DECIPHER) per chromosome')

total_cnvs <- tmp_distribution %>% pull(n_cnvs) %>% sum()

p_distribution2 <- tmp_distribution %>% 
  mutate(perc = n_cnvs / total_cnvs) %>%
  group_by(chrom) %>%
  arrange(desc(perc)) %>% 
  ungroup() %>%
  ggplot(aes(reorder(Name, -perc), perc)) +
    geom_col(aes(fill = chrom), color = 'black', show.legend = FALSE) +
  scale_y_continuous(label = percent) +
  facet_wrap(~ chrom, scales = 'free') +
    theme_minimal() +
  ggtitle('Distribution of CNVs (DECIPHER) across genome') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))


p_distribution3 <- tmp_distribution %>% 
  mutate(perc = n_cnvs / total_cnvs) %>%
  group_by(chrom) %>%
   mutate(perc = n_cnvs / sum(n_cnvs)) %>%
  arrange(desc(perc)) %>% 
  ungroup() %>%
  ggplot(aes(reorder(Name, -perc), perc)) +
  geom_col(aes(fill = chrom), color = 'black', show.legend = FALSE) +
  scale_y_continuous(label = percent) +
  facet_wrap(~ chrom, scales = 'free') +
  theme_minimal() +
  ggtitle('Distribution of CNVs (DECIPHER) across genome (percentages by chromosome)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

p_distribution4 <- tmp_distribution %>% 
  ggplot(aes(reorder(Name, -n_cnvs), n_cnvs)) +
  geom_col(aes(fill = chrom), color = 'black', show.legend = FALSE) +
  facet_wrap(~ chrom, scales = 'free') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

p_distribution5 <- tmp_distribution %>%
  left_join(coord_chrom_hg19 %>% rename(length_chrom = length), by = 'chrom') %>%
  group_by(chrom) %>%
  mutate(total_cnvs = sum(n_cnvs)) %>%
  ungroup() %>%
  select(chrom, length_chrom, total_cnvs) %>% distinct() %>%
  ggplot(aes(length_chrom, total_cnvs)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  geom_label(aes(label = chrom)) +
  labs(x = 'Chromosome length', y = 'Number of CNVs') +
  theme_minimal()
  


ggsave(glue("cnvscore_results/plot_model/count_chromosome.png"), p_distribution1, width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/plot_model/count_per_chromosome.png"), p_distribution4, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("cnvscore_results/plot_model/perc_across_genome.png"), p_distribution2, width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/plot_model/perc_across_chromosome.png"), p_distribution3, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("cnvscore_results/plot_model/ratio_length_n_cnvs.png"), p_distribution5, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

# ------------------------------------------------------------------------------
# RANDOM GENERATION - PARALLELIZATION TEST
# ------------------------------------------------------------------------------

gen_random(input_check_cnv$chrom[1], 
           input_check_cnv$start[1], 
           input_check_cnv$end[1],
           input_check_cnv$id_tmp[1],
           n_rep = 9,
           input_n_try = 400
)

# ------------------------------------------------------------------------------
# RANDOM GENERATION - PARALLELIZATION
# ------------------------------------------------------------------------------

plan("multiprocess", workers = 10)


tic()

output_list <- future_pmap(list(input_check_cnv$chrom, 
                                input_check_cnv$start, 
                                input_check_cnv$end,
                                input_check_cnv$id_tmp, 
                                n_rep = 9,
                                input_n_try = 400), 
                                gen_random, 
                                .progress = TRUE)

random_cnvs <- bind_rows(lapply(output_list, as.data.frame.list)) %>% as_tibble()

toc()

random_cnvs %>% count(chrom) %>% filter(chrom == 'error') # 124 CNVs error

# ------------------------------------------------------------------------------
# RANDOM GENERATION - REMOVING CNVS THAT DO NOT FIT THE RULES
# ------------------------------------------------------------------------------

remove_cnvs_error <- random_cnvs %>% 
  filter(chrom == 'error') %>% 
  pull(id_tmp) %>%
  as.character()

input_check_cnv <- input_check_cnv %>% filter(! id_tmp %in% remove_cnvs_error)

random_cnvs <- random_cnvs %>% 
  filter(chrom != 'error') %>%
  mutate(inheritance = NA, variant_class = NA )

# ------------------------------------------------------------------------------
# RANDOM GENERATION - BIND DECIPHER AND RANDOM CNVs
# ------------------------------------------------------------------------------

input_check_cnv <- input_check_cnv %>% 
  mutate(id_tmp = as.character(id_tmp)) %>% 
  bind_rows(random_cnvs %>% mutate(id_tmp = as.character(id_tmp)))



# ------------------------------------------------------------------------------
# ANNOTATION - FUNCTION 
# ------------------------------------------------------------------------------


check_cnv <- function(input_id, input_clinical, input_variant, input_inheritance,
                      input_chrom, input_start, input_end) {
  
  
  # input_id <- input_check_cnv %>% slice(1) %>% pull(id_tmp)
  # input_clinical <- input_check_cnv %>% slice(1) %>% pull(pathogenicity)
  # input_variant <- input_check_cnv %>% slice(1) %>% pull(variant_class)
  # input_inheritance <- input_check_cnv %>% slice(1) %>% pull(inheritance)
  # input_chrom <- input_check_cnv %>% slice(1) %>% pull(chrom)
  # input_start <- input_check_cnv %>% slice(1) %>% pull(start)
  # input_end <- input_check_cnv %>% slice(1) %>% pull(end)
  
  
  # input_id <- training_tbl %>% filter(id == 1251) %>% pull(id)
  # input_clinical <- training_tbl %>% filter(id == 1251)  %>% pull(clinical)
  # input_variant <- training_tbl %>% filter(id == 1251)  %>% pull(type_variant)
  # input_inheritance <- training_tbl %>% filter(id == 1251)  %>% pull(type_inheritance)
  # input_chrom <- training_tbl %>% filter(id == 1251)  %>% pull(chrom)
  # input_start <- training_tbl %>% filter(id == 1251)  %>% pull(start)
  # input_end <- training_tbl %>% filter(id == 1251)  %>% pull(end)
  
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
    filter(chrom == chrom_tmp) %>%
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
  
  
  # Nº proteins / Nº high-confidence protein-protein interactions
  
  
  n_interactions <- interactions_db %>% 
    filter(from %in% vector_genes & to %in% vector_genes) %>% nrow()
  
  result_density_interactions <- if_else(n_interactions == 0, 0, ( length(vector_genes) / n_interactions))
  
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
    'density_int' = result_density_interactions,
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
# ANNOTATION - PARALLELIZATION TEST
# ------------------------------------------------------------------------------

check_cnv(input_check_cnv$id_tmp[1],
          input_check_cnv$pathogenicity[1],
          input_check_cnv$variant_class[1],
          input_check_cnv$inheritance[1],
          input_check_cnv$chrom[1],
          input_check_cnv$start[1],
          input_check_cnv$end[1]
)

# ------------------------------------------------------------------------------
# ANNOTATION - PARALLELIZATION
# ------------------------------------------------------------------------------

plan("multiprocess", workers = 40)


tic()

output_list <- future_pmap(list(input_check_cnv$id_tmp, 
                                input_check_cnv$pathogenicity, 
                                input_check_cnv$variant_class,
                                input_check_cnv$inheritance,
                                input_check_cnv$chrom, 
                                input_check_cnv$start, 
                                input_check_cnv$end), 
                           check_cnv, .progress = TRUE)

output_df <- bind_rows(lapply(output_list, as.data.frame.list)) %>% as_tibble()

toc()

# ------------------------------------------------------------------------------
# RMARKDOWN VISUALIZATION
# ------------------------------------------------------------------------------

output_df %>% write_tsv('output_df.tsv')

rmarkdown::render("cnvscore_results/output_annotation_cnvs.Rmd", 'html_document')


# ------------------------------------------------------------------------------
# MODEL FORMULE
# ------------------------------------------------------------------------------

formule_models <- clinical ~  n_cnv_syndromes + n_genes_hpo + max_ccr + max_hi + max_pli + embryo_mouse +
  essent_cl + max_cpg + n_target_drugs + n_prot_complex + patho_cnv + nonpatho_cnv + n_systems + density_int + 
  pubmed_del


# ------------------------------------------------------------------------------
# RANDOM GENERATION - REMOVE CNVs WITH NO 10 GROUP CNVs
# ------------------------------------------------------------------------------

model_df <- output_df %>% 
  mutate(clinical = as.character(clinical)) %>%
  mutate(clinical = if_else(is.na(clinical), 'unlabeled', clinical)) %>%
  separate(id, sep = '_', c('id', 'random_n'))

remove_ids <- model_df %>% group_by(id) %>% count() %>% filter(n != 10) %>% pull(id)

length(remove_ids) == 0
# model_df <- model_df %>% filter(! id %in% remove_ids)


# ------------------------------------------------------------------------------
#  RANDOM GENERATION - LABEL PROPERLY PATHOGENIC/NON-PATHOGENIC UNLABELED CNVs
# ------------------------------------------------------------------------------

id_pathogenics <- model_df %>% filter(clinical == 'Pathogenic') %>% pull(id)
id_benigns <- model_df %>% filter(clinical == 'Benign') %>% pull(id)


model_df <- model_df %>% 
  mutate(clinical = if_else(id %in% id_pathogenics & clinical == 'unlabeled', paste(clinical, 'patho', sep = '_'), clinical)) %>%
  mutate(clinical = if_else(id %in% id_benigns & clinical == 'unlabeled', paste(clinical, 'benign', sep = '_'), clinical)) %>%
  mutate(clinical = if_else(clinical == 'unlabeled', 'unlabeled_likely_benign', clinical))

# ------------------------------------------------------------------------------
#  RANDOM GENERATION - SELECTION PATHOGENIC vs. UNLABELED OR BENIGN vs. UNLABELED
# ------------------------------------------------------------------------------

# target class

target_class <- 'decipher'

# benign CNVs
first_split <- initial_split(data = model_df %>%
                               mutate(clinical = if_else(clinical == 'Likely benign', 'Benign', clinical)) %>%
                               filter(clinical %in% c('Benign', 'unlabeled_benign')) %>%
                               mutate(clinical = as.factor(clinical)), strata = clinical)
# pathogenic CNVs
first_split <- initial_split(data = model_df %>% 
                               filter(clinical %in% c('Pathogenic', 'unlabeled_patho')) %>%
                               mutate(clinical = as.factor(clinical)), strata = clinical)

training_tbl <- first_split %>% training()
test_tbl <- first_split %>% testing()


# ------------------------------------------------------------------------------
# SET TARGET CLASS NAME
# ------------------------------------------------------------------------------

target_class <- '.pred_pathogenic'

# ------------------------------------------------------------------------------
# TRAINING AND TEST DATASETS
# ------------------------------------------------------------------------------

first_split_deletion <- initial_split(output_df %>% filter(type_variant == 'deletion'), strata = clinical)

training_tbl_deletion <- first_split_deletion %>% training()
test_tbl_deletion <- first_split_deletion %>% testing()

first_split_duplication <- initial_split(output_df %>% filter(type_variant == 'duplication'), strata = clinical)

training_tbl_duplication <- first_split_duplication %>% training()
test_tbl_duplication <- first_split_duplication %>% testing()

# ------------------------------------------------------------------------------
# CREATION OF STRATIFICATION VARIABLE (POSITIVE/UNLABELED AND CHROMOSOME)
# ------------------------------------------------------------------------------

# training_tbl <- training_tbl %>%
#   mutate(strata_variable = paste(chrom, clinical, sep = '_'))



# ------------------------------------------------------------------------------
# TRAINING DATA - GENERATION OF CROSS-VALIDATION
# ------------------------------------------------------------------------------

training_cv_deletion <- vfold_cv(training_tbl_deletion, prop = 0.9, strata = clinical)
training_cv_duplication <- vfold_cv(training_tbl_duplication, prop = 0.9, strata = clinical)


# ------------------------------------------------------------------------------
# TRAINING DATA - TRAINING MODELS
# ------------------------------------------------------------------------------


logistic_cv_deletion <- fit_resamples(logistic_reg(mode = 'classification') %>% 
                               set_engine(engine = 'glm'),
                             formule_models,
                             control = control_resamples(save_pred = TRUE),
                             training_cv_deletion)


xgboost_cv_deletion <- fit_resamples(boost_tree(mode = 'classification') %>% 
                              set_engine(engine = 'xgboost'),
                            formule_models,
                            control = control_resamples(save_pred = TRUE),
                            training_cv_deletion)

forest_cv_deletion <- fit_resamples(rand_forest(mode = 'classification') %>% 
                              set_engine(engine = 'ranger'),
                            formule_models,
                            control = control_resamples(save_pred = TRUE),
                            training_cv_deletion)

rulefit_cv_deletion <- rulefit_cv(training_cv_deletion)



logistic_cv_duplication <- fit_resamples(logistic_reg(mode = 'classification') %>% 
                                           set_engine(engine = 'glm'),
                                         formule_models,
                                         control = control_resamples(save_pred = TRUE),
                                         training_cv_duplication)

xgboost_cv_duplication <- fit_resamples(boost_tree(mode = 'classification') %>% 
                                          set_engine(engine = 'xgboost'),
                                        formule_models,
                                        control = control_resamples(save_pred = TRUE),
                                        training_cv_duplication)

forest_cv_duplication <- fit_resamples(rand_forest(mode = 'classification') %>% 
                                          set_engine(engine = 'ranger'),
                                        formule_models,
                                        control = control_resamples(save_pred = TRUE),
                                        training_cv_duplication)

rulefit_cv_duplication <- rulefit_cv(training_cv_duplication)

# ------------------------------------------------------------------------------
# TRAINING DATA - 10CV - ROC AND PR CURVES
# ------------------------------------------------------------------------------

logistic_roc_pr_deletion <- evaluate_model(logistic_cv_deletion, model_name = 'Logistic regression', cv = TRUE, cum_dist = FALSE)
xgboost_roc_pr_deletion <- evaluate_model(xgboost_cv_deletion, model_name = 'XGboost', cv = TRUE, cum_dist = FALSE)
randomf_roc_pr_deletion <- evaluate_model(forest_cv_deletion, model_name = 'Random forest', cv = TRUE, cum_dist = FALSE)
rulefit_roc_pr_deletion <- from_cv_rulefit(rulefit_cv_deletion, model_name = 'RuleFit') 

logistic_roc_pr_deletion[[1]] +  logistic_roc_pr_deletion[[2]]
xgboost_roc_pr_deletion[[1]] +  xgboost_roc_pr_deletion[[2]]
randomf_roc_pr_deletion[[1]] +  randomf_roc_pr_deletion[[2]]
rulefit_roc_pr_deletion[[1]] +  rulefit_roc_pr_deletion[[2]]

ggsave(glue("cnvscore_results/evaluation_model/deletion/logistic_{target_class}_10cv.png"), logistic_roc_pr_deletion[[1]] +  logistic_roc_pr_deletion[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/evaluation_model/deletion/xgboost_{target_class}_10cv.png"), xgboost_roc_pr_deletion[[1]] +  xgboost_roc_pr_deletion[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/evaluation_model/deletion/forest_{target_class}_10cv.png"), randomf_roc_pr_deletion[[1]] +  randomf_roc_pr_deletion[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/evaluation_model/deletion/rulefit_{target_class}_10cv.png"), rulefit_roc_pr_deletion[[1]] +  rulefit_roc_pr_deletion[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')


logistic_roc_pr_duplication <- evaluate_model(logistic_cv_duplication, model_name = 'Logistic regression', cv = TRUE, cum_dist = FALSE)
xgboost_roc_pr_duplication <- evaluate_model(xgboost_cv_duplication, model_name = 'XGboost', cv = TRUE, cum_dist = FALSE)
randomf_roc_pr_duplication <- evaluate_model(forest_cv_duplication, model_name = 'Random forest', cv = TRUE, cum_dist = FALSE)
rulefit_roc_pr_duplication <- from_cv_rulefit(rulefit_cv_duplication, model_name = 'RuleFit') 

logistic_roc_pr_duplication[[1]] +  logistic_roc_pr_duplication[[2]]
xgboost_roc_pr_duplication[[1]] +  xgboost_roc_pr_duplication[[2]]
randomf_roc_pr_duplication[[1]] +  randomf_roc_pr_duplication[[2]]
rulefit_roc_pr_duplication[[1]] +  rulefit_roc_pr_duplication[[2]]

ggsave(glue("cnvscore_results/evaluation_model/duplication/logistic_{target_class}_10cv.png"), logistic_roc_pr_duplication[[1]] +  logistic_roc_pr_duplication[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/evaluation_model/duplication/xgboost_{target_class}_10cv.png"), xgboost_roc_pr_duplication[[1]] +  xgboost_roc_pr_duplication[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/evaluation_model/duplication/forest_{target_class}_10cv.png"), randomf_roc_pr_duplication[[1]] +  randomf_roc_pr_duplication[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/evaluation_model/duplication/rulefit_{target_class}_10cv.png"), rulefit_roc_pr_duplication[[1]] +  rulefit_roc_pr_duplication[[2]], width = 17, height = 9.6, dpi = 300, units = "in", device='png')


# ------------------------------------------------------------------------------
# TRAINING DATA - 10CV - CUMULATIVE FREQUENCY
# ------------------------------------------------------------------------------

select_filter <- if_else(str_detect(target_class, 'pathogenic'), 'pathogenic', 'benign')

logistic_cv_deletion_cumfreq <- evaluate_model(logistic_cv_deletion, model_name = 'Logistic regression', cv = TRUE, cum_dist = TRUE)
forest_cv_deletion_cumfreq <- evaluate_model(forest_cv_deletion, model_name = 'Random forest', cv = TRUE, cum_dist = TRUE)
xgboost_cv_deletion_cumfreq <- evaluate_model(xgboost_cv_deletion, model_name = 'XGboost', cv = TRUE, cum_dist = TRUE)
rulefit_cv_deletion_cumfreq <- evaluate_model(rulefit_cv_deletion, model_name = 'RuleFit', rulefit = TRUE, cv = TRUE, cum_dist = TRUE)

logistic_cv_duplication_cumfreq <- evaluate_model(logistic_cv_duplication, model_name = 'Logistic regression', cv = TRUE, cum_dist = TRUE)
forest_cv_duplication_cumfreq <- evaluate_model(forest_cv_duplication, model_name = 'Random forest', cv = TRUE, cum_dist = TRUE)
xgboost_cv_duplication_cumfreq <- evaluate_model(xgboost_cv_duplication, model_name = 'XGboost', cv = TRUE, cum_dist = TRUE)
rulefit_cv_duplication_cumfreq <- evaluate_model(rulefit_cv_duplication, model_name = 'RuleFit', rulefit = TRUE, cv = TRUE, cum_dist = TRUE)


# ------------------------------------------------------------------------------
# TRAINING DATA - 10CV - CUMULATIVE FREQUENCY
# ------------------------------------------------------------------------------

# random_line <- tibble(id = 'Random', cum_perc = seq(0, 100, 10), model = 'Random')

perc_10_cv_total_deletion <- bind_rows(logistic_cv_deletion_cumfreq,
                              forest_cv_deletion_cumfreq,
                              xgboost_cv_deletion_cumfreq,
                              rulefit_cv_deletion_cumfreq) %>%
                                ggplot(aes(p_patho, cum_perc)) +
                                geom_path(aes(p_patho, cum_perc, group = id, color = model), show.legend = FALSE) +
                                geom_point(aes(p_patho, cum_perc)) +
                                scale_x_continuous(breaks = scales::pretty_breaks(10)) +
                                scale_y_continuous(labels = scales::percent) +
                                facet_wrap(~ model) +
                                theme_bw() +
                                labs(x = 'Percentile', y = 'Cumulative Frequency')



ggsave(glue("cnvscore_results/plot_results/10cv_cum_freq_{target_class}_all.png"), perc_10_cv_total_deletion, width = 17, height = 9.6, dpi = 300, units = "in", device='png')



perc_10_cv_total_duplication <- bind_rows(logistic_cv_duplication_cumfreq,
                                       forest_cv_duplication_cumfreq,
                                       xgboost_cv_duplication_cumfreq,
                                       rulefit_cv_duplication_cumfreq) %>%
  ggplot(aes(p_patho, cum_perc)) +
  geom_path(aes(p_patho, cum_perc, group = id, color = model), show.legend = FALSE) +
  geom_point(aes(p_patho, cum_perc)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = 'Percentile', y = 'Cumulative Frequency')

ggsave(glue("cnvscore_results/plot_results/10cv_cum_freq_{target_class}_all.png"), perc_10_cv_total_duplication, width = 17, height = 9.6, dpi = 300, units = "in", device='png')



# ------------------------------------------------------------------------------
# TEST DATA - MODEL TRAINING
# ------------------------------------------------------------------------------

logistic_model_deletion <- logistic_reg() %>%
  set_mode('classification') %>%
  set_engine('glm') %>%
  fit(formule_models, data = training_tbl_deletion)

xgboost_model_deletion <- boost_tree() %>%
  set_mode('classification') %>%
  set_engine('xgboost') %>%
  fit(formule_models, data = training_tbl_deletion)

forest_model_deletion <- rand_forest() %>%
  set_mode('classification') %>%
  set_engine('ranger') %>%
  fit(formule_models, data = training_tbl_deletion)


rulefit_model_deletion <- xrftest::xrf(formule_models,
                      data = training_tbl_deletion,
                      xgb_control = list(nrounds = 20, scale_pos_weight = 9),
                      family = "binomial")


bayesian_model_deletion <- rstanarm::stan_glm(formula = clinical ~ .,
                   family = 'binomial',
                   data = rulefit_model_deletion$full_data %>%
                     select_if(~  sum(.x == 1) > 2) %>%
                     select(-chrom),
                   cores = 4,
                   iter = 2000,
                   chains = 4,
                   algorithm = 'sampling', # variational inference algorithms
                   QR = FALSE,
                   prior = hs(),
                   prior_intercept = student_t())

  
logistic_model_duplication <- logistic_reg() %>%
  set_mode('classification') %>%
  set_engine('glm') %>%
  fit(formule_models, data = training_tbl_duplication)

xgboost_model_duplication <- boost_tree() %>%
  set_mode('classification') %>%
  set_engine('xgboost') %>%
  fit(formule_models, data = training_tbl_duplication)

forest_model_duplication <- rand_forest() %>%
  set_mode('classification') %>%
  set_engine('ranger') %>%
  fit(formule_models, data = training_tbl_duplication)


rulefit_model_duplication <- xrftest::xrf(formule_models,
                                       data = training_tbl_duplication,
                                       xgb_control = list(nrounds = 20, scale_pos_weight = 9),
                                       family = "binomial")


bayesian_model_duplication <- rstanarm::stan_glm(formula = formule_models,
                                              family = 'binomial',
                                              data = rulefit_model_duplication$full_data,
                                              cores = 2,
                                              iter = 1000,
                                              chains = 2,
                                              algorithm = 'sampling', # variational inference algorithms
                                              QR = TRUE,
                                              prior = hs(),
                                              prior_intercept = student_t())

# ------------------------------------------------------------------------------
# TEST DATA - CUMULATIVE FREQUENCY
# ------------------------------------------------------------------------------


prev_logistic_deletion <- evaluate_model(logistic_model_deletion, test_object = test_tbl_deletion, model_name = 'Logistic regression', cv = FALSE, 
                                cum_dist = TRUE)
prev_xgboost_deletion <- evaluate_model(xgboost_model_deletion, test_object = test_tbl_deletion, model_name = 'XGBoost', cv = FALSE, 
                                cum_dist = TRUE)
prev_forest_deletion <- evaluate_model(forest_model, test_object = test_tbl_deletion, model_name = 'Random Forest', cv = FALSE, 
                                cum_dist = TRUE)

prev_rulefit_deletion <- evaluate_model(rulefit_model_deletion, test_object = test_tbl_deletion, rulefit = TRUE, model_name = 'RuleFit', cv = FALSE, 
                              cum_dist = TRUE)

prev_bayesian_deletion <- evaluate_model(bayesian_model_deletion, test_object = test_tbl_deletion, bay = TRUE, model_name = 'Bayesian', cv = FALSE, 
                               cum_dist = TRUE)


prev_logistic_duplication <- evaluate_model(logistic_model_duplication, test_object = test_tbl_duplication, model_name = 'Logistic regression', cv = FALSE, 
                                         cum_dist = TRUE)
prev_xgboost_duplication <- evaluate_model(xgboost_model_duplication, test_object = test_tbl_duplication, model_name = 'XGBoost', cv = FALSE, 
                                        cum_dist = TRUE)
prev_forest_duplication <- evaluate_model(forest_model_duplication, test_object = test_tbl_duplication, model_name = 'Random Forest', cv = FALSE, 
                                       cum_dist = TRUE)

prev_rulefit_duplication <- evaluate_model(rulefit_model_duplication, test_object = test_tbl_duplication, rulefit = TRUE, model_name = 'RuleFit', cv = FALSE, 
                                        cum_dist = TRUE)

prev_bayesian_duplication <- evaluate_model(bayesian_model_duplication, test_object = test_tbl_duplication, bay = TRUE, model_name = 'Bayesian', cv = FALSE, 
                                         cum_dist = TRUE)



# ------------------------------------------------------------------------------
# TEST DATA - PLOT CUMULATIVE FREQUENCY
# ------------------------------------------------------------------------------
 
 target_n_deletion <- test_tbl_deletion %>% count(clinical) %>% filter(clinical == select_filter) %>% pull(n)
 target_n_duplication <- test_tbl_duplication %>% count(clinical) %>% filter(clinical == select_filter) %>% pull(n)

 
 cumulative_test_deletion <- cum_plot(by_percentage = TRUE, 
                select_filter = select_filter, 
                target_n = target_n_deletion,
                total_n = test_tbl_deletion %>% nrow(),
                prev_logistic_deletion, prev_xgboost_deletion,
                prev_forest_deletion, prev_rulefit_deletion, prev_bayesian_deletion)
 
 cumulative_test_duplication <- cum_plot(by_percentage = TRUE, 
                                      select_filter = select_filter, 
                                      target_n = target_n_duplication,
                                      total_n = test_tbl_duplication %>% nrow(),
                                      prev_logistic_deletion, 
                                      prev_xgboost_duplication, prev_forest_duplication, 
                                      prev_rulefit_duplication, prev_bayesian_duplication)

ggsave("cnvscore_results/plot_model/test_cum_frequency.png", cumulative_test, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

# ------------------------------------------------------------------------------
# TEST DATA - ROC AND PR CURVES
# ------------------------------------------------------------------------------

logistic_result_deletion <- evaluate_model(logistic_model_deletion, test_object = test_tbl_deletion, model_name = 'Logistic regression', cv = FALSE, 
                                cum_dist = FALSE)
xgboost_result_deletion <- evaluate_model(xgboost_model_deletion, test_object = test_tbl_deletion, model_name = 'XGBoost', cv = FALSE, 
                               cum_dist = FALSE)
forest_result_deletion <- evaluate_model(forest_model_deletion, test_object = test_tbl_deletion, model_name = 'Random Forest', cv = FALSE, 
                              cum_dist = FALSE)

rulefit_result_deletion <- evaluate_model(rulefit_model_deletion, test_object = test_tbl_deletion, rulefit = TRUE, model_name = 'RuleFit', cv = FALSE, 
                                cum_dist = FALSE)

bayesian_result_deletion <- evaluate_model(bayesian_model_deletion, test_object = test_tbl_deletion, bay = TRUE, model_name = 'Bayesian', cv = FALSE, 
                                 cum_dist = FALSE)



roc_auc_deletion <- bind_rows(logistic_result_deletion[[1]], xgboost_result_deletion[[1]],
                              forest_result_deletion[[1]], rulefit_result_deletion[[1]], bayesian_result_deletion[[1]]) %>%
  ggplot(aes(1-specificity, sensitivity)) +
  geom_path(aes(group = model, color = model), size = 2.5, show.legend = FALSE) +
  theme_minimal() +
  theme(plot.title = element_text(size=10, face="bold"),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_text(size=19, face="bold"),
        axis.text.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14, face="bold")) +
  ggtitle(glue('random_forest = ', {forest_result_deletion[[3]]}, ' - ', 
               'logistic = ', {logistic_result_deletion[[3]]}, ' - ', 
               'xgboost = ', {xgboost_result_deletion[[3]]}, ' - ',
               'rulefit = ', {rulefit_result_deletion[[3]]}, ' - ',
               'bayesian = ', {bayesian_result_deletion[[3]]}))


pr_auc_deletion <-  bind_rows(logistic_result_deletion[[2]], xgboost_result_deletion[[2]],
                              forest_result_deletion[[2]], rulefit_result_deletion[[2]], bayesian_result_deletion[[2]]) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = model, color = model), size = 2.5, show.legend = TRUE) +
  theme_minimal() +
  theme(plot.title = element_text(size=10, face="bold"),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_text(size=19, face="bold"),
        axis.text.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14, face="bold")) + 
  ggtitle(glue('random_forest = ', {forest_result_deletion[[4]]}, ' - ',
               'logistic = ', {logistic_result_deletion[[4]]}, ' - ',
               'xgboost = ', {xgboost_result_deletion[[4]]}, ' - ',
               'rulefit = ', {rulefit_result_deletion[[4]]}, ' - ',
               'bayesian = ', {bayesian_result_deletion[[4]]}
  ))

ggsave("cnvscore_results/plot_model/test_roc_pr_auc.png", roc_auc_deletion + pr_auc_deletion, width = 17, height = 9.6, dpi = 300, units = "in", device='png')



logistic_result_duplication <- evaluate_model(logistic_model_duplication, test_object = test_tbl_duplication, model_name = 'Logistic regression', cv = FALSE, 
                                           cum_dist = FALSE)
xgboost_result_duplication <- evaluate_model(xgboost_model_duplication, test_object = test_tbl_duplication, model_name = 'XGBoost', cv = FALSE, 
                                          cum_dist = FALSE)
forest_result_duplication <- evaluate_model(forest_model_duplication, test_object = test_tbl_duplication, model_name = 'Random Forest', cv = FALSE, 
                                         cum_dist = FALSE)

rulefit_result_duplication <- evaluate_model(rulefit_model_duplication, test_object = test_tbl_duplication, rulefit = TRUE, model_name = 'RuleFit', cv = FALSE, 
                                          cum_dist = FALSE)

bayesian_result_duplication <- evaluate_model(bayesian_model_duplication, test_object = test_tbl, bay = TRUE, model_name = 'Bayesian', cv = FALSE, 
                                           cum_dist = FALSE)



roc_auc_duplication <- bind_rows(logistic_result_duplication[[1]], xgboost_result_duplication[[1]],
                              forest_result_duplication[[1]], rulefit_result_duplication[[1]], bayesian_result_duplication[[1]]) %>%
  ggplot(aes(1-specificity, sensitivity)) +
  geom_path(aes(group = model, color = model), size = 2.5, show.legend = FALSE) +
  theme_minimal() +
  theme(plot.title = element_text(size=10, face="bold"),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_text(size=19, face="bold"),
        axis.text.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14, face="bold")) +
  ggtitle(glue('random_forest = ', {forest_result_duplication[[3]]}, ' - ', 
               'logistic = ', {logistic_result_duplication[[3]]}, ' - ', 
               'xgboost = ', {xgboost_result_duplication[[3]]}, ' - ',
               'rulefit = ', {rulefit_result_duplication[[3]]}, ' - ',
               'bayesian = ', {bayesian_result_duplication[[3]]}))


pr_auc_duplication <-  bind_rows(logistic_result_duplication[[2]], xgboost_result_duplication[[2]],
                              forest_result_duplication[[2]], rulefit_result_duplication[[2]], bayesian_result_duplication[[2]]) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = model, color = model), size = 2.5, show.legend = TRUE) +
  theme_minimal() +
  theme(plot.title = element_text(size=10, face="bold"),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_text(size=19, face="bold"),
        axis.text.x = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14, face="bold")) + 
  ggtitle(glue('random_forest = ', {forest_result_duplication[[4]]}, ' - ',
               'logistic = ', {logistic_result_duplication[[4]]}, ' - ',
               'xgboost = ', {xgboost_result_duplication[[4]]}, ' - ',
               'rulefit = ', {rulefit_result_duplication[[4]]}, ' - ',
               'bayesian = ', {bayesian_result_duplication[[4]]}
  ))

ggsave("cnvscore_results/plot_model/test_roc_pr_auc.png", roc_auc_duplication + pr_auc_duplication, width = 17, height = 9.6, dpi = 300, units = "in", device='png')




# ------------------------------------------------------------------------------
# INDEPENDENT DATASET (CNVs - ClinVar variants)
# ------------------------------------------------------------------------------

url <- 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz'
download.file(url, destfile = basename(url))

clinvar_tab_raw <- read_tsv('variant_summary.txt.gz')

clinvar_cnvs_hg37 <- clinvar_tab_raw %>% 
  filter(Type %in% c('copy number gain', 'copy number loss', 'Deletion', 'Duplication')) %>% 
  filter(ClinicalSignificance %in% c('Benign', 'Pathogenic')) %>% 
  filter(Assembly %in% c('GRCh37')) %>% 
  filter(!str_detect(Origin, 'somatic')) %>% 
  filter(str_detect(ReviewStatus, 'submitter')) %>% 
  select(Chromosome, Start, Stop, Type, ClinicalSignificance) %>% 
  rename(chrom = Chromosome, start = Start, end = Stop, variant_class = Type, 
         clinical = ClinicalSignificance) %>%
  # mutate(chrom = paste0('chr', chrom)) %>%
  mutate(id_tmp = paste0(row_number(), '_clinvar')) %>%
  mutate(inheritance = 'no_somatic') %>%
  distinct(chrom, start, end, .keep_all = TRUE) %>%
  mutate(length_cnv = end - start + 1) %>% 
  filter(length_cnv >= 50) %>%
  select(-length_cnv) %>%
  mutate(clinical = tolower(clinical)) %>% 
  mutate(variant_class = tolower(variant_class))

file.remove('variant_summary.txt.gz')

keep_ids_clinvar <- reciprocal_overlap(bind_rows(clinvar_cnvs_hg37 %>% select(id_tmp, chrom, start, end),
                                                 input_check_cnv %>% select(id_tmp, chrom, start, end)))

keep_ids <- reciprocal_overlap(clinvar_cnvs_hg37)

input_check_cnv <- input_check_cnv %>% filter(id_tmp %in% keep_ids)



plan("multiprocess", workers = 40)

check_cnv(clinvar_cnvs_hg37$id_tmp[1],
          clinvar_cnvs_hg37$clinical[1],
          clinvar_cnvs_hg37$variant_class[1],
          clinvar_cnvs_hg37$inheritance[1],
          clinvar_cnvs_hg37$chrom[1],
          clinvar_cnvs_hg37$start[1],
          clinvar_cnvs_hg37$end[1]
)

tic()

output_list <- future_pmap(list(clinvar_cnvs_hg37$id_tmp, 
                                clinvar_cnvs_hg37$clinical, 
                                clinvar_cnvs_hg37$variant_class,
                                clinvar_cnvs_hg37$inheritance,
                                clinvar_cnvs_hg37$chrom, 
                                clinvar_cnvs_hg37$start, 
                                clinvar_cnvs_hg37$end), 
                           check_cnv, .progress = TRUE)

clinvar_ind <- bind_rows(lapply(output_list, as.data.frame.list)) %>% as_tibble()

toc()

clinvar_ind_deletion <- clinvar_ind %>% filter(type_variant == 'deletion')
clinvar_ind_duplication <- clinvar_ind %>% filter(type_variant == 'duplication')

# ------------------------------------------------------------------------------
# INDEPENDENT - ROC CURVES
# ------------------------------------------------------------------------------

logistic_clinvar_deletion <- evaluate_model(logistic_model_deletion, test_object = clinvar_ind_deletion, model_name = 'Logistic regression', cv = FALSE, 
               cum_dist = FALSE)

xgboost_clinvar_deletion <- evaluate_model(xgboost_model_deletion, test_object = clinvar_ind_deletion, model_name = 'XGBoost', cv = FALSE, 
               cum_dist = FALSE)

forest_clinvar_deletion <- evaluate_model(forest_model_deletion, test_object = clinvar_ind_deletion, model_name = 'Random Forest', cv = FALSE, 
               cum_dist = FALSE)

rulefit_clinvar_deletion <- evaluate_model(rulefit_model_deletion, rulefit = TRUE, test_object = clinvar_ind_deletion, model_name = 'RuleFit', cv = FALSE, 
                                 cum_dist = FALSE)

bayesian_clinvar_deletion <- evaluate_model(bayesian_model_deletion, bay = TRUE, test_object = clinvar_ind_deletion, model_name = 'Bayesian', cv = FALSE, 
                                  cum_dist = FALSE)

clinvar_roc_deletion <- bind_rows(logistic_clinvar_deletion[[1]], xgboost_clinvar_deletion[[1]], 
                         forest_clinvar_deletion[[1]], rulefit_clinvar_deletion[[1]], bayesian_clinvar_deletion[[1]]) %>%
 
  ggplot(aes(1 - specificity, sensitivity)) +
    geom_path(aes(color = model)) +
    theme_bw() +
  ggtitle(glue('random_forest = ', {round(forest_clinvar_deletion[[3]], 2)}, ' - ', 
               'logistic = ', {round(logistic_clinvar_deletion[[3]], 2)}, ' - ', 
               'xgboost = ', {round(xgboost_clinvar_deletion[[3]], 2)}, ' - ', 
               'rulefit = ', {round(rulefit_clinvar_deletion[[3]], 2)}, ' - ',
               'bayesian = ', {round(bayesian_clinvar_deletion[[3]], 2)}))

clinvar_roc_pr_deletion <- bind_rows(logistic_clinvar_deletion[[2]], 
                            xgboost_clinvar_deletion[[2]], 
                            forest_clinvar_deletion[[2]], 
                            rulefit_clinvar_deletion[[2]], 
                            bayesian_clinvar_deletion[[2]])  %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(color = model)) +
  theme_minimal() +
  ggtitle(glue('random_forest = ', {round(forest_clinvar_deletion[[4]], 2)}, ' - ', 
               'logistic = ', {round(logistic_clinvar_deletion[[4]], 2)}, ' - ', 
               'xgboost = ', {round(xgboost_clinvar_deletion[[4]], 2)}, ' - ', 
               'rulefit = ', {round(rulefit_clinvar_deletion[[4]], 2)}, ' - ',
               'bayesian = ', {round(bayesian_clinvar_deletion[[4]], 2)}))


ggsave("cnvscore_results/evaluation_model/deletion/independent_roc_auc.png", clinvar_roc_deletion + clinvar_roc_pr_deletion, width = 17, height = 9.6, dpi = 300, units = "in", device='png')


logistic_clinvar_duplication <- evaluate_model(logistic_model_duplication, test_object = clinvar_ind_duplication, model_name = 'Logistic regression', cv = FALSE, 
                                            cum_dist = FALSE)

xgboost_clinvar_duplication <- evaluate_model(xgboost_model_duplication, test_object = clinvar_ind_duplication, model_name = 'XGBoost', cv = FALSE, 
                                           cum_dist = FALSE)

forest_clinvar_duplication <- evaluate_model(forest_model_duplication, test_object = clinvar_ind_duplication, model_name = 'Random Forest', cv = FALSE, 
                                          cum_dist = FALSE)

rulefit_clinvar_duplication <- evaluate_model(rulefit_model_duplication, rulefit = TRUE, test_object = clinvar_ind_duplication, model_name = 'RuleFit', cv = FALSE, 
                                           cum_dist = FALSE)

bayesian_clinvar_duplication <- evaluate_model(bayesian_model_duplication, bay = TRUE, test_object = clinvar_ind_duplication, model_name = 'Bayesian', cv = FALSE, 
                                            cum_dist = FALSE)

clinvar_roc_duplication <- bind_rows(logistic_clinvar_duplication[[1]], xgboost_clinvar_duplication[[1]], 
                                  forest_clinvar_duplication[[1]], rulefit_clinvar_duplication[[1]], bayesian_clinvar_duplication[[1]]) %>%
  
  ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(color = model)) +
  theme_bw() +
  ggtitle(glue('random_forest = ', {round(forest_clinvar_duplication[[3]], 2)}, ' - ', 
               'logistic = ', {round(logistic_clinvar_duplication[[3]], 2)}, ' - ', 
               'xgboost = ', {round(xgboost_clinvar_duplication[[3]], 2)}, ' - ', 
               'rulefit = ', {round(rulefit_clinvar_duplication[[3]], 2)}, ' - ',
               'bayesian = ', {round(bayesian_clinvar_duplication[[3]], 2)}))

clinvar_roc_pr_duplication <- bind_rows(logistic_clinvar_duplication[[2]], 
                                     xgboost_clinvar_duplication[[2]], 
                                     forest_clinvar_duplication[[2]], 
                                     rulefit_clinvar_duplication[[2]], 
                                     bayesian_clinvar_duplication[[2]])  %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(color = model)) +
  theme_minimal() +
  ggtitle(glue('random_forest = ', {round(forest_clinvar_duplication[[4]], 2)}, ' - ', 
               'logistic = ', {round(logistic_clinvar_duplication[[4]], 2)}, ' - ', 
               'xgboost = ', {round(xgboost_clinvar_duplication[[4]], 2)}, ' - ', 
               'rulefit = ', {round(rulefit_clinvar_duplication[[4]], 2)}, ' - ',
               'bayesian = ', {round(bayesian_clinvar_duplication[[4]], 2)}))


ggsave("cnvscore_results/evaluation_model/duplication/independent_roc_auc.png", clinvar_roc_duplication + clinvar_roc_pr_duplication, width = 17, height = 9.6, dpi = 300, units = "in", device='png')


# ------------------------------------------------------------------------------
# INDEPENDENT - PLOT CUMULATIVE FREQUENCY
# ------------------------------------------------------------------------------

select_filter <- if_else(str_detect(target_class, 'pathogenic'), 'pathogenic', 'benign')

target_n_deletion <- clinvar_ind_deletion %>% count(clinical) %>% filter(clinical == select_filter) %>% pull(n)


logistic_clinvar_cum_dist_deletion <- evaluate_model(logistic_model_deletion, test_object = clinvar_ind_deletion, model_name = 'Logistic regression', cv = FALSE, 
                                cum_dist = TRUE)
xgboost_clinvar_cum_dist_deletion <- evaluate_model(xgboost_model_deletion, test_object = clinvar_ind_deletion, model_name = 'XGBoost', cv = FALSE, 
                               cum_dist = TRUE)
forest_clinvar_cum_dist_deletion <- evaluate_model(forest_model_deletion, test_object = clinvar_ind_deletion, model_name = 'Random Forest', cv = FALSE, 
                              cum_dist = TRUE)

rulefit_clinvar_cum_dist_deletion <- evaluate_model(rulefit_model_deletion, rulefit = TRUE, test_object = clinvar_ind_deletion, model_name = 'RuleFit', cv = FALSE, 
                                          cum_dist = TRUE)

bayesian_clinvar_cum_dist_deletion <- evaluate_model(bayesian_model_deletion, bay = TRUE, test_object = clinvar_ind_deletion, model_name = 'Bayesian', cv = FALSE, 
                                           cum_dist = TRUE)

cumulative_test_deletion <- cum_plot(by_percentage = TRUE, 
                            select_filter = select_filter, 
                            target_n = target_n_deletion,
                            total_n = clinvar_ind_deletion %>% nrow(),
                            logistic_clinvar_cum_dist_deletion, 
                            xgboost_clinvar_cum_dist_deletion, 
                            forest_clinvar_cum_dist_deletion, 
                            rulefit_clinvar_cum_dist_deletion,
                            bayesian_clinvar_cum_dist_deletion)



ggsave("cnvscore_results/evaluation_model/deletion/independent_cum_frequency.png", cumulative_test_deletion, width = 17, height = 9.6, dpi = 300, units = "in", device='png')


target_n_duplication <- clinvar_ind_duplication %>% count(clinical) %>% filter(clinical == select_filter) %>% pull(n)


logistic_clinvar_cum_dist_duplication <- evaluate_model(logistic_model_duplication, test_object = clinvar_ind_duplication, model_name = 'Logistic regression', cv = FALSE, 
                                                     cum_dist = TRUE)
xgboost_clinvar_cum_dist_duplication <- evaluate_model(xgboost_model_duplication, test_object = clinvar_ind_duplication, model_name = 'XGBoost', cv = FALSE, 
                                                    cum_dist = TRUE)
forest_clinvar_cum_dist_duplication <- evaluate_model(forest_model_duplication, test_object = clinvar_ind_duplication, model_name = 'Random Forest', cv = FALSE, 
                                                   cum_dist = TRUE)

rulefit_clinvar_cum_dist_duplication <- evaluate_model(rulefit_model_duplication, rulefit = TRUE, test_object = clinvar_ind_duplication, model_name = 'RuleFit', cv = FALSE, 
                                                    cum_dist = TRUE)

bayesian_clinvar_cum_dist_duplication <- evaluate_model(bayesian_model_duplication, bay = TRUE, test_object = clinvar_ind_duplication, model_name = 'Bayesian', cv = FALSE, 
                                                     cum_dist = TRUE)

cumulative_test_duplication <- cum_plot(by_percentage = TRUE, 
                                     select_filter = select_filter, 
                                     target_n = target_n_duplication,
                                     total_n = clinvar_ind_duplication %>% nrow(),
                                     logistic_clinvar_cum_dist_duplication, 
                                     xgboost_clinvar_cum_dist_duplication, 
                                     forest_clinvar_cum_dist_duplication, 
                                     rulefit_clinvar_cum_dist_duplication,
                                     bayesian_clinvar_cum_dist_duplication)


ggsave("cnvscore_results/evaluation_model/duplication/independent_cum_frequency.png", cumulative_test_duplication, width = 17, height = 9.6, dpi = 300, units = "in", device='png')



# ------------------------------------------------------------------------------
# HYPERPARAMETER TUNING
# ------------------------------------------------------------------------------

# penalty = 0.003
# tree_depth = 8

library(doParallel)



cl <- makePSOCKcluster(40)
registerDoParallel(cl)


xgboost_mod <-
  boost_tree(learn_rate = tune(), trees = tune()) %>%
  set_mode("classification") %>%
  set_engine("xgboost")

evaluation_strategy <- metric_set(yardstick::roc_auc)
resampling_strategy <- vfold_cv(training_tbl, prop = 0.9, strata = clinical)
control_strategy <- control_grid(verbose = TRUE, allow_par = TRUE)

formula_res <-
  xgboost_mod %>% 
  tune_grid(
    formule_models,
    metrics = evaluation_strategy,
    control = control_strategy,
    resamples = resampling_strategy,
    grid = 100
  )


formula_res %>%
  collect_metrics() %>%
  ggplot(aes(factor(trees), factor(tree_depth))) +
  geom_tile(aes(fill = mean)) +
  scale_fill_viridis_c()

# 2 - grid
tree_grid <- grid_regular(tree_depth(), penalty(), levels = 10)

# 1 - model 
test_model <- boost_tree() %>%
  set_mode('classification') %>%
  set_engine('xgboost')

library("doFuture")


registerDoFuture()
cl <- makeCluster(40)
plan(future::cluster, workers = cl)


a <- tune_grid(formule_models, test_model, resamples = training_cv, grid = tree_grid)


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
# END
# ------------------------------------------------------------------------------



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



# ------------------------------------------------------------------------------
# RULEFIT - EXPLORATION
# ------------------------------------------------------------------------------

coef(bayesian_model_deletion)
# COEFFICIENTS
# lambda.1se or lambda.min
rules_tbl <- coef(rulefit_model_deletion, lambda = 'lambda.min') %>% 
  as_tibble() %>%
  rename(coefficient = coefficient_lambda.min) %>%
  filter(coefficient != 0)

# SUMMARY CV
print(rulefit_model_deletion$glm$model)
# Df -> nº of nonzero coefficients
print(rulefit_model_deletion$glm$model$glmnet.fit)

# CROSS-VALIDATION GRID LAMBDA PLOT
plot(rulefit_model_deletion$glm$model)


rules_tbl <- rules_tbl %>%
  rowwise() %>%
  mutate(tmp_column = get_rules_information(training_tbl, rule)) %>%
  separate(col = tmp_column, into = c('support', 'risk'), sep = ' ') %>%
  mutate(support = as.integer(support), risk = as.numeric(risk)) %>%
  relocate(term, rule, coefficient, support, risk) %>%
  arrange(desc(risk))

rules_tbl %>%
  mutate(color_support = if_else(support <= 5, 'yes', 'no')) %>%
  ggplot(aes(coefficient, risk)) +
  geom_point(aes(fill = color_support), shape = 21, color = 'black') + 
  geom_hline(yintercept = 0.5, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") + 
  theme_minimal()


rules_tbl %>%
  select(rule) %>%
  separate_rows(rule, sep = ' & ') %>%
  mutate(rule = str_remove(rule, '\\<.*')) %>%
  mutate(rule = str_remove(rule, '\\>.*')) %>%
  count(rule) %>%
  na.omit() %>%
  ggplot(aes(reorder(rule, n), n)) +
    geom_col(aes(fill = n), color = 'black') +
    coord_flip() +
  scale_fill_viridis_c() +
  theme_minimal()




generate_rules <- function(model, input_obs) {
  
  # model <- rulefit_model_deletion
  # input_obs <- clinvar_ind_deletion
  model <- coef(model) %>% as_tibble() %>% filter(str_detect(term, '^[r]'))
  
  tmp_result <- tibble()

  for (i in 1:nrow(input_obs)) {
    
    print(glue('{i}', {nrow(input_obs)}))
    
  obs <- input_obs[i,]

  tmp_obs <- input_obs[i,] %>%
    bind_cols(
      model %>% 
      rowwise() %>%
      mutate(yes_no = if_else(nrow(obs %>% filter_(rule)) == 1, 1, 0)) %>%
        ungroup() %>%
        select(term, yes_no) %>%
      pivot_wider(values_from = yes_no, names_from = term))
    
  tmp_result <- tmp_result %>% bind_rows(tmp_obs)
  }

  return(tmp_result)

}

tic('Creating rules')
clinvar_ind_deletion_rules <- generate_rules(rulefit_model_deletion, clinvar_ind_deletion %>% sample_n(300))
toc()


enrich_prediction <- function(input_model, input_test, model_name) {
  
  # input_model <- bayesian_model_deletion
  # input_test <- clinvar_ind_deletion_rules
  
  tmp_result <- posterior_linpred(input_model, newdata = input_test, transform = TRUE) %>%
    as_tibble()
  
  # point estimate
  tmp_point_estimate <- tmp_result %>% map_dbl(~ median(.x)) %>% enframe(name = NULL) %>%
    rename(pred_target = value) %>%
    mutate(pred_target = 1 - pred_target) 
  
  # uncertainty
  tmp_uncertainty <- tmp_result %>% map_dbl(~ sd(.x)) %>% enframe(name = NULL) %>%
    rename(sd = value)
  
  # applicability (larger distance -> lower rank)
  tmp_applicability <- score(applicability_pca, clinvar_ind_deletion_rules) %>%
    select(distance_pctl) %>%
    mutate(distance_pctl = 100 - distance_pctl)
  
    final_result <- input_test %>% 
      select(clinical) %>%
      bind_cols(tmp_point_estimate,
              tmp_uncertainty,
              tmp_applicability) 
  
  
  return(final_result)
}


p_corr1 <- enrich_prediction(bayesian_model_deletion, clinvar_ind_deletion_rules) %>% 
  select(-clinical) %>% 
  correlate(method = 'spearman') %>%
  stretch(remove.dups = TRUE) %>%
  na.omit() %>%
  mutate(id = paste(x, y, sep =' - ')) %>%
  mutate(positive = if_else(r >= 0, 'yes', 'no'), show.legend = FALSE) %>%
  ggplot(aes(id, r)) + 
  geom_col(aes(fill = positive), color = 'black') +
  coord_flip() +
  labs(x = 'Association', y = 'Spearman correlation', title = 'Correlation') +
  theme_minimal()

p_corr2 <- enrich_prediction(bayesian_model_deletion, clinvar_ind_deletion_rules) %>% 
  ggplot(aes(sd, distance_pctl)) +
    geom_point() +
  theme_minimal()

p_corr3 <- enrich_prediction(bayesian_model_deletion, clinvar_ind_deletion_rules) %>% 
  ggplot(aes(sd, pred_target)) +
  geom_point() +
  theme_minimal()

p_corr4 <- enrich_prediction(bayesian_model_deletion, clinvar_ind_deletion_rules) %>% 
  ggplot(aes(distance_pctl, pred_target)) +
  geom_point() +
  theme_minimal()


ggsave(glue("cnvscore_results/plot_model/bayesian_predictions_cor1"), p_corr1, width = 17, height = 9.6, dpi = 300, units = "in", device='png')
ggsave(glue("cnvscore_results/plot_model/bayesian_predictions_cor2"), p_corr2 + p_corr3 + p_corr4, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

# ------------------------------------------------------------------------------
# RULEFIT - APPLICABILITY DOMAIN
# # https://www.marlycormar.com/presentations/R-Pharma-2019/presentation.html#1
# distance between each principal component and its mean.
# ------------------------------------------------------------------------------
library(recipes)

recipe_applicable <-
  recipe( ~ ., data = training_tbl %>% select(any_of(all.vars(formule_models)[-1]))) %>%
  step_dummy(all_nominal()) %>%
  # Remove variables that have the same value for every data point.
  step_zv(all_predictors()) %>%
  # Transform variables to be distributed as Gaussian-like as possible.
  step_YeoJohnson(all_numeric()) %>%
  # Normalize numeric data to have a mean of zero and
  # standard deviation of one.
  step_normalize(all_numeric())

applicability_pca <- apd_pca(recipe_applicable, data = training_tbl_deletion)

pca_score <- score(applicability_pca, clinvar_ind_deletion_rules)

  
pca_output <- test_tbl %>% 
  # bind_cols(pca_score %>% select(distance_pctl) %>% rename(applicability = distance_pctl)) %>%
  select(any_of(all.vars(formule_models)[-1])) %>%
    prcomp(scale = TRUE, center = TRUE)


pca_output$x %>%
  as_tibble() %>%
  select(PC1, PC2) %>%
  bind_cols(pca_score %>% select(distance_pctl) %>% rename(applicability = distance_pctl)) %>%
  mutate(applicability = if_else(applicability >= 95, 'No', 'Yes')) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = factor(applicability)), shape = 21, color = 'black') +
  labs(fill = 'Applicability?') +
  theme_minimal()

# ------------------------------------------------------------------------------
# RULEFIT - UNCERTAINTY
# laplacian problem -> in order to have a lot of zeros in the signal, 
# you are also forcing the non-zero elements to be very small
# BAYESIAN APPROACH (laplacian -> horsehoe -> finnish horseshoe)
# ------------------------------------------------------------------------------
# plot(bayesian_model, "areas", prob = 0.95, prob_outer = 1)

library(rstan)
library(rstanarm)

# nº iterations per chain (iter)
# nº chains (chains)
# algorithm (algorithm = "sampling", "optimizing", "meanfield", "fullrank")
# QR
# posterior distribution -> point estimate (mean, median)





bayesian_model <- rstanarm::stan_glm(formula = formule_models,
                                       family = 'binomial',
                                       data = rulefit_model_deletion$full_data,
                                       cores = 4,
                                       iter = 2000,
                                       chains = 4,
                                       algorithm = 'sampling', # variational inference algorithms
                                       QR = TRUE,
                                       prior = hs(),
                                       prior_intercept = hs())

library(bayesplot)
print(bayesian_model)
summary(bayesian_model)

as_tibble(bayesian_model_deletion) %>%
  mcmc_areas(prob = 0.8)

 prev_bayesian <- posterior_linpred(bayesian_model_deletion, newdata = clinvar_ind_deletion_rules, transform = TRUE) %>%
    as_tibble() %>%
    map_dbl(~ mean(.x)) %>% # mean
    as_tibble() %>%
    bind_cols(clinvar_ind_deletion_rules %>% select(clinical)) %>%
    rename(pred_target = value) %>%
    mutate(pred_target = 1 - pred_target)
    
  bayesian_roc <- prev_bayesian %>%
    roc_curve(clinical, pred_target) %>%
    mutate(model = 'bayesian')
  
  bayesian_pr <- prev_bayesian %>%
    pr_curve(clinical, pred_target) %>%
    mutate(model = 'bayesian')
  
  bayesian_auc <- prev_bayesian %>%
    roc_auc(clinical, pred_target) %>%
    mutate(model = 'bayesian') %>% 
    pull(.estimate) %>% 
    round(2)
  
  bayesian_pr_auc <- prev_bayesian %>%
    pr_auc(clinical, pred_target) %>%
    mutate(model = 'bayesian') %>% 
    pull(.estimate) %>% 
    round(2)

  
 bayesian_roc %>%
    ggplot(aes(1-specificity, sensitivity)) +
    geom_path(aes(group = model, color = model), size = 2.5, show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size=10, face="bold"),
          axis.title.x = element_text(size=15, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    ggtitle(glue('bayesian = ', {bayesian_auc})) +
    theme(legend.position="bottom")
  
  
  pr_auc <-   bayesian_pr %>%
    bind_rows(logistic_pr, forest_pr, rulefit_pr, xgboost_pr) %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = model, color = model), size = 2.5) +
    theme_bw() +
    theme(plot.title = element_text(size=10, face="bold"),
          axis.title.x = element_text(size=15, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold"),
          legend.text = element_text(size=19),
          legend.title = element_text(size = 15)) +
    ggtitle(glue('random_forest = ', {forest_pr_auc}, ' - ', 'logistic = ', {logistic_pr_auc}, ' - ',  'xgboost = ', {xgboost_pr_auc},
                 ' - ', 'rulefit = ', {rulefit_pr_auc},' - ', 'bayesian = ', {bayesian_pr_auc})) +
    theme(legend.position="bottom")

  
  roc_auc + pr_auc
  
  
# ------------------------------------------------------------------------------
# API CNVSCORE
# https://medium.com/@JB_Pleynet/how-to-do-an-efficient-r-api-81e168562731
# ------------------------------------------------------------------------------
  
library(plumber)
  

# ------------------------------------------------------------------------------
# FUNCTION - LENGTH DISTRIBUTION - PLOT
# ------------------------------------------------------------------------------

plot_length_distribution <- function(data) {
  
  p_distribution0 <- data %>%
    # mutate(source = if_else(source == 'decipher', paste(source, pathogenicity), source)) %>%
    ggplot(aes(length_cnv)) +
    geom_histogram(aes(fill = source),color = 'black', alpha = 0.4) +
    scale_x_log10() +
    facet_wrap(~ source) +
    theme_minimal()
  
  
  p_distribution0_2 <- data %>%
    # mutate(source = if_else(source == 'decipher', paste(source, pathogenicity), source)) %>%
    group_by(source) %>%
    mutate(p_length_cnv = ntile(length_cnv, 100)) %>%
    group_by(source, p_length_cnv) %>%
    mutate(median_length_cnv = median(length_cnv)) %>%
    ungroup() %>%
    ggplot(aes(p_length_cnv, median_length_cnv)) +
    geom_point(aes(fill = source), shape = 21) +
    scale_y_log10() +
    theme_minimal()
  
  
  p_distribution0_3 <- data %>%
    # mutate(source = if_else(source == 'decipher', paste(source, pathogenicity), source)) %>%
    ggplot(aes(length_cnv)) +
    geom_density(aes(fill = source), alpha = 0.4) +
    scale_x_log10() +
    facet_wrap(~ source) +
    theme_minimal()
  
  p_distribution0_4 <- data %>%
    ggplot(aes(length_cnv)) +
      geom_density(aes(fill = clinical), alpha = 0.4) +
      scale_x_log10() + 
      theme_minimal()
  
  p_distribution0 + p_distribution0_2 + p_distribution0_3 + p_distribution0_4 + plot_layout(nrow = 2)
  
}

# ------------------------------------------------------------------------------
# FUNCTION - CUMULATIVE PLOT
# ------------------------------------------------------------------------------

cum_plot <- function(by_percentage = TRUE, target_n, total_n,  select_filter, ...) {
  
  input_total <- bind_rows(...)
  
  random_perc <- (target_n /  total_n) / target_n
  random_tbl <- tibble('model' = 'random', cum_line = rep(random_perc, total_n))
  random_tbl <- random_tbl %>% 
    mutate(perc_line = cumsum(cum_line)) %>%
    mutate(number_cnv = row_number())
    
  
  
  # target_n <- input_model %>%
  #   count(clinical) %>%
  #   filter(clinical == select_filter) %>% pull(n)

  # 
  #  a <- logistic_rs %>% unnest(.predictions) %>% select(id, .pred_Pathogenic, clinical) %>% 
  #    mutate(model = 'pepito')  %>% rename(pred_target = target_class)
  # 
  # 
  # a %>%
  #   group_by(id, model) %>%
  #   arrange(desc(pred_target)) %>%
  #   mutate(cum_line1 = if_else(clinical == select_filter, 1, 0)) %>%
  #   mutate(cum_line = cumsum(cum_line1)) %>%
  #   mutate(perc_line = cum_line / target_n) %>%
  #   mutate(number_cnv = row_number()) %>%
  #   ggplot(aes(number_cnv, perc_line)) +
  #   geom_path(aes(group = id, color = id), size = 2) +
  #   scale_y_continuous(label = scales::percent) +
  #   facet_wrap(~ model) +
  #   xlab('Total number of CNVs') +
  #   ylab(glue('Cumulative frequency')) +
  #   ggtitle(glue('Cumulative frequency {select_filter} CNVs')) +
  #   theme_minimal()
    
  
  if (by_percentage) {
    
    input_total %>%
    group_by(model) %>%
    arrange(desc(pred_target)) %>%
    mutate(cum_line1 = if_else(clinical == select_filter, 1, 0)) %>%
    mutate(cum_line = cumsum(cum_line1)) %>%
    mutate(perc_line = cum_line / target_n) %>%
    mutate(number_cnv = row_number()) %>%
    ggplot(aes(number_cnv, perc_line)) +
    geom_path(aes(group = model, color = model), size = 2) +
    geom_path(data = random_tbl, aes(number_cnv, perc_line, group = model, color = model), size = 2, linetype = 'dashed') +
    scale_y_continuous(label = scales::percent) +
    xlab('Total number of CNVs') +
    ylab(glue('Cumulative frequency')) +
    ggtitle(glue('Cumulative frequency {select_filter} CNVs')) +
    theme_minimal()
    
  } else {
    
    input_total %>%
      group_by(model) %>%
      arrange(desc(pred_target)) %>%
      mutate(cum_line1 = if_else(clinical == select_filter, 1, 0)) %>%
      mutate(cum_line = cumsum(cum_line1)) %>%
      mutate(number_cnv = row_number()) %>%
      ggplot(aes(number_cnv, cum_line)) +
      geom_path(aes(group = model, color = model), size = 2) +
      xlab('Total number of CNVs') +
      ylab(glue('Cumulative frequency')) +
      ggtitle(glue('Cumulative frequency {select_filter} CNVs')) +
      geom_hline(yintercept = target_n, linetype = 'dashed') +
      theme_minimal()

    
  }

}

# ------------------------------------------------------------------------------
# FUNCTION - RECIPROCAL OVERLAP
# ------------------------------------------------------------------------------

reciprocal_overlap <- function(input_check_cnv) {
  
  print(glue('Number of input rows: {nrow(input_check_cnv)}'))

  vector_id_to_keep <- c()
  
  for (i in 1:nrow(input_check_cnv)) {
    
    print(glue('{i}/{nrow(input_check_cnv)}'))
    
  
  tmp_a <- input_check_cnv %>% slice(i) %>% select(chrom, start, end, id_tmp)
  tmp_a_length = tmp_a %>% mutate(length_cnv = end - start + 1) %>% pull(length_cnv)
  tmp_a_id_tmp = tmp_a %>% pull(id_tmp)
  tmp_b <- input_check_cnv %>% slice(-i) %>% select(chrom, start, end, id_tmp)
  
  
  id_overlap_a <- tmp_a %>% bed_intersect(tmp_b) %>% 
    mutate(overlap = .overlap / tmp_a_length) %>%
    filter(overlap >= 0.9) %>%
    pull(id_tmp.y)
  
  
  id_overlap_b <- tmp_b %>% bed_intersect(tmp_a) %>%
    mutate(length_cnv = end.x - start.x + 1) %>%
    mutate(overlap = .overlap / length_cnv) %>%
                 filter(overlap >= 0.9) %>%
                 pull(id_tmp.x)
  
  both_a_b <- id_overlap_a[id_overlap_a %in% id_overlap_b]
  
  if (length(both_a_b) == 0) {
    
    vector_id_to_keep <- c(vector_id_to_keep, tmp_a_id_tmp)
    
  } else {
    
    take_smallest_cnv <- input_check_cnv %>%
      filter(id_tmp %in% c(tmp_a_id_tmp, both_a_b)) %>%
      mutate(length_cnv = end - start + 1) %>%
      arrange(length_cnv) %>% 
      slice(1) %>%
      pull(id_tmp)
    
    vector_id_to_keep <- c(vector_id_to_keep, take_smallest_cnv) 
  
   }  
  }
  vector_id_to_keep <- vector_id_to_keep %>% unique()
  print(glue('Number of input rows: {nrow(input_check_cnv)}, \n after filtering: {length(vector_id_to_keep)}'))
  return(vector_id_to_keep)
}


# ------------------------------------------------------------------------------
#  FUNCTION - GENERATE RANDOM CNVs FOR EACH DECIPHER CNV
# ------------------------------------------------------------------------------


gen_random <- function(input_chrom, input_start, input_end, input_tmp_id, n_rep = 9, input_n_try) {
  
  result_df <- tibble()
  
  tmp_chrom <- input_chrom
  tmp_start <- input_start
  tmp_end <- input_end
  tmp_id <- input_tmp_id
  
  tmp_length <- tmp_end - tmp_start + 1
  
  target_cnv <- tibble('chrom' = tmp_chrom,'start' = tmp_start, 'end' = tmp_end)
  
  
  tmp_cyto <- bed_intersect(coord_cytobands, target_cnv) %>% 
    pull(Name.x)
  
  
  tmp_n_genes <- hgcn_genes %>% bed_intersect(target_cnv) %>% nrow()
  
  result_tmp_cnv <- tibble()
  
  n_try <- 0
  
  while (nrow(result_tmp_cnv) < n_rep) {
    
    n_try <- n_try + 1
    
    if (n_try < input_n_try) {
      
      # print(n_try)
      
      
      selected_cytoband <- coord_cytobands %>% filter(Name %in% tmp_cyto & chrom == tmp_chrom)
      random_from <- selected_cytoband %>% slice_head() %>% pull(start)
      random_to <-  selected_cytoband %>% slice_tail() %>% pull(end) # - tmp_length
      random_start <- sample(random_from:random_to, 1)
      random_end <- random_start + tmp_length
      
      random_tbl <- tibble('chrom' = tmp_chrom, 
                           'start' = random_start, 
                           'end' = random_end)
      
      # Filter 1: Remove random CNVs mapping the original CNV
      overlap_original_cnv <- random_tbl %>% bed_intersect(target_cnv)
      if (nrow(overlap_original_cnv) > 0) next
      
      # Filter 2: eliminate cnvs with no genes if the original cnv has genes
      random_genes <- hgcn_genes %>% bed_intersect(random_tbl) %>% nrow()
      if (tmp_n_genes > 0 & random_genes == 0) next 
      
      result_tmp_cnv <- result_tmp_cnv %>% bind_rows(random_tbl)
      
    } else {
      
      result_tmp_cnv <- tibble('chrom' = 'error', 'start' = 0, 'end' = 0)
      break
    }
    
  }
  result_df <- result_df %>% bind_rows(result_tmp_cnv %>% mutate('id_tmp' = tmp_id))
  # }
  
  return(result_df)
  
}

# ------------------------------------------------------------------------------
# FUNCTION - PLOT CROSS-VALIDATION RESULTS (ROC-AUC + PR-AUC)
# ------------------------------------------------------------------------------

  
evaluate_model <- function(model, model_name, 
                           test_object = NULL, 
                           rulefit = FALSE,
                           bay = FALSE,
                           cv = FALSE, 
                           cum_dist = FALSE) {
  

  if (isTRUE(cv) & isFALSE(cum_dist)) {
    
  # ROC AUC
  tmp_mean <- model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    roc_auc(clinical, target_class) %>% 
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <- model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    roc_auc(clinical, target_class) %>% 
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p1 <- model %>% 
    unnest(.predictions) %>%
    # filter(id == 'Fold01') %>%
    group_by(id) %>% 
    # ggplot(aes(.pred_Pathogenic)) + geom_density(aes(fill = clinical), alpha = 0.3) + facet_wrap(~ id, scales = 'free')
    roc_curve(clinical, target_class) %>% 
    ggplot(aes(1-specificity, sensitivity)) +
    geom_path(aes(group = id, color = id),  show.legend = T) +
    theme_bw() +
    theme(plot.title = element_text(size=19, face="bold"),
          axis.title.x = element_text(size=19, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    
    labs(title = glue('{model_name}  (Mean AUC: {tmp_mean} ± {tmp_err})'))
  
  # PR AUC
  tmp_mean <- model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    pr_auc(clinical, target_class) %>% 
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <-model %>% 
    unnest(.predictions) %>% 
    group_by(id) %>% 
    pr_auc(clinical, target_class) %>% 
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p2 <- model %>% 
    unnest(.predictions) %>%
    group_by(id) %>%
    pr_curve(clinical, target_class) %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = id, color = id), show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size=19, face="bold"),
          axis.title.x = element_text(size=19, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    labs(title = glue('{model_name} (Mean AUCpr: {tmp_mean} ± {tmp_err})'))
  
  return(list(p1, p2))
  
  
  } else if (isTRUE(cv) & isTRUE(cum_dist)) {
    
    
    if (isFALSE(rulefit)) {
    perc_10_cv_tmp <- model %>% 
      unnest(.predictions) %>%
      group_by(id) %>%
      rename(p_patho = target_class) %>%
      mutate(p_patho = ntile(-p_patho, 10)) %>% 
      select(id, p_patho, clinical) %>%
      group_by(id, p_patho) %>%
      count(clinical) %>% 
      filter(clinical == select_filter) %>%
      group_by(id) %>%
      mutate(percentage = n / sum(n)) %>%
      group_by(id) %>%
      mutate(cum_perc = cumsum(percentage)) %>%
      mutate(model = model_name)
    
    } else {
      
      perc_10_cv_tmp <- model %>% 
        select(id, prob_predicted, clinical) %>%
        group_by(id) %>%
        mutate(prob_predicted = ntile(-prob_predicted, 10)) %>% 
        select(id, prob_predicted, clinical) %>%
        group_by(id, prob_predicted) %>%
        count(clinical) %>% 
        filter(clinical == select_filter) %>%
        group_by(id) %>%
        mutate(percentage = n / sum(n)) %>%
        mutate(cum_perc = cumsum(percentage)) %>%
        mutate(model = model_name) %>%
        rename(p_patho = prob_predicted)
    }
    
    return(perc_10_cv_tmp)
    
  } else if (isFALSE(cv) & isTRUE(cum_dist)) {
    
    if (isFALSE(rulefit) & isFALSE(bay)) {
      
    tmp_object <- predict(model, test_object, type = 'prob') %>%
    rename(pred_target = target_class) %>%
    select(pred_target) %>%
    bind_cols(test_object %>% select(clinical)) %>% mutate(model = model_name)
    
    } else if (isTRUE(bay)) {
      
      
      tmp_object <-  posterior_linpred(model, newdata = test_object, transform = TRUE) %>%
        as_tibble() %>%
        map_dbl(~ median(.x)) %>% # mean, median, map_estimate
        as_tibble() %>%
        bind_cols(test_object %>% select(clinical)) %>%
        rename(pred_target = value) %>%
        mutate(pred_target = 1 - pred_target) %>%
        mutate(model = model_name)
      
    
    } else {
      
      tmp_object <-  tibble(target_class = 1 - as.vector(predict(model, test_object, type = 'response'))) %>%
        rename(pred_target = target_class) %>%
        select(pred_target) %>%
        bind_cols(test_object %>% select(clinical))  %>% 
        mutate(model = model_name)
      
      
    }
  
    return(tmp_object)
    
  } else if (isFALSE(cv) & isFALSE(cum_dist)) {
    
    if (isFALSE(rulefit) & isFALSE(bay)) {
      
    tmp_object <- predict(model, test_object, type = 'prob') %>%
      rename(pred_target = target_class) %>%
      select(pred_target) %>%
      bind_cols(test_object %>% select(clinical)) %>% mutate(model = model_name)
    
    } else if (isTRUE(bay)) {
      
      tmp_object <-  posterior_linpred(model, newdata = test_object, transform = TRUE) %>%
        as_tibble() %>%
        map_dbl(~ median(.x)) %>% # mean, median, map_estimate
        as_tibble() %>%
        bind_cols(test_object %>% select(clinical)) %>%
        rename(pred_target = value) %>%
        mutate(pred_target = 1 - pred_target) %>%
        mutate(model = model_name)

      
    } else {
      
      tmp_object <-  tibble(target_class = 1 - as.vector(predict(model, test_object, type = 'response'))) %>%
        rename(pred_target = target_class) %>%
        select(pred_target) %>%
        bind_cols(test_object %>% select(clinical))  %>% 
        mutate(model = model_name)
    }
    
    tmp_roc <- tmp_object %>%
      roc_curve(clinical, pred_target) %>%
      mutate(model = model_name)
  
    tmp_pr <- tmp_object %>%
      pr_curve(clinical, pred_target) %>%
      mutate(model = model_name)
  
    tmp_auc <- tmp_object %>%
      roc_auc(clinical, pred_target) %>%
      mutate(model = model_name) %>% 
      pull(.estimate) %>% 
      round(2)
  
    tmp_pr_auc <- tmp_object %>%
      pr_auc(clinical, pred_target) %>%
      mutate(model = model_name) %>% 
      pull(.estimate) %>% 
      round(2)
      
    
    
    return (list(tmp_roc, tmp_pr, tmp_auc, tmp_pr_auc))
    
  }
  
}



# ------------------------------------------------------------------------------
# FUNCTION (ONLY FOR RULEFIT) - PLOT CROSS-VALIDATION RESULTS (ROC-AUC + PR-AUC)
# ------------------------------------------------------------------------------


from_cv_rulefit <- function(model, model_name = 'RuleFit') {
  
  # model <- rulefit_rs
  
  model <- model %>% mutate(prob_predicted = 1 - prob_predicted)
  
  # ROC AUC
  tmp_mean <- model %>% 
    group_by(id) %>% 
    roc_auc(clinical, prob_predicted) %>% 
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <- model %>% 
    group_by(id) %>% 
    roc_auc(clinical, prob_predicted) %>% 
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p1 <- model %>% 
    group_by(id) %>% 
    roc_curve(clinical, prob_predicted) %>% 
    ggplot(aes(1-specificity, sensitivity)) +
    geom_path(aes(group = id, color = id),  show.legend = T) +
    theme_bw() +
    theme(plot.title = element_text(size=19, face="bold"),
          axis.title.x = element_text(size=19, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    
    labs(title = glue('{model_name}  (Mean AUC: {tmp_mean} ± {tmp_err})'))
  
  # PR AUC
  tmp_mean <- model %>% 
    group_by(id) %>% 
    pr_auc(clinical, prob_predicted) %>% 
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <- model %>% 
    group_by(id) %>% 
    pr_auc(clinical, prob_predicted) %>% 
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p2 <- model %>% 
    group_by(id) %>%
    pr_curve(clinical, prob_predicted) %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = id, color = id), show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size=19, face="bold"),
          axis.title.x = element_text(size=19, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    labs(title = glue('{model_name} (Mean AUCpr: {tmp_mean} ± {tmp_err})'))
  
  return(list(p1, p2))
  
}

# ------------------------------------------------------------------------------
# FUNCTION - BAYESIAN RULEFIT - CROSS-VALIDATION
# ------------------------------------------------------------------------------


bay_rulefit_cv <- function(cv_object) {
  
  bayesian_result <- tibble()
  
  for (i in 1:10) {
    
    print(glue('Rulefit - split nº {i}/10 '))
    
    
    tmp_training <- cv_object$splits[[i]] %>% training()
    tmp_testing <- cv_object$splits[[i]] %>% testing()
    tmp_name_fold <- cv_object %>% slice(i) %>% pull(id)
    
    
    data_pre_lasso <- xrftest::xrf(formule_models,
                                   data = tmp_training,
                                   # xgb_control = list(nrounds = 100, max_depth = 5, min_child_weight = 3),
                                   family = "binomial")
    
    print(glue('Bayesian - split nº {i}/10 '))
    
    
    # stan(file = 'stan_model.stan' , iter= 2000, data = data_pre_lasso$full_data,
    #      chains= 4, seed=194838)
    
    # special_columns <-    data_pre_lasso %>%
    #   coef(s = "lambda.min") %>%
    #   as_tibble() %>% filter(coefficient_lambda.min >= 1 | coefficient_lambda.min <= -1) %>%
    #   na.omit() %>%
    #   pull(term)
    
    bayesian_model <- rstanarm::stan_glm(formule_models,
                                         family = 'binomial',
                                         data = data_pre_lasso$full_data,
                                         # data = data_pre_lasso$full_data[,colnames(data_pre_lasso$full_data) %in% special_columns],
                                         # cores = 4,
                                         iter = 2000,
                                         # chains = 4,
                                         algorithm = 'meanfield', # variational inference algorithms
                                         QR = TRUE,
                                         prior = normal())
    # QR = TRUE,
    # prior = laplace())
    
    pred_bayesian  <- posterior_linpred(bayesian_model, newdata = tmp_testing, transform = TRUE) %>%
      as_tibble() %>%
      map_dbl(~ mean(.x))
    
    tmp_probs <- pred_bayesian %>% as_tibble() %>% rename(prob_predicted = value)
    
    tmp_result <- tmp_probs %>% bind_cols(tmp_testing %>% select(clinical)) %>% mutate(id = tmp_name_fold)
    
    if (i == 1) {
      
      bayesian_result <- tmp_result
      
    } else {
      
      bayesian_result <- bayesian_result %>% bind_rows(tmp_result)
      
    }
    
  }
  return(bayesian_result)
}

# ------------------------------------------------------------------------------
# FUNCTION - RULEFIT MODEL - CROSS-VALIDATION
# ------------------------------------------------------------------------------


rulefit_cv <- function(cv_object) {
  
  tbl_result <- tibble()
  
  for (i in 1:10) {
    print(glue('Split nº {i}/10'))
    
    tmp_training <- cv_object$splits[[i]] %>% training()
    tmp_testing <- cv_object$splits[[i]] %>% testing()
    tmp_name_fold <- cv_object %>% slice(i) %>% pull(id)
    
    
    tmp_xrf <- xrf(formule_models,
                   data = tmp_training,
                   xgb_control = list(nrounds = 20, scale_pos_weight = 9),
                   family = "binomial")
    
    tmp_probs <- predict(tmp_xrf, tmp_testing) %>% as_tibble() %>% rename(prob_predicted = `1`)
    
    tmp_result <- tmp_probs %>% bind_cols(tmp_testing %>% select(clinical)) %>% mutate(id = tmp_name_fold)
    
    if (i == 1) {
      
      tbl_result <- tmp_result
    } else {
      tbl_result <- tbl_result %>% bind_rows(tmp_result)
      
      
    }
  }
  return(tbl_result)
}


# ------------------------------------------------------------------------------
# FUNCTION - EXTRACT SUPPORT AND RISK
# ------------------------------------------------------------------------------

get_rules_information <- function(data, rule) {
  
  
  tmp_filter <-   data %>%
    filter_(rule)
  
  result_support <- tmp_filter %>% nrow()
  result_risk <- tmp_filter %>%
    filter(clinical == 'pathogenic') %>%
    nrow() / result_support
  
  return(paste(result_support, result_risk))
  
}
  
# ------------------------------------------------------------------------------
# GET PREDICTIONS FROM STRUCTURE
# ------------------------------------------------------------------------------
  
  get_structure <- function(input_df, input_type = c('DEL', 'DUP')) {
    
    # Setting liftover
    library(rtracklayer)
    library(liftOver)

    from_hg19_to_hg38 = import.chain('/data-cbl/liftover/hg19ToHg38.over.chain')
    
    input_df_liftover <- clinvar_ind_deletion %>%
      mutate(tmp_id = row_number()) %>%
      select(chrom, start, end, tmp_id) %>% 
      GRanges() 
    
    seqlevelsStyle(input_df_liftover) = "UCSC"
    
    # WHY ALMOST X2
    input_df_liftover <- input_df_liftover %>%
      liftOver(from_hg19_to_hg38) %>%
      as_tibble()

    main_path <- '/data-cbl/frequena_data/rival_cnvscore/structure/StrVCTVRE-v.1.6/'
    
    use_python('/home/frequena/.conda/envs/py3.6/bin/python', required = TRUE)
    use_condaenv(condaenv = "py3.6", conda = "/usr/local/miniconda3/condabin/conda", required = TRUE)
    
    
    input_df <- input_df_liftover %>% 
      mutate(start = start - 1) %>%
      mutate(type = input_type)
      # 1-based -> 0-based
    
    write_tsv(input_df, paste0(main_path, 'test_frequena/test.bed'),
              col_names = FALSE)
    
    setwd(main_path)
    
    system(glue('python {main_path}StrVCTVRE.py \\
           -i {main_path}test_frequena/test.bed \\
           -o {main_path}test_frequena/test_annotated.bed \\
           -p {main_path}data/hg38.phyloP100way.bw \\
           -f bed')) 
    
    setwd('/data-cbl/frequena_data/cnvxplorer')
    
    output_df <- read_tsv('/data-cbl/frequena_data/rival_cnvscore/structure/StrVCTVRE-v.1.6/test_frequena/test_annotated.bed',
                          col_names = c('chrom', 'start', 'end', 'variant_type', 'prob')) %>%
      mutate(start = start + 1) # 0-based -> 1-based
    
    return(output_df)
    
  }
  

