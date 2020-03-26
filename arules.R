library(arules)
library(arulesCBA)




to_cba <- first_df %>% 
  # mutate(clinical = if_else(clinical == 'Likely benign', 'Benign', clinical)) %>%
  filter(clinical %in% c('Pathogenic', 'Benign', 'Likely benign', 'Likely pathogenic')) %>%
  mutate(clinical = as.factor(clinical)) %>%
  mutate(type_variant = if_else(type_variant == 'Deletion', 1, 0)) %>%
  mutate(type_inheritance = if_else(type_inheritance == 'De novo constitutive', 1, 0)) 
  
to_cba[,-c(1:4)] <- map_df(to_cba[,-c(1:4)], function(x) if_else(x > 0, 1, 0))

model_cba <- arulesCBA::CBA(clinical ~ n_cnv_syndromes + type_variant + type_inheritance + patho_cnv + nonpatho_cnv + disease_genes + disease_variants, 
                       data = to_cba, 
               support = 0.05, confidence = 0.8)


inspect(rules(model_cba))

result <- predict(model_cba , to_cba[to_cba$clinical == 'Benign',][1:300,])


apriori_test <- apriori(arules_sparse, parameter = list(support = 0.1, confidence = 0.5, minlen = 2))
inspect(head(apriori_test))


result_df[,-c(1:2)] <- map_df(result_df[,-c(1:2)], function(x) if_else(x > 0, 1, 0))

