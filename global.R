source('functions.R')

# system('gunzip -c local_data.RData.gz > local_data.RData')
load('local_data.RData')

ridges_home <- cnv_df %>%
  mutate(source =
           case_when(
             source == 'dgv' ~ 'DGV',
             source == 'decipher' ~ 'DECIPHER',
             source == 'gnomad_v2.1' ~ 'gnomAD v2.1',
             source == 'decipher_control' ~ 'DECIPHER Control'
           )) %>%
  ggplot(aes(length_cnv, y = source)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = source), alpha = 0.6, show.legend = FALSE, size = 1.25) +
  # geom_vline(aes(xintercept = size_cnv_query), linetype = 2, color = 'red', size = 1.5) +
  scale_x_log10() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis_d() +
  # scale_fill_manual(values = c('#CD5C5C','#32CD32', '#32CD32')) +
  xlab('log10(CNVs size)') +
  ylab('Database') +
  theme_ridges()



theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size=22))
  
}