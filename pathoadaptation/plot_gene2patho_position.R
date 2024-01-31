library(tidyverse)
library(ggplot2)

# ranking of mutations in pathoadaptive trajectoies 
patho_positions = read_tsv("gene2position.txt")

# pathoadaptive genes
sig_poisson <- read_tsv("sig_poisson_test_extended.txt") %>% select(PAO1, padj)

patho_positions_norm <- patho_positions %>% group_by(gene) %>% add_tally() %>% filter(position <= 40) %>%  group_by(position, gene) %>% count() %>% ungroup() %>% complete(position, gene, fill = list(n = 0)) %>% group_by(position) %>% mutate(position_sum = sum(n), n = n/position_sum)

patho_positions_norm <- sig_poisson %>% slice(1:40) %>% inner_join(patho_positions_norm %>% separate(gene, c("PAO1", "gene_name"), remove = F, sep = "_") )  %>% select(-PAO1, -gene_name, -padj)

patho_positions_norm %>% inner_join(patho_clust_cut) %>% ggplot(aes(position, n, color = as.factor(cluster_id))) + geom_point() + facet_wrap(~fct_reorder(gene, cluster_id)) + geom_smooth()
ggsave("patho_density_selected.pdf", width = 15, height = 15)
