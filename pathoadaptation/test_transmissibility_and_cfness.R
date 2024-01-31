library(tidyverse)

# pathoadaptive genes
comb <- read_tsv("sig_poisson_test_extended.txt")

# all mutations without recombination
v <-  read_tsv("all_mutations_wo_recomb.txt") %>%
    filter(impact != "LOW") %>% 
    filter(impact != "MODIFIER")

# sample meta data
meta <- read_tsv("samples_qc_filtered_and_annotated_v2.txt")

# node to infection type/transmitted
node2type <- read_tsv("node2type.txt") 
v <- v %>% rename(level7000 =st) %>% left_join(node2type) %>% rename(sample_name = node) %>% left_join(meta %>% select(-infection_type))

#transmissible CF vs transmissible non-CF
# CF
v_inf <- v %>% filter(!is.na(infection_type))
cf <- v_inf   %>% filter( infection_type == "cystic fibrosis") 
# non-CF
non_cf <-  v_inf   %>% filter(infection_type != "cystic fibrosis", !is.na(infection_type), infection_type != "multiple")



cf <- cf  %>% group_by(PAO1, gene_name) %>% count()  %>% arrange(-n) 
non_cf <- non_cf  %>% group_by(PAO1, gene_name) %>% count()  %>% arrange(-n) 
non_cf_count <- non_cf  %>% inner_join(comb %>% filter(padj < 0.05)) %>% ungroup() %>%  summarize(sum(n))
cf_count <- cf  %>% right_join(comb %>% filter(padj < 0.05)) %>% replace_na(list(n = 0)) %>%  ungroup() %>%  summarize(sum(n))
cf_vs_other <- cf %>% rename(cf = n) %>% full_join(non_cf %>% rename(non_cf = n)) %>%replace_na(list(cf = 0, non_cf = 0)) %>%  mutate(cf_ratio = cf/non_cf) %>% inner_join(comb) %>% filter(padj < 0.05)  %>% arrange(padj) %>% select(PAO1, gene_name, cf_ratio, cf, non_cf)
cf_vs_other <- cf_vs_other %>% mutate(cf_all = cf_count$`sum(n)`, non_cf_all = non_cf_count$`sum(n)`)
cf_vs_other_f_test <- cf_vs_other %>% rowwise() %>% mutate(pval = fisher.test(matrix(c(cf, non_cf, cf_all - cf, non_cf_all - non_cf), nrow = 2), alternative = "two.sided")$p.value, odds_ratio = fisher.test(matrix(c(cf, non_cf, cf_all - cf, non_cf_all - non_cf), nrow = 2), alternative = "two.sided")$estimate) 
cf_vs_other_f_test$padj <- p.adjust(cf_vs_other_f_test$pval, method = "BH")
cf_vs_other_f_test %>% arrange(padj) %>% write_tsv("cf_vs_non_cf.txt")

# transmissible vs non-transmissible
# transmitted
cf <- v   %>% filter( is_cf_transmissible| is_non_cf_transmissible) 
# not yet transmitted
non_cf <-  v   %>% filter(is_within_patient)

cf <- cf  %>% group_by(PAO1, gene_name) %>% count()  %>% arrange(-n) 
non_cf <- non_cf  %>% group_by(PAO1, gene_name) %>% count()  %>% arrange(-n) 
non_cf_count <- non_cf  %>% inner_join(comb %>% filter(padj < 0.05)) %>% ungroup() %>%  summarize(sum(n))
cf_count <- cf  %>% right_join(comb %>% filter(padj < 0.05)) %>% replace_na(list(n = 0)) %>%  ungroup() %>%  summarize(sum(n))
cf_vs_other <- cf %>% rename(cf = n) %>% full_join(non_cf %>% rename(non_cf = n)) %>%replace_na(list(cf = 0, non_cf = 0)) %>%  mutate(cf_ratio = cf/non_cf) %>% inner_join(comb) %>% filter(padj < 0.05)  %>% arrange(padj) %>% select(PAO1, gene_name, cf_ratio, cf, non_cf)
cf_vs_other <- cf_vs_other %>% mutate(cf_all = cf_count$`sum(n)`, non_cf_all = non_cf_count$`sum(n)`)
cf_vs_other_f_test <- cf_vs_other %>% rowwise() %>% mutate(pval = fisher.test(matrix(c(cf, non_cf, cf_all - cf, non_cf_all - non_cf), nrow = 2), alternative = "two.sided")$p.value, odds_ratio = fisher.test(matrix(c(cf, non_cf, cf_all - cf, non_cf_all - non_cf), nrow = 2), alternative = "two.sided")$estimate) 
cf_vs_other_f_test$padj <- p.adjust(cf_vs_other_f_test$pval, method = "BH")
cf_vs_other_f_test %>% arrange(padj) %>% write_tsv("transmitted_vs_non_transmitted.txt")
