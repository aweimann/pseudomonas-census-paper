library(tidyverse)
type2color <-setNames(c("#fb8072", "#80b1d3"), c("cf", "non_cf")) 

scale_fill_type <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = type2color, 
        ... 
    )   
}
hc_clust <- read_tsv("hc_clustering_upgma.txt")
mlst <- read_tsv("samples_qc_filtered_and_annotated_per_patient_v2.txt")  %>% filter(origin == "human") 
cf_counts <- hc_clust %>% inner_join(mlst) %>% group_by(level7000, majority_ST, infection_type) %>% filter(!is.na(infection_type), !is.na(collection_date), infection_type == "cystic fibrosis") %>% count()
mlst  <- mlst %>% group_by(majority_ST) %>% add_tally() %>% filter(n > 30, !is.na(collection_date)) %>% select(-n) %>% mutate(infection_type_broad = ifelse(infection_type != "cystic fibrosis", "non-cf", "cf"))  #%>%

clone_counts <-  hc_clust %>% group_by(level7000, majority_ST)%>% inner_join(mlst)%>%  count()
cf_counts <- cf_counts%>% rename(cf_count = n) %>% right_join(clone_counts)  %>% mutate(cf_count = ifelse(is.na(cf_count), 0, cf_count)) %>% mutate(cf_proportion = cf_count/n) 
cf_non_cf_counts <-cf_counts %>% filter(n >= 20) %>%  mutate(non_cf = n - cf_count, cf = cf_count) 
cf_non_cf_counts %>% pivot_longer(cf:non_cf) %>% 
    ggplot(aes(x = fct_reorder(majority_ST, cf_proportion), y = value, fill = name)) + 
    geom_bar(stat = "identity", position = "fill") +
    labs(fill = "") + 
    xlab("Majority ST") + 
    ylab("Proportion of CF/Non-CF patients") + 
    scale_fill_type( labels = c("CF", "Non-CF"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("cf_per_clone.pdf")

#Turton et al data
cf_counts_turton <- tibble(st = c(17, 27, 395, 111, 146, 155, 253), cf = c(88, 61, 29, 0, 124, 18, 39), non_cf = c(43, 10, 23, 36, 0, 10, 84))
cf_counts_turton %>% mutate(cf_proportion = cf/non_cf) %>% pivot_longer(cf:non_cf, names_to = "Infection type") %>% ggplot(aes(x = fct_reorder(as.factor(st), cf_proportion), y = value, fill = `Infection type`)) + geom_bar(stat = "identity", position = "fill") + xlab("Majority ST") + 
    scale_fill_type()
ggsave("cf_per_clone_prop_selected_turton.pdf")




mlst <- read_tsv("samples_qc_filtered_and_annotated_per_patient_v2.txt")  %>% filter(origin == "human") 
cf_counts <- hc_clust %>% inner_join(mlst) %>% group_by(level7000, majority_ST, infection_type) %>% filter(!is.na(infection_type), !is.na(collection_date), infection_type == "cystic fibrosis") %>% count()
#mlst  <- mlst %>% mutate(infection_type_broad = ifelse(infection_type != "cystic fibrosis", "non-cf", "cf"))  #%>%
mlst <- mlst  %>%  
    mutate(infection_type = case_when(
        infection_type == "cystic fibrosis" ~ "cystic fibrosis",
        infection_type == "respiratory tract" ~ "non-CF lung",
        infection_type == "bronchiectasis" ~ "non-CF lung",
        infection_type == "pneumonia" ~ "non-CF lung",
        infection_type == "copd" ~ "non-CF lung",
        infection_type == "pcd" ~ "non-CF lung",
        infection_type == "epyema" ~ "non-CF lung",
        infection_type == "chronic obstructive bronchopneumopath" ~ "non-CF lung",
        is.na(infection_type) ~ "unknown",
        TRUE  ~ "non-CF other")) 
mlst_sel <- mlst %>% filter(infection_type != "unknown", infection_type != "non-CF other") 



# scatter plots

# non-cf lung v non-cf other
cf_counts <- hc_clust %>% inner_join(mlst) %>% group_by(level7000, majority_ST, infection_type) %>% filter(!is.na(infection_type), !is.na(collection_date), infection_type == "cystic fibrosis") %>% count()
mlst_sel  <- mlst %>% group_by(majority_ST) %>% add_tally() %>% filter(n > 30, !is.na(collection_date), infection_type != "non-CF lung", infection_type != "unknown") %>% select(-n) %>% mutate(infection_type_broad = ifelse(infection_type != "cystic fibrosis", "non-cf", "cf"))  #%>%





clone_counts <-  hc_clust %>% group_by(level7000, majority_ST)%>% inner_join(mlst_sel)%>%  count()
cf_counts <- cf_counts%>% rename(cf_count = n) %>% right_join(clone_counts)  %>% mutate(cf_count = ifelse(is.na(cf_count), 0, cf_count)) %>% mutate(cf_proportion = cf_count/n) 
cf_non_cf_other_counts <-cf_counts  %>%  mutate(non_cf = n - cf_count, cf = cf_count) 


infection_counts <- mlst %>% group_by(majority_ST) %>% add_tally() %>%  filter(n > 30, !is.na(collection_date), infection_type != "unknown") %>% select(-n) %>% add_tally()%>% rename(total = n) 
infection_counts <- infection_counts %>% group_by(majority_ST, infection_type, total) %>% count() %>% inner_join(cf_non_cf_counts %>% ungroup() %>% select(majority_ST, cf_proportion))  %>% mutate(n = n/total)
infection_counts %>% ggplot(aes(fct_reorder(majority_ST, cf_proportion), n, fill = infection_type)) + geom_bar(stat = "identity", position = "dodge")
ggsave("cf_proportions_w_non_cf_lung.pdf")

cf_proportions <- cf_non_cf_lung_counts %>% ungroup() %>% rename(cf_proportion_non_cf_lung = cf_proportion) %>% select(majority_ST, cf_proportion_non_cf_lung)  %>% inner_join(cf_non_cf_other_counts) %>% select(majority_ST, cf_proportion, cf_proportion_non_cf_lung) %>% rename(cf_proportion_non_cf_other = cf_proportion)
non_cf_lung_proportions <- non_cf_non_cf_lung_counts %>% ungroup() %>% rename(non_cf_lung_proportion_non_cf_other = cf_proportion) %>% select(majority_ST, non_cf_lung_proportion_non_cf_other)

cf_proportions <- cf_proportions %>% inner_join(non_cf_lung_proportions)


summary(lm(cf_proportion_non_cf_other ~ cf_proportion_non_cf_lung, cf_proportions)) 

cf_proportions %>% ggplot(aes(cf_proportion_non_cf_other, cf_proportion_non_cf_lung))  + geom_point()+  geom_smooth(method = "lm") + xlab("CF/non-CF other - ratio") + ylab("CF/non-CF lung - ratio")
ggsave("cf_vs_non_cf_other_scatter.pdf")


# turton 
cf_proportions <- cf_counts_turton %>% mutate(majority_ST =  as.character(st)) %>% mutate(cf_proportion_turton = cf/(cf+non_cf)) %>% select(cf_proportion_turton, majority_ST) %>%  inner_join(cf_non_cf_counts)
summary(lm(cf_proportion ~ cf_proportion_turton, cf_proportions))
cf_proportions %>% ggplot(aes(cf_proportion, cf_proportion_turton)) + geom_point() + geom_smooth(method = 'lm')
ggsave("cf_turton_vs_cf_our_scatter.pdf")


