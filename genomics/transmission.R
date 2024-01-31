library(tidyverse)
t <- read_tsv("dist_aggr.txt") %>%
    rename(snp_distance = dist)
mlst <- read_tsv("samples_qc_filtered_and_annotated_v2.txt") %>% 
    mutate(infection_type = ifelse(cohort == "sickkids_sterile_sites", "acute", infection_type)) %>%
    mutate(infection_type = ifelse(cohort == "telemed", "cystic fibrosis", infection_type)) %>%
    select(sample_name, infection_type, isolation_country, isolation_city, patient, majority_ST, level7000, origin) %>% 
    filter(isolation_city != "East Midlands" | !str_detect(patient, "LiP")) 
# seed for subsampling
set.seed(1) 
mlst  <- mlst %>% filter(origin == "human") %>% mutate(infection_type_broad = ifelse(infection_type != "cystic fibrosis", "non-cf", "cf"))    %>%
# comment out for downsample/subsample non-cf
# filter(!is.na(infection_type_broad)) %>%group_by(majority_ST, patient) %>% sample_n(1) %>%  group_by(infection_type_broad)  %>% sample_n(740)
t_thresholded <- t %>% filter(snp_distance < 100) # , sample_name < sample_name_2, sample_name < sample_name_2)
t_thresholded_annotated <- t_thresholded %>% inner_join(mlst) %>% rename(sample_name_3 = sample_name) %>% rename(sample_name = sample_name_2, sample_name2 = sample_name_3) %>% inner_join(mlst, by = c("sample_name"))  %>% arrange(sample_name, sample_name2) 
# remove duplicate links between patients due to several patient isolates
duplicate_links <- bind_rows(t_thresholded_annotated26 %>%  select(patient.x, patient.y), t_thresholded_annotated26 %>% mutate(patient = patient.x, patient.x = patient.y, patient.y = patient) %>% select(patient.x, patient.y)) %>% group_by(patient.x, patient.y) %>% count() %>% filter(n > 1, patient.x > patient.y)
t_thresholded_annotated <- t_thresholded_annotated %>% anti_join(duplicate_links)
#per patient
t_thresholded_annotated <-t_thresholded_annotated %>% group_by(patient.x, patient.y) %>% filter(snp_distance == min(snp_distance)) %>% slice(1)  %>% filter(patient.x != patient.y)
# by sequence type
# sts <- mlst %>% group_by(majority_ST) %>% count() %>% arrange(-n) %>% ungroup() %>% slice(1:21) 
# sts <- sts$majority_ST
 sts <- ""
# comment out for down/subsampled analysis
#sts <- "subsampled" 
for (st in sts){
    link_stats <- matrix(nrow = 100, ncol = 7)
    for (i in 0:50)
    {
        t_thresholded_annotated_sel <- t_thresholded_annotated #%>% filter(majority_ST.x == st, majority_ST.y == st)
        links <- t_thresholded_annotated_sel  %>% 
            filter(snp_distance < !!i, isolation_country.x != isolation_country.y)
        links_same_country <- t_thresholded_annotated_sel %>% filter(snp_distance < !!i, isolation_country.x == isolation_country.y, isolation_city.x != isolation_city.y)
        links_local <- t_thresholded_annotated_sel %>% filter(snp_distance < !!i, isolation_country.x == isolation_country.y, isolation_city.x == isolation_city.y)
        links_cf_cf <- t_thresholded_annotated_sel %>% filter(snp_distance < !!i, infection_type.x == "cystic fibrosis", infection_type.y == "cystic fibrosis")
        links_non_cf_cf <- t_thresholded_annotated_sel %>% filter(snp_distance < !!i, infection_type.x == "cystic fibrosis", infection_type.y != "cystic fibrosis")
        links_cf_non_cf <- t_thresholded_annotated_sel %>% filter(snp_distance < !!i, infection_type.x != "cystic fibrosis", infection_type.y == "cystic fibrosis")
        links_non_cf_non_cf <- t_thresholded_annotated_sel %>% filter(snp_distance < !!i, infection_type.x != "cystic fibrosis", infection_type.y != "cystic fibrosis")
        link_stats[i, 1] = i
        link_stats[i, 2] = nrow(links)
        link_stats[i, 3] = nrow(links_same_country)
        link_stats[i, 4] = nrow(links_local)
        link_stats[i, 5] = nrow(links_cf_cf)
        link_stats[i, 6] = nrow(links_non_cf_cf) + nrow(links_cf_non_cf)
        link_stats[i, 7] = nrow(links_non_cf_non_cf)
    }

colnames(link_stats) <- c("snp_threshold", "different_country", "same_country", "local", "cf_cf", "non_cf_to_cf", "non_cf_to_non_cf")
# out_mod <- "cf_vs_other_lung_"

# plotting
theme_set(theme_classic())
theme_set(theme_classic(base_size = 18))
#7B8475
transmission2color <-setNames(c("#E05A5B", "#F2B76F", "#546DF7"), c("cf_cf", "non_cf_to_cf", "non_cf_to_non_cf")) 
scale_fill_transmission <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = transmission2color, 
        ... 
    )   
}
out_mod <- ""
link_stats <- as_tibble(link_stats)
link_stats %>% pivot_longer(different_country:local) %>% 
    ggplot(aes(snp_threshold, value,  fill = name)) +
    geom_bar(stat = 'identity', position = "stack")
ggsave(str_c("transmission_links_stratified_", out_mod,  st, ".png"))
link_stats %>% pivot_longer(different_country:local) %>% 
    ggplot(aes(snp_threshold, value,  fill = name)) + geom_bar(stat = 'identity', position = "fill") 
ggsave(str_c("transmission_links_stratified_proportions_", out_mod,st, ".png"))
link_stats %>% pivot_longer(cf_cf:non_cf_to_non_cf) %>% 
    ggplot(aes(snp_threshold, value,  fill = name)) +
    xlab("SNP pairwise transmission threshold")+
    ylab("Number of transmission links")+
    scale_fill_transmission(labels = c("CF <-> CF", "CF <-> Non-CF", "Non-CF <-> Non-CF"), guide = 'none') +
    geom_bar(stat = 'identity', position = "fill") 
ggsave(str_c("cf_transmission_links_stratified_proportions_", out_mod,st, ".pdf"))
link_stats %>% pivot_longer(cf_cf:non_cf_to_non_cf) %>% 
    ggplot(aes(snp_threshold, value,  fill = name)) +
    ylab("Number of transmission links")+
    xlab("SNP pairwise transmission threshold")+
    scale_fill_transmission(labels = c("CF <-> CF", "CF <-> Non-CF", "Non-CF <-> Non-CF"), name = "Transmission") +
    geom_bar(stat = 'identity', position = "stack") +
   theme(legend.position = c(0.25, 0.8))
ggsave(str_c("cf_transmission_links_stratified_", out_mod,st, ".pdf"))

}
t_thresholded_annotated26 <- t_thresholded_annotated %>% filter(snp_distance <= 26) 
mlst_patient <- read_tsv("../collect_meta/samples_qc_filtered_and_annotated_per_patient_v2.txt")
mlst_patient %>% filter((patient %in% t_thresholded_annotated26$patient.x) | (patient %in% t_thresholded_annotated26$patient.y)) %>% unite(patient, patient, majority_ST, remove = F) %>%  write_tsv("meta_threshold26_updated.txt")
t_thresholded_annotated <- t_thresholded_annotated %>% filter(patient.x != patient.y) %>%  select(patient.x, patient.y, majority_ST.x, majority_ST.y, infection_type.x, infection_type.y,  snp_distance)
t_thresholded_annotated %>% filter(snp_distance <= 26) %>% mutate(majority_ST = majority_ST.x) %>% unite(patient.y, patient.y, majority_ST, remove = F) %>% unite(patient.x, patient.x, majority_ST) %>%    select(-majority_ST.x, -majority_ST.y) %>% write_tsv("pairsnp_thresholded26_annotated_updated.txt") 
