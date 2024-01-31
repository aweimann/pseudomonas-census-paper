library(tidyverse)
theme_set(theme_classic(base_size = 18))
# pathoadaptation per sample
patho_table <- read.table("patho_per_sample.txt", sep = "\t", header = 1, row.names = 1, comment.char = "")
# meta data for samples
meta <- read_tsv("samples_qc_filtered_and_annotated_v2.txt")
# UMAP coords
umap_tibble <- read_csv("umap_coords_jaccard_per_patient.txt") 
umap_tibble <- meta %>% select(-ST)%>% inner_join(umap_tibble) 
umap_tibble <-  as_tibble(patho_table, rownames = "sample_name") %>% inner_join(umap_tibble)
umap_tibble <- inner_join(umap_tibble, as_tibble(apply(patho_table, 1, sum), rownames = "sample_name")) %>% rename(no_patho = value)
# omit samples without any pathoadaptive mutations
umap_tibble <- filter(umap_tibble, no_patho > 0)
# define infection types
umap_tibble <- umap_tibble  %>%  group_by(infection_type) %>% add_tally()  %>% 
    mutate(infection_type = case_when(
        infection_type == "cystic fibrosis" ~ "cystic fibrosis",
        infection_type == "respiratory tract" ~ "non-CF lung",
        infection_type == "bronchiectasis" ~ "bronchiectasis",
        infection_type == "pneumonia" ~ "non-CF lung",
        infection_type == "copd" ~ "non-CF lung",
        infection_type == "pcd" ~ "non-CF lung",
        infection_type == "epyema" ~ "non-CF lung",
        infection_type == "chronic obstructive bronchopneumopathy" ~ "non-CF lung",
        is.na(infection_type) ~ "unknown",
        TRUE  ~ "non-CF other")) 
# UMAP infection type centroids
infection_means <- umap_tibble %>% filter(infection_type != "unknown") %>% filter(!is.na(infection_type)) %>%  group_by(infection_type) %>% summarise(umap.1 = mean(umap.1), umap.2 = mean(umap.2)) 
# infection type annotated UMAP
umap_tibble  %>%  ggplot() + geom_point(aes(umap.1, umap.2, fill = as.factor(infection_type),  color = as.factor(infection_type)), show.legend = T, alpha = 0.5, shape = 21)  +
    geom_point(data = infection_means, size = 5, aes(umap.1, umap.2, fill = infection_type), alpha = 0.75, shape = 21, color = 'black') +
    labs(fill = "") + 
    guides(colour = "none")
ggsave("umap_patho_infection_type_per_patient.pdf", width = 9)

# burden annotated UMAP
umap_tibble %>% ggplot() + geom_point(aes(umap.1, umap.2, color = no_patho), show.legend = T, alpha = 0.25) +
    scale_color_viridis_c(option = "turbo") 
ggsave("umap_patho_patho_burden_per_patient.pdf", width = 9)


