library(tidyverse)
library(broom)


combined <- read_tsv("measurements_all_vf.txt")

combined %>% ggplot(aes(cf_proportion, halo_width)) + geom_point() + facet_wrap(~assay, scales = "free_y")
ggsave("scatter_all_strains_cf_proportions.pdf")

combined_wide <- combined %>% pivot_wider(names_from = assay, values_from = halo_width)
combined_wide <- combined_wide %>% mutate(siderophore = ifelse(is.na(siderophore), mean(siderophore, na.rm = T), siderophore))
combined_feats <- as.matrix(combined_wide %>%ungroup() %>%select(caseinase, gelatinase, swim, siderophore, biofilm, twitch))
rownames(combined_feats) <- combined_wide$sample_name
feats_pca <- prcomp(combined_feats, scale = T, center = T)
loadings <-as_tibble(feats_pca$rotation, rownames = "assay")
pca_coords <- as_tibble(feats_pca$x)
pca_coords$sample_name <- rownames(feats_pca$x)
combined_wide <- combined_wide %>% inner_join(pca_coords)
combined_wide %>% ggplot(aes(PC1, PC2, color = cf_proportion)) + geom_point() +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1 * 2),
     yend = (PC2 * 2)), arrow = arrow(length = unit(1/2, "picas")),
     color = "black") +
  geom_point(size = 1) +
  annotate("text", x = (loadings$PC1*2), y = (loadings$PC2*2),
     label = loadings$assay) +
     xlab("PC1 (48.7% variance explained)") + ylab("PC2 (18.8% variance explained)")
ggsave("pca_all_strains_w_cf_proportions.pdf")

# test association of phenotype with CF proportions
combined %>% 
    ungroup()  %>% 
    nest_by(assay) %>% 
    mutate(mod = list(lm(cf_proportion ~ halo_width, data = data))) %>% 
    reframe(tidy(mod)) %>% 
    filter(term == "halo_width") %>% 
    arrange(p.value) %>% 
    mutate(padj = p.adjust(p.value, method = "BH"))

