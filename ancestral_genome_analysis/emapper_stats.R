library(tidyverse)
library(forcats)
library(aggregation)
library(ggpmisc)
theme_set(theme_classic(base_size = 18))
meta <- read_tsv("samples_qc_filtered_and_annotated_v2.txt")
# read in CF proportions
cf_abund <- read_tsv("cf_per_clone_sts21.txt") %>% mutate(majority_ST = as.character(majority_ST)) %>% select(-n)
# read emapper annotations
clone2cog <- read_tsv("emapper_clone2function.txt")
# ignore unassigned function
clone2cog <- clone2cog%>% filter(COG_category != "S")
# add sample info
clone2cog <- clone2cog %>% inner_join(meta %>% rename(clone = sample_name))
# add CF info
clone_counts <- clone_counts %>% inner_join(cf_abund)
# counts by event type 
cog_counts <- clone2cog %>% group_by(majority_ST, COG_category, event_type) %>% count() %>% ungroup() %>% complete(majority_ST, COG_category, event_type, fill = list(n = 0)) %>%    inner_join(clone_counts) %>% mutate(prop = n/clone_count)
# stratify sporadic/epidemic
clone2cog <- clone2cog %>% inner_join(meta %>% rename(clone = sample_name)) %>% mutate(majority_ST = ifelse(!(majority_ST %in% cf_abund$majority_ST), "sporadic", "epidemic")) 
clone_counts <- clone2cog %>% group_by(majority_ST) %>% count() %>% rename(clone_count = n) 
cog_counts <- clone2cog %>% group_by(majority_ST, COG_category, event_type) %>% count() %>% inner_join(clone_counts) %>% mutate(prop = n/clone_count)

# fisher exact test  
fisher_res <- cog_counts %>% select(-prop, -clone_count) %>% pivot_wider(values_from = n, names_from = majority_ST) %>% filter(event_type == "gain") %>% ungroup() %>% mutate(total_epidemic = sum(epidemic), total_sporadic = sum(sporadic)) %>% rowwise() %>% mutate(pval = fisher.test(x = matrix(c(epidemic, total_epidemic-epidemic, sporadic, total_sporadic - sporadic), nrow = 2))$p.value) %>% select(-event_type)
#correct p-value
cog_descs <- read_tsv("cog_categories.txt")
fisher_res <- fisher_res %>% left_join(cog_descs)
fisher_res$padj <- p.adjust(fisher_res$pval)
fisher_res <- fisher_res %>% mutate(sig_level = case_when(
	padj < 0.001 ~ "***",
	padj < 0.01 ~ "**",
	padj < 0.05 ~ "*",
    padj >= 0.05 ~ "NS"))

segments_df <- expand.grid(COG_category = unique(fisher_res$COG_category),
                           x = 'epidemic',
                           xend = 'sporadic',
                           y = .24,
                           yend = .24)

#renormalise cog proportions to only gain events
cog_counts <- cog_counts %>% group_by(majority_ST) %>% filter(event_type == "gain") %>% mutate(prop = prop/sum(prop))
#add COG long descriptions
cog_combined <- segments_df %>% left_join(cog_descs) %>% inner_join(cog_counts) %>% inner_join(fisher_res) %>% filter(padj < 0.05)
cog_combined  %>% 
    #filter(event_type == "gain", COG_description %in% c("Defense mechanisms", "Inorganic ions", "Transcription")) %>%
    ggplot(.) + 
    geom_bar(aes(x = majority_ST, y = prop, fill = majority_ST), stat = "identity") +  
    geom_segment( aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(aes(label = sig_level, y = 0.27, x = as.factor('epidemic') ), nudge_x = 0.5) +
    ylab("Proportion of acquired annotated genes") + xlab("") +
    facet_wrap(~ COG_description)  + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(fill = "") +
    scale_fill_manual(values = c("#b3e2cd", "#cbd5e8"), guide = "none") 

ggsave("COG_barplot_props_epidemic_vs_sporadic_sel.pdf", width = 9)

# CF proportions vs COG categories
epidemic2cog <- clone2cog %>% 
	 filter(majority_ST == "epidemic") %>% select(-majority_ST) %>%
	filter(event_type == "gain") %>% select(-event_type) 
epidemic2cog_count <-epidemic2cog %>%  group_by(level7000) %>% count() %>% arrange(-n)  %>% rename(clone_count = n)
epidemic2cog <-epidemic2cog  %>% 
     mutate(COG_category = as.factor(COG_category), level7000 = as.factor(level7000)) %>% 
 	group_by(level7000, COG_category, .drop = F) %>% count() %>% ungroup()  %>%
	mutate(level7000 = as.double(as.character(level7000)))
                                        
epidemic2cog <- epidemic2cog %>% inner_join(cf_abund) %>% inner_join(epidemic2cog_count) %>% mutate(n = n)
epidemic2cog <- epidemic2cog %>% inner_join(cog_descs)

epidemic2cog %>% 
    ggplot(aes(x = cf_proportion, y = n/clone_count)) + 
    geom_smooth(method='lm') + 
    geom_point(stat = "identity") + 
    facet_wrap(~COG_description, scales = 'free_y')  + 
    theme(axis.text.x = element_text(angle = 90))  +
    theme(strip.text = element_text(size = 8)) +
   stat_poly_eq(mapping = use_label(c("R2", "P")), p.digits = 2) +
   ylab("Proportion of acquired genes")
ggsave("COG_scatter_props.pdf", width = 12)

