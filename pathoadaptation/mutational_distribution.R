
library(tidyverse)
library("scales")
library(ggrepel)
library(forcats)
library(fuzzyjoin)
library(stringr)
library(ggtree)
library(ape)

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

#manhattan plot function
manhattan <- function(g, by){
    
    
g %>% mutate(label = ifelse(padj < 0.05, gene_name, NA)) %>% ggplot(aes(position, pval,  label = label)) +#, color = !!sym(by))) +
    geom_point(size = 2) +
    guides(size = "none")+
    geom_hline(aes(yintercept = (g %>% filter(padj < 0.05) %>% ungroup() %>% summarize(max(pval)))[[1]])) + 
    theme_light(base_size = 8)+
    geom_text_repel(label.size = 4) +
    scale_y_continuous(trans = reverselog_trans(10))+
    ylab(paste("p-value")) 

#https://slowkow.com/notes/ggplot2-qqplot/
}

#qqplot function
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  print(sort(ps))
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], "p-value"))
  log10Po <- expression(paste("Observed -log"[10], "p-value"))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 1) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    geom_hline(yintercept = 2.657577, size = 0.1) +
    theme_light(base_size = 26)+
    xlab(log10Pe) +
    ylab(log10Po)
}

meta <- read_tsv("../collect_meta/samples_st_leq2_stats.txt")
#sts <- c("182", "596", "6", "252", "595", "3", "37", "45", "5", "7")
#sts <- c("6", "3", "11", "80")
#sts <- c("252")
#sts <- c("37", "35", "40", "252")
#sts <- c("7", "45", "76")
sts <- meta$level7000
sts <- sts[!(sts %in% c(221))]
st_df <- tibble(st = sts)
st_df <- mutate(st_df, mutations = str_c("~/aw27/all_paerug/v2_all_studies/muttui/ST/v2/", st, "/variant_effect_predictions.txt")) %>%
    mutate(st_df, indels = str_c("~/aw27/all_paerug/v2_all_studies/parsimony/ST/indels/indels_reconstructed/", st, "/events.txt")) %>%
    mutate(st_df, recombination = str_c("~/aw27/all_paerug/v2_all_studies/gubbins/ST/parsed/", st, ".recombination_prediction.txt")) %>%
    mutate(st_df, recombination_pos = str_c("~/aw27/all_paerug/v2_all_studies/gubbins/ST/parsed/", st, ".recombination_pos.txt"))

v <-  st_df %>% select(st, mutations) %>% mutate(mutations = map(mutations, ~ read_tsv(., col_types = cols(.default = col_character(), pos = col_integer())))) %>% unnest(cols = c(mutations))  %>% rename(ref = upstream_allele, alt = downstream_allele)

#when looking at intergenic regions only consider upstream gene variants
#v <- v %>% filter(mutation_type != "downstream_gene_variant")

annot <- read_tsv("~/aw27/public_within_patient_cf_paerug/results/mutational_distribution/dtu/Pseudomonas_aeruginosa_PAO1_107.tsv")
annot <- annot %>% rename(PAO1 = `Locus Tag`)
v <- v %>% rename(PAO1 = locus_tag) %>% inner_join(annot) %>% rename(gene_name = `Gene Name`) 
#clade 178 in ST235 has quite a few homoplasic SNPs found in other ST235 lineages

#load effects for cohorts
effects_dtu <- read_tsv("~/aw27/public_within_patient_cf_paerug/results/snpeff/dtu/pseudo_effects_indels.txt")
effects_sickkids <- read_tsv("~/aw27/public_within_patient_cf_paerug/results/snpeff/sickkids/pseudo_effect_indels.txt")
effects_hzi <- read_tsv("~/aw27/public_across_patient_paerug/results/snpeff/hzi_amr/pseudo_effect_indels.txt")
effects_bacteremia <- read_tsv("~/aw27/public_across_patient_paerug/results/snpeff/bacteremia/pseudo_effect_indels.txt")
effects_all_pseudomonas <- read_tsv("~/aw27/all_paerug/v2_all_studies/variant_calling/pseudo_effect_indels.txt")
effect_indels <- effects_dtu %>% full_join(effects_sickkids) %>% full_join(effects_hzi) %>% full_join(effects_bacteremia) %>% full_join(effects_all_pseudomonas)

#read in indels
indels <- st_df %>% select(st, indels) %>%mutate(indels = map(indels, ~ read_tsv(., col_types = cols(.default = col_character()) )))  %>% unnest(cols = c(indels))
indels <- indels %>% separate(variant_id, c("pos", "ref", "alt"), sep = "_") %>%
    mutate(ref = toupper(ref), alt = toupper(alt), pos = as.double(pos))
indels <- indels %>% group_by(pos, ref, alt) %>% add_tally() %>% filter(n == 1)
indels <- indels %>% inner_join(effect_indels) %>% select(- parent_node, -node_state, -parent_node_state, -n) %>% rename(mutation_type = eff_type)
indels <- indels %>% inner_join(annot)

v <- bind_rows(v, indels)
#v <- v %>% filter(impact != "MODIFIER")

#remove recombination sites
#read in recombination position and interval 
recombination <- st_df %>% select(st, recombination) %>%mutate(recombination = map(recombination, ~ read_tsv(., col_types = cols(.default = col_character(), start = col_integer(), stop = col_integer()))))  %>% unnest(cols = c(recombination)) %>% rename(end = stop)
gubbins_embl <- st_df %>% select(st, recombination_pos) %>%mutate(recombination_pos = map(recombination_pos, ~ read_tsv(., col_types = cols(.default = col_character(), pos = col_integer()))))  %>% unnest(cols = c(recombination_pos)) 
recomb_pos <- gubbins_embl %>% mutate(start = pos, end = pos + 1)  %>% genome_join(recombination,  by = c("st", "start", "end")) %>% filter(node.x == node.y) %>% select(-node.y, - start.y, -st.y) %>% rename(st = st.x) %>% select(pos, "st")
#remove recombination from multi codon substitutions and standard substitutions separately as MCS positions might not agree with Gubbins positions
v_mcs <- v  %>% filter(!is.na(multi_codon_substitution)) %>% mutate(start = pos, end =  pos + 2) %>%  genome_anti_join(recomb_pos %>% mutate(start = pos, end = pos + 1), by = c("st", "start", "end"))  %>% select(-start, -end)
v <- v %>% filter(is.na(multi_codon_substitution)) %>% anti_join(recomb_pos)
v <- bind_rows(v, v_mcs)
v <- v %>% group_by(PAO1) %>% mutate(gene_name = ifelse(any(!is.na(gene_name)), na.omit(gene_name)[1], gene_name))

write_tsv(v, "all_mutations_wo_recomb.txt")



#stratify by synonymous vs. non-synonymous variants 
mod <-  v  %>% group_by(PAO1, gene_name) %>% filter(impact == "LOW") %>% count() %>% arrange(-n) 
high <-  v  %>% group_by(PAO1, gene_name) %>%   filter(impact != "LOW") %>% count() %>% arrange(-n) 

#count global number of mutations
comb <- mod %>% full_join(high, by = c("PAO1", "gene_name")) %>% mutate(n.x = ifelse(is.na(n.x), 0, n.x), n.y = ifelse(is.na(n.y), 0, n.y))
comb <- comb %>% mutate(dn_ds = n.y/n.x)  %>% mutate(dn_ds = ifelse(is.infinite(dn_ds), n.y, dn_ds ))  %>% arrange(-dn_ds) 
#stratify between point mutations and structural variation
mutation_type <- v %>% mutate(mutation_type = "point_mutation") %>%   group_by(PAO1, gene_name, mutation_type)  %>%  count()  %>% spread(mutation_type, n, fill = 0)
#stratify by impact 
impact <- v %>%   group_by(PAO1, gene_name, impact)  %>% count() %>% spread(impact, n, fill = 0)
#non-synonymous variants per ST
per_patient <- v  %>%  filter(impact != "LOW") %>%  group_by(PAO1, gene_name, st) %>% count() %>% spread(st, n, fill = 0)
#count number of STs mutations in a particular gene are found 
per_st <- v  %>%  filter(impact != "LOW") %>%  group_by(PAO1, gene_name, st) %>% count()  %>%  group_by(PAO1, gene_name) %>% count() %>% rename(no_sts = n)
comb <- comb %>% inner_join(impact) %>% inner_join(mutation_type)
uq_variants <-  v %>% group_by(PAO1, gene_name, pos, ref, alt) %>% filter(impact != "LOW") %>% 
    summarize(unique_variants = 1) %>% group_by(PAO1, gene_name) %>% 
    summarize(unique_variants = sum(unique_variants))
ref_genome_coverage <- read_tsv("reference_genes_coverage.txt")
#use fixed subset of genes to calculate the total length of coding DNA
#total_length <- (mutate(annot, gene_length = End - Start)  %>% ungroup() %>% slice_sample(n = 5300) %>%  summarize(total_length = sum(gene_length)))$total_length
annot <- ref_genome_coverage %>% rename(PAO1 = PA_ORF) %>% right_join(annot) %>% mutate(proportion = ifelse(is.na(proportion), 1, proportion))
#subtract numebr of overlapping base pairs (9972) 
total_length <- (mutate(annot, gene_length = End - Start)  %>% 
                 ungroup() %>% 
                 summarize(total_length = sum(gene_length * proportion) - 9972 ))$total_length
comb <- comb %>% inner_join(uq_variants) %>%
    inner_join(annot %>% select(-gene_name)) %>% 
    mutate(gene_length = End - Start, position = Start) %>% 
    select(PAO1, gene_name, n.x, n.y, dn_ds, gene_length,  position, HIGH:point_mutation, unique_variants, proportion) %>% 
    mutate(d.x_mod = n.y * 1000/gene_length) %>% arrange(-d.x_mod) 
no_mutations <- comb %>% ungroup() %>% summarize(sum(n.y)) %>% as_vector()
out = ""

comb <- comb   %>%  rowwise() %>% mutate(pval = poisson.test(n.y,r=no_mutations*(gene_length /(total_length)), alternative = "greater")['p.value'][[1]])
comb$padj <- p.adjust(comb$pval, method = "BH")
comb <- comb %>% inner_join(per_st)
ggplot(comb, aes(no_sts, padj, color = gene_length)) + geom_point(size = 1) + scale_y_log10() + scale_color_continuous(trans = 'log2') + geom_hline(yintercept = 0.05)
ggsave("padj_vs_no_sts.png")



#write out all patho-adaptive mutations per branch and cluster type
v %>% filter(impact != "LOW", PAO1 %in% (comb %>% filter(padj < 0.05))$PAO1) %>% select(node, PAO1, gene_name, st) %>% write_tsv("all_mutations_per_node.txt")
comb %>%    
    ungroup() %>% 
    slice(1:30)  %>%
    mutate(`non-unique variants` = n.y - unique_variants) %>% 
    gather(key, value, `non-unique variants`, unique_variants) %>% 
    mutate(PAO1_w_gene = ifelse(is.na(gene_name), PAO1, paste(PAO1, gene_name, sep = "_"))) %>% 
    ggplot(aes(fct_reorder(PAO1_w_gene, padj), value, fill = key))  +
    geom_bar(stat = "identity") + 
    coord_flip() +
    labs(x = "locus tag/gene name", fill = "") + 
    theme(axis.text.y = element_text(size = 10))
ggsave(str_c(out,"barplot_unique_mutations.png"))
#extended prefix means using Panaroo of ~800 isolates to infer total length of coding DNA rather than using subset of  5000 genes
comb <- comb  %>% arrange(padj) %>%  write_tsv(str_c(out,"sig_poisson_test_extended.txt"))

manhattan(comb)
ggsave(str_c(out,"manhattan_poisson_test_extended.png"), width = 18)
gg_qqplot(comb$pval)
ggsave(str_c(out,"qqplot_poisson_test_extended.png"))
