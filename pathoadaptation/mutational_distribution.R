
library(tidyverse)


v <- read_tsv("all_mutations_wo_recomb.txt")
v <- v %>%  filter(mutation_type != "upstream_gene_variant" & mutation_type != "downstream_gene_variant")
annot <- read_tsv("Pseudomonas_aeruginosa_PAO1_107.tsv")
annot <- annot %>% rename(PAO1 = `Locus Tag`)



# stratify by synonymous vs. non-synonymous variants 
mod <-  v  %>% group_by(PAO1, gene_name) %>% filter(impact == "LOW") %>% count() %>% arrange(-n) 
high <-  v  %>% group_by(PAO1, gene_name) %>%   filter(impact != "LOW") %>% count() %>% arrange(-n) 

# count global number of mutations
comb <- mod %>% full_join(high, by = c("PAO1", "gene_name")) %>% mutate(n.x = ifelse(is.na(n.x), 0, n.x), n.y = ifelse(is.na(n.y), 0, n.y))
comb <- comb %>% mutate(dn_ds = n.y/n.x)  %>% mutate(dn_ds = ifelse(is.infinite(dn_ds), n.y, dn_ds ))  %>% arrange(-dn_ds) 
#stratify between point mutations and structural variation
mutation_type <- v %>% mutate(mutation_type = ifelse(str_length(ref) == str_length(alt), "point_mutations", "structural_variants")) %>%   group_by(PAO1, gene_name, mutation_type)  %>%  count()  %>% spread(mutation_type, n, fill = 0)
# stratify by impact 
impact <- v %>%   group_by(PAO1, gene_name, impact)  %>% count() %>% spread(impact, n, fill = 0)
# non-synonymous variants per ST
per_patient <- v  %>%  filter(impact != "LOW") %>%  group_by(PAO1, gene_name, st) %>% count() %>% spread(st, n, fill = 0)
# count number of STs mutations in a particular gene are found 
per_st <- v  %>%  filter(impact != "LOW") %>%  group_by(PAO1, gene_name, st) %>% count()  %>%  group_by(PAO1, gene_name) %>% count() %>% rename(no_sts = n)
comb <- comb %>% inner_join(impact) %>% inner_join(mutation_type)
uq_variants <-  v %>% group_by(PAO1, gene_name, pos, ref, alt) %>% filter(impact != "LOW") %>% 
summarize(unique_variants = 1) %>% group_by(PAO1, gene_name) %>% 
summarize(unique_variants = sum(unique_variants))
ref_genome_coverage <- read_tsv("reference_genes_coverage.txt")
annot <- ref_genome_coverage %>% rename(PAO1 = PA_ORF) %>% right_join(annot) %>% mutate(proportion = ifelse(is.na(proportion), 1, proportion))

# subtract numebr of overlapping base pairs (9972) 
total_length <- (mutate(annot, gene_length = End - Start)  %>% 
             ungroup() %>% 
             summarize(total_length = sum(gene_length * proportion) - 9972 ))$total_length

comb <- comb %>% inner_join(uq_variants) %>%
inner_join(annot %>% select(-gene_name)) %>% 
mutate(gene_length = End - Start, position = Start) %>% 
select(PAO1, gene_name, n.x, n.y, dn_ds, gene_length,  position, HIGH:structural_variants, unique_variants, proportion) %>% 
mutate(d.x_mod = n.y * 1000/gene_length) %>% arrange(-d.x_mod) 
no_mutations <- comb %>% ungroup() %>% summarize(sum(n.y)) %>% as_vector()
out = ""

comb <- comb   %>%  rowwise() %>% mutate(pval = poisson.test(n.y,r=no_mutations*(gene_length /(total_length)), alternative = "greater")['p.value'][[1]])
comb$padj <- p.adjust(comb$pval, method = "BH")
comb <- comb %>% inner_join(per_st)

comb <- comb  %>% arrange(padj) %>%  write_tsv(str_c(out,"sig_poisson_test_extended.txt"))

