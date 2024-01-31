#manhattan plot function
library(tidyverse)
library("scales")
library(ggrepel)
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-200, Inf))
}

scientific_10 <- function(x) {
      ifelse(x == 0, "0",
      parse(text=gsub("e[+]", " %*% 10^", scales::scientific_format()(x))))
}

theme_set(theme_classic(base_size = 18))
manhattan <- function(g, by){
    
    
g  %>% arrange(desc(pval)) %>% mutate(label = ifelse(padj < 1e-50, label, NA)) %>% ggplot(aes(position, pval,  label = label, size = is_sig, color = is_sig)) +#, color = !!sym(by))) +
    geom_hline(aes(yintercept = (g %>% filter(padj < 0.05) %>% ungroup() %>% summarize(max(pval)))[[1]])) + 
    geom_point() +
    geom_text_repel(size = 6) +
    scale_x_continuous(label = scientific_10) +
    scale_y_continuous(trans = reverselog_trans(10),
    	breaks = trans_breaks("log10", function(x) 10^x),
    	labels = trans_format("log10", math_format(.x)))+
    scale_colour_manual(values = c("grey", "black"), guide = 'none') +
    scale_size_manual(guide = 'none',  values = c(sig = 3, `non-sig` = 1), 1 ) +
    ylab(expression(log[10]~group("(", p-value, ")")))+
    xlab("Genomic position")

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
    #geom_hline(yintercept = 2.657577, size = 0.1) +
    xlab(log10Pe) +
    ylab(log10Po)
}

out = ""
comb <- read_tsv("sig_poisson_test_extended.txt")
comb <- comb %>% mutate(label = ifelse(!is.na(gene_name), gene_name, PAO1))
comb <- comb %>% mutate(is_sig = ifelse(padj < 0.05, "sig", "non-sig"))

manhattan(comb)
ggsave(str_c(out,"manhattan_poisson_test_extended.pdf"), width = 18)
gg_qqplot(comb$pval)
ggsave(str_c(out,"qqplot_poisson_test_extended.pdf"))
