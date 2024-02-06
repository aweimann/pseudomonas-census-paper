library(optparse)
library(BiocParallel)
library(DESeq2)
library(EnhancedVolcano)
library(magrittr)
library(plyr)
library(RColorBrewer)
library(tidyverse)
library(umap)

option_list <- list(
    make_option(c("-r", "--raw_counts"),
        type = "character", default = NULL,
        help = "table of raw read counts per gene",
        metavar = "character"
    ),
    make_option(c("-n", "--norm_counts"),
        type = "character", default = NULL,
        help = "table of normalised read counts per gene",
        metavar = "character"
    ),
    make_option(c("-m", "--metadata"),
        type = "character", default = NULL,
        help = "sample metadata",
        metavar = "character"
    ),
    make_option(c("-g", "--g2r"),
        type = "character", default = NULL,
        help = "table mapping Panaroo gene IDs to reference strain gene IDs",
        metavar = "character"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "directory for results",
        metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

raw_f <- opt$raw_counts
norm_f <- opt$norm_counts
meta_f <- opt$metadata
g2r_f <- opt$g2r
outdir <- opt$outdir


dir.create(file.path(outdir, "plots"))


############################################
## Read in and pre-process data
############################################
metadata <- read_tsv(meta_f, show_col_types = FALSE)
g2r_data <- read_tsv(g2r_f, show_col_types = FALSE)
raw_counts <- read.table(raw_f, header = TRUE, sep = "\t", row.names = 1)
norm_counts <- read.table(norm_f, header = TRUE, sep = "\t", row.names = 1)

g2r_data_sub <- g2r_data[!(
    duplicated(g2r_data$Panaroo) | duplicated(g2r_data$Panaroo, fromLast = TRUE)
), ]


############################################
## UMAP plot
############################################
norm_counts_scaled <- t(scale(t((norm_counts)), scale = FALSE, center = TRUE))

umap_res <- umap(t(norm_counts_scaled))

umap_coords <- as.data.frame(umap_res$layout)
colnames(umap_coords) <- paste0("UMAP_", c(1, 2))
umap_coords$sample_name <- rownames(umap_coords)
umap_coords <- left_join(umap_coords, metadata)

umap_plot <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
        size = 4, shape = 21, colour = "black",
        aes_string(fill = "cf_proportion")
    ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    scale_fill_gradient2(
        "CF proportion",
        limits = c(
            min(umap_coords[["cf_proportion"]]),
            max(umap_coords[["cf_proportion"]])
        ),
        low = "#80B1D3", mid = "#FDB462", high = "#FB8072",
        midpoint = (max(umap_coords[["cf_proportion"]]) + min(
            umap_coords[["cf_proportion"]]
        )) / 2,
        guide = "colourbar", aesthetics = "fill"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")


ggsave(
    umap_plot,
    file = file.path(outdir, "plots", "umap_CF_assoc.pdf"),
    device = "pdf", units = "in",
    width = 7, height = 6
)


############################################
## Assess strain clustering by CF proportion
############################################
pca_res <- prcomp(t(norm_counts_scaled))

kmeans_result <- kmeans(pca_res$x, centers = 20, nstart = 100)
cluster_assignments <- kmeans_result$cluster

cluster_df <- data.frame(cluster_assignments)
cluster_df$sample_name <- rownames(cluster_df)
cluster_df <- left_join(cluster_df, metadata)

summary_clusters <- cluster_df %>%
    group_by(cluster_assignments) %>%
    summarise(
        mean_cf_proportion = mean(cf_proportion),
        cf_proportion_sd = sd(cf_proportion),
        n = n()
    ) %>%
    arrange(cf_proportion_sd)

cat("Mean cluster size: ", mean(summary_clusters$n), "\n")

cat(
    "Mean SD of CF proportion per cluster: ",
    mean(summary_clusters$cf_proportion_sd, na.rm = TRUE), "\n"
)

randomised_sd <- sapply(1:10000, function(x) {
    cluster_df2 <- cluster_df
    cluster_df2$cf_proportion <- sample(cluster_df2$cf_proportion)
    summary_clusters <- cluster_df2 %>%
        group_by(cluster_assignments) %>%
        summarise(
            mean_cf_proportion = mean(cf_proportion),
            cf_proportion_sd = sd(cf_proportion),
            n = n()
        ) %>%
        arrange(cf_proportion_sd)
    mean(summary_clusters$cf_proportion_sd, na.rm = TRUE)
})
cat("Empirical P-value for CF diversity per cluster: ", sum(
    randomised_sd < mean(summary_clusters$cf_proportion_sd,
        na.rm = TRUE
    )
) / length(randomised_sd), "\n")


############################################
## DESeq2 model
############################################
dds <- DESeqDataSetFromMatrix(
    countData = round(raw_counts),
    colData = metadata,
    design = ~cf_proportion
)

dds <- DESeq(dds, test = "Wald", parallel = TRUE)
deseq_res <- results(dds, name = "cf_proportion", parallel = TRUE)

# add PAO1 locus tags and export results
deseq_res$Panaroo <- rownames(deseq_res)
deseq_res <- left_join(as_tibble(deseq_res), g2r_data_sub, by = "Panaroo")
deseq_res <- deseq_res[c(
    "PA_ORF", "PA_name", "Panaroo", "baseMean", "log2FoldChange",
    "lfcSE", "stat", "pvalue", "padj"
)]

write.table(
    deseq_res,
    file = file.path(outdir, "CF_assoc_DESeq2_res.tsv"),
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE
)


############################################
## Volcano plot
############################################
# For gene labels: use gene name if available; otherwise use PAO1 locus tag for
#   PAO1 genes; otherwise use Panaroo ID (for non-PAO1 genes)
deseq_res$PA_name[is.na(deseq_res$PA_name)] <- deseq_res$PA_ORF[
    is.na(deseq_res$PA_name)
]
deseq_res$PA_name[is.na(deseq_res$PA_name)] <- deseq_res$Panaroo[
    is.na(deseq_res$PA_name)
]

# set axis limits
xmin <- round_any(
    min(subset(deseq_res, padj < 0.05)$log2FoldChange), 0.1, floor
)
xmax <- round_any(
    max(subset(deseq_res, padj < 0.05)$log2FoldChange), 0.1, ceiling
)

if (abs(xmax) > abs(xmin)) {
    xmin <- xmax * -1
} else if (abs(xmin) > abs(xmax)) {
    xmax <- xmin * -1
}

ymin <- 0
ymax <- max(-log10(subset(deseq_res, padj < 0.05)$padj))

# show gene labels for genes with large log2FC or low padj
show_labs <- subset(deseq_res, padj < 0.05 & (
    (abs(log2FoldChange) > xmax * 0.6) | (abs(-log10(padj)) > ymax * 0.6)
))$PA_name

set1pal <- brewer.pal(9, "Set1")
keyvals <- ifelse(
    deseq_res$log2FoldChange < (0.5 * -1) & deseq_res$padj < 0.05,
    set1pal[3],
    ifelse(
        deseq_res$log2FoldChange > 0.5 & deseq_res$padj < 0.05,
        set1pal[1], "grey45"
    )
)
keyvals[is.na(keyvals)] <- "grey45"
names(keyvals)[keyvals == set1pal[1]] <- "up"
names(keyvals)[keyvals == "grey45"] <- "NS"
names(keyvals)[keyvals == set1pal[3]] <- "down"

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
    deseq_res,
    lab = deseq_res[["PA_name"]],
    selectLab = show_labs, x = "log2FoldChange", y = "padj",
    pointSize = 3.0, labSize = 4.0, pCutoff = 0.05, FCcutoff = 0.5,
    colCustom = keyvals, drawConnectors = TRUE, max.overlaps = 20,
    arrowheads = FALSE, min.segment.length = 0.5,
    title = "", subtitle = "", xlim = c(xmin, xmax), ylim = c(ymin, ymax)
) +
    theme(
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 18, hjust = 0.5)
    ) +
    ylab(bquote(~ -Log[10] ~ adjusted ~ italic(P)))

ggsave(
    volcano_plot,
    file = file.path(outdir, "plots", "volcano_CF_assoc.pdf"),
    device = "pdf", units = "in",
    width = 8, height = 7
)


############################################
## Plot counts for top CF-associated genes
############################################
top_genes_positive <- subset(
    deseq_res,
    (PA_ORF %in% show_labs | PA_name %in% show_labs | Panaroo %in% show_labs) &
        log2FoldChange > 0
)$Panaroo

top_genes_negative <- subset(
    deseq_res,
    (PA_ORF %in% show_labs | PA_name %in% show_labs | Panaroo %in% show_labs) &
        log2FoldChange < 0
)$Panaroo


dir.create(file.path(outdir, "plots", "positive_gene_assoc"))
dir.create(file.path(outdir, "plots", "negative_gene_assoc"))

lapply(top_genes_positive, function(gene) {
    norm_counts <- plotCounts(
        dds,
        gene = gene, intgroup = "cf_proportion", returnData = TRUE,
        normalized = TRUE, transform = TRUE
    )
    ref_dat <- subset(g2r_data_sub, Panaroo == gene)
    if (gene %in% g2r_data_sub$Panaroo) {
        if (is.na(pull(ref_dat, PA_name))) {
            id <- pull(ref_dat, PA_ORF)
        } else {
            id <- pull(ref_dat, PA_name)
        }
    } else {
        id <- gene
    }
    p1 <- ggplot(norm_counts, aes(x = cf_proportion, y = count)) +
        geom_point(size = 4, shape = 21, colour = "black", fill = "blue") +
        xlab("CF proportion") +
        ylab("Normalised count") +
        ggtitle(id) +
        theme_bw(base_size = 12) +
        geom_smooth(method = "glm.nb") +
        theme(legend.position = "bottom")
    ggsave(
        p1,
        file = file.path(
            outdir,
            "plots", "positive_gene_assoc", paste0(gene, ".pdf")
        ),
        device = "pdf", units = "in",
        width = 7, height = 6
    )
})

lapply(top_genes_negative, function(gene) {
    norm_counts <- plotCounts(
        dds,
        gene = gene, intgroup = "cf_proportion", returnData = TRUE,
        normalized = TRUE, transform = TRUE
    )
    ref_dat <- subset(g2r_data_sub, Panaroo == gene)
    if (gene %in% g2r_data_sub$Panaroo) {
        if (is.na(pull(ref_dat, PA_name))) {
            id <- pull(ref_dat, PA_ORF)
        } else {
            id <- pull(ref_dat, PA_name)
        }
    } else {
        id <- gene
    }
    p1 <- ggplot(norm_counts, aes(x = cf_proportion, y = count)) +
        geom_point(size = 4, shape = 21, colour = "black", fill = "blue") +
        xlab("CF proportion") +
        ylab("Normalised count") +
        ggtitle(id) +
        theme_bw(base_size = 12) +
        geom_smooth(method = "glm.nb") +
        theme(legend.position = "bottom")
    ggsave(
        p1,
        file = file.path(
            outdir,
            "plots", "negative_gene_assoc", paste0(gene, ".pdf")
        ),
        device = "pdf", units = "in",
        width = 7, height = 6
    )
})
