#Plots BSPs for each of the 7 DCCs
#Run on the 7 DCC BSPs from calculateBayesianSkyline.py

args <- commandArgs(TRUE)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)

library(ggplot2)
library(gridExtra)
library(LaplacesDemon)


#meta <- read_tsv("../real_root_age.txt")

meta <- read_tsv("real_root_age.txt")
#clones <- c(179, 252, 274, 308, 348, 358,  381, 395, 446 )

clones <- meta$st 
#Import the DCC BSPs
#dccBSPs <- lapply(clones, function(clone){read.table(paste("../bayesian_skyline/", clone, "_combined.txt", sep = ""), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)})
dccBSPs <- lapply(clones, function(clone){read.table(paste("bayesian_skyline_v2/", clone, "_combined.txt", sep = ""), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)})
dccsFiltered <- lapply(dccBSPs, function(x){x[, 2:ncol(x)]})
clone_names <- as.list(clones)

get_dcc_stats <- function(dcc_filtered_list, probs = c(0.05, 0.95)){
        dcc_filtered <- dcc_filtered_list[[1]]
        clone_name <- dcc_filtered_list[[2]]
    dcc <- data.frame(matrix(0, ncol = 5, nrow = length(dcc_filtered[1,])))
    names(dcc) <- c("Date", "Median", "HPD.min", "HPD.max", "st")
    dcc[, "st"] <- clone_name
    for (year in 1:length(dcc_filtered[1, ])) {
      dcc[year,1] <- as.numeric(names(dcc_filtered)[year])
      dcc[year,2] <- median(dcc_filtered[,year])
      yearInterval <- quantile(dcc_filtered[, year], probs = probs, na.rm = T)
      #yearInterval <- p.interval(dcc_filtered[,year], prob = 0.95)
      #dcc[year,2] <- yearInterval[2] 
      dcc[year,3] <- yearInterval[1]
      dcc[year,4] <- yearInterval[2]}
      dcc
}
dcc_stats <- map(transpose(list(dccsFiltered, clone_names)), get_dcc_stats, probs = c(0.05, 0.95))
#Plot the BSPs

for (n in 2:ncol(dccBSPs[[1]])){
    
    print(sum(dccBSPs[[1]][, n] != dccBSPs[[1]][, n - 1]))


}

plot_dcc <- function(dcc_list, is_log = T)
{
    dcc <- dcc_list[[1]]
    clone_name <- dcc_list[[2]]
    pdf(paste(clone_name, ".pdf", sep = ""))
    dcc_plot <- ggplot(dcc) + theme_bw() + geom_line(aes(x = Date, y = Median)) + geom_ribbon(aes(x = Date, ymin = HPD.min, ymax = HPD.max), colour = "grey70", alpha = 0.5) + ggtitle(clone_name) + theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5))  + 
        #xlim(1800, 2005) + 
        xlim(1800, 2005) 
    if(is_log){
        dcc_plot <- dcc_plot + 
            scale_y_log10(
                              breaks = scales::trans_breaks("log2", function(x) 2^x),
                                 labels = scales::trans_format("log2", scales::math_format(2^.x))
                               ) 
            #scale_y_log10(limits = c(0.001, 100000))
    }
    print(dcc_plot)
    dev.off()
    dcc_plot
}

pdf("dccs_combined_log.pdf")
dcc_plots_log <- map(transpose(list(dcc_stats, clone_names)),  plot_dcc, is_log = T)
combinedPlot <- grid.arrange(grobs = dcc_plots_log, nrow = 4, left = "Relative genetic diversity", bottom = "Date")
print(combinedPlot)
dev.off()

pdf("dccs_combined.pdf")
dcc_plots <- map(transpose(list(dcc_stats, clone_names)),  plot_dcc, is_log = F)
combinedPlot <- grid.arrange(grobs = dcc_plots, nrow = 3, left = "Relative genetic diversity", bottom = "Date")

print(combinedPlot)
dev.off()

