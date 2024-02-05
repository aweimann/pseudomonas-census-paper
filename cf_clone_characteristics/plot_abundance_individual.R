library(tidyverse)

abds_norm <- read_tsv("abundances_individual_1h_normalised.txt")

abds_norm %>% filter(time != 1) %>% ggplot(aes(y = mean_theta_norm, x = as.factor(time), color = as.factor(majority_ST))) + geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.7)) + geom_point(position=position_jitterdodge(jitter.width = 0.5,dodge.width = 0.7)) + facet_wrap(~CFTR)

ggsave("abundance_over_time_individual_bp.pdf")

#t-tests
st235_4_WT <- abds_norm %>% filter(time == 4, CFTR == "WT", majority_ST == "235")
st27_4_WT <- abds_norm %>% filter(time == 4, CFTR == "WT", majority_ST == "27")
st111_4_WT <- abds_norm %>% filter(time == 4, CFTR == "WT", majority_ST == "111")
st235_4_CFTR <- abds_norm %>% filter(time == 4, CFTR == "CF", majority_ST == "235")
st27_4_CFTR <- abds_norm %>% filter(time == 4, CFTR == "CF", majority_ST == "27")
st111_4_CFTR <- abds_norm %>% filter(time == 4, CFTR == "CF", majority_ST == "111")
st235_2_WT <- abds_norm %>% filter(time == 2, CFTR == "WT", majority_ST == "235")
st27_2_WT <- abds_norm %>% filter(time == 2, CFTR == "WT", majority_ST == "27")
st111_2_WT <- abds_norm %>% filter(time == 2, CFTR == "WT", majority_ST == "111")
st235_2_CFTR <- abds_norm %>% filter(time == 2, CFTR == "CF", majority_ST == "235")
st27_2_CFTR <- abds_norm %>% filter(time == 2, CFTR == "CF", majority_ST == "27")
st111_2_CFTR <- abds_norm %>% filter(time == 2, CFTR == "CF", majority_ST == "111")

# 2h
# WT
t.test(x = st27_2_WT$mean_theta_norm,  y = c(st111_2_WT$mean_theta_norm, st235_2_WT$mean_theta_norm), var.equal = F)
t.test(x = st27_2_WT$mean_theta_norm,  y = st235_2_WT$mean_theta_norm, var.equal = F)
# 4h
# CFTR
t.test(x = st235_4_CFTR$mean_theta_norm, y = st27_4_CFTR$mean_theta_norm, var.equal = F)
t.test(x = c(st111_4_CFTR$mean_theta_norm, st235_4_CFTR$mean_theta_norm),  y = st27_4_CFTR$mean_theta_norm, var.equal = F)
# WT
t.test(x = c(st111_4_WT$mean_theta_norm, st235_4_WT$mean_theta_norm), y = st27_4_WT$mean_theta_norm, var.equal = F)
t.test(x = st235_4_WT$mean_theta_norm, y = st27_4_WT$mean_theta_norm, var.equal = F)


