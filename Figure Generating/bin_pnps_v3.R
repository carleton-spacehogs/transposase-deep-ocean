setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

selected <- c("percent_trans", "log_percent_trans", "median_bin_pnps", "depth", 
              "size_fraction", "complete genome size (Mbp)")

all_x <- bin_taxon[,selected] %>% filter(!is.na(depth)) 
all_y <- malaspina_bins[,selected] %>% filter(size_fraction != "error")
all_x <- filter_outliers(all_x, "median_bin_pnps")
all_y <- filter_outliers(all_y, "median_bin_pnps")
all_p <- filter_outliers(rbind(all_x, all_y), "percent_trans")

# all_p$genome_size_plot <- ifelse(all_p$`complete genome size (Mbp)` > 7, 7, all_p$`complete genome size (Mbp)`)

# used in individual_metagenome_pnps.R
# scale <- c(0, 0.2, 0.4, 0.6)
p2 <- all_p %>% 
  ggplot(aes(y = percent_trans, x = median_bin_pnps)) +
  facet_wrap(~depth, ncol = 2, scales = "free") +
  # scale_x_continuous(breaks = scale, limits=c(0, 0.65)) +
  xlab("Median MAG pN/pS") +
  ylab("% of transposase ORFs in MAG") + 
  geom_smooth(se = FALSE,method = lm) +
  scale_color_manual(values=c('green','orange', "red"))+
  geom_point(aes(color = size_fraction)) +
  labs(color='Lifestyle') +
  theme_classic()


pnps_depth.lm <- lm(percent_trans~depth*median_bin_pnps, all_p)
summary(pnps_depth.lm)


# median_bin_pnps is decently normally distributed
hist(filter_outliers(all_x,"median_bin_pnps")$median_bin_pnps, 
     xlab = "median bin pnps",
     main = "Tara Oceans bins median bin pnps")
hist(filter_outliers(all_y,"median_bin_pnps")$median_bin_pnps, 
     xlab = "median bin pnps",
     main = "Malaspina bins median bin pnps")
hist(all_p$median_bin_pnps)

