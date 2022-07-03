setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins, low_trans, high_trans) %=% init_bins()

selected <- c("percent_trans", "log_percent_trans", "median_bin_pnps", "depth", 
              "size_fraction", "complete genome size (Mbp)")

all_x <- bin_taxon[,selected] %>% filter(!is.na(depth)) 
all_y <- malaspina_bins[,selected] %>% filter(size_fraction != "error")

all_p <- rbind(all_x, all_y) %>%
  mutate(log_bin_pnps = log10(median_bin_pnps)) %>%
  # -2 is for looks nice on graph
  mutate(log_trans_graph = ifelse(percent_trans == 0, -2, log10(percent_trans))) %>%
  # -5 is an arbitrary small number that's close to 0
  mutate(log_trans_reg = ifelse(percent_trans == 0, -5, log10(percent_trans)))

# used in individual_metagenome_pnps.R
trans_percent <- c(0.01, 0.03, 0.09, 0.27, 0.81)
trans_scale <- log10(trans_percent)
trans_percent <- as.character(trans_percent)
trans_percent[1] = "< 0.01"
bin_pnps <- c(0.033, 0.1, 0.3, 0.9, 2.7)
pnps_scale <- log10(bin_pnps)

p2 <- all_p %>% 
  ggplot(aes(y = log_trans_graph, x = log_bin_pnps)) +
  facet_wrap(~depth, ncol = 2) + # , scales = "free"
  scale_y_continuous(breaks = trans_scale, labels = trans_percent) + 
  scale_x_continuous(breaks = pnps_scale, labels = bin_pnps) + 
  xlab("Median MAG pN/pS") +
  ylab("% of transposase ORFs in MAG") + 
  geom_smooth(se = FALSE,method = lm) +
  scale_color_manual(values=c('green','gray','orange'))+
  geom_point(aes(color = size_fraction)) +
  labs(color='Lifestyle') +
  theme_classic()

pnps_depth.lm <- lm(log_trans_reg~depth*log_bin_pnps, all_p)
summary(pnps_depth.lm)

# log(median_bin_pnps) itself does not have a significant 
# correlation with the log(transposase abundance in bins),
# even though the interaction terms are significant.


# median_bin_pnps is decently normally distributed
hist(filter_outliers(all_x,"median_bin_pnps")$median_bin_pnps, 
     xlab = "median bin pnps",
     main = "Tara Oceans bins median bin pnps")
hist(filter_outliers(all_y,"median_bin_pnps")$median_bin_pnps, 
     xlab = "median bin pnps",
     main = "Malaspina bins median bin pnps")
hist(all_p$median_bin_pnps)

all_x <- filter_outliers(all_x, "median_bin_pnps")
all_y <- filter_outliers(all_y, "median_bin_pnps")
all_p <- filter_outliers(rbind(all_x, all_y), "percent_trans")
