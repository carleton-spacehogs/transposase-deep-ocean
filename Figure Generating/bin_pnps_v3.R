setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

selected <- c("percent_trans", "median_bin_pnps", "depth", "size_fraction")
# select("median_bin_pnps", "percent_trans", "log_percent_trans", "complete genome size (Mbp)", 
# "depth", "Class", "percent_biofilm", "log_percent_biofilm", "graphing_log_trans", "class_trans")

all_x <- bin_taxon[,selected] %>% filter(!is.na(depth)) 
all_y <- malaspina_bins[,selected] %>% filter(size_fraction != "error")
all_p <- filter_outliers(rbind(all_x, all_y), "percent_trans")


# used in individual_metagenome_pnps.R
scale <- c(0, 0.2, 0.4, 0.6)
p2 <- all_p %>% 
  ggplot(aes(y = percent_trans, x = median_bin_pnps)) +
  facet_wrap(~depth, ncol = 2) +
  scale_x_continuous(breaks = scale, limits=c(0, 0.65)) +
  xlab("Median MAG pN/pS") +
  ylab("% of transposase ORFs in MAG") + 
  geom_smooth(se = FALSE,method = lm) +
  scale_color_manual(values=c('green','orange', "red"))+
  geom_point(aes(color = size_fraction), alpha = 0.5) +
  labs(color='Lifestyle') +
  theme_classic()








all_p <- all_p %>% 
  filter(!is.na(median_bin_pnps)) %>%
  filter(median_bin_pnps < 0.5)

all_p %>% filter(!is.na(median_bin_pnps)) %>% with(table(depth))


# median_bin_pnps is decently normally distributed
hist(all_p$median_bin_pnps)
hist(filter_outliers(all_y,"median_bin_pnps")$median_bin_pnps, 
     xlab = "median bin pnps",
     main = "Malaspina bins median bin pnps")
hist(filter_outliers(all_x,"median_bin_pnps")$median_bin_pnps, 
     xlab = "median bin pnps",
     main = "Tara Oceans bins median bin pnps")

old.lm <- lm(log_percent_trans~depth*Class +log_percent_biofilm 
              +`complete genome size (Mbp)` + I(`complete genome size (Mbp)`^2), all_p)

pnps.lm <- update(old.lm, .~. +median_bin_pnps)
anova(pnps.lm)
depth.lm <- lm(log_percent_trans~depth, all_p)
pnps2.lm <- update(depth.lm, .~. +median_bin_pnps)
anova(pnps2.lm)






all_x <- filter_outliers(all_x, "median_bin_pnps")
color2 <- ifelse(all_x$depth=="MES","blue4",ifelse(all_x$depth=="DCM", "steelblue", "sky blue"))
scatterplot3d(all_x[,c("median_bin_pnps", "graphing_log_trans", "depth")], 
              angle = 23, pch = 20, # pch = all$shape3D, # type = "h", lty.hplot = 2,
              lab.z = 2,
              lab = c(4, 4, 1),
              y.ticklabs = c("<0.01%", "0.03%", "0.1%", "0.32%", "1%"),
              # scale.y = 0.5,
              z.ticklabs = c("     MES", "DCM", "SRF"), 
              grid=TRUE, box=TRUE, color=color2)


  

