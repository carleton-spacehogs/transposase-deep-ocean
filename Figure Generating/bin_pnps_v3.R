setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

all_x <- bin_taxon %>% # filter_outliers(bin_taxon, "percent_trans") %>% 
  select("median_bin_pnps", "percent_trans", "log_percent_trans", "complete genome size (Mbp)",
         "depth", "Class", "percent_biofilm", "log_percent_biofilm", "graphing_log_trans",
         "class_trans") %>% filter(!is.na(depth)) 

all_y <- malaspina_bins %>% 
  select("median_bin_pnps", "percent_trans", "log_percent_trans", "complete genome size (Mbp)",
         "depth", "Class", "percent_biofilm", "log_percent_biofilm", "graphing_log_trans",
         "class_trans")
all_p <- filter_outliers(rbind(all_x, all_y), "percent_trans")

scale <- c(0.05, 0.15, 0.25, 0.35, 0.45)
# log_scale <- c("0.05", "0.15", "0.25", "0.35")
# supplementary 5
all_p %>% 
  ggplot(aes(y = percent_trans, x = median_bin_pnps)) +
  facet_wrap(~depth, ncol = 2) +
  scale_x_continuous(breaks = scale, limits=c(0, 0.6)) +
  xlab("Median MAG pN/pS") +
  ylab("% of transposase ORFs\n(among all ORFs in a MAG)") + 
  geom_smooth(se = FALSE,method = lm) +
  geom_point(aes(color = class_trans), alpha = 0.5) +
  scale_color_manual(labels = c("High %-transposase (\u0251- \u03b2- \u03b3- proteobact., and Actinobact.)", 
                                "Normal transposase abundance", 
                                "Low %-transposase (Flavobact., Acidimicrobidae, etc.)"),
                     values = c('orange','gray', "green"))+
  labs(color='Taxonomical Class of:') +
  theme_classic() +
  theme(legend.position="bottom",
        legend.direction="vertical")

ggsave("supplementary_bin_pnps_trans.png", plot = last_plot())

# all_p <- filter_outliers(rbind(all_x, all_y), "median_bin_pnps")
# all_p <- rbind(filter_outliers(all_x,"median_bin_pnps"), filter_outliers(all_y,"median_bin_pnps"))

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


  

