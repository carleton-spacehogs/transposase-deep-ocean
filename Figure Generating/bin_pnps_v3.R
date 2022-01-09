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
all_p <- filter_outliers(rbind(all_x, all_y), "median_bin_pnps")
all_p <- filter_outliers(all_p, "percent_trans")

# median_bin_pnps is decently normally distributed
hist(all_p$median_bin_pnps)

old.lm <- lm(log_percent_trans~depth*Class +log_percent_biofilm 
              +`complete genome size (Mbp)` + I(`complete genome size (Mbp)`^2), all_p)

pnps.lm <- update(old.lm, .~. +median_bin_pnps)
anova(pnps.lm)
depth.lm <- lm(log_percent_trans~depth, all_p)
pnps2.lm <- update(depth.lm, .~. +median_bin_pnps)
anova(pnps2.lm)

scale <- c(0.05, 0.15, 0.25, 0.35)
# log_scale <- c("0.05", "0.15", "0.25", "0.35")
p3 <- all_p %>% 
  ggplot(aes(y = percent_trans, x = median_bin_pnps)) +
  facet_wrap(~depth, ncol = 1) +
  scale_x_continuous(breaks = scale) +
  xlab("Median MAG pN/pS") +
  ylab("%-transposase ORF") + 
  geom_smooth(se = FALSE,method = lm) +
  geom_jitter(aes(color = class_trans), alpha = 0.5, height = 0.05) +
  scale_color_manual(labels = c("High %-transposase", "Normal %-transposase", "Low %-transposase"),
                     values = c('orange','gray', "green"))+
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y=element_blank())
# used in bin_taxon_genomesize_graph.R




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


  

