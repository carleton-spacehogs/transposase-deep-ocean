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

# median_bin_pnps is decently normally distributed
hist(all_p$median_bin_pnps)

p3 <- filter_outliers(all_p, "percent_trans") %>% 
  ggplot(aes(y = percent_trans, x = median_bin_pnps)) +
  facet_wrap(~depth, ncol = 4) +
  # scale_y_continuous(breaks = scale, labels = log_scale) +
  xlab("Median MAG pN/pS") +
  ylab("%-transposase ORF") + 
  geom_jitter(aes(color = class_trans), alpha = 0.5, height = 0.05) +
  scale_color_manual(labels = c("High %-transposase", "Normal %-transposase", "Low %-transposase"),
                     values = c('orange','gray', "green"))+
  theme_classic() +
  theme(legend.position = "none")
# used in bin_taxon_genomesize_graph.R

scale <- c(-2.5, -2, -1.5, -1, -0.5, 0)
log_scale <- c("no hits", "0.01%", "0.03%", "0.10%", "0.33%", "1.0%")


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


  

