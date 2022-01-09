setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

# median_bin_pnps is normally distributed
hist(filter_outliers(bin_taxon, "median_bin_pnps")$median_bin_pnps)

all <- filter_outliers(bin_taxon, "percent_trans") %>% 
  select("median_bin_pnps", "percent_trans", "log_percent_trans", "complete genome size (Mbp)",
         "depth", "Class", "percent_biofilm", "log_percent_biofilm")
all <- filter_outliers(all, "median_bin_pnps")
all$graphing_log_trans <- ifelse(all$log_percent_trans < -9, -2, all$log_percent_trans)
all$graphing_log_biofilm <- ifelse(all$log_percent_biofilm < -9, -2, all$log_percent_biofilm)

all$depth <- fct_rev(all$depth)
# all$shape3D <- ifelse(all$percent_trans %in% c(out), "`", "*")
color2 <- ifelse(all$depth=="MES","blue4",ifelse(all$depth=="DCM", "steelblue", "sky blue"))

scatterplot3d(all[,c("median_bin_pnps", "graphing_log_trans", "depth")], 
              angle = 23, pch = 20, # pch = all$shape3D, # type = "h", lty.hplot = 2,
              lab.z = 2,
              lab = c(4, 4, 1),
              y.ticklabs = c("<0.01%", "0.03%", "0.1%", "0.32%", "1%"),
              # scale.y = 0.5,
              z.ticklabs = c("     MES", "DCM", "SRF"), 
              grid=TRUE, box=TRUE, color=color2)

scatterplot3d(all[,c("graphing_log_biofilm", "graphing_log_trans", "depth")], 
              angle = 23, pch = 20, # pch = all$shape3D, # type = "h", lty.hplot = 2,
              lab.z = 2,
              lab = c(4, 4, 1),
              y.ticklabs = c("<0.01%", "0.03%", "0.1%", "0.32%", "1%"),
              x.ticklabs = c("<0.01%", "0.03%", "0.1%", "0.32%", "1%     "),
              z.ticklabs = c("     MES", "DCM", "SRF"), 
              grid=TRUE, box=TRUE, color=color2)


depth.lm <- lm(log_percent_trans~depth, data=all)
summary(depth.lm)

taxon_depth.lm <- lm(log_percent_trans~depth+Class, data=all)
anova(taxon_depth.lm)

biofilm_depth.lm <- lm(log_percent_trans~depth+log_percent_biofilm, all)
anova(biofilm_depth.lm)

gsize_depth.lm <- lm(log_percent_trans~depth+`complete genome size (Mbp)`, all)
anova(gsize_depth.lm)

all2.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm, all)
anova(all2.lm)
summary(all2.lm)

all2.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm, all)
all_interact.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm+depth*Class, all)

anova(all_interact.lm)
get_r(all_interact.lm)

get_r(all2.lm)
get_r(depth.lm)
get_r(taxon_depth.lm)
get_r(biofilm_depth.lm)

summary(lm(log_percent_trans~log_percent_biofilm, data=all%>%filter(depth=="SRF")))
summary(lm(log_percent_trans~log_percent_biofilm, data=all))

pnps_depth.lm <- lm(log_percent_trans~depth+median_bin_pnps, data=all)
anova(pnps_depth.lm) 
summary(lm(log_percent_trans~median_bin_pnps, data=all))
summary(pnps_depth.lm)

all.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm+median_bin_pnps, all)
anova(all.lm)
summary(all.lm)


split_layer_taxon <- function(taxon_name){
  t_name <- taxon_name
  strain <- bin_taxon %>% filter(Class == t_name)
  strain <- strain %>% mutate(color2 = ifelse(depth=="MES","blue4",ifelse(depth=="DCM", "steelblue", "sky blue")))
  strain$depth <- factor(strain$depth, levels = c("MES","DCM","SRF"))
  jpeg(paste(t_name,'.jpg', sep = ""))
  scatterplot3d(strain[,c("percent_biofilm", "percent_trans", "depth")], 
                angle = 35, pch = 20, type = "h", lty.hplot = 2,
                lab.z = 2,
                lab = c(4, 4, 1),
                scale.y = 1,
                z.ticklabs = c("     MES", "DCM", "SRF"),
                grid=TRUE, box=FALSE, color=strain$color2)
  dev.off()
}

for (t_name in c("Acidimicrobidae", "Flavobacteria", "Gammaproteobacteria",
                 "Betaproteobacteria", "Deltaproteobacteria", "OM190")){
  split_layer_taxon(t_name)
}


############ trash

Yrandom <- 10*rnorm(1000)
Xrandom <- 10*rnorm(1000)
anova(lm(Yrandom~Xrandom))

library("gg3D")
# install.packages("scatterplot3d") # Install
 # load

scatterplot3d(to_graph1[-1,c(3,13,2)], angle = 60,lty.hplot = 2, pch = 16,
              type="h", grid=TRUE, box=FALSE, color = to_graph1[-1,]$color2)


depth.lm <- lm(percent_trans~depth, data=strain)
summary(depth.lm)

biofilm_depth.lm <- lm(percent_trans~depth+percent_biofilm, data=strain)
anova(biofilm_depth.lm)

pnps_depth.lm <- lm(percent_trans~depth+median_bin_pnps, data=strain)
anova(pnps_depth.lm)

taxon_depth.lm <- lm(percent_trans~depth+Class, data=strain)
anova(taxon_depth.lm)

all.lm <- lm(percent_trans~depth+Class+percent_biofilm+median_bin_pnps, data=strain)
anova(all.lm)
summary(all.lm)

strain <- bin_taxon 
out <- boxplot(percent_trans~depth, data=strain)$out
# out <- boxplot(strain$percent_trans)$out
out_ind <- which(strain$percent_trans %in% c(out))
strain$color2 <- ifelse(strain$depth=="MES","blue4",ifelse(strain$depth=="DCM", "steelblue", "sky blue"))
strain$shape3D <- ifelse(strain$percent_trans %in% c(out), "`", "*")
strain$depth <- factor(strain$depth, levels = c("MES","DCM","SRF"))
scatterplot3d(strain[,c("percent_biofilm", "percent_trans", "depth")], 
              angle = 35, pch = strain$shape3D, # type = "h", lty.hplot = 2,
              lab.z = 2,
              lab = c(4, 4, 1),
              scale.y = 1,
              z.ticklabs = c("     MES", "DCM", "SRF"), 
              grid=TRUE, box=FALSE, color=strain$color2)
legend("right", legend = c("outliers", "non-outliers"), pch = c("`", "*"))

pnps_upper <- min(boxplot(bin_taxon$median_bin_pnps)$out)

trans_upper <- min(boxplot(bin_taxon$percent_trans)$out)
# out <- boxplot(percent_trans~depth, data=strain)$out
# out_ind <- which(strain$percent_trans %in% c(out))

strain <- bin_taxon %>% 
  filter(log_percent_trans > -10) %>%
  filter(log_percent_biofilm > -10) %>%
  filter(log_median_bin_pnps > -10)
# filter(median_bin_pnps < pnps_upper) %>% 
# filter(strain$percent_trans < trans_upper)
strain$color2 <- ifelse(strain$depth=="MES","blue4",ifelse(strain$depth=="DCM", "steelblue", "sky blue"))
# strain$shape3D <- ifelse(strain$percent_trans %in% c(out), "`", "*")
strain$depth <- factor(strain$depth, levels = c("MES","DCM","SRF"))
# scatterplot3d(strain[,c("median_bin_pnps", "percent_trans", "depth")], 
# scatterplot3d(strain[,c("log_percent_biofilm", "log_percent_trans", "depth")],
scatterplot3d(strain[,c("log_median_bin_pnps", "log_percent_trans", "depth")],
              angle = 45, pch = 20, # pch = strain$shape3D, # type = "h", lty.hplot = 2,
              lab.z = 2,
              lab = c(4, 4, 1),
              scale.y = 0.5,
              z.ticklabs = c("     MES", "DCM", "SRF"), 
              grid=TRUE, box=TRUE, color=strain$color2)
# legend("right", legend = c("outliers", "non-outliers"), pch = c("`", "*"))

depth.lm <- lm(log_percent_trans~depth, data=strain)
summary(depth.lm)

biofilm_depth.lm <- lm(log_percent_trans~depth+log_percent_biofilm, data=strain)
anova(biofilm_depth.lm)

pnps_depth.lm <- lm(log_percent_trans~depth+log_median_bin_pnps, data=strain)
anova(pnps_depth.lm)

taxon_depth.lm <- lm(log_percent_trans~depth+Class, data=strain)
anova(taxon_depth.lm)

all.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm+median_bin_pnps, data=strain)

# dont_show <- c("SAR202-2", "novelClass_E", "Phycisphaerae", "Opitutae", "OM190", 
#               "Gemmatimonadetes", "Planctomycetia","Others Or Unknown", "Flavo.")
# low_trans <- c("Acidi.", "Verruco.")

# "Acidimicrobidae", "Flavobacteria", "Gammaproteobacteria"
# "sky blue", "steelblue", "blue4"