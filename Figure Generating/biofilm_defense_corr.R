setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()
biofilm_scale2 <- seq(-3.25, -4.0, by = -0.125)
biofilm_percent_scale2 <- c("0.056","0.042","0.032","0.024","0.018", "0.013", "0.010") # 10**(seq(-3.25, -4.0, by = -0.125))
defense_scale <- c(-2.4, -2.2, -2.0, -1.8)
defense_percent_scale <- c("0.40%", "0.63%", "1.00%", "1.58%")

mala_cov$size_fraction <- ifelse(mala_cov$lower_filter_size == "0.2", "planktonic", "particle")
DNA_tara$size_fraction <- ifelse(DNA_tara$upper_size_dna == "1.6", "planktonic", "particle")
bd_sel_col <- c("log_dna_biofilm", "log_dna_trans", "log_dna_defense", "Ocean_DNA", "Layer_DNA", "size_fraction")
bd_to_graph <- rbind(mala_cov[,bd_sel_col], DNA_tara[,bd_sel_col])
bd_to_graph$Layer_DNA <- factor(bd_to_graph$Layer_DNA, levels = c("SRF","DCM","MES","Malaspina"))

bd_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  filter(log_dna_biofilm > -4) %>% # filter outlier
  ggplot(aes(x=log_dna_biofilm, y=log_dna_defense)) +
  #scale_y_continuous(breaks = defense_scale, labels = defense_percent_scale, limits=c(-2.4, -1.8)) + 
  facet_wrap(~Layer_DNA) + # scales = "free"
  scale_x_continuous(breaks = biofilm_scale, labels = biofilm_percent_scale) + # outlier
  geom_point() + # aes(color = size_fraction)
  theme_classic() +
  ylab("% Defense Mechanisms DNA reads") +
  xlab("% DNA reads mapped to biofilm") +
  geom_smooth(method = "lm", se = F)  +
  theme() # legend.position = "none"

defense_biofilm.lm <- lm(log_dna_defense~Layer_DNA + log_dna_biofilm, bt_to_graph)
summary(defense_biofilm.lm)