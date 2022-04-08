setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()
# biofilm_scale2 <- seq(-3.25, -4.0, by = -0.25)

biofilm_scale_bd <- c("0.0004","0.0002","0.0001") 
biofilm_scale2 <- log10(as.numeric(biofilm_scale_bd))
biofilm_percent_scale2 <- c("0.04%","0.02%","0.01%") 
defense_scale <- c(-2.2, -2.0, -1.8)
defense_percent_scale <- c("0.63%", "1.00%", "1.58%")

mala_cov$size_fraction <- ifelse(mala_cov$lower_filter_size == "0.2", "planktonic", "particle")
DNA_tara$size_fraction <- ifelse(DNA_tara$upper_size_dna == "1.6", "planktonic", "particle")
bd_sel_col <- c("log_dna_biofilm", "log_dna_trans", "log_dna_defense", 
                "log_dna_toxin", "Ocean_DNA", "Layer_DNA", "size_fraction")
bd_to_graph <- rbind(mala_cov[,bd_sel_col], DNA_tara[,bd_sel_col])
bd_to_graph$Layer_DNA <- factor(bd_to_graph$Layer_DNA, levels = c("SRF","DCM","MES","BAT"))

p_bd <- bd_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=log_dna_biofilm, y=log_dna_defense)) +
  facet_wrap(~Layer_DNA, ncol = 1) + # scales = "free"
  scale_y_continuous(breaks = defense_scale, labels = defense_percent_scale, limits=c(-2.3, -1.75)) + 
  scale_x_continuous(breaks = biofilm_scale2, labels = biofilm_percent_scale2, limits=c(-4, -3.25)) + # outlier
  geom_point() + # aes(color = size_fraction)
  theme_classic() +
  ylab("% DNA reads mapped to defense mechanism ORFs") +
  xlab("% DNA reads mapped to\nbiofilm ORFs") +
  geom_smooth(method = "lm", se = F)  +
  theme() # legend.position = "none"

defense_biofilm.lm <- lm(log_dna_defense~Layer_DNA + log_dna_biofilm, bd_to_graph)
summary(lm(log_dna_defense~Layer_DNA, bd_to_graph))
summary(defense_biofilm.lm)


bd_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=log_dna_biofilm, y=log_dna_toxin)) +
  facet_wrap(~Layer_DNA, ncol = 1) + #,  scales = "free"
  geom_point() + #aes(color = size_fraction)
  theme_classic() +
  ylab("% Defense Mechanisms DNA reads") +
  xlab("% DNA reads mapped to biofilm") +
  geom_smooth(method = "lm", se = F)  +
  theme() # legend.position = "none"

toxin_biofilm.lm <- lm(log_dna_toxin~Layer_DNA + log_dna_biofilm, bd_to_graph)
summary(toxin_biofilm.lm)
summary(lm(log_dna_toxin~Layer_DNA, bd_to_graph))
anova(toxin_biofilm.lm)
