setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()
# biofilm_scale2 <- seq(-3.25, -4.0, by = -0.25)

biofilm_scale_bd <- c("0.0004","0.0002","0.0001") 
biofilm_scale2 <- log10(as.numeric(biofilm_scale_bd))
biofilm_percent_scale2 <- c("0.04%","0.02%","0.01%")
defense_percent_scale <- c(0.5, 1, 2)
defense_scale <- log10(defense_percent_scale/100)

CAZ_percent_scale <- c(0.8, 0.4, 0.2, 0.1)
CAZ_scale <- log10(CAZ_percent_scale/100)

mala_cov$size_fraction <- ifelse(mala_cov$lower_filter_size == "0.2", "planktonic", "particle")
DNA_tara$size_fraction <- ifelse(DNA_tara$upper_size_dna == "1.6", "planktonic", "particle")
bd_sel_col <- c("log_dna_biofilm", "log_dna_trans", "log_dna_defense",
                "log_dna_sect_CAZ", "percent_sect_CAZ", 
                "percent_sect_pep","DNA_sect_pep",
                "log_dna_sect_pep","avg_percent",
                "Ocean_DNA", "Layer_DNA", "size_fraction")
bd_to_graph <- rbind(mala_cov[,bd_sel_col], DNA_tara[,bd_sel_col])
bd_to_graph$Layer_DNA <- factor(bd_to_graph$Layer_DNA, levels = c("SRF","DCM","MES","BAT"))

bd_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=log_dna_sect_CAZ, y=log_dna_defense)) +
  facet_wrap(~Layer_DNA, ncol = 2) + # scales = "free"
  scale_y_continuous(breaks = defense_scale, labels = defense_percent_scale, limits=c(-2.3, -1.65)) + 
  scale_x_continuous(breaks = CAZ_scale, labels = CAZ_percent_scale) +
  geom_point(aes(color=size_fraction)) +
  theme_classic() +
  ylab("% DNA reads mapped to defense mechanism ORFs") +
  xlab("% DNA reads mapped to secretory CAZyme") +
  geom_smooth(method = "lm", se = F)  +
  theme()

bd_to_graph %>%
  filter(!Layer_DNA %in% c("MIX",NA)) %>%
  # filter(!(Layer_DNA == "BAT" & size_fraction == "particle")) %>%
  ggplot(aes(x=log_dna_sect_CAZ, y=log_dna_trans, color=Layer_DNA)) +
  facet_wrap(~Layer_DNA, ncol = 2, scales = "free") + # 
  # scale_y_continuous(breaks = defense_scale, labels = defense_percent_scale) + # , limits=c(-2.3, -1.65) 
  geom_point() +
  theme_classic() +
  ylab("% DNA reads mapped to transposase") +
  # xlab("secretory CAZenzyme among all CAZenzyme (%)") +
  xlab("% DNA reads mapped to CAZenzyme") +
  geom_smooth(method = "lm", se = F)  +
  theme()


bd_to_graph %>%
  filter(!Layer_DNA %in% c("MIX",NA)) %>%
  ggplot(aes(x=percent_sect_CAZ, y=log_dna_defense)) +
  facet_wrap(~Layer_DNA, ncol = 2) + # scales = "free"
  scale_y_continuous(breaks = defense_scale, labels = defense_percent_scale) + # , limits=c(-2.3, -1.65) 
  geom_point(aes(color=size_fraction)) +
  theme_classic() +
  ylab("% DNA reads mapped to defense mechanism ORFs") +
  # xlab("abundance of secretory CAZyme in a metagenome") +
  xlab("secretory CAZenzyme among all CAZyme (%)") +
  geom_smooth(method = "lm", se = F)  +
  theme()

bd_to_graph %>%
  ggplot(aes(x=DNA_sect_pep, y=percent_sect_pep))+
  geom_point(aes(color=size_fraction, shape = Layer_DNA)) +
  xlab("abundance of secretory peptidase in a metagenome") +
  ylab("abundance of secretory peptidase in all peptidase")

bd_to_graph %>%
  ggplot(aes(x=log_dna_sect_CAZ, y=percent_sect_CAZ))+
  geom_point(aes(color=size_fraction, shape = Layer_DNA)) +
  xlab("abundance of secretory CAZyme in a metagenome") +
  ylab("abundance of secretory CAZyme in all CAZyme")

plot(DNA_sect_pep ~ percent_sect_pep, bd_to_graph)

# used in gen_figure5_integron_finder_category_v2.R
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
