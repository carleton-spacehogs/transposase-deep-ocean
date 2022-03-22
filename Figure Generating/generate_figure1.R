setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()

scale <- c(-2, -2.5, -3, -3.5, -4, -4.5)
# log_scale <- round(10^scale, digits = 5)
percent_scale <- c("1.00%", "0.316%", "0.100%", "0.032%", "0.010%", "0.003%")

mala_cov$size_fraction <- ifelse(mala_cov$lower_filter_size == "0.2", "planktonic", "particle-\nassociated")
DNA_tara$size_fraction <- ifelse(DNA_tara$upper_size_dna == "1.6", "planktonic", "particle-\nassociated")

sel_col <- c("Layer_DNA","log_dna_trans","size_fraction","Depth")
to_graph <- rbind(mala_cov[,sel_col], DNA_tara[,sel_col]%>%filter(Layer_DNA != "MIX"))
to_graph$Layer_DNA <- factor(to_graph$Layer_DNA, levels = c("SRF","DCM","MES","BAT"))

counts = to_graph %>% group_by(size_fraction, Layer_DNA) %>% tally
# generate figure 1
to_graph %>% 
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=fct_rev(Layer_DNA), y=log_dna_trans, fill=fct_rev(size_fraction))) + 
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -4.5, ymax = -2), fill = "#57C1FF", alpha = 0.1) +
  geom_rect(aes(xmin = 1.5, xmax = 4.5, ymin = -4.5, ymax = -2), fill = "#A9E6FF", alpha = 0.1) +
  geom_boxplot(varwidth = TRUE, outlier.alpha = 0.5) +
  theme_classic() +
  ylab("% reads mapped to transposase ORFs") +
  xlab("Depth") +
  scale_y_continuous(breaks = scale, labels = percent_scale, limits=c(-4.5, -2)) +
  geom_text(data=counts, aes(label=n, y=-2.33), position=position_dodge(0.6)) +
  scale_x_discrete(labels=c("SRF" = "SRF\n(depth < 10 m)", 
                            "DCM" = "DCM\n(17 - 120 m)",
                            "MES" = "MES\n(250 - 1000 m)",
                            "BAT" = "BAT\n(2400 - 4000 m)")) + 
  coord_flip() +
  scale_fill_manual(values=c("green","orange"))+
  guides(fill = guide_legend(title = "Size fraction", reverse = TRUE)) +
  annotate(geom="text", x=2.25, y=-4.1, 
           label="italic(Tara)~Oceans", parse=TRUE, color="blue") +
  annotate(geom="text", x=2, y=-4.1, 
           label="metagenomes", color="blue") +
  annotate(geom="text", x=1.15, y=-4.1, 
           label="Malaspina", color="dark blue") +
  annotate(geom="text", x=0.9, y=-4.1, 
           label="metagenomes", color="dark blue") +
  theme(legend.background = element_rect(fill="transparent",color="transparent"))

ggsave("F1_depth_size_fraction_trans_coor.png", plot = last_plot())

DNA_trans_depth.lm <- lm(log_dna_trans~Layer_DNA, to_graph)
DNA_trans_with_size.lm <- lm(log_dna_trans~Layer_DNA + size_fraction, to_graph)
anova(DNA_trans_depth.lm, DNA_trans_with_size.lm)
summary(DNA_trans_with_size.lm)
anova(DNA_trans_with_size.lm)


sel_col2 <- c("log_dna_trans","is_MES", "Ocean_DNA")
supplement_graph <- rbind(mala_cov[,sel_col2], DNA_tara[,sel_col2]%>%filter(is_MES != "MIX"))

for (ocean in unique(supplement_graph$Ocean_DNA)){
  print(ocean)
  MES_rows <- DNA_tara %>% filter(Ocean_DNA %in% ocean) %>% filter(is_MES == "MES, Malaspina")
  if (nrow(MES_rows) < 5){
    print("less than 5 MES, skip")
    next
  }
  t_test <- t.test(log_dna_trans~is_MES, data = supplement_graph %>% 
                     filter(Ocean_DNA %in% ocean) %>% 
                     filter(is_MES != "MIX"))
  print(t_test$p.value)
}


# supplement 1
supplement_graph <- supplement_graph %>% # these oceans don't have enough samples for boxplots
  filter(!Ocean_DNA %in% c("Southern", "Red Sea", "Mediterranean")) 
counts = supplement_graph %>% group_by(Ocean_DNA, is_MES) %>% tally

s11 <- supplement_graph%>%
  ggplot(aes(x=fct_rev(Ocean_DNA), y=log_dna_trans, fill=is_MES)) +
  theme_classic() +
  geom_boxplot(outlier.color = "gray") + 
  # stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0)) +
  geom_text(data=counts, aes(label=n, y=-2.3), position=position_dodge(0.8)) +
  xlab("") + 
  labs(fill = "Depth") +
  scale_y_continuous(breaks = scale, labels = percent_scale) +
  ylab("% DNA reads mapped to transposase") +
  scale_fill_manual(breaks=c("SRF, DCM", "MES, BAT"), values=c('azure','blue')) +
  coord_flip()

# supplementary 2 
colors <- c("light blue","sky blue", "steelblue","blue")
color_breaks <- c('SRF','DCM','MES','BAT')
depth_scale <- c(5, 10, 40, 100, 400, 1000, 4000)
d_scale <- log10(depth_scale)
to_graph$log_depth <- log10(to_graph$Depth)

s12 <- to_graph %>% 
  ggplot(aes(x = log_depth, y = log_dna_trans)) +
  theme_classic() +
  scale_y_continuous(breaks = scale, labels = percent_scale, limits=c(-4.5, -2)) + 
  scale_x_continuous(breaks = d_scale, labels = depth_scale) + 
  labs(x="Depth (m)", y="% DNA reads mapped to transposase", colour="Depth") +
  geom_jitter(aes(color = Layer_DNA), width = 0.02) +
  geom_smooth(method = "lm", se = F, color = "orange") +
  scale_color_manual(breaks=color_breaks, 
                     values=colors) 
  # theme(legend.position = "none")
ggarrange(s12, s11, labels = c("A", "B"), 
          ncol = 2, nrow = 1, widths = c(0.48, 0.55))

ggsave("S1_transposase-depth_each-ocean.png", plot = last_plot())





# Not relavant to published graph from this point
# exploratory graph for biofilm
DNA_tara %>% 
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=Layer_DNA, y=log_dna_biofilm))  +
  geom_boxplot() +
  stat_compare_means(comparisons = depth_comparison, method = "t.test") + 
  stat_summary(fun.data = boxplot.give.n, geom = "text", 
               position=position_nudge(x = 0, y = 0.02)) +
  coord_flip()


sel_col3 <- c("log_dna_trans","Layer_DNA", "Ocean_DNA")
individual_depth <- rbind(mala_cov[,sel_col3], DNA_tara[,sel_col3]%>%filter(Layer_DNA != "MIX"))
individual_depth$Layer_DNA <- factor(individual_depth$Layer_DNA, levels = c("SRF", "DCM", "MES", "Malaspina"))
individual_depth %>% 
  filter(!Ocean_DNA %in% c("Southern", "Red Sea", "Mediterranean")) %>%
  ggplot(aes(x=fct_rev(Ocean_DNA), y=log_dna_trans, fill=Layer_DNA)) +
  theme_classic() +
  geom_boxplot(outlier.color = "gray") + 
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0)) +
  xlab("") + 
  labs(fill = "Depth") +
  scale_y_continuous(breaks = scale, 
                     labels = percent_scale) +
  ylab("% DNA reads mapped to transposase") +
  # scale_fill_manual(breaks=c("SRF, DCM", "MES, Malaspina"), values=c('azure','blue')) +
  coord_flip()


