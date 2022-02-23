setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()

scale <- c(-2, -2.5, -3, -3.5, -4, -4.5)
# log_scale <- round(10^scale, digits = 5)
percent_scale <- c("1.00%", "0.316%", "0.100%", "0.032%", "0.010%", "0.003%")
colors <- c("sky blue", "steelblue","gray","blue")
color_breaks <- c('SRF','DCM','MIX','MES')

# a fake upper size to conform with tara Ocean data
mala_cov$upper_size_dna <- ifelse(mala_cov$size == "0.2", "1.6", "3")
sel_col <- c("Layer_DNA","log_dna_trans","upper_size_dna")
to_graph <- rbind(mala_cov[,sel_col],
                  DNA_tara[,sel_col]%>%filter(Layer_DNA != "MIX"))
# repaired the fake upper size fraction
to_graph$size_fraction <- ifelse(to_graph$upper_size_dna == "1.6", "planktonic", "particle-associated")
to_graph$Layer_DNA <- factor(to_graph$Layer_DNA, levels = c("SRF","DCM","MES","Malaspina"))

counts = to_graph %>% group_by(size_fraction, Layer_DNA) %>% tally

to_graph %>% 
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=Layer_DNA, y=log_dna_trans, fill=size_fraction)) + 
  geom_boxplot(varwidth = TRUE, outlier.alpha = 0.5) +
  theme_classic() +
  ylab("% reads mapped to transposase ORFs") +
  xlab("Depth") +
  scale_y_continuous(breaks = scale, labels = percent_scale, limits=c(-4.5, -2)) +
  geom_text(data=counts, aes(label=n, y=-2.27), position=position_dodge(0.6)) +
  scale_x_discrete(labels=c("SRF" = "SRF\n(depth < 10 m)", 
                            "DCM" = "DCM\n(17 - 120 m)",
                            "MES" = "MES\n(250 - 1000 m)",
                            "Malaspina" = "Malaspina\n(2400 - 4000 m)")) +
  guides(fill = guide_legend(title = "Metagenome size fraction:")) + # reverse = TRUE
  theme(legend.position = c(0.73, 0.18),
        # axis.title.x=element_blank(),
        legend.background = element_rect(fill="transparent",color="transparent"))

ggsave("depth_size_fraction_trans_coor2.png", plot = last_plot())



for (ocean in unique(DNA_tara$Ocean_DNA)){
  print(ocean)
  MES_rows <- DNA_tara %>% filter(Ocean_DNA %in% ocean) %>% filter(is_MES == "MES")
  if (nrow(MES_rows) < 5){
    print("less than 5 MES, skip")
    next
  }
  t_test <- t.test(log_dna_trans~is_MES, data = DNA_tara %>% 
                     filter(Ocean_DNA %in% ocean) %>% 
                     filter(Layer_DNA %in% c("MES", "SRF", "DCM")))
  print(t_test$p.value)
}

# supplement 1
trans_scale <- c(-2, -2.5, -3, -3.5, -4, -4.5)
trans_percent_scale <- c("1.00%", "0.316%", "0.100%", "0.032%", "0.010%", "0.003%")
DNA_tara %>% 
  filter(!is.na(is_MES)) %>%
  filter(!Ocean_DNA %in% c("Southern", "Red Sea", "Mediterranean")) %>%
  ggplot(aes(x=fct_rev(Ocean_DNA), y=log_dna_trans, fill=is_MES)) +
  theme_classic() +
  geom_boxplot(outlier.color = "gray") + 
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0)) +
  xlab("") + 
  labs(fill = "Depth") +
  scale_y_continuous(breaks = trans_scale, 
                     labels = trans_percent_scale, 
                     limits=c(-4.5, -2)) +
  ylab("% DNA reads mapped to transposase") +
  scale_fill_manual(breaks=c("SRF, DCM", "MES"), values=c('azure','blue')) +
  coord_flip()

oceans <- DNA_tara %>% filter(!Ocean_DNA %in% c("Southern", "Red Sea", "Mediterranean")) %>% 
  select("Ocean_DNA") %>% unique() %>% unlist()



DNA_tara %>% 
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=Layer_DNA, y=log_dna_biofilm))  +
  geom_boxplot() +
  stat_compare_means(comparisons = depth_comparison, method = "t.test") + 
  stat_summary(fun.data = boxplot.give.n, geom = "text", 
               position=position_nudge(x = 0, y = 0.02)) +
  coord_flip()



DNA_tara %>% 
  filter(Layer_DNA != "MIX") %>%
  filter(!Ocean_DNA %in% c("Southern", "Red Sea", "Mediterranean")) %>% # they have too little data points/no MES
  ggplot(aes(x=fct_rev(Ocean_DNA), y=DNA_Transposase, fill=fct_rev(Layer_DNA))) +
  geom_boxplot(outlier.color = "gray") + 
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0)) +
  xlab("") + labs(fill = "Depth") +
  scale_y_continuous(trans='log10')+ ylab("Proportion of transposase reads") +
  scale_fill_manual(breaks=c('SRF','DCM','MES'), values=c('azure','sky blue','blue')) +
  coord_flip() + theme(legend.position = "none")



my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(aes(color=data$Layer_DNA)) + geom_smooth(method=lm, se = F, ...) +
    stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label..)), 
               parse = TRUE)
  p
}

DNA_tara %>%
  filter(!DNA_TA %in% boxplot(DNA_tara$DNA_TA)$out) %>%
  select(c("log_dna_biofilm", "log_dna_trans", "DNA_TA")) %>%
  ggpairs(ggplot2::aes(), lower = list(continuous = my_fn))

DNA_tara %>% 
  select(c("log_dna_biofilm", "log_dna_trans", "DNA_Defense", "Layer_DNA")) %>%
  ggpairs(ggplot2::aes(), lower = list(continuous = my_fn))

hist(DNA_tara$DNA_Defense)

summary(lm(log_dna_trans~Ocean_short, DNA_tara))


