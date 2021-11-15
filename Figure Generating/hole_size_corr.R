setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()

# for making figure 1, the 2 size filter, see biofilm_transposase_corr

DNA_tara %>% 
  # filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=upper_size_dna, y=DNA_Defense)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_jitter(aes(color=Layer_DNA)) +
  labs(color = "Depth")+
  xlab("size fraction (\u03BCm)") +
  scale_x_discrete(labels=c("1.6" = "0.22-1.6", "3" = "0.22-3.0")) +
  stat_summary(fun.data = boxplot.give.n, geom = "text") + 
  # scale_y_continuous(breaks = scale, labels = log_scale) +
  #stat_pvalue_manual(stat.test, label = "p.adj") + 
  stat_compare_means(comparisons = list( c("1.6", "3") ))
  scale_color_manual(breaks=c('SRF','DCM','MES', "MIX"), 
                     values=c("sky blue", "steelblue", "blue", "gray"),
                     labels = c("SRF (5 or 9 m)", "DCM (17-188m)", "MIX (25-200m)", "MES (250-1000m)"))
  
DNA_tara %>% 
  #filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x = Depth, y = DNA_Defense)) +
  theme_classic() +
  geom_point(aes(color = Layer_DNA)) +
  geom_smooth(method = "lm", se = F)
