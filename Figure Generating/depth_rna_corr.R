setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()

RNA_tara %>% 
  filter(Layer_RNA != "MIX") %>%
  filter(!Ocean_RNA %in% c("Southern Ocean", "Red Sea", "Mediterranean Sea")) %>% # they have too little data points/no MES
  ggplot(aes(x=Ocean_RNA, y=log_rna_trans, fill=Layer_RNA)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  geom_boxplot() + coord_flip()

RNA_tara %>% 
  filter(Layer_RNA != "MIX") %>%
  ggplot(aes(x=Layer_RNA, y=log_rna_biofilm)) +
  geom_violin() +
  geom_jitter(aes(color=Ocean_RNA)) + xlab("") +
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  stat_compare_means(comparisons = depth_comparison)

  

RNA_tara %>% 
  # filter(Layer_RNA != "MIX") %>%
  filter(!Ocean_RNA %in% c("Southern Ocean", "Red Sea", "Mediterranean Sea")) %>% # they have too little data points/no MES
  ggplot(aes(x=Ocean_RNA, y=log_rna_toxin, fill=Layer_RNA)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  geom_boxplot() + coord_flip()
