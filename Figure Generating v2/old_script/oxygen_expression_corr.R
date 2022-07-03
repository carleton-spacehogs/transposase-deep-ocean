setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()

DNA_RNA_tara %>% 
  #filter(upper_size_rna == 3.0) %>%
  filter(Ocean_RNA != "Arctic Ocean") %>%
  ggplot(aes(x = Oxygen_RNA, y = trans_exp_rate)) + 
  facet_wrap(~Layer_RNA) + 
  geom_point(aes(color = Layer_RNA)) +
  geom_smooth(method = "lm", se = F) 

DNA_RNA_tara %>% 
  filter(upper_size_rna == 3.0) %>%
  filter(Ocean_RNA != "Arctic Ocean") %>%
  ggplot(aes(x = Oxygen_RNA, y = trans_exp_rate)) + 
  geom_point(aes(color = Layer_RNA)) +
  geom_smooth(method = "lm", se = F) 

DNA_RNA_tara %>% 
  # filter(upper_size_rna == 3.0) %>%
  filter(Ocean_RNA != "Arctic Ocean") %>%
  ggplot(aes(x = Oxygen_RNA, y = biofilm_exp_rate)) + 
  facet_wrap(~Layer_RNA) + 
  geom_point(aes(color = Ocean_RNA)) +
  geom_smooth(method = "lm", se = F) 

DNA_RNA_tara %>% 
  #filter(upper_size_rna == 3.0) %>%
   # filter(Ocean_RNA != "Arctic Ocean") %>%
  ggplot(aes(x = Oxygen_RNA, y = toxin_exp_rate)) + 
  facet_wrap(~Layer_RNA) + 
  geom_point(aes(color = Layer_RNA)) +
  geom_smooth(method = "lm", se = F)


