setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()

# supplement 1
RNA_tara %>% 
  # filter(upper_size_rna == 3.0) %>%
  filter(Ocean_RNA != "Arctic Ocean") %>%
  ggplot(aes(x = Oxygen_RNA, y = log_rna_trans)) + 
  # facet_wrap(~Layer_RNA) + 
  geom_point(aes(color = Layer_RNA)) +
  geom_smooth(method = "lm", se = F) 

# sup 4
DNA_RNA_tara %>%
  #filter(Ocean_RNA != "Arctic Ocean") %>%
  ggplot(aes(x = Oxygen_RNA, y = log_rna_biofilm)) + 
  facet_wrap(~Layer_RNA) + 
  geom_point(aes(color = Ocean_RNA)) +
  geom_smooth(method = "lm", se = F) 


RNA_tara %>% 
  # filter(upper_size_rna == 3.0) %>%
  # filter(Ocean_RNA != "Arctic Ocean") %>%
  ggplot(aes(x = Oxygen_RNA, y = log_rna_trans)) + 
  # facet_wrap(~Layer_RNA) + 
  geom_point(aes(color = Layer_RNA)) +
  geom_smooth(method = "lm", se = F) 



