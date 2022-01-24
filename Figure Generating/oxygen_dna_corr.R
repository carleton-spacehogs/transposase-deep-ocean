setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()


trans_scale <- c(-2, -2.5, -3, -3.5, -4, -4.5)
trans_percent_scale <- c("1.00%", "0.316%", "0.100%", "0.032%", "0.010%", "0.003%")
p1 <- DNA_tara %>% 
  ggplot(aes(x = Oxygen_DNA, y = log_dna_trans)) + 
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) +
  ylab("% transposase reads") +
  xlab("Oxygen Concentration (\U00B5M)")+
  facet_wrap(~Layer_DNA) + 
  theme_classic() +
  geom_point(aes(color = Layer_DNA)) +
  geom_smooth(method = "lm", se = F) +
  theme(legend.position = "none")

p2 <- DNA_tara %>% 
  ggplot(aes(x = Oxygen_DNA, y = log_dna_trans)) + 
  geom_point(aes()) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) +
  ylab("% reads mapped to transposase ORFs") +
  xlab("Oxygen Concentration (\U00B5M)")+
  labs(color="Dpeth")+
  theme_classic() +
  geom_point(aes(color = Layer_DNA)) +
  geom_smooth(method = "lm", se = F) 

# supplementary 2
ggarrange(p2, p1, labels = c("   A", "B"), ncol = 2, nrow = 1, widths = c(0.5, 0.5))

with_oxygen <- lm(log_dna_trans~Layer_DNA+Oxygen_DNA, DNA_tara)
anova(with_oxygen)
