setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
options(scipen=10000)
DNA_RNA_tara <- DNA_RNA_tara %>% filter(Layer_DNA != "MIX")

pTrans <- DNA_RNA_tara %>%
  ggplot(aes(x=Layer_DNA, y=trans_exp_rate, fill = mean_dna_trans)) +
  geom_boxplot() + 
  labs(y = "ORF Expression Ratio") + 
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank()) +
  coord_cartesian(ylim = c(0, 4)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = -0.1)) +
  scale_fill_gradient(limits = c(0, max(DNA_RNA_tara$mean_dna_trans)), 
                      low="white",
                      high="dodgerblue3", 
                      name = "Mean DNA\nabundance of\ntransposase\n",
                      breaks = c(0,
                                 # max(DNA_RNA_tara$mean_dna_trans)/2,
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "SRF", "mean_dna_trans"][[1,1]], 
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "DCM", "mean_dna_trans"][[1,1]],
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "MES", "mean_dna_trans"][[1,1]]), 
                      labels = c("0%", 
                                 # paste(round(max(DNA_RNA_tara$mean_dna_trans)/2*100, 3), "%", sep = ""),
                                 paste("SRF ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "SRF", "mean_dna_trans"][[1,1]]*100, 3), "%", sep = ""),
                                 paste("DCM ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "DCM", "mean_dna_trans"][[1,1]]*100, 3), "%", sep = ""), 
                                 paste("MES ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "MES", "mean_dna_trans"][[1,1]]*100, 3), "%", sep = ""))) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8, 
                                ticks.linewidth = 3, ticks.colour = "black")) 


pBiofilm <- DNA_RNA_tara %>% 
  ggplot(aes(x=Layer_DNA, y=biofilm_exp_rate, fill = mean_dna_biofilm)) +
  geom_boxplot() + 
  labs(y = "ORF Expression Ratio\n(% RNA reads / % DNA reads)") + 
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank()) +
  coord_cartesian(ylim = c(0, 4)) +
  # stat_summary(fun.data = boxplot.give.n, geom = "text") +
  scale_fill_gradient(limits = c(0, max(DNA_RNA_tara$mean_dna_biofilm)), 
                      low="white", high="chartreuse4",
                      # name = "mean % of\nBiofilm DNA\nreads", 
                      name = "Mean DNA\nabundance of\nbiofilm\n",
                      breaks = c(0,
                                 # max(DNA_RNA_tara$mean_dna_biofilm)/2,
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "SRF", "mean_dna_biofilm"][[1,1]], 
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "DCM", "mean_dna_biofilm"][[1,1]],
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "MES", "mean_dna_biofilm"][[1,1]]), 
                      labels = c("0%", 
                                 # paste("SRF ", round(max(DNA_RNA_tara$mean_dna_biofilm)/2*100, 3), "%", sep = ""),
                                 paste("SRF ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "SRF", "mean_dna_biofilm"][[1,1]]*100, 3), "%", sep = ""),
                                 paste("DCM ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "DCM", "mean_dna_biofilm"][[1,1]]*100, 3), "%", sep = ""), 
                                 paste("MES ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "MES", "mean_dna_biofilm"][[1,1]]*100, 3), "%", sep = ""))) + 
  stat_compare_means(comparisons = list(c("DCM","MES"), c("SRF","MES")),
                     label = "p.signif",
                     label.y = c(3.3, 3.7)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8, 
                                ticks.linewidth = 3, ticks.colour = "black")) 

pDefense <- DNA_RNA_tara %>% 
  ggplot(aes(x=Layer_DNA, y=defense_exp_rate, fill = mean_dna_defense)) +
  geom_boxplot() + 
  labs(y = "ORF Expression Ratio") +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank()) +
  stat_compare_means(comparisons = list(c("DCM","MES"), c("SRF","MES")),
                     label.y = c(1.1, 1.25),
                     label = "p.signif") +
  scale_fill_gradient(limits = c(0, max(DNA_RNA_tara$mean_dna_defense)), 
                      low="white", high="red", 
                      # name = "mean % of\nDefense DNA\nreads", 
                      name = "Mean DNA\nabundance of\ndefense\nmechanisms\n",
                      breaks = c(0,
                                 # max(DNA_RNA_tara$mean_dna_defense)/2,
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "SRF", "mean_dna_defense"][[1,1]], 
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "DCM", "mean_dna_defense"][[1,1]],
                                 DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "MES", "mean_dna_defense"][[1,1]]), 
                      labels = c("0%",
                                 # paste(round(max(DNA_RNA_tara$mean_dna_defense)/2*100, 3), "%", sep = ""),
                                 "", 
                                 "",
                                 paste("SRF ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "SRF", "mean_dna_defense"][[1,1]]*100, 3), "%\n", 
                                       "DCM ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "DCM", "mean_dna_defense"][[1,1]]*100, 3), "%\n",
                                       "MES ", round(DNA_RNA_tara[DNA_RNA_tara$Layer_DNA == "MES", "mean_dna_defense"][[1,1]]*100, 3), "%",
                                       sep = ""))) + 
  coord_cartesian(ylim = c(0, 1.3)) + 
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8, 
                                ticks.linewidth = 3, ticks.colour = "black")) 
#stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0))


#im_A <- ggplot() + 
#  background_image(readPNG("blank.pNG")) 
#theme(plot.margin = margin(r=1, unit = "cm", t=1, l=1, b=1))

ggarrange(
  ggarrange(pCassette, pDefense, labels = c("A", "C"), widths = c(0.55, 0.47), ncol =2),
  ggarrange(pBiofilm, pTrans, labels = c("B", "D"), widths = c(0.55, 0.47), ncol =2),
  nrow = 2
)

ggsave("expression_and_cassette_pnps.png", plot = last_plot())

