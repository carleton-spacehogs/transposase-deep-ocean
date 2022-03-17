setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()
options(scipen=10000)
DNA_RNA_tara <- DNA_RNA_tara %>% filter(Layer_DNA %in% c("SRF", "DCM", "MES"))

get_depth_median <- function(depth, gene) {
  sel_depth <- c("Layer_DNA", paste("DNA_", gene, sep = ""))
  all_depth <- rbind(mala_cov[,sel_depth], DNA_tara[,sel_depth])
  gene_depth <- all_depth[all_depth[,"Layer_DNA"] == depth,2]
  return(median(gene_depth[[1]]))
}

get_depth_median2 <- function(depth, gene) {
  gene_depth <- DNA_RNA_tara[DNA_RNA_tara[,"Layer_DNA"] == depth, paste("DNA_", gene, sep = "")]
  return(median(gene_depth[[1]]))
}

get_depth_percent <- function(depth, num){paste(depth, " ", round(num*100, 3), "%", sep = "")}

SRF_trans_median <- get_depth_median("SRF", "Transposase")
DCM_trans_median <- get_depth_median("DCM", "Transposase")
MES_trans_median <- get_depth_median("MES", "Transposase")
BAT_trans_median <- get_depth_median("BAT", "Transposase")
trans_array <- c(SRF_trans_median, DCM_trans_median, MES_trans_median, BAT_trans_median)

pTrans <- DNA_RNA_tara %>%
  ggplot(aes(x=Layer_DNA, y=trans_exp_rate)) +
  geom_boxplot(aes(fill = median_dna_trans)) + 
  labs(y = "ORF Expression Ratio") + 
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank()) +
  coord_cartesian(ylim = c(0, 4)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = -0.1)) +
  scale_fill_gradient(limits = c(0, max(trans_array)), 
                      low="white",
                      high="dodgerblue3", 
                      name = "Median DNA\nabundance of\ntransposase\n",
                      breaks = c(0, trans_array), 
                      labels = c("0%", 
                                 paste("SRF  ", round(SRF_trans_median*100, 3), "%", sep = ""),
                                 paste("DCM ", round(DCM_trans_median*100, 3), "%\n", sep = ""),
                                 get_depth_percent("MES", MES_trans_median),
                                 get_depth_percent("BAT", BAT_trans_median))) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8, 
                                ticks.linewidth = 3, ticks.colour = "black")) 


SRF_bio_median <- get_depth_median2("SRF", "Biofilm")
DCM_bio_median <- get_depth_median2("DCM", "Biofilm")
MES_bio_median <- get_depth_median2("MES", "Biofilm")
BAT_bio_median <- get_depth_median("BAT", "Biofilm")
bio_array <- c(SRF_bio_median, DCM_bio_median, MES_bio_median, BAT_bio_median)

pBiofilm <- DNA_RNA_tara %>% 
  ggplot(aes(x=Layer_DNA, y=biofilm_exp_rate, fill = median_dna_biofilm)) +
  geom_boxplot() + 
  labs(y = "ORF Expression Ratio\n(% RNA reads / % DNA reads)") + 
  theme_classic() + 
  theme(axis.title.x=element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(ylim = c(0, 4)) +
  # stat_summary(fun.data = boxplot.give.n, geom = "text") +
  scale_fill_gradient(limits = c(0, max(bio_array)), # bio_array is based on DNA only, DNA_RNA_tara has different data
                      low="white", high="chartreuse4",
                      name = "Median DNA\nabundance of\nbiofilm\n",
                      breaks = c(0, bio_array), 
                      labels = c("0%","","",get_depth_percent("MES ", MES_bio_median),
                                 paste("SRF  ", round(SRF_bio_median*100, 3), "%\n", 
                                       "DCM ", round(DCM_bio_median*100, 3), "%\n",
                                       "BAT   ", round(BAT_bio_median*100, 3), "%",
                                       sep = ""))) + 
  stat_compare_means(comparisons = list(c("DCM","MES"), c("SRF","MES")),
                     label = "p.signif",
                     label.y = c(3.3, 3.7)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8, 
                                ticks.linewidth = 3, ticks.colour = "black")) 

SRF_def_median <- get_depth_median("SRF", "Defense")
DCM_def_median <- get_depth_median("DCM", "Defense")
MES_def_median <- get_depth_median("MES", "Defense")
BAT_def_median <- get_depth_median("BAT", "Defense")

pDefense <- DNA_RNA_tara %>% 
  ggplot(aes(x=Layer_DNA, y=defense_exp_rate, fill = median_dna_defense)) +
  geom_boxplot() + 
  labs(y = "ORF Expression Ratio") +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank()) +
  stat_compare_means(comparisons = list(c("DCM","MES"), c("SRF","MES")),
                     label.y = c(1.1, 1.25),
                     label = "p.signif") +
  scale_fill_gradient(limits = c(0, max(SRF_def_median, DCM_def_median, DCM_def_median, BAT_def_median)), 
                      low="white", high="red", 
                      name = "Median DNA\nabundance of\ndefense\nmechanisms\n",
                      breaks = c(0,SRF_def_median, DCM_def_median, DCM_def_median, BAT_def_median),
                      labels = c("0%",
                                 "", 
                                 "",
                                 paste("SRF  ", round(SRF_def_median*100, 3), "%\n", 
                                       "DCM ", round(DCM_def_median*100, 3), "%\n",
                                       "MES  ", round(MES_def_median*100, 3), "%",
                                       sep = ""),
                                 get_depth_percent("BAT", BAT_def_median))
                      ) + 
  coord_cartesian(ylim = c(0, 1.3)) + 
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8, 
                                ticks.linewidth = 3, ticks.colour = "black")) 
#stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0))


#im_A <- ggplot() + 
#  background_image(readPNG("blank.pNG")) 
#theme(plot.margin = margin(r=1, unit = "cm", t=1, l=1, b=1))

# from metagenome_pnps_v3.R -> pCassette
ggarrange(
  ggarrange(pCassette, pDefense, labels = c("A", "C"), widths = c(0.5, 0.47), ncol =2),
  ggarrange(pBiofilm, pTrans, labels = c("B", "D"), widths = c(0.5, 0.47), ncol =2),
  nrow = 2
)

ggsave("F6_expression_and_cassette_pnps.png", plot = last_plot())

