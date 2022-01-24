setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
SRF_s <- read_csv("SRF_0.22-1.6_summary.csv")
DCM_s <- read_csv("DCM_0.22-1.6_summary.csv")
MES_s <- read_csv("MES_0.22-1.6_summary.csv")
SRF_b <- read_csv("SRF_0.22-3_summary.csv")
DCM_b <- read_csv("DCM_0.22-3_summary.csv")
MES_b <- read_csv("MES_0.22-3_summary.csv")

all <- rbind(SRF_s, DCM_s, MES_s, SRF_b, DCM_b, MES_b)

all %>%
  filter(COG_function %in% c("total", ""))
  ggplot(aes(x=`pnps-median`, fill = depth)) +
  geom_boxplot(outlier.colour = NA, notch = TRUE) +
  scale_x_continuous(breaks = scale1, labels = log_scale1, 
                     name = "pN/pS of metagenomic ORFs", limits = c(-1.5, 1)) +
  theme_classic() +
  facet_wrap(~depth, nrow = 4, strip.position = "left") +
  guides(fill=guide_legend(title="ORF Type:",reverse = TRUE)) +
  scale_fill_manual(labels = c("Transposase","Non-Transposase"),
                    values=c("orange","gray")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = c(0.85, 0.47),
        plot.margin=unit(c(0.2,0.2,0.2,0.5), "cm"),
        legend.background = element_rect(fill = "transparent", 
                                         color = "transparent"))