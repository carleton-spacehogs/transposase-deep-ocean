setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(pn_ps_metagenome) %=% init_individual_metagenomes()

to_graph <- pn_ps_metagenome %>%
  mutate(is_trans = ifelse(gene_type == "transposase", 'y', 'n'))

# run after running individual_metagenome_pnps.R
counts1 <- to_graph %>% count(depth, size_fraction)

scale <- c(-2, -1.5, -1, -0.5, 0, 0.5)
log_scale <- c("< 0.01", "0.032", "0.1", "0.316", "1", "3.162")

to_graph %>%
  ggplot(aes(x=fct_rev(depth), y=log_pnps, fill = size_fraction)) +
  geom_boxplot(outlier.shape=NA) +
  geom_text(data=counts1, aes(label=n, y=0.85), position=position_dodge(0.85)) +
  scale_fill_manual(labels = c("particle-associated","planktonic"),
                    values=c("orange","green"))+
  ylab(expression(paste(italic("pN/pS"), " ratio"))) +
  scale_y_continuous(breaks = scale, limits = c(-2, 1),
                     labels = log_scale) + #, limits = c(-2, 0.8)
  theme_classic() +
  xlab("Depth") + 
  coord_flip() +
  guides(fill = guide_legend(title = "Size fraction", reverse = TRUE)) +
  theme(# legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill="transparent",color="transparent"))

ggsave("S5_metagenome-pnps-size-fraction.png", plot = last_plot())
ggsave("S5_metagenome-pnps-size-fraction.svg", plot = last_plot())


bat_particle <- pn_ps_metagenome %>% filter(depth == "BAT") %>% filter(size_fraction == "b")
median(bat_particle$pnps, na.rm = TRUE)

bat_plankton <- pn_ps_metagenome %>% filter(depth == "BAT") %>% filter(size_fraction == "s")
median(bat_plankton$pnps, na.rm = TRUE)
