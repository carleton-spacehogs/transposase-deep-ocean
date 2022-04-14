setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()
g(pn_ps_metagenome) %=% init_individual_metagenomes()

scale1 <- c(-2, -1.5, -1, -0.5, 0, 0.60206)
log_scale1 <- c("<0.01", "<0.032", "0.10", "0.32", "1", ">4")

to_graph <- pn_ps_metagenome %>%
  filter(log_pnps < 1.4) %>%
  filter(log_pnps > -3) %>%
  mutate(is_trans = ifelse(gene_type == "transposase", 'y', 'n'))

counts <- to_graph %>% count(depth, is_trans)
# counts <- mutate(counts, n = signif(n, digits = 4))

p1<- to_graph %>%
  ggplot(aes(x=fct_rev(depth), y=log_pnps, fill = is_trans)) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_continuous(breaks = scale1, labels = log_scale1, 
                     name = expression(paste(italic("pN/pS"), " of metagenomic ORFs")), 
                     limits = c(-1.5, 1.4)) +
  theme_classic() +
  # facet_wrap(~depth, nrow = 4, strip.position = "left") +
  guides(fill=guide_legend(title="ORF Type:", reverse = TRUE)) +
  scale_fill_manual(labels = c("Non-transposase", "Transposase"),
                     values= c("gray","orange")) +
  geom_text(data=counts, aes(label=n, y=1.25), position=position_dodge(0.8)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = c(0.72, 0.58),
        plot.margin=unit(c(0.2,0.2,0.2,0.5), "cm"),
        legend.background = element_rect(fill = "transparent", 
                                         color = "transparent"))

# p2 comes from bin_pnps_v3.R

ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(0.5,0.5),labels = c("A", "B"))
ggsave("F4_pnps_depth_bin_pnps.png", plot = last_plot())



# last supplementary figure
bin_taxon %>%
  ggplot(aes(x=median_bin_pnps, y=`complete genome size (Mbp)`)) +
  # geom_point(aes(color=depth)) +
  geom_smooth(method = lm) + 
  theme_classic()+
  labs(color = "Depth", x ="Median pN/pS of each MAG")+
  scale_x_continuous(trans = 'log10')

# statistic analysis
# no relation between bin pnps and transposase
bin_taxon$log_bin_pnps <- log(bin_taxon$median_bin_pnps)
summary(lm(`transposase gene calls in genome (%)`~log_bin_pnps, data = bin_taxon))

# deep ocean bins dont have higher pnps
bin_taxon <- bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")) %>%
  mutate(is_MES = ifelse(depth == "MES", "MES", "SRF&DCM"))
t.test(log_bin_pnps ~ is_MES, data = bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")))



# trash

reg <- pn_ps_metagenome %>% 
  filter(pnps < 1) %>%
  filter(!is.na(pnps)) %>%
  filter(gene_type %in% c("defense", "n")) %>% 
  filter(depth != "Malaspina") %>%
  mutate(depth2 = ifelse(depth == "MES", "MES", "SRF&DCM"))

reg$gene_type <- factor(reg$gene_type, levels = c("n", "defense"))
depth_defense.lm <- lm(pnps~depth2*gene_type, data = reg) 
summary(depth_defense.lm)

mean((filter(reg, gene_type == "defense"))$pnps)
mean((filter(reg, gene_type == "n"))$pnps)
filter(reg, gene_type == "n") %>% nrow()

