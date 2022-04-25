setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(pn_ps_metagenome) %=% init_individual_metagenomes()

scale1 <- c(-2, -1.5, -1, -0.5, 0, 0.60206)
log_scale1 <- c("<0.01", "<0.032", "0.10", "0.32", "1", ">4")

to_graph <- pn_ps_metagenome %>%
  filter(log_pnps < 1.4) %>%
  filter(log_pnps > -3) %>%
  mutate(is_trans = ifelse(gene_type == "transposase", 'y', 'n'))

counts <- to_graph %>% count(depth, is_trans) %>% as.data.frame()
# counts$n <- formatC(counts$n, format = "e", digits = 2)

counts[1,3] <- "2.2e+06"
counts[3,3] <- "2.8e+06"
counts[5,3] <- "1.5e+06"
counts[7,3] <- "3.7e+05"

p1<- to_graph %>%
  ggplot(aes(x=fct_rev(depth), y=log_pnps, fill = is_trans)) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_continuous(breaks = scale1, labels = log_scale1, 
                     name = expression(paste(italic("pN/pS"), " of metagenomic ORFs")), 
                     limits = c(-1.5, 1.25)) +
  theme_classic() +
  # facet_wrap(~depth, nrow = 4, strip.position = "left") +
  guides(fill=guide_legend(title="", reverse = TRUE)) +
  scale_fill_manual(labels = c("Non-transposase", "Transposase"),
                     values= c("gray","steelblue")) +
  geom_text(data=counts, aes(label=n, y=1.05), position=position_dodge(0.8)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = c(0.7, 0.82),
        plot.margin=unit(c(0.2,0.2,0.2,0.5), "cm"),
        legend.background = element_rect(fill = "transparent", 
                                         color = "transparent"))

# p2 comes from bin_pnps_v3.R

ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(0.45,0.5),labels = c("A", "B"))
ggsave("F4_pnps_depth_bin_pnps.png", 
       plot = last_plot(),
       height = 4,
       width = 10)
ggsave("F4_pnps_depth_bin_pnps.pdf", 
       plot = last_plot(), 
       height = 4,
       width = 10)



