setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(Tara_pnps_integron, pn_ps_bins, Tara_pnps_non_cassette, 
  malaspina_pnps_non_cassette, malaspina_integron_all) %=% init_integron()

# statistic supporting the section:
# "Signatures of selection in defense mechanism ORFs within cassette sequences"

Tara_pnps_integron %>% filter(!is.na(pnps)) %>% nrow() # 1286
malaspina_integron_all %>% filter(!is.na(pnps)) %>% nrow() # 447

cols <- c("gene_type", "pnps", "log_pnps")
n_cassette <- "other\ncassettes"
ord <- c("defense","normal",n_cassette)
pn_ps_merged <- rbind(Tara_pnps_non_cassette[,cols],
                      Tara_pnps_integron[,cols],
                      malaspina_pnps_non_cassette[cols],
                      malaspina_integron_all[,cols])

pn_ps_merged <- pn_ps_merged %>%
  mutate(gene_type = ifelse(gene_type %in% c("non_defense","no_call"), 
                            n_cassette, gene_type)) %>%
  mutate(gene_type = factor(gene_type, levels = ord))

pnps_less4 <- pn_ps_merged %>% filter(pnps < 4)

stat.pnps <- pnps_less4 %>% pairwise_t_test(log_pnps ~ gene_type)
stat.pnps.g <- stat.pnps[c(2,3),]
stat.pnps.g <- stat.pnps.g %>% add_xy_position()
stat.pnps.g$y.position <- c(0.73, 1.09)

pnps_scale <- c(-2, -1.5, -1, -0.5, 0, 0.5)
log_pnps_scale <- c("< 0.01    ", "0.032", "0.1", "0.316", "1", "3.162")

pn_ps_merged %>%
  filter(log_pnps < 0.603) %>%
  ggplot(aes(x = gene_type, y = log_pnps)) +
  geom_boxplot(outlier.color = "gray", aes(fill=gene_type)) + #notch = TRUE,
  ylab(expression(paste(italic("pN/pS"), " ratio"))) +
  scale_y_continuous(breaks = pnps_scale, 
                     labels = log_pnps_scale, 
                     limits = c(-2, 1.1)) + # 
  stat_summary(fun.data = boxplot.give.nr, 
               geom = "text", position=position_nudge(x = 0.2, y = 0.8)) +
  stat_pvalue_manual(stat.pnps.g, 
                     tip.length = 0.01,
                     label = "p.signif") + 
  scale_x_discrete(labels=c("defense\nmechanism\ncassettes",
                            "non-cassettes",
                            "non-defense\ncassettes")) +
  theme_classic() +
  theme(axis.title.y = element_blank(), legend.position="none") + 
  scale_fill_manual(values=c("red", "transparent","yellow")) +
  coord_flip()

ggsave("AGU_cassette_pnps.png", plot = last_plot(),
       height = 2, width = 4.5)



