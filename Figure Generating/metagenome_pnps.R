setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(pn_ps_integron, pn_ps_bins, pn_ps_total, tara_integron_summary, malaspina_integron_summary, tara_integron_all, malaspina_integron_all) %=% init_integron()

pn_ps_merged <- rbind(select(pn_ps_total, gene_type, pnps, log_pnps), 
                      select(pn_ps_integron, pnps, gene_type, log_pnps))
pn_ps_merged$gene_type <- factor(pn_ps_merged$gene_type, levels = c(
  "transposase", "all cassette genes", # transposase and all cassette genes is not currently used
  "normal", "defense", "non_defense", "no_call"))

scale <- c(-3, -2, -1, 0, 0.60, 1)
log_scale <- c("lower\nthan\n0.001", "0.01", "0.1", "1", "max \npnps:4", "10")
pB <- pn_ps_merged %>%
  filter(pnps < 4) %>%
  filter(gene_type %in% c("normal","non_defense", "defense", "no_call")) %>%
  ggplot(aes(x = gene_type, y = log_pnps, fill=gene_type)) +
  geom_boxplot(outlier.color = "gray") +
  ylab("pN/pS ratio") +
  scale_y_continuous(breaks = scale, labels = log_scale, limits = c(-3, 1)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.17)) +
  stat_compare_means(comparisons =list(c("non_defense", "defense"),
                                       c("no_call","non_defense")),
                     label.y = c(0.6, 0.8),
                     method = "t.test",
                     tip.length = 0.005) + 
  scale_x_discrete(labels=c("normal (not\ntranposase\nor cassette)", 
                            "Defense\nmechanism\ncassette ORFs",
                            "cassette with\nnon-defense\nCOG-function",
                            "cassette\nwithout\nCOG calls")) +
  theme_classic() +
  theme(axis.title.x = element_blank(), legend.position="none") + 
  scale_fill_manual(values=c("transparent", "red", "orange", "yellow"))

# pB used in conjunction with pC in bin_lifestyle.R


pn_ps_merged %>% 
  filter(pnps < 10) %>%
  ggplot(aes(x = gene_type, y = pnps)) +
  #geom_violin() + 
  #ylim(0,0.95) +
  geom_boxplot(outlier.colour=NA) + 
  coord_cartesian(ylim = c(0, 1.05)) + ylab("pN/pS") + xlab ("") +
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  stat_compare_means(comparisons = list( c("all cassette genes", "normal"), 
                                         c("transposase", "normal"),
                                         c("all cassette genes", "non-defense cassette gene calls"),
                                         c("non-defense cassette gene calls", "defense mech")),
                     label.y = c(1.05, 0.8, 0.9, 0.75),
                     tip.length = 0.003)


p1 <- pn_ps_total %>% 
  filter(pnps < 10) %>%
  ggplot(aes(x = gene_type, y = pnps)) +
  #geom_violin() + 
  #ylim(0,0.95) +
  geom_boxplot(outlier.colour=NA) + 
  coord_cartesian(ylim = c(0, 1.05)) + ylab("pN/pS") + xlab ("") +
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  stat_compare_means(comparisons = list( c("all cassette genes", "normal"), c("transposase", "normal") ),
                     label.y = c(1.05, 0.8),
                     tip.length = 0.003)

ggarrange(p1, p2, ncol = 2, nrow = 1)

pn_ps_integron <- pn_ps_integron %>% filter(!is.na(category)) %>%
  mutate(is_metabolism = ifelse(category %in% c("Lipid transport and metabolism","Nucleotide transport and metabolism"), "lipid", "non-lipid"))
# mutate(is_metabolism = ifelse(category %like% "metabolism", "metabolism", "non-metabolism"))

t.test(log_pnps~ is_metabolism, data = pn_ps_integron %>% 
         filter(pnps < 10) %>% filter(!is.na(category)))


