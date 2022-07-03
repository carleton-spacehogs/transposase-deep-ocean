setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(Tara_pnps_integron, pn_ps_bins, Tara_pnps_non_cassette, 
  malaspina_pnps_non_cassette, malaspina_integron_all) %=% init_integron()

# statistic supporting the section:
# "Signatures of selection in defense mechanism ORFs within cassette sequences"

Tara_pnps_integron %>% filter(!is.na(pnps)) %>% nrow() # 1286
malaspina_integron_all %>% filter(!is.na(pnps)) %>% nrow() # 447

cols <- c("gene_type", "pnps", "log_pnps")
ord <- c("defense","normal","non_defense","no_call")
pn_ps_merged <- rbind(Tara_pnps_non_cassette[,cols],
                      Tara_pnps_integron[,cols],
                      malaspina_pnps_non_cassette[cols],
                      malaspina_integron_all[,cols]) %>%
  mutate(gene_type = factor(gene_type, levels = ord))

pnps_cassette <- pn_ps_merged %>% 
  filter(gene_type != "normal") %>% 
  filter(!is.na(pnps)) %>% 
  filter(pnps < 100) # take out INf

quantile(pnps_cassette$pnps, probs = c(0.97), na.rm = TRUE)
# 97% : 4.361391 
# upper cutoff: pnps = 4

pnps_less4 <- pn_ps_merged %>% filter(pnps < 4)

pnps_less4 %>%
  group_by(gene_type) %>%
  summarise(count = n(), median = median(pnps))

stat.pnps <- pnps_less4 %>% pairwise_t_test(log_pnps ~ gene_type)
stat.pnps.g <- stat.pnps[c(2,6),]
stat.pnps.g <- stat.pnps.g %>% add_xy_position()

pnps_scale <- c(-2, -1.5, -1, -0.5, 0, 0.5)
log_pnps_scale <- c("< 0.01", "0.032", "0.1", "0.316", "1", "3.162")

pn_ps_merged %>%
  filter(log_pnps < 0.603) %>%
  ggplot(aes(x = gene_type, y = log_pnps)) +
  geom_boxplot(outlier.color = "gray", aes(fill=gene_type)) + #notch = TRUE,
  ylab(expression(paste(italic("pN/pS"), " ratio"))) +
  scale_y_continuous(breaks = pnps_scale, 
                     labels = log_pnps_scale, 
                     limits = c(-2, 0.8)) + # 
  stat_summary(fun.data = boxplot.give.nr, 
               geom = "text", position=position_nudge(x = 0.2, y = 0.6)) +
  stat_pvalue_manual(stat.pnps.g, 
                     tip.length = 0.01,
                     label = "p.signif") + 
  # annotate(geom="text", x=3.2, y=0.68, label="***") + 
  scale_x_discrete(labels=c("defense\nmechanism\ncassettes",
                            "non-cassettes",
                            "\ncassettes of\nnon-defense\nfunction",
                            "cassettes\nwithout\nCOG-calls")) +
  theme_classic() +
  theme(axis.title.y = element_blank(), legend.position="none") + 
  scale_fill_manual(values=c("red", "transparent", "orange", "yellow")) +
  coord_flip()

ggsave("F6_cassette_pnps.png", plot = last_plot(),
       height = 2.8, width = 5)

ggsave("F6_cassette_pnps.pdf", plot = last_plot(),
       height = 2.8, width = 5)



x_lipid <- pn_ps_integron %>% filter(grepl("Lipid|Carbohydrate|Amino",category) )
y_lipid <- malaspina_integron_all %>% filter(COG20_CATEGORY == "I")
cassette_lipid <- c(x_lipid$pnps, y_lipid$pnps)

x_amino <- pn_ps_integron %>% filter(grepl("Amino",category) )
y_amino <- malaspina_integron_all %>% filter(COG20_CATEGORY == "E")
cassette_amino <- c(x_amino$pnps, y_amino$pnps)

median(cassette_lipid, na.rm = TRUE)
median(cassette_amino, na.rm = TRUE)

x_metabolism <- pn_ps_integron %>% filter(grepl("Inorganic|Lipid|Carbohydrate|Amino|Secondary metabolites",category) )
y_metabolism <- malaspina_integron_all %>% filter(grepl("P|E|G|I|Q",COG20_CATEGORY))
cassette_metabolism <- c(x_metabolism$pnps, y_metabolism$pnps)
median(cassette_metabolism, na.rm = TRUE)

# Small malaspina cassette categories: (only 3 and 1)
# lipid (accession: I) pnps: 0.327, 1.809, 1.778
# amino acid (accession: E) pnps: 0.775


