setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(pn_ps_integron, pn_ps_bins, pn_ps_total, tara_integron_summary, 
  malaspina_integron_summary, malaspina_total, malaspina_integron_all) %=% init_integron()

scale <- c(-2, -1.5, -1, -0.5, 0, 0.5)
log_scale <- c("less\nthan\n0.01", "0.032", "0.1", "0.316", "1", "3.162")
cols <- c("gene_type", "pnps", "log_pnps")
pn_ps_merged <- rbind(pn_ps_total[,cols] %>% 
                        filter(gene_type == "normal"), # non-transposase, non-cassette 
                      malaspina_integron_all[,cols],
                      malaspina_total[,cols],
                      pn_ps_integron[,cols])

pn_ps_merged$gene_type<-factor(pn_ps_merged$gene_type,levels=c("normal","defense","non_defense", "no_call"))

pB <- pn_ps_merged %>%
  ggplot(aes(x = gene_type, y = log_pnps, fill=gene_type)) +
  geom_boxplot(outlier.color = "gray", notch = TRUE) +
  ylab("pN/pS ratio") +
  scale_y_continuous(breaks = scale, labels = log_scale, limits = c(-2, 0.603)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.17)) +
  scale_x_discrete(labels=c("Non-tranposase\nnon-cassette\nORFs", 
                           "Defense\nmechanism\ncassette",
                           "Cassette with\nnon-defense\nCOG-function",
                           "Cassette\nwithout\nCOG-calls")) +
  theme_classic() +
  theme(axis.title.x = element_blank(), legend.position="none") + 
  scale_fill_manual(values=c("transparent", "red", "orange", "yellow"))

# pB used in depth_expression_corr.R

has_COG <- pn_ps_integron %>% filter(gene_type %in% c("non_defense", "defense")) %>% filter(pnps < 4)
t.test(log_pnps~gene_type, data=has_COG)

has_COG1 <- pn_ps_merged %>% filter(gene_type %in% c("non_defense", "defense")) %>% filter(pnps < 4)
t.test(pnps~gene_type, data=has_COG1)
boxplot(pnps~gene_type, data=has_COG1)

pn_ps_integron %>% filter(!is.na(pnps) & pnps< 4) %>% nrow() # 1208
malaspina_integron_all %>% filter(!is.na(pnps) & pnps< 4) %>% nrow() # 443

pn_ps_total[,cols] %>% filter(gene_type == "normal") %>% filter(pnps < 4) %>% nrow()
malaspina_total %>% filter(pnps < 4) %>% nrow()

pnps_cassette<- (pn_ps_merged %>% filter(gene_type != "normal") %>% filter(pnps < 100))$pnps
quantile(pnps_cassette, probs = c(0.97), na.rm = TRUE)


