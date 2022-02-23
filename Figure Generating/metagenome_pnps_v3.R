setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(pn_ps_integron, pn_ps_bins, pn_ps_total, malaspina_total, malaspina_integron_all) %=% init_integron()

scale <- c(-2, -1.5, -1, -0.5, 0, 0.5)
log_scale <- c("less\nthan\n0.01", "0.032", "0.1", "0.316", "1", "3.162")
cols <- c("gene_type", "pnps", "log_pnps")
pn_ps_merged <- rbind(pn_ps_total[,cols] %>% 
                        filter(gene_type == "normal"), # non-transposase, non-cassette 
                      malaspina_integron_all[,cols],
                      malaspina_total[,cols],
                      pn_ps_integron[,cols])

pn_ps_merged$gene_type<-factor(pn_ps_merged$gene_type,levels=c("normal","defense","non_defense", "no_call"))

pCassette <- pn_ps_merged %>%
  ggplot(aes(x = gene_type, y = log_pnps, fill=gene_type)) +
  geom_boxplot(outlier.color = "gray", notch = TRUE) +
  ylab("pN/pS ratio") +
  scale_y_continuous(breaks = scale, labels = log_scale, limits = c(-2, 0.603)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.17)) +
  scale_x_discrete(labels=c("Non-\ntransposon\nORFs", 
                           "Defense\nmechanism \ncassette",
                           "Cassette with\nnon-defense\nCOG-function",
                           "Cassette\nwithout\nCOG-calls")) +
  theme_classic() +
  theme(axis.title.x = element_blank(), legend.position="none") + 
  scale_fill_manual(values=c("transparent", "red", "orange", "yellow"))

# pCassette used in depth_expression_corr.R

# statistic supporting the section:
# "Signatures of selection in defense mechanism ORFs within cassette sequences"

pn_ps_integron %>% filter(!is.na(pnps)) %>% nrow() # 1286
malaspina_integron_all %>% filter(!is.na(pnps)) %>% nrow() # 447

pnps_cassette<- (pn_ps_merged %>% filter(gene_type != "normal") %>% filter(!is.na(pnps)))$pnps
pnps_cassette<- (pn_ps_merged %>% filter(gene_type != "normal") %>% filter(pnps < 100))$pnps
# if don't filter inf: 95% quantile: 7.1, ln = 1.96
# if filter out inf: 95% quantile: 2.6, ln = 0.9555
quantile(pnps_cassette, probs = c(0.935), na.rm = TRUE) 
# e^((1.96+0.955)/2) = 4.3 -> upper cutoff: pnps = 4

pn_ps_total[,cols] %>% filter(gene_type == "normal") %>% filter(pnps < 4) %>% nrow()
malaspina_total %>% filter(pnps < 4) %>% nrow()

with_pnps <- pn_ps_merged %>% filter(gene_type != "normal") %>% filter(!is.na(pnps)) # 1753
with_pnps_fun <- with_pnps %>% filter(gene_type != "no_call") # 617

pnps_less4 <- with_pnps %>% filter(pnps < 4) # 1651
func_less4 <- pnps_less4 %>% filter(gene_type != "no_call")
no_func_less4 <- pnps_less4 %>% filter(gene_type == "no_call")
defense_less4 <- pnps_less4 %>% filter(gene_type == "defense")
median(pnps_less4$pnps)
median(func_less4$pnps)
median(no_func_less4$pnps)
median(defense_less4$pnps)
t.test(func_less4$pnps, no_func_less4$pnps)

t.test(pnps~gene_type, data=func_less4)
boxplot(pnps~gene_type, data=has_COG1)


background_pnps <- pn_ps_merged %>% filter(gene_type == "normal") %>% filter(pnps < 4)
median(background_pnps$pnps)

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


