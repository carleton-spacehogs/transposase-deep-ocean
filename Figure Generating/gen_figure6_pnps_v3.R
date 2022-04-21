setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(Tara_pnps_integron, pn_ps_bins, Tara_pnps_non_cassette, 
  malaspina_pnps_non_cassette, malaspina_integron_all) %=% init_integron()

# statistic supporting the section:
# "Signatures of selection in defense mechanism ORFs within cassette sequences"

pn_ps_integron %>% filter(!is.na(pnps)) %>% nrow() # 1286
malaspina_integron_all %>% filter(!is.na(pnps)) %>% nrow() # 447

pnps_cassette <- pn_ps_merged %>% 
  filter(gene_type != "normal") %>% 
  filter(!is.na(pnps)) %>% 
  filter(pnps < 100) # take out INf

quantile(pnps_cassette$pnps, probs = c(0.97), na.rm = TRUE)
# 97% 
# 4.361391 
# upper cutoff: pnps = 4

pn_ps_total[,cols] %>% filter(gene_type == "normal") %>% filter(pnps < 4) %>% nrow()
malaspina_total %>% filter(pnps < 4) %>% nrow()

with_pnps <- pn_ps_merged %>% filter(gene_type != "normal") %>% filter(!is.na(pnps)) # 1753
with_pnps_fun <- with_pnps %>% filter(gene_type != "no_call") # 617

# if I don't take out Inf
quantile(with_pnps_fun$pnps, c(.95, .97, .99))
# 95%      97%      99% 
# 3.182308 5.508244      Inf 

# if I take out Inf
pnps_no_inf <- with_pnps_fun %>% filter(pnps < 100)
quantile(pnps_no_inf$pnps, c(.95, .97, .99))
# 95%      97%      99% 
# 2.429197 3.312591 6.150663 

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


pnps_scale <- c(-2, -1.5, -1, -0.5, 0, 0.5)
log_pnps_scale <- c("< 0.01", "0.032", "0.1", "0.316", "1", "3.162")
cols <- c("gene_type", "pnps", "log_pnps")
pn_ps_merged <- rbind(Tara_pnps_non_cassette[,cols],
                      Tara_pnps_integron[,cols],
                      malaspina_pnps_non_cassette[cols],
                      malaspina_integron_all[,cols])

pn_ps_merged$gene_type<-factor(pn_ps_merged$gene_type,
                               levels=c("defense","normal","non_defense","no_call"))

df_p_val <- data.frame(
  group1 = "non_defense",
  group2 = "defense",
  p = "", # see the stat calculation above 0.0005
  y.position = 0.7
)

pn_ps_merged %>%
  filter(log_pnps < 0.603) %>%
  ggplot(aes(x = gene_type, y = log_pnps)) +
  geom_boxplot(outlier.color = "gray", aes(fill=gene_type)) + #notch = TRUE,
  ylab(expression(paste(italic("pN/pS"), " ratio"))) +
  scale_y_continuous(breaks = pnps_scale, 
                     labels = log_pnps_scale, limits = c(-2, 0.7)) + # 
  stat_summary(fun.data = boxplot.give.nr, 
               geom = "text", position=position_nudge(x = 0.2, y = 0.6)) +
  stat_pvalue_manual(df_p_val, tip.length = 0.01) + 
  annotate(geom="text", x=3.2, y=0.68, label="***") + 
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





