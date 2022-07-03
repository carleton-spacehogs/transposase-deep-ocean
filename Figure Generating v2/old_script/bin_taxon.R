setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

trans_stats_taxon <- bin_taxon %>%
  group_by(`Class with more than 10 MAGs`) %>%
  summarise(mean_trans_percent = round(mean(`transposase gene calls in genome (%)`, na.rm = TRUE), digit = 3), 
            mean_biofilm_percent = round(mean(`biofilm gene calls in genome (%)`, na.rm = TRUE), digit = 3), 
            mean_trans_count = round(mean(Transposases, na.rm = TRUE), digit =2), 
            mean_biofilm_count = mean(biofilm_count, na.rm = TRUE), 
            mean_total_ORF = round(mean(`Total ORFs`), digit = 1),
            mean_complete_total_ORF = round(mean(`complete_ORF_count`), digit = 1),
            mean_bin_size = round(mean(`bin_size(Mbp)`), digits = 3),
            `# bins total` = n()) %>% 
  filter(`# bins total` > 5)

depth_stats_taxon <- bin_taxon %>% 
  with(table(`Class with more than 10 MAGs`, depth)) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column(var = "Class with more than 10 MAGs")

taxon_transCount_depth <- merge(trans_stats_taxon, depth_stats_taxon, by="Class with more than 10 MAGs")
taxon_transCount_depth <- taxon_transCount_depth[order(-taxon_transCount_depth$mean_trans_percent),]

summary(lm(mean_trans_percent~mean_biofilm_percent, data = trans_stats_taxon))
summary(lm(`Percent Transposases`~Class, data = bin_taxon))

write.csv(taxon_transCount_depth,"taxon_transposases_summary.csv", row.names = FALSE)

to_graph1 <- trans_stats_taxon %>% 
  # mutate_at("Class with more than 10 MAGs", str_replace, "proteo", "pro.") %>% 
  mutate_at("Class with more than 10 MAGs", str_replace, "bacteria", ".") %>%
  mutate_at("Class with more than 10 MAGs", str_replace, "microbidae", ".") %>%
  mutate_at("Class with more than 10 MAGs", str_replace, "microbiae", ".") %>%
  mutate_at("Class with more than 10 MAGs", str_replace, "Alpha", "\u0251-") %>%
  mutate_at("Class with more than 10 MAGs", str_replace, "Beta", "\u03b2-") %>%
  mutate_at("Class with more than 10 MAGs", str_replace, "Gamma", "\u03b3-") %>%
  mutate_at("Class with more than 10 MAGs", str_replace, "Delta", "\u03b4-")

`%ni%` <- Negate(`%in%`)
dont_show <- c("SAR202-2", "novelClass_E", "Phycisphaerae", "Opitutae", "OM190", 
               "Gemmatimonadetes", "Planctomycetia","Others Or Unknown", "Flavo.")
low_trans <- c("Acidi.", "Verruco.")
graph_df <- to_graph1 %>% filter(`Class with more than 10 MAGs` %ni% dont_show)
to_graph1 <- to_graph1 %>%
  mutate(color = ifelse(`Class with more than 10 MAGs` %in% dont_show, "don't show", ifelse(`Class with more than 10 MAGs` %in% low_trans, "low trans", "high trans")))

to_graph1  %>%
  ggplot(aes(y = mean_trans_percent, x = mean_biofilm_percent)) +
  geom_point(aes(color = color), size = 3) + 
  # geom_smooth(method=lm, se = F)  +
  xlab("Mean % of biofilm-associated ORFs in a Class") +
  ylab("Mean % of transposase ORFs") +
  xlim(-0.03,0.35) + ylim(-0.02,0.33) +
  scale_color_manual(breaks=c("don't show", "low trans", "high trans"), 
                     values=c("gray", "sky blue","deepskyblue3")) +
  theme_classic()+
  theme(legend.position = "none") +
  geom_text_repel(aes(label=`Class with more than 10 MAGs`), 
                  data = graph_df, box.padding = 0.6) 

  # now use order, instead of classes, not included in the paper
trans_stats_taxon <- bin_taxon %>%
  group_by(Order) %>%
  summarise(`mean (% transposase in this Order)` = round(mean(`transposase gene calls in genome (%)`, na.rm = TRUE), digit = 3), 
            `mean (% biofilm in this Order's bin genome)` = round(mean(`biofilm gene calls in genome (%)`, na.rm = TRUE), digit = 3), 
            `# bins total` = n()) %>% 
  filter(`# bins total` > 5)

trans_stats_taxon %>% 
  ggplot(aes(y = `mean (% transposase in this Order)`, x = `mean (% biofilm in this Order's bin genome)`)) +
  geom_point() + geom_smooth(method=lm, se = F)  +
  xlim(0,0.3) + ylim(0,0.25) +
  geom_text(aes(label=Order), check_overlap = TRUE, hjust=-0.1)
summary(lm(`mean (% transposase in this Order)`~`mean (% biofilm in this Order's bin genome)`, data = trans_stats_taxon))

taxon_biofilm.lm <- lm(`transposase gene calls in genome (%)`~`biofilm gene calls in genome (%)` + Order, data = bin_taxon)
taxon_biofilm_depth.lm <- lm(`transposase gene calls in genome (%)`~`biofilm gene calls in genome (%)` + Order + depth, data = bin_taxon)
anova(taxon_biofilm.lm, taxon_biofilm_depth.lm)
summary(taxon_biofilm_depth.lm)

colnames(bin_taxon)
mean(filter(bin_taxon, Class=="Betaproteobacteria")$`Total ORFs`)
filter(bin_taxon, Class=="Betaproteobacteria")%>%select(`bin_size(Mbp)`)
#bin_size(Mbp)
