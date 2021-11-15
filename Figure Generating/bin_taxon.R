setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

# mean(bin_taxon$`transposase gene calls in genome (%)`)

trans_stats_taxon <- bin_taxon %>%
  group_by(`Class with more than 10 MAGs`) %>%
  summarise(mean_trans_percent = round(mean(`transposase gene calls in genome (%)`, na.rm = TRUE), digit = 3), 
            mean_biofilm_percent = round(mean(`biofilm gene calls in genome (%)`, na.rm = TRUE), digit = 3), 
            mean_trans_count = round(mean(Transposases, na.rm = TRUE), digit =2), 
            mean_total_ORF = round(mean(`Total ORFs`), digit = 1),
            mean_complete_total_ORF = round(mean(`complete_ORF_count`), digit = 1),
            mean_bin_size = round(mean(`bin_size(Mbp)`), digits = 3),
            `# bins total` = n()) %>% 
  filter(`# bins total` > 5)

depth_stats_taxon <- bin_taxon %>% 
  with(table(`Class with more than 10 MAGs`, depth)) %>% 
  as.data.frame.matrix() %>% 
  select(-c("mixed", "epi")) %>%
  rownames_to_column(var = "Class with more than 10 MAGs")

taxon_transCount_depth <- merge(trans_stats_taxon, depth_stats_taxon, by="Class with more than 10 MAGs")
taxon_transCount_depth <- taxon_transCount_depth[order(-taxon_transCount_depth$mean_trans_percent),]

summary(lm(mean_trans_percent~mean_biofilm_percent, data = trans_stats_taxon))
summary(lm(`Percent Transposases`~Class, data = bin_taxon))

write.csv(taxon_transCount_depth,"taxon_transposases_summary.csv", row.names = FALSE)

highlight_df <- trans_stats_taxon %>% filter(`Class with more than 10 MAGs`=="Others Or Unknown")

trans_stats_taxon %>% 
  ggplot(aes(y = mean_trans_percent, x = mean_biofilm_percent)) +
  geom_point() + 
  geom_point(data=highlight_df, 
             aes(y = mean_trans_percent, x = mean_biofilm_percent), 
             color='red',
             size=3)+
  annotate(geom="text",x = 0.065, y = 0.045,label = "Others/Unknown", color = "red") +
  geom_smooth(method=lm, se = F)  +
  xlab("Average % of biofilm-associated ORFs in a taxonomic class") +
  ylab("Average % of transposase ORFs") +
  xlim(0,0.3) + ylim(0,0.225) +
  theme_classic()+
  geom_text(aes(label=`Class with more than 10 MAGs`), check_overlap = TRUE, hjust=-0.1)


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
