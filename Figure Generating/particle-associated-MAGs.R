setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins, low_trans, high_trans) %=% init_bins()
quantitative_particle_association = init_MAGs_CAZenzyme_peptide()
quantitative_MAG_eukaryote = init_MAGs_eukaryote()

selected <- c("bin","percent_trans", "log_percent_trans", "median_bin_pnps", "depth", 
              "size_fraction", "complete genome size (Mbp)", "Class", "Genus","Order",
              "Total ORFs")

all_x <- bin_taxon[selected] %>% filter(!is.na(depth)) 
all_y <- malaspina_bins[selected] %>% filter(size_fraction != "error")

all_p <- rbind(all_x, all_y) # %>%
#   mutate(log_bin_pnps = log10(median_bin_pnps)) %>%
#   mutate(log_trans_graph = ifelse(percent_trans == 0, -2, log10(percent_trans))) %>%
#   mutate(log_trans_reg = ifelse(percent_trans == 0, -5, log10(percent_trans)))

MAGs_p = merge(all_p, quantitative_particle_association, by = "bin")

p_euk = merge(MAGs_p, quantitative_MAG_eukaryote, bu = "bin")
p_euk$abundance_sect_CAZ = p_euk$signal_CAZ_count/p_euk$`Total ORFs`
p_euk$abundance_sect_pep = p_euk$signal_pep_count/p_euk$`Total ORFs`
p_euk$avg_sect_abundance = (p_euk$abundance_sect_CAZ + p_euk$abundance_sect_pep)/2

plot(`complete genome size (Mbp)` ~ sqrt(percent_sect_CAZ), data = p_euk)

plot(log(`complete genome size (Mbp)`) ~ sqrt(percent_sect_CAZ), 
     data = filter(p_euk, eukarya_prop < 0.05))

get_r(lm(log(`complete genome size (Mbp)`) ~ log(abundance_sect_CAZ + 0.0005),
         data = filter(p_euk, eukarya_prop < 0.05)))
summary(lm(log(`complete genome size (Mbp)`) ~ Class + depth + log(abundance_sect_CAZ + 0.0005),
           data = filter(p_euk, eukarya_prop < 0.05)))

tmp <- all_p %>% group_by(Order, depth) %>%
  count() %>% filter(n > 5)

p_euk %>%
  filter(eukarya_prop < 0.03) %>%
  # filter(abundance_sect_CAZ > 0) %>%
  filter(!grepl("MED",bin)) %>%
  ggplot(aes(y = log10(`complete genome size (Mbp)`), x = log(abundance_sect_CAZ + 0.0005))) +
  geom_point() +
  facet_wrap(~depth, ncol = 2) + 
  geom_smooth(method = lm, se = TRUE) +
  xlab("log proportion of secretory CAZyme ORFs in MAGs")

p_euk %>%
  filter(eukarya_prop < 0.03 & percent_sect_CAZ != 100) %>%
  ggplot(aes(x = sqrt(percent_sect_CAZ), y = log(abundance_sect_CAZ))) + 
  geom_point() +
  labs(y="log(prop secretory CAZyme in TOTAL MAG ORFs)",
       x="Sqrt(percent secretory CAZyme among all CAZyme ORFs in a MAG)")

p_euk %>%
  filter(eukarya_prop < 0.03) %>%
  ggplot(aes(x = sqrt(percent_sect_pep), y = log(abundance_sect_pep))) + 
  geom_point() +
  labs(y="log(prop secretory peptidase in TOTAL MAG ORFs)",
       x="Sqrt(percent secretory peptidase among all peptidase ORFs in a MAG)")


p_euk %>%
  ggplot(aes(x = log(abundance_sect_CAZ + 0.0005), y = size_fraction)) +
  xlab("log(prop secretory CAZyme in TOTAL MAG ORFs)") +
  geom_boxplot()




