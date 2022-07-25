setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins, low_trans, high_trans) %=% init_bins()

selected <- c("bin","percent_trans", "log_percent_trans", "depth",
              "Class with more than 10 MAGs", "Class_25_MAGs", "Completeness",
              "complete genome size (Mbp)", "Class", "Genus","Order", "Total ORFs")
bin_taxon$Completeness = bin_taxon$`Percent Complete`

all_x <- bin_taxon[selected] %>% filter(!is.na(depth)) 
all_y <- malaspina_bins[selected]
all_p <- rbind(all_x, all_y)

MAGs_p = merge_MAGs_eukaryote(merge_MAGs_CAZyme_peptide(all_p))
MAGs_p = filter(MAGs_p, eukarya_prop < 0.05)

colourCount = length(unique(MAGs_p$Class_25_MAGs))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
MAGs_p %>%
  # filter(log_abun_sect_CAZ > -8) %>%
  ggplot(aes(x=log_abun_sect_CAZ, y=log(`complete genome size (Mbp)`))) +
  facet_wrap(~depth) +
  geom_smooth(method = "lm", se = T) +
  # scale_color_manual(values = getPalette(colourCount)) +
  scale_color_viridis(option = "A") +
  geom_point(aes(color=log_percent_trans), alpha = 0.5) +
  theme_classic()

regress_MAGs = MAGs_p %>%
  filter(!is.na(Class)) %>%
  filter(log_abun_sect_CAZ > -8) %>%
  group_by(depth, Class) %>%
  mutate(n = n()) %>%
  filter(n > 11) %>%
  select(c("bin", "Class", "depth"))

MAGs_p2 = MAGs_p %>%
  filter(Class %in% regress_MAGs$Class)

MAGs_p %>%
  ggplot(aes(x=log_abun_sect_CAZ, y=log(`complete genome size (Mbp)`))) +
  facet_grid(depth~Class_25_MAGs) +
  geom_smooth(method = "lm", se = T) +
  scale_color_viridis(option = "A") +
  geom_point(alpha = 0.3, aes(color = log_percent_trans))

ggsave("per-taxon-CAZ-genomesize.pdf", plot = last_plot(),
       height = 7, width = 18)

# genome size (ORFs) compared to depth
MAGs_p %>%
  ggplot(aes(y=fct_rev(depth), x=log(`complete genome size (Mbp)`))) +
  facet_wrap(~Class, scale = "free_x") +
  # xlim(-1, max(MAGs_p$`complete genome size (Mbp)`)) +
  geom_boxplot()

pre_CAZ = lm(log(`complete genome size (Mbp)`) ~ Class + depth, data = MAGs_p)
post_CAZ = lm(log(`complete genome size (Mbp)`) ~ Class + depth + 
                log_abun_sect_CAZ, data = MAGs_p)
anova(pre_CAZ, post_CAZ)

get_r(pre_CAZ)
get_r(post_CAZ)
# > get_r(pre_CAZ)
# [1] 0.4694659
# > get_r(post_CAZ)
# [1] 0.5516766

# > get_r(pre_CAZ)
# [1] 0.4694659
# > get_r(post_CAZ)
# [1] 0.5307099

# > get_r(pre_CAZ)
# [1] 0.4694659
# > get_r(post_CAZ)
# [1] 0.5220544

trans1.lm = lm(log_percent_trans ~ depth + log_abun_sect_CAZ + 
                log(`complete genome size (Mbp)`), data = MAGs_p)
trans2.lm = lm(log_percent_trans ~ depth + log_abun_sect_CAZ + 
                log(`complete genome size (Mbp)`) + Class, data = MAGs_p)
summary(trans.lm)
get_r(trans1.lm)
get_r(trans2.lm)



plot(`complete genome size (Mbp)` ~ log(abundance_sect_CAZ), data = p_euk)

plot(log(`complete genome size (Mbp)`) ~ sqrt(percent_sect_CAZ), 
     data = filter(p_euk, eukarya_prop < 0.05))

get_r(lm(log(`complete genome size (Mbp)`) ~ log(abundance_sect_CAZ + 0.0005),
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




