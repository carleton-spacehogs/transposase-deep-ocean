setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins, low_trans, high_trans) %=% init_bins()
quantitative_particle_association = init_MAGs_CAZenzyme_peptide()
quantitative_MAG_eukaryote = init_MAGs_eukaryote()

selected <- c("bin","percent_trans", "log_percent_trans", "median_bin_pnps", "depth", 
              "size_fraction", "complete genome size (Mbp)", "Class", "Genus","Order")

all_x <- bin_taxon[,selected] %>% filter(!is.na(depth)) 
all_y <- malaspina_bins[,selected] %>% filter(size_fraction != "error")

all_p <- rbind(all_x, all_y) %>%
  mutate(log_bin_pnps = log10(median_bin_pnps)) %>%
  # -2 is for looks nice on graph
  mutate(log_trans_graph = ifelse(percent_trans == 0, -2, log10(percent_trans))) %>%
  # -5 is an arbitrary small number that's close to 0
  mutate(log_trans_reg = ifelse(percent_trans == 0, -5, log10(percent_trans)))

all_p = merge(all_p, quantitative_particle_association, by = "bin")

p_euk = merge(all_p, quantitative_MAG_eukaryote, bu = "bin")

plot(`complete genome size (Mbp)` ~ percent_CAZ, data = all_p)

plot(log(`complete genome size (Mbp)`) ~ sqrt(avg_percent), data = p_euk)

plot(log(`complete genome size (Mbp)`) ~ sqrt(avg_percent),
     data = filter(p_euk, eukarya_prop < 0.01))

plot(log(`complete genome size (Mbp)`) ~ sqrt(percent_CAZ), data = all_p)

get_r(lm(log(`complete genome size (Mbp)`) ~ Class, data = all_p))
summary(lm(log(`complete genome size (Mbp)`) ~ Class + depth + sqrt(avg_percent),
           data = filter(p_euk, eukarya_prop < 0.05)))

tmp <- all_p %>% group_by(Order, depth) %>%
  count() %>% filter(n > 5)

all_p %>%
  filter(Order == "Flavobacteriales") %>%
  ggplot(aes(y = `complete genome size (Mbp)`, x = avg_percent)) +
  geom_point() +
  facet_wrap(~depth, ncol = 2) + 
  geom_smooth(se = TRUE,method = lm)

all_p %>%
  filter(Order == "Acidimicrobiales") %>%
  ggplot(aes(y = `complete genome size (Mbp)`, x = percent_CAZ)) +
  geom_point() +
  facet_wrap(~depth, ncol = 2) + 
  geom_smooth(se = TRUE,method = lm)


all_p %>%
  ggplot(aes(y = percent_CAZ, x = size_fraction)) +
  geom_boxplot()

all_p %>%
  filter(Class == "Alphaproteobacteria") %>%
  ggplot(aes(y = avg_percent, x = depth)) +
  geom_boxplot()




