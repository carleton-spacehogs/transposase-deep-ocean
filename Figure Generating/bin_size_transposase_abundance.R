setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()
transposase_in_bins <- init_transposase_in_bins()

# bin pnps and transposases => nothing correlative (use the old, a few fasta map all-bins-concat approch)

lm(`transposase gene calls in genome (%)`~ `complete genome size (Mbp)`*depth + I(`complete genome size (Mbp)`^2), 
   data = bin_taxon) -> transposase_size_depth.lm
pred <- bin_taxon %>% select(`complete genome size (Mbp)`, depth)
bin_taxon$predicted_transposase_percent <- predict(transposase_size_depth.lm, pred)

outlier_cut <- min(boxplot(bin_taxon$`transposase gene calls in genome (%)`)$out)
p1 <- bin_taxon %>% 
  filter(`transposase gene calls in genome (%)` < outlier_cut) %>%
  filter(`complete genome size (Mbp)` > 2.4) %>%
  ggplot(aes(y = `transposase gene calls in genome (%)`, 
             x = `complete genome size (Mbp)`)) +
  facet_wrap(~depth) +
  ylab("% of transposase ORFs in genome") +
  geom_point(aes(color = `Class with more than 10 MAGs`), alpha = 0.6) +
  geom_smooth(se=F) +
  guides(col = guide_legend(nrow = 9)) +
  labs(color='Class with \u2265 10 MAGs') +
  theme_classic()
  
p2 <- bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")) %>%
  group_by(`Class with more than 10 MAGs`) %>%
  mutate(depth = fct_rev(depth)) %>% 
  ggplot(aes(x = depth, y = `complete genome size (Mbp)`)) +
  ylim(0,12.5) +
  stat_compare_means(comparisons = list(c("DCM","MES"), c("SRF","MES")), 
                     label = "p.signif", hide.ns = TRUE, 
                     method = "t.test", 
                     label.y = c(11,12)) +
  geom_jitter(aes(color = `Class with more than 10 MAGs`), alpha = 0.6) +
  geom_violin(alpha = 0.5) + 
  stat_summary(fun.data = boxplot.give.n, 
               geom = "text", 
               position=position_nudge(x = 0, y = 0)) +
  coord_flip() +
  theme_classic() + 
  theme(legend.position = "none")

p3 <- bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")) %>%
  group_by(`Class with more than 10 MAGs`) %>%
  mutate(depth = fct_rev(depth)) %>% 
  ggplot(aes(x = depth, y = `transposase gene calls in genome (%)`)) +
  ylab("% of transposase ORFs in genome")+
  stat_compare_means(comparisons = list(c("DCM","MES"), c("SRF","MES")), 
                     label = "p.signif", method = "t.test", 
                     tip.length = 0.01,
                     label.y = c(0.75, 0.85)) +
  # stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  xlab("") + 
  coord_flip(ylim = c(0, 0.85)) + 
  theme_classic() + 
  theme(legend.position = "none")

ggarrange(ggarrange(p2, p3, ncol = 2, labels = c("A", "B")), # Second row with box and dot plots
          p1,
          nrow = 2, 
          heights = c(2/5, 3/5),
          labels = c("", "  C"))

bin_taxon %>% 
  filter(depth %in% c("MES")) %>% 
  ggplot(aes(x = `complete genome size (Mbp)`, y = median_bin_pnps)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F)

# summary(lm(mean_bin_pnps~`complete genome size (Mbp)`, data = bin_taxon))
summary(lm(median_bin_pnps~`complete genome size (Mbp)`, 
           data = bin_taxon))


#  'viral' (<0.22 ??m), 'girus' (0.22-0.8 ??m), 'bacterial' (0.22-1.6 ??m), and 'protistan' (0.8-5.0 ??m),
bin_taxon %>%
  filter(size_fraction %in% c("girus", "prot")) %>%
  ggplot(aes(x = size_fraction, y = `transposase gene calls in genome (%)`)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0)) +
  geom_boxplot(alpha = 0.5)  + coord_flip(ylim = c(0, 1))

bin_taxon %>%
  filter(size_fraction %in% c("bact", "girus", "prot")) %>%
  ggplot(aes(x = size_fraction, y = `biofilm gene calls in genome (%)`)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.45, y = 0)) +
  geom_boxplot(alpha = 0.5)  + coord_flip()

# useless ?

tara_and_malaspina_complete_size <- bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")) %>%
  mutate(is_MES = ifelse(depth == "MES", "MES", "not_MES")) %>%
  select(is_MES, `complete genome size (Mbp)`) %>%
  rbind(malaspina_bins %>% select(is_MES, `complete genome size (Mbp)`))

tara_and_malaspina_complete_size %>%
  group_by(is_MES) %>%
  summarise(`mean bin genome size` = mean(`complete genome size (Mbp)`), 
            `median genome size` = median(`complete genome size (Mbp)`), 
            sd = sd(`complete genome size (Mbp)`),
            `# bins total` = n())

bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")) %>%
  group_by(depth) %>%
  summarise(`mean bin genome size` = mean(`complete genome size (Mbp)`), 
            `median genome size` = median(`complete genome size (Mbp)`), 
            `mean percentage` = mean(`transposase gene calls in genome (%)`),
            sd = sd(`complete genome size (Mbp)`),
            `# bins total` = n())
