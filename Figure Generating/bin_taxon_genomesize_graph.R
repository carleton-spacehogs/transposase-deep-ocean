setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()
transposase_in_bins <- init_transposase_in_bins()

graph_taxon <- c("Acidimicrobidae", "Deltaproteobacteria", "Gammaproteobacteria")

x_tax <- bin_taxon[,c("Class", "biofilm_count", "log_percent_trans", 
                      "depth", "percent_trans", "is_biofilm")] %>% 
  filter(Class %in% graph_taxon) %>% filter(!is.na(depth))

y_tax <- malaspina_bins[,c("Class", "biofilm_count", "log_percent_trans", 
                          "depth", "percent_trans", "is_biofilm")]%>% 
  filter(Class %in% graph_taxon)

g_tax <- filter_outliers(rbind(x_tax, y_tax), "percent_trans")
# g_tax <- rbind(x_tax, y_tax)
# g_tax <- filter_outliers(x_tax, "percent_trans")
g_tax$Class <- factor(g_tax$Class, levels = graph_taxon)

g_tax$g_depth <- gsub('Deep Malaspina', 'Deep\nMalaspina', g_tax$depth)
g_tax$g_depth <- factor(g_tax$g_depth, levels = c("SRF", "DCM", "MES", 'Deep\nMalaspina'))
g_tax$is_biofilm <- gsub('absent', 'Biofilm\nabsent', g_tax$is_biofilm)
g_tax$is_biofilm <- gsub('present', 'Biofilm\n     present', g_tax$is_biofilm)
p_depth <- g_tax %>%
  ggplot(aes(x=percent_trans, y = fct_rev(g_depth))) +
  facet_wrap(~Class, ncol = 4) +
  geom_boxplot(outlier.alpha = 0.3) + 
  ylab("depth")+
  stat_summary(fun.data = boxplot.give.n, geom = "text",
               position=position_nudge(x = 0.1, y = 0)) + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank())

p_biofilm <- g_tax %>% 
  ggplot(aes(y = percent_trans, x = fct_rev(is_biofilm))) +
  facet_wrap(~Class, ncol = 4) +
  geom_boxplot(outlier.alpha = 0.3) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.1)) + 
  ylab("% of transposase ORFs in MAG") +
  # xlab("biofilm ORF")+
  coord_flip() + 
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title.y=element_blank())

ggarrange(p_depth, p_biofilm, ncol = 1, 
          nrow = 2, labels = c("A", "B"),
          heights = c(4/7,3/7))
ggsave("some_taxa_trans.png", plot = last_plot())






all_tara <- filter_outliers(bin_taxon, "percent_trans") %>% 
  select("complete genome size (Mbp)", "percent_trans", "log_percent_trans", 
         "depth", "Class", "class_trans", "graphing_log_trans") %>%
  filter(!is.na(depth))

all_mala <- malaspina_bins[,c("complete genome size (Mbp)", "percent_trans", "log_percent_trans", 
                           "depth", "Class", "class_trans", "graphing_log_trans")]

all <- filter_outliers(rbind(all_tara, all_mala),"percent_trans") %>%
  mutate(size_bin=cut(`complete genome size (Mbp)`, breaks = c(0,2,3,4,5,20)))
# mutate(size_bin = cut_interval(`complete genome size (Mbp)`, 5))

all$g_depth <- gsub('Deep Malaspina', 'Deep\nMalaspina', all$depth)
all$g_depth <- factor(all$g_depth, levels = c("SRF", "DCM", "MES", 'Deep\nMalaspina'))

all$size_bin <- as.character(all$size_bin)
all$size_bin[all$size_bin == '(0,2]'] <- '0-2'
all$size_bin[all$size_bin == '(2,3]'] <- '2-3'
all$size_bin[all$size_bin == '(3,4]'] <- '3-4'
all$size_bin[all$size_bin == '(4,5]'] <- '4-5'
all$size_bin[all$size_bin == '(5,20]'] <- '>5'
all$size_bin <- factor(all$size_bin, levels = c(">5","4-5","3-4","2-3","0-2"))

p1<-ggplot(all, aes(x=fct_rev(size_bin), y=percent_trans)) + 
  geom_jitter(aes(color = class_trans), alpha = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_color_manual(labels = c("High %-transposase", "Normal %-transposase", "Low %-transposase"),
                     values = c('orange','gray', "green"))+
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.03)) + 
  facet_wrap(~depth, ncol= 4) +
  # labs(color='Taxonomical Classes of:') +
  #coord_flip()+
  ylab("%-transposase ORF") +
  xlab("complete genome size (Mbp) of MAGs") +
  theme_classic() +
  theme(legend.position = "none")

# filter_outliers(all, "complete genome size (Mbp)")
p2<-all %>% 
  mutate(depth = fct_rev(g_depth)) %>% 
  ggplot(aes(x = depth, y = `complete genome size (Mbp)`)) +
  ylim(0,9) +
  stat_compare_means(comparisons = list(c("DCM","MES"), c("SRF","MES")), 
                     label = "p.signif", hide.ns = TRUE, method = "t.test", 
                     tip.length = 0.01, label.y = c(8,9)) +
  scale_color_manual(labels = c("High %-transposase (\u0251- \u03b2- \u03b3- proteobact.,\n and Actinobacteria)\n", 
                                "Normal transposase abundance\n", 
                                "Low %-transposase\n(Verrucomicrobia, Acidimicrobidae)"),
                     values = c('orange','gray', "green"))+
  labs(color='Taxonomical Class of:') +
  geom_jitter(aes(color = class_trans), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.2)) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y=element_blank())

# p3 comes from bin_pnps_v3.R

ggarrange(p1,p2,p3,
          nrow = 3, 
          # heights = c(3/5,2/5,3/5),
          labels = c("A", "B", "C"))

ggsave("trans_genome_size.png", plot = last_plot())








