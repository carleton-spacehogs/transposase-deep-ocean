setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()
transposase_in_bins <- init_transposase_in_bins()

# bin pnps and transposases => nothing correlative (use the old, a few fasta map all-bins-concat approch)
low_trans <- c("Verrucomicrobia", "Acidimicrobidae")
high_trans <- c("Alphaproteobacteria", "Betaproteobacteria", "Deltaproteobacteria", 
                "Gammaproteobacteria", "Sphingobacteria", "Actinobacteria")

all <- filter_outliers(bin_taxon, "percent_trans") %>% 
  select("complete genome size (Mbp)", "percent_trans", "log_percent_trans", 
         "Class with more than 10 MAGs", "depth", "Class", "percent_biofilm", 
         "log_percent_biofilm") %>%
  filter(!is.na(depth))
all <- all%>%  #filter_outliers(all,"complete genome size (Mbp)") %>%
  mutate(size_bin=cut(`complete genome size (Mbp)`, breaks = c(0,2,3,4,5,20)))
  # mutate(size_bin = cut_interval(`complete genome size (Mbp)`, 5))
all$graphing_log_trans <- ifelse(all$log_percent_trans < -9, -2, all$log_percent_trans)
all$graphing_log_biofilm <- ifelse(all$log_percent_biofilm < -9, -2, all$log_percent_biofilm)
all$class_trans <- ifelse(all$Class %in% low_trans, "low",
                          ifelse(all$Class %in% high_trans, "high", "normal"))

all$size_bin <- as.character(all$size_bin)
all$size_bin[all$size_bin == '(0,2]'] <- '0-2'
all$size_bin[all$size_bin == '(2,3]'] <- '2-3'
all$size_bin[all$size_bin == '(3,4]'] <- '3-4'
all$size_bin[all$size_bin == '(4,5]'] <- '4-5'
all$size_bin[all$size_bin == '(5,20]'] <- '>5'
all$size_bin <- factor(all$size_bin, levels = c(">5","4-5","3-4","2-3","0-2"))

ggplot(all, aes(x=fct_rev(size_bin), y=percent_trans)) + 
  geom_jitter(aes(color = class_trans), alpha = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_color_manual(labels = c("High %-transposase", "Low %-transposase", "Normal %-transposase"),
                     values = c('orange', "green",'gray'))+
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.03)) + 
  facet_wrap(~depth, ncol= 1) +
  labs(color='Taxonomical Classes of:') +
  #coord_flip()+
  ylab("% of transposase ORFs in genome") +
  xlab("MAG genome size (Mbp)") +
  theme_classic()


# complete genome size (Mbp) is approximately normally distributed 
hist(all$`complete genome size (Mbp)`)

old.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm, all)
linear_bin_size.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm+`complete genome size (Mbp)`, all)
quadra_bin_size.lm <- lm(log_percent_trans~depth+Class+log_percent_biofilm
                         +`complete genome size (Mbp)`+ I(`complete genome size (Mbp)`^2), all)

anova(quadra_bin_size.lm)
get_r(linear_bin_size.lm)
get_r(quadra_bin_size.lm)

all %>% ggplot(aes(y = percent_trans, x = `complete genome size (Mbp)`)) +
  facet_wrap(~depth) +
  ylab("% of transposase ORFs in genome") +
  geom_point(aes(color = class_trans), alpha = 0.7) +
  geom_smooth(se=F) + 
  scale_color_manual(labels = 
    c("High %-transposase",
      "Normal %-transposase", 
      "Low %-transposase"),
    values = c('darkorange','gold', "khaki")) +
  guides(col = guide_legend(nrow = 9)) +
  labs(color='Taxonomical Classes of:') +
  theme_classic()












graph_taxon <- as.factor(c("Acidimicrobidae", "Flavobacteria", "Gammaproteobacteria"))
# graph_taxon <- factor(graph_taxon, levels= c("Acidimicrobidae", "Flavobacteria","Alphaproteobacteria", "Gammaproteobacteria"))

g_tax <- bin_taxon[,c("Class", "biofilm_count", "log_percent_trans", "depth", "percent_trans")] %>% 
  filter(Class %in% graph_taxon) %>% filter(!is.na(depth)) %>%
  mutate(depth = fct_rev(depth)) 

lm(percent_trans~ `complete genome size (Mbp)`*depth + I(`complete genome size (Mbp)`^2), 
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
