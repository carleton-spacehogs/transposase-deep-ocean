setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()

# from bin_taxon_genomesize_statistics.R
low_trans <- c("Flavobacteria","Acidimicrobidae","novelClass_E",
               "Gemmatimonadetes","SAR202-2","Marinisomatia")
high_trans <- c("Alphaproteobacteria","Gammaproteobacteria",
                "Betaproteobacteria","Actinobacteria")

normal_taxa <- setdiff(setdiff(big_taxa, less_taxa), more_taxa)
normal_taxa_rows <- rbind(bin_taxon%>%select("Class"), malaspina_bins%>%select("Class")) %>%
  filter(Class %in% normal_taxa) 
as.data.frame(table(normal_taxa_rows[, "Class"]))
table(all$Class, all$size_fraction)
# Class Freq
#         Bacteroidia   30
#     Dehalococcoidia   27
# Deltaproteobacteria   35
#               OM190   11
#            Opitutae   16
#       Phycisphaerae   24
#      Planctomycetia   32
#     Sphingobacteria   14
#     Verrucomicrobia   43 # use this because it has the most MAGs

graph_taxon <- c("Flavobacteria", "Verrucomicrobia", "Gammaproteobacteria")
graph_taxon1 <- c("Flavobacteria", "Verrucomicrobia", "\u03b3-proteobacteria")

col_list <- c("complete genome size (Mbp)", "percent_trans", "log_percent_trans", 
              "percent_biofilm", "is_biofilm", "size_fraction",
              "depth", "Class", "class_trans", "graphing_log_trans", "graphing_log_biofilm")

x_tax <- bin_taxon[,col_list] %>% 
  filter(Class %in% graph_taxon) %>% filter(!is.na(depth))

y_tax <- malaspina_bins[,col_list]%>% 
  filter(Class %in% graph_taxon)

g_tax <- rbind(x_tax, y_tax)

g_tax$g_depth <- gsub('BAT', '         BAT', g_tax$depth)
g_tax$g_depth <- factor(g_tax$g_depth, levels = c("SRF", "DCM", "MES", '         BAT'))
g_tax$is_biofilm <- gsub('absent', 'Biofilm\nabsent', g_tax$is_biofilm)
g_tax$is_biofilm <- gsub('present', 'Biofilm\n     present', g_tax$is_biofilm)
g_tax$Class <- gsub('Gammaproteobacteria', '\u03b3-proteobacteria', g_tax$Class)
g_tax$Class <- factor(g_tax$Class, levels = graph_taxon1)

tax_depth <- g_tax %>%
  ggplot(aes(x=percent_trans, y = fct_rev(g_depth), fill = class_trans)) +
  scale_fill_manual(values = c('orange','gray', "green"))+
  facet_wrap(~Class, ncol = 3) +
  geom_boxplot(outlier.alpha = 0.3,outlier.shape = 21) + 
  ylab("depth")+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits=c(0, 0.8)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.1, y = 0)) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

ls_order <- c("planktonic", "mixed", "particle", "all")
g_tax_label <- c("High %-transposase classes (\u0251- \u03b2- \u03b3- proteobact. and Actinobact.)", 
                "Normal transposase abundance classes", 
                "Low %-transposase classes (Flavobact., Acidimicrobidae, etc.)")
g_tax_col <- c('orange','gray', "green")

# lab_text <- expression("Taxonomical "*italic("Class")*" of")

tax_lifestyle <- g_tax %>%
  ggplot(aes(x = percent_trans, y = size_fraction, fill = class_trans)) +
  facet_wrap(~Class, ncol = 3) +
  geom_boxplot(outlier.alpha = 0.3,outlier.shape = 21) +
  scale_fill_manual(guide = guide_legend(reverse = TRUE),
    labels = g_tax_label, values = g_tax_col)+
  # labs(fill=lab_text) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0.1, y = 0)) + 
  xlab("% of transposase ORFs in MAGs") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits=c(0, 0.8)) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.justification='left',
        legend.margin = margin(0),
        legend.position="bottom",
        legend.direction="vertical")


all_tara <- bin_taxon[,col_list] %>%
  filter(!is.na(depth))

all_mala <- malaspina_bins[,col_list] %>%
  filter(size_fraction != "error")

all <- rbind(all_tara, all_mala) %>%
  mutate(biofilm_bin=cut(percent_biofilm, breaks = c(-1,0.001,0.2,5)))%>%
  mutate(size_bin=cut(`complete genome size (Mbp)`, breaks = c(0,2,3,4,5,20)))

fact <- c("None", '<0-0.2%', '>0.2%')
all$biofilm_bin <- as.character(all$biofilm_bin)
all$biofilm_bin[all$biofilm_bin == '(-1,0.001]'] <- fact[1]
all$biofilm_bin[all$biofilm_bin == '(0.001,0.2]'] <- fact[2]
all$biofilm_bin[all$biofilm_bin == '(0.2,5]'] <- fact[3]
all$biofilm_bin <- factor(all$biofilm_bin, levels = fact)

all$size_bin <- as.character(all$size_bin)
all$size_bin[all$size_bin == '(0,2]'] <- '0-2'
all$size_bin[all$size_bin == '(2,3]'] <- '2-3'
all$size_bin[all$size_bin == '(3,4]'] <- '3-4'
all$size_bin[all$size_bin == '(4,5]'] <- '4-5'
all$size_bin[all$size_bin == '(5,20]'] <- '>5'
all$size_bin <- factor(all$size_bin, levels = c(">5","4-5","3-4","2-3","0-2"))

all$g_depth <- factor(all$depth, levels = c("SRF", "DCM", "MES", 'BAT'))


all <- all %>%
  group_by(depth) %>%
  mutate(countEachDepth= sum(n())) %>%
  group_by(size_fraction, .add=TRUE) %>%
  mutate(lifestyle_prop_to_depth=n()/countEachDepth)

particle_gsize <- all %>% 
  mutate(size_fraction = factor(size_fraction, levels = ls_order)) %>% 
  ggplot(aes(x=size_fraction, y = `complete genome size (Mbp)`)) +
  facet_wrap(~depth, ncol= 1) +
  xlab("MAG lifestyle") +
  ylab("MAG complete genome size (Mbp)") +
  coord_cartesian(ylim =c(0.7, 7.7)) +
  geom_boxplot(outlier.alpha = 0.1, varwidth = TRUE, aes(weight=sqrt(all$lifestyle_prop_to_depth)))+
  stat_summary(fun.data = boxplot.give.n, geom = "text",
               position=position_nudge(x = 0.3, y = 1.5)) +
  theme_classic() +
  theme()# strip.text.x = element_blank()

p_depth_gsize<-all %>% 
  mutate(depth = fct_rev(g_depth)) %>% 
  ggplot(aes(x = depth, y = `complete genome size (Mbp)`)) +
  ylim(0,9) +
  stat_compare_means(comparisons = 
                     list(c("BAT","DCM"),c("BAT","SRF")), 
                     label = "p.signif", hide.ns = TRUE, method = "t.test", 
                     tip.length = 0.01, label.y = c(7.2,7.8)) +
  scale_color_manual(labels = c("High %-transposase", "Normal %-transposase", "Low %-transposase"),
                     values = c('orange','gray', "green"))+
  geom_jitter(aes(color = class_trans), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) +
  ylab("MAG estimated complete genome size (Mbp)") +
  # stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.4)) +
  scale_x_discrete(labels=c("SRF" = "SRF (n = 339)", 
                            "DCM" = "DCM (n = 346)", 
                            "MES" = "MES (n = 497)",
                            "BAT" = "BAT (n = 255)")) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        legend.position="none")

size_trans<- ggplot(all, aes(x=fct_rev(size_bin), y=percent_trans)) + 
  geom_jitter(aes(color = class_trans), alpha = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_color_manual(labels = c("High %-transposase", "Normal %-transposase", "Low %-transposase"),
                     values = c('orange','gray', "green"))+
  stat_summary(fun.data = boxplot.give.n, geom = "text", 
               position=position_nudge(x = 0, y = 0.3)) + 
  facet_wrap(~depth, ncol= 1) +
  ylim(c(0,0.8)) +
  ylab("%-transposase ORF in MAGs") +
  xlab("complete genome size") +
  theme_classic() +
  theme(legend.position = "none")

pcom <- ggarrange(p_depth_gsize, tax_depth, tax_lifestyle,
                  nrow = 3, 
                  heights = c(4/14,3.5/14,5.7/14),
                  labels = c("A", "D", "E"))

ggarrange(pcom,particle_gsize,size_trans,
          ncol = 3, 
          widths = c(5/10,2.2/10,2.4/10),
          labels = c("", "B","C"))

ggsave("F3_MAG_transposase_genome_size_varwidth1.png", plot = last_plot())

library(Cairo)
cairo_pdf(filename = "F3_MAG_transposase_genome_size.pdf",
          width = 9,
          height = 6.5)

ggarrange(pcom,particle_gsize,size_trans,
          ncol = 3, 
          widths = c(5/10,2.2/10,2.4/10),
          labels = c("", "B","C"))
dev.off()
