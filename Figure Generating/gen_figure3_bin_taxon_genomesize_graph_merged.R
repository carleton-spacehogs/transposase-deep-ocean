setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins, low_trans, high_trans) %=% init_bins()

col_list <- c("complete genome size (Mbp)", "percent_trans", "log_percent_trans", 
              "percent_biofilm", "is_biofilm", "size_fraction",
              "depth", "Class", "class_trans", "graphing_log_trans", "graphing_log_biofilm")

all <- rbind(bin_taxon[,col_list], malaspina_bins[,col_list]) %>%
  mutate(biofilm_bin=cut(percent_biofilm, breaks = c(-1,0.001,0.2,5)))%>%
  mutate(size_bin=cut(`complete genome size (Mbp)`, breaks = c(0,2,3,4,5,20)))

all_taxa <- unique(all$Class)
normal_taxa <- setdiff(setdiff(all_taxa, low_trans), high_trans)
normal_taxa_rows <- all %>% filter(Class %in% normal_taxa)
count_nor_class <- as.data.frame.matrix(table(normal_taxa_rows$Class, normal_taxa_rows$depth))
count_nor_class[order(-count_nor_class$SRF),]
#                        SRF DCM MES BAT
# Opitutae                7   4   0   0
# Actinobacteria          5   3  11   9 
# Planctomycetia          5  11  10   0 # has the least worst distribution for each depth 
# Phycisphaerae           4  10   3   1
# Kiritimatiellae         3   0   0   1

high_taxa_rows <- all %>% filter(Class %in% high_trans)
table(high_taxa_rows$Class, high_taxa_rows$depth)
#                     SRF DCM MES BAT
# Alphaproteobacteria  68  46 120  60
# Betaproteobacteria    1   0  10   0
# Deltaproteobacteria   9   4  15   0
# Gammaproteobacteria  52  45 140  64


low_taxa_rows <- all %>% filter(Class %in% low_trans)
as.data.frame(table(low_taxa_rows[, "Class"]))
#             Var1 Freq
#  Acidimicrobidae   52
#    Flavobacteria  155 # use this
# Gemmatimonadetes   24

all <- all %>%
  filter(!is.na(depth))

graph_taxon <- c("Flavobacteria", "Planctomycetia", "Alphaproteobacteria")
graph_taxon1 <- c("Flavobacteria", "Planctomycetia", "\u0251-proteobacteria")

x_tax <- bin_taxon[,col_list] %>% 
  filter(Class %in% graph_taxon) %>% filter(!is.na(depth))

y_tax <- malaspina_bins[,col_list]%>% 
  filter(Class %in% graph_taxon)

g_tax <- rbind(x_tax, y_tax)

g_tax$g_depth <- gsub('BAT', '         BAT', g_tax$depth)
g_tax$g_depth <- factor(g_tax$g_depth, levels = c("SRF", "DCM", "MES", '         BAT'))
g_tax$is_biofilm <- gsub('absent', 'Biofilm\nabsent', g_tax$is_biofilm)
g_tax$is_biofilm <- gsub('present', 'Biofilm\n     present', g_tax$is_biofilm)
g_tax$Class <- gsub('Alphaproteobacteria', '\u0251-proteobacteria', g_tax$Class)
g_tax$Class <- factor(g_tax$Class, levels = graph_taxon1)

# high, mid, low
color_code <- c("steelblue","cyan","yellow")

tax_depth <- g_tax %>%
  ggplot(aes(x=percent_trans, y = fct_rev(g_depth), fill = class_trans)) +
  scale_fill_manual(values = color_code)+
  facet_wrap(~Class, ncol = 3) +
  geom_boxplot(outlier.alpha = 0.3,
               outlier.shape = 21,
               outlier.color = NA) + 
  ylab("depth")+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), limits=c(0, 0.8)) +
  stat_summary(fun.data = boxplot.give.n, 
               geom = "text", 
               position=position_nudge(x = 0.4, y = 0)) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

ls_order <- c("planktonic", "mixed", "particle", "all")
g_tax_label <- c("High %-transposase classes (\u0251- \u03b2- \u03b3- \u03b4- proteobact.)", 
                "Normal transposase abundance classes", 
                "Low %-transposase classes (Flavobact., Acidimicrobidae, etc.)")
g_tax_col <- color_code

# lab_text <- expression("Taxonomical "*italic("Class")*" of")

tax_lifestyle <- g_tax %>%
  ggplot(aes(x = percent_trans, y = fct_rev(size_fraction), fill = class_trans)) +
  facet_wrap(~Class, ncol = 3) +
  geom_boxplot(outlier.alpha = 0.3,
               outlier.shape = 21,
               outlier.color = NA) +
  scale_fill_manual(guide = guide_legend(reverse = TRUE),
    labels = g_tax_label, values = g_tax_col)+
  # labs(fill=lab_text) +
  stat_summary(fun.data = boxplot.give.n, 
               geom = "text", 
               position=position_nudge(x = 0.4, y = 0)) + 
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
                     values = color_code)+
  geom_jitter(aes(color = class_trans), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, 
               outlier.color = NA,
               alpha = 0.1) + # shape = NA, alpha = 0.1
  ylab("MAG estimated complete genome size (Mbp)") +
  # stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.4)) +
  scale_x_discrete(labels=c("SRF" = "SRF (n = 300)", 
                            "DCM" = "DCM (n = 303)", 
                            "MES" = "MES (n = 483)",
                            "BAT" = "BAT (n = 255)")) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        legend.position="none")

size_trans<- ggplot(all, aes(x=fct_rev(size_bin), y=percent_trans)) + 
  # geom_jitter(aes(color = class_trans), alpha = 0.8) +
  geom_boxplot(outlier.alpha = 0.2) + # .shape = NA, alpha = 0.2
  scale_color_manual(labels = c("High %-transposase",
                                "Normal %-transposase",
                                "Low %-transposase"),
                     values = color_code)+
  stat_summary(fun.data = boxplot.give.n, geom = "text", 
               position=position_nudge(x = 0, y = 0.3)) + 
  facet_wrap(~depth, ncol= 1) +
  coord_cartesian(ylim = c(0,0.8)) +
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

ggsave("F3_MAG_transposase_genome_size.png",
       plot = last_plot(),
       width = 9,
       height = 6.5)

library(Cairo)
cairo_pdf(filename = "F3_MAG_transposase_genome_size.pdf",
          width = 9,
          height = 6.5)

ggarrange(pcom,particle_gsize,size_trans,
          ncol = 3, 
          widths = c(5/10,2.2/10,2.4/10),
          labels = c("", "B","C"))
dev.off()
