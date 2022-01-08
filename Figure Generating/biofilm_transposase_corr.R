setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()

trans_scale <- c(-2, -2.5, -3, -3.5, -4, -4.5)
trans_percent_scale <- c("1.00%", "0.316%", "0.100%", "0.032%", "0.010%", "0.003%")
biofilm_scale <- c(-3.4, -3.6, -3.8, -4.0)
biofilm_percent_scale <- c("0.040%", "0.025%", "0.016%", "0.01%")

p1 <- DNA_tara %>%
  filter(DNA_Biofilm > 0) %>%
  filter(DNA_Transposase > 0) %>%
  ggplot(aes(x=log_dna_biofilm, y=log_dna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale, limits=c(-4.5, -2)) + 
  scale_x_continuous(breaks = biofilm_scale, labels = biofilm_percent_scale) +
  geom_point(aes(color = Ocean_short)) + 
  theme_classic() +
  ylab("% transposase DNA reads") +
  xlab("% DNA reads mapped to biofilm") +
  geom_smooth(method = "lm", se = F)  +
  theme(legend.position = "none")

to_graph <- DNA_tara[,c("log_dna_biofilm", "Layer_DNA", "log_dna_trans")]%>%filter(Layer_DNA !="MIX")
to_graph$Layer_DNA <- factor(to_graph$Layer_DNA, levels = c("SRF", "DCM", "MES"))

ocean_names <- levels(DNA_tara$Ocean_short)
nine_colors <- scales::hue_pal()(9) # scales::show_col(scales::hue_pal()(9))
all_oceans <- DNA_tara$Ocean_short
color_DNA <- case_when(all_oceans==ocean_names[1]~nine_colors[1],
                       all_oceans==ocean_names[2]~nine_colors[2],
                       all_oceans==ocean_names[3]~nine_colors[3],
                       all_oceans==ocean_names[4]~nine_colors[4],
                       all_oceans==ocean_names[5]~nine_colors[5],
                       all_oceans==ocean_names[6]~nine_colors[6],
                       all_oceans==ocean_names[7]~nine_colors[7],
                       all_oceans==ocean_names[8]~nine_colors[8],
                       all_oceans==ocean_names[9]~nine_colors[9])

#color_DNA = ifelse(to_graph$Layer_DNA==ocean_names[1],nine_colors[1],
#                   ifelse(to_graph$Layer_DNA=="DCM", "steelblue", "sky blue"))

scatterplot3d(to_graph, color=color_DNA,
              angle = -30, pch = 20, #type = "h", lty.hplot = 2,
              grid=TRUE, box=FALSE)

biofilm_RNA_scale <- c(-4.5, -4, -3.5, -3, -2.5, -2)
biofilm_RNA_percent_scale <- c("0.003%", "0.01%", "0.032%", "0.10%", "0.316%", "1.00%")

p2 <- RNA_tara %>% 
  filter(Layer_RNA != "MIX") %>%
  mutate(Layer_RNA = fct_rev(Layer_RNA)) %>% 
  ggplot(aes(x=Layer_RNA, y=log_rna_biofilm)) +
  geom_violin() +
  geom_jitter(aes(color=Ocean)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  stat_compare_means(comparisons = depth_comparison, label = "p.signif") + 
  theme_classic() + 
  scale_y_continuous(breaks = biofilm_RNA_scale, 
                     labels = biofilm_RNA_percent_scale,
                     guide = guide_axis(n.dodge = 2)) +
  ylab("% RNA transcript mapped to biofilm") +
  theme(axis.title.y = element_blank()) +
  coord_flip()
  


# ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 2, nrow = 2, widths = c(0.4, 0.6))
ggarrange(p2, p1, labels = c("A", "B"), ncol = 2, nrow = 1, widths = c(0.55, 0.45))
ggarrange(p3, p4, labels = c("A", "B"), ncol = 2, nrow = 1)



# ggarrange(p3, p1, p2, labels = c("A", "B", "C"), widths = c(0.4, 0.3, 0.3), ncol = 3, nrow = 1)

depth_biofilm.lm <- lm(log_dna_trans~Layer_DNA+log_dna_biofilm, DNA_tara)
anova(depth_biofilm.lm)
get_r(depth_biofilm.lm) - get_r(lm(log_dna_trans~Layer_DNA, DNA_tara))

get_r(lm(log_dna_trans~log_dna_biofilm, DNA_tara))
summary(lm(log_rna_trans~log_rna_biofilm, RNA_tara))
summary(lm(log_dna_trans~log_rna_biofilm, DNA_RNA_tara))

# statistic learned from http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r
var.test(log_dna_trans~upper_size_dna, DNA_tara) # equal variance, F-test to test for homogeneity in variances
t.test(log_dna_trans~upper_size_dna, DNA_tara, var.equal = TRUE)
with(DNA_tara, shapiro.test(log_dna_trans[upper_size_dna == "1.6"]))# small p-value, not normal
with(DNA_tara, shapiro.test(log_dna_trans[upper_size_dna == "3"])) # assume normality

# supplement figure 3
# biofilm RNA + DNA scale, transposase DNA scale defined above 

trans_RNA_scale <- c(-4.5, -4, -3.5, -3)
trans_RNA_percent_scale <- c("0.003%", "0.01%", "0.032%", "0.10%")

p1 <- RNA_tara %>%
  ggplot(aes(x=log_rna_biofilm, y=log_rna_trans)) +
  theme_classic() +
  scale_y_continuous(breaks = trans_RNA_scale, labels = trans_RNA_percent_scale) + 
  scale_x_continuous(breaks = biofilm_RNA_scale, labels = biofilm_RNA_percent_scale) +
  geom_point(aes(color = Ocean)) +
  xlab("% RNA transcript mapped to biofilm") +
  ylab("% transposase RNA transcript") +
  geom_smooth(method = "lm", se = F) 

p2 <- DNA_RNA_tara %>%
  ggplot(aes(x=log_rna_biofilm, y=log_dna_trans)) +
  geom_point(aes(color = Ocean_short)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) + 
  scale_x_continuous(breaks = biofilm_RNA_scale, labels = biofilm_RNA_percent_scale) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("% transposase DNA reads") +
  xlab("% RNA transcript mapped to biofilm") +
  geom_smooth(method = "lm", se = F)

ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1, widths = c(0.57, 0.43))




stat.test1 <- compare_means(
  log_rna_biofilm ~ Layer_RNA, data = RNA_tara %>%
    filter(Layer_RNA %in% c("SRF", "DCM", "MES")),
  method = "t.test"
)
stat.test1 <- stat.test1 %>% mutate(y.position = c(0.0009, 0.0011, 0.0015))
#stat.test1$.y. <- c("RNA_Biofilm", "RNA_Biofilm", "RNA_Biofilm")


# trashed
pn <-RNA_tara %>% 
  filter(Layer_RNA != "MIX") %>%
  ggplot(aes(x=Layer_RNA, y=log_rna_biofilm)) +
  geom_violin() +
  geom_jitter(aes(color=Ocean)) + xlab("") +
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  stat_compare_means(comparisons = depth_comparison) + 
  theme_bw()

