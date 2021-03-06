setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()

biofilm_RNA_scale <- c(-4.5, -4, -3.5, -3, -2.5, -2)
biofilm_RNA_percent_scale <- c("0.003", "0.01", "0.03", "0.1", "0.32", "1.00")
trans_RNA_scale <- c(-4.5, -4, -3.5, -3)
trans_RNA_percent_scale <- c("0.003", "0.010", "0.032", "0.100")

trans_scale <- c(-2.5, -3, -3.5, -4, -4.5)
trans_percent_scale <- c("0.316", "0.100", "0.032", "0.01", "0.003")
biofilm_percent_scale <- c(0.04, 0.02, 0.01)
biofilm_scale <- log10(biofilm_percent_scale/100)
CAZ_percent_scale <- c(0.8, 0.4, 0.2, 0.1)
CAZ_scale <- log10(CAZ_percent_scale/100)
mala_cov$size_fraction <- ifelse(mala_cov$lower_filter_size == "0.2", "planktonic", "particle")
DNA_tara$size_fraction <- ifelse(DNA_tara$upper_size_dna == "1.6", "planktonic", "particle")

bt_sel_col <- c("log_dna_trans", "log_dna_biofilm", "log_dna_defense",
                "percent_sect_CAZ","DNA_sect_CAZ","log_dna_sect_CAZ",
                "percent_sect_pep","DNA_sect_pep","log_dna_sect_pep",
                "Ocean_DNA", "Layer_DNA", "Oxygen_DNA", "size_fraction")
bt_to_graph <- rbind(mala_cov[,bt_sel_col], DNA_tara[,bt_sel_col])
bt_to_graph$Layer_DNA <- factor(bt_to_graph$Layer_DNA, levels = c("SRF","DCM","MES","BAT"))

RNA_tara %>%
  filter(Layer_RNA != "MIX") %>%
  ggplot(aes(x=percent_sect_CAZ, y=log_rna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) +
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Layer_RNA) +
  geom_point(aes(color = Ocean)) + 
  theme_classic() +
  ylab("% RNA reads mapped to transposase") +
  xlab("Percentage of secretory CAZenzyme") +
  geom_smooth(method = "lm", se = F)  +
  theme(legend.position = "none")

# RNA_tara$log_rna_sect_CAZ = log(RNA_tara$RNA_sect_CAZ)

RNA_tara %>%
  filter(Layer_RNA != "MIX") %>%
  ggplot(aes(x=log_rna_sect_CAZ, y=log_rna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) + 
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Layer_RNA) +
  geom_point(aes(color = Ocean)) + 
  theme_classic() +
  ylab("% RNA reads mapped to transposase") +
  xlab("(log 10 scale) % RNA reads mapped to secretory CAZenzyme") +
  geom_smooth(method = "lm", se = F)  +
  theme(legend.position = "none")

bt_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=percent_sect_CAZ, y=log_dna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) +
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Layer_DNA) + # , scales = "free"
  geom_point(aes(color = size_fraction)) + 
  theme_classic() +
  ylab("% DNA reads mapped to transposase") +
  xlab("secretory CAZenzyme among all CAZenzyme (%)") +
  geom_smooth(method = "lm", se = F)  +
  theme() # legend.position = "none"

bt_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=log_dna_sect_CAZ, y=log_dna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) +
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Layer_DNA) + # , scales = "free"
  geom_point(aes(color = size_fraction)) + 
  theme_classic() +
  ylab("% DNA reads mapped to transposase") +
  xlab("secretory CAZenzyme among all CAZenzyme (%)") +
  geom_smooth(method = "lm", se = F)  +
  theme() # legend.position = "none"

bt_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=percent_sect_pep, y=log_dna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) +
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Layer_DNA) + # , scales = "free"
  geom_point(aes(color = size_fraction)) + 
  theme_classic() +
  ylab("% DNA reads mapped to transposase") +
  xlab("secretory peptidase among all peptidase (%)") +
  geom_smooth(method = "lm", se = F)  +
  theme() # legend.position = "none"

bt_to_graph$avg_percent = (bt_to_graph$percent_sect_CAZ + bt_to_graph$percent_sect_pep)/2

bt_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=avg_percent, y=log_dna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) +
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Layer_DNA) + # , scales = "free"
  geom_point(aes(color = size_fraction)) + 
  theme_classic() +
  ylab("% DNA reads mapped to transposase") +
  xlab("average secretory (%)") +
  geom_smooth(method = "lm", se = F)  +
  theme() # legend.position = "none"




bt_DNA <- bt_to_graph %>%
  filter(Layer_DNA != "MIX") %>%
  ggplot(aes(x=log_dna_biofilm, y=log_dna_trans)) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) + 
  scale_x_continuous(breaks = biofilm_scale, labels = biofilm_percent_scale, limits=c(-4, -3.3)) +
  # scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale, limits=c(-4.5, -2)) + 
  # scale_x_continuous(breaks = biofilm_scale, labels = biofilm_percent_scale) +
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Layer_DNA) + # , scales = "free"
  geom_point(aes(color = Ocean_DNA)) + 
  theme_classic() +
  ylab("% DNA reads mapped to transposase") +
  xlab("% DNA reads mapped to biofilm") +
  geom_smooth(method = "lm", se = F)  +
  theme(legend.position = "none") # legend.position = "none"

ann_text <- data.frame(log_rna_biofilm = -3.83, log_rna_trans = -3.6,
            label = "No metatranscriptome\nsamples available",
            Layer_RNA = factor("BAT",levels = c("SRF", "DCM", "MES", "BAT")))

RNA_tara$Layer_RNA <- factor(RNA_tara$Layer_RNA, levels = c("SRF", "DCM", "MES", "BAT"))

bt_RNA <- RNA_tara %>% #
  filter(!is.na(Layer_RNA)) %>%
  ggplot(aes(x=log_rna_biofilm, y=log_rna_trans)) +
  theme_classic() +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale) + 
  scale_x_continuous(breaks = biofilm_RNA_scale, labels = biofilm_RNA_percent_scale) +
  facet_wrap(~Layer_RNA, nrow = 2, drop = FALSE) +
  geom_point(aes(color = Ocean)) + 
  scale_color_brewer(palette="Set1")+
  xlab("% RNA transcript mapped to biofilm") +
  ylab("% RNA transcript mapped to transposase") +
  geom_smooth(method = "lm", se = F) +
  geom_text(data = ann_text,label=ann_text$label)

ggarrange(bt_DNA, bt_RNA, labels = c("A", "B"), widths = c(0.42,0.6), ncol =2)

ggsave("F2_biofilm_transposase_correlation_RNA_DNA.png", plot = last_plot()) 
ggsave("F2_biofilm_transposase_correlation_RNA_DNA.pdf", plot = last_plot()) 

RNA_tara$Layer_RNA <- factor(RNA_tara$Layer_RNA, levels = c("SRF", "DCM", "MES", "BAT"))

biofilm_RNA_percent_scale2 <- c(0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
biofilm_RNA_scale2 <- log10(biofilm_RNA_percent_scale2/100)
RNA_tara %>%
  filter(Layer_RNA != "MIX") %>%
  mutate(Layer_RNA = fct_rev(Layer_RNA)) %>% 
  ggplot(aes(x=Layer_RNA, y=log_rna_biofilm))+ 
  scale_x_discrete(labels=c("SRF" = "SRF (n = 103)", 
                            "DCM" = "DCM (n = 51)", 
                            "MES" = "MES (n = 31)")) +
  geom_boxplot(outlier.colour = NA) + # geom_violin() +
  geom_jitter(aes(color=Ocean)) +
  # stat_summary(fun.data = boxplot.give.n, geom = "text") +
  stat_compare_means(label = "p.signif",
                     method = "t.test",
    comparisons = depth_comparison ) +
  scale_color_brewer(palette="Set1")+
  theme_classic() + 
  scale_y_continuous(breaks = biofilm_RNA_scale2, 
                     labels = biofilm_RNA_percent_scale2) + # limits=c(-4.5, -2)
  ylab("% RNA transcript mapped to biofilm") + 
  xlab("Depth") +
  # scale_x_discrete(drop=FALSE) + # position = "right"
  # annotate(geom="text", label="No metatranscriptome\nsamples available", y=-3.6, x=1) +
  coord_flip()

ggsave("S3_biofilm_depth_RNA.png", plot = last_plot())
ggsave("S3_biofilm_depth_RNA.svg", plot = last_plot())

summary(lm(log_dna_trans~Layer_DNA, DNA_tara))
trans_biofilm.lm <- lm(log_dna_trans~Layer_DNA + log_dna_biofilm, DNA_tara)
summary(trans_biofilm.lm)


summary(lm(log_rna_trans~Layer_RNA, RNA_tara))
RNA_trans_biofilm.lm <- lm(log_rna_trans~Layer_RNA + log_rna_biofilm, RNA_tara)
summary(RNA_trans_biofilm.lm)
anova(RNA_trans_biofilm.lm)

b <- RNA_tara[,c("log_rna_biofilm", "Layer_RNA")] %>% 
  mutate(is_MES = ifelse(Layer_RNA == "MES", "Y", "N"))
b.test <- t.test(log_rna_biofilm~is_MES, a)

1/10**b.test$conf.int[1:2]
# 95% CI: 3.126015 1.700280


RNA_tara %>%
  group_by(Layer_RNA) %>%
  summarise_at(vars(RNA_Biofilm), list(mean = mean, median = median))


DNA_DCM_trans <- DNA_tara %>%
  filter(Layer_DNA == "DCM") %>%
  select(c("DNA_Transposase"))

RNA_DCM_trans <- RNA_tara %>%
  filter(Layer_RNA == "DCM") %>%
  select(c("RNA_Trans"))

median(DNA_DCM_trans[[1]])
median(RNA_DCM_trans[[1]])

# trash

p2 <- RNA_tara %>% 
  filter(Layer_RNA != "MIX") %>%
  mutate(Layer_RNA = fct_rev(Layer_RNA)) %>% 
  ggplot(aes(x=log_rna_biofilm, y=Layer_RNA)) +
  geom_violin() +
  labs(title="    Biofilm transcript by depth",
       y = "% RNA transcript mapped to biofilm") +
  geom_jitter(aes(color=Ocean)) +
  scale_color_brewer(palette="Set1")+
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  theme_classic() + 
  scale_x_continuous(breaks = biofilm_RNA_scale, 
                     labels = biofilm_RNA_percent_scale) +
  scale_y_discrete(position = "right") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.position = "none")

print(p1)
print(p2, vp = viewport(width = 0.37, height = 0.45, x = 0.62, y = 0.262))
# ggsave doesn't work for print out graphs,
# save it manually as F2_biofilm_transposase_RNA.pdf


p1 <- bt_to_graph %>%
  ggplot(aes(x=log_dna_biofilm, y=log_dna_trans)) +
  facet_wrap(~Layer_DNA) +
  scale_y_continuous(breaks = trans_scale, labels = trans_percent_scale, limits=c(-4.5, -2)) + 
  scale_x_continuous(breaks = biofilm_scale, labels = biofilm_percent_scale) +
  geom_point(aes(color = Ocean_DNA)) + 
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


depth_loc_biofilm.lm <- lm(log_dna_trans~Layer_DNA+Ocean_short+log_dna_biofilm, DNA_tara)
anova(depth_loc_biofilm.lm)
get_r(depth_loc_biofilm.lm) - get_r(lm(log_dna_trans~Layer_DNA+Ocean_short, DNA_tara))
get_r(depth_loc_biofilm.lm)


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
  filter(!is.na(Depth_RNA)) %>%
  ggplot(aes(x=log_rna_biofilm, y=log_rna_trans)) +
  theme_classic() +
  scale_y_continuous(breaks = trans_RNA_scale, labels = trans_RNA_percent_scale) + 
  scale_x_continuous(breaks = biofilm_RNA_scale, labels = biofilm_RNA_percent_scale) +
  facet_wrap(~Layer_RNA) +
  geom_point(aes(color = Ocean)) +
  xlab("% RNA transcript mapped to biofilm") +
  ylab("% transposase RNA transcript") +
  geom_smooth(method = "lm", se = F) 



p2 <- DNA_RNA_tara %>%
  ggplot(aes(x=log_rna_biofilm, y=log_dna_trans)) +
  facet_wrap(~Layer_RNA) +
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

