setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(pn_ps_integron, pn_ps_bins, pn_ps_total, tara_integron_summary, malaspina_integron_summary, tara_integron_all, malaspina_integron_all) %=% init_integron()
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()
g(trans_in_tara_bins, trans_in_mala_bins) %=% init_transposase_in_bins()

scale1 <- c(-2, -1.5, -1, -0.5, 0, 0.5)
log_scale1 <- c("<0.01", "0.03", "0.10", "0.32", "1", "3.16")
graph_col <- c("bin","log_median_bin_pnps", "is_tara")
rbind(malaspina_bins %>% select(graph_col),
      bin_taxon %>% filter(is_tara != "unsure") %>% select(graph_col)) %>% 
  ggplot(aes(x= is_tara, y = log_median_bin_pnps, fill=is_tara)) +
  geom_boxplot(alpha = 0.5, width=.5, outlier.colour = NA) +
  theme_classic() +
  scale_y_continuous(breaks = scale1, labels = log_scale1, name = "pN/pS", limits = c(-1.5,0)) +
  stat_compare_means(comparisons = list( c("SRF&DCM", "MES"), c("Deep_Malaspina", "MES")),
                    tip.length = 0.01, label.y = c(-0.1,0.7)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = -0.7)) +
  # scale_fill_manual(breaks=c("SRF&DCM","MES","Deep_Malaspina"), values=c("sky blue", "steelblue", "blue4"))+ 
  theme(legend.position="none") +
  coord_flip() # 
  
  

tara_bin_trans <- trans_in_tara_bins %>% select(c("pNpS"))%>% filter(pNpS < 10)
tara_bin_trans$is_MES <- "Tara_trans"

mala_bin_trans <- trans_in_mala_bins %>% select(c("pNpS")) %>% filter(pNpS < 10)
mala_bin_trans$is_MES <- "Mala_trans"

colnames(mala_bin_trans) <- c("median_bin_pnps", "is_MES")
colnames(tara_bin_trans) <- c("median_bin_pnps", "is_MES")

malaspina <- malaspina_bins %>% select(c("median_bin_pnps", "is_MES")) %>% filter(median_bin_pnps < 10)
to_graph <- rbind(mala_bin_trans,
                  tara_bin_trans,
                  malaspina,
                  bin_taxon %>% 
                    filter(depth %in% c("SRF", "DCM", "MES")) %>% 
                    select(c("is_MES", "median_bin_pnps")))
to_graph$is_MES <- factor(to_graph$is_MES, levels = c("Mala_trans","Deep_Malaspina",
  "MES","SRF&DCM","Tara_trans"))

to_graph$log_pnps <- log10(to_graph$median_bin_pnps)
tara_median <- log10(0.1066)
mala_median <- log10(0.455)

scale1 <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1)
log_scale1 <- c("<0.01", "0.03", "0.10", "0.32", "1", "3.16", "10")
p1 <- to_graph %>% 
  ggplot(aes(x = is_MES, y = log_pnps, fill=is_MES)) +
  scale_x_discrete(labels=c("Tara_trans" = "Transposases\nin MAGs", 
                            "SRF&DCM" = "SRF&DCM\nMAG medians",
                            "MES" = "MES\nMAG medians",
                            "Mala_trans" = "Transposases\nin MAGs",
                            "Deep_Malaspina" = "Malaspina\nMAG medians")) +
  geom_segment(aes(x=0.5,xend=2.4,y=mala_median,yend=mala_median),size=1, linetype=1, color = "blue") +
  geom_segment(aes(x=2.5,xend=5.5,y=tara_median,yend=tara_median),size=1, linetype=1, color = "red") +
  annotate(geom="label",x = 2.5, y = mala_median,label.size = NA, label = "Malaspina median: 0.455", fill="white", color = "blue") +
  annotate(geom="label",x = 5.5, y = tara_median,label.size = NA, label = "Tara Oceans median: 0.1066", fill="white", color = "red") +
  theme_classic() +
  geom_boxplot(alpha = 0.5, width=.5, outlier.colour = NA) + 
  stat_compare_means(comparisons = list( c("SRF&DCM", "MES"), 
                                         c("Deep_Malaspina", "MES"),
                                         c("MES", "Tara_trans"),
                                         c("Deep_Malaspina", "Mala_trans")),
                     label = "p.signif", hide.ns = TRUE, tip.length = 0.01,
                     label.y = c(-0.1,0.7, 0.9, 0.95)) +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = -0.7)) +
  scale_y_continuous(breaks = scale1, labels = log_scale1, name = "pN/pS", limits = c(-2, 1)) +
  xlab("Malaspina          Tara Oceans")+
  scale_fill_manual(breaks=c("Tara_trans", "SRF&DCM","MES","Deep_Malaspina", "Mala_trans"), 
                     values=c("orange", "sky blue", "steelblue", "blue4", "orange"))+ 
  theme(legend.position="none")+
  coord_flip() # 

to_graph1 <- bin_taxon %>%
  select("transposase gene calls in genome (%)", "median_bin_pnps", "depth") %>%
  mutate(`transposase gene calls in genome (%)`, `transposase gene calls in genome (%)` = ifelse(
    `transposase gene calls in genome (%)` == 0, 0.01, `transposase gene calls in genome (%)`))
to_graph1$`transposase gene calls in genome (%)` = log10(to_graph1$`transposase gene calls in genome (%)`)

scale <- c(-2, -1.5, -1, -0.5, 0)
log_scale <- c("0%\nno hits", "0.03%", "0.10%", "0.33%", "1.0%")
p2 <- to_graph1 %>% 
  ggplot(aes(x = median_bin_pnps, y = `transposase gene calls in genome (%)`)) + 
  geom_point(alpha=0.3, stroke = 0, color = "red") +
  theme_classic() +
  scale_x_log10() +
  scale_y_continuous(breaks = scale, labels = log_scale) +
  xlab("Median MAG pN/pS (Tara Oceans only)") +
  ylab("% of transposase ORF in MAG")
  # geom_point(aes(color = depth)) +
  # labs(color = "Depth")+
  # annotate(geom="text", x=0.8, y=-1.75, label="Tara Oceans\nMAGs only",color="red") +
  # scale_color_manual(breaks=c("SRF","DCM","MES"), values=c("sky blue", "steelblue", "blue"))

ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(0.6,0.4),labels = c("A", "  B"))

# last supplementary figure
bin_taxon %>%
  ggplot(aes(x=median_bin_pnps, y=`complete genome size (Mbp)`)) +
  # geom_point(aes(color=depth)) +
  geom_smooth(method = lm) + 
  theme_classic()+
  labs(color = "Depth", x ="Median pN/pS of each MAG")+
  scale_x_continuous(trans = 'log10')

# statistic analysis
# no relation between bin pnps and transposase
summary(lm(`transposase gene calls in genome (%)`~log_bin_pnps, data = bin_taxon))

# deep ocean bins have higher pnps
t.test(log_bin_pnps ~ is_MES, data = bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")))
t.test(log_pnps ~ is_MES, data = to_graph %>% filter(is_MES %in% c("MES", "Deep_Malaspina")))


# transposases in bins have higher pnps than bin medians
t.test(
  (transposase_in_bins %>% filter(pNpS < 10))$pNpS,
  (bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")))$log_bin_pnps)

# deep ocean bins have higher pnps
summary(lm(log_bin_pnps ~ `complete genome size (Mbp)`, data = bin_taxon))

bin_taxon %>% filter(`complete genome size (Mbp)` < 1) %>% 
  select("ocean", "depth", "complete genome size (Mbp)")

bin_taxon %>% filter(!is.na(median_bin_pnps)) %>% count()
