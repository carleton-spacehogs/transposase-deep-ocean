setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()
g(pn_ps_metagenome) %=% init_individual_metagenomes()

scale1 <- c(-2, -1.5, -1, -0.5, 0, 0.60206)
log_scale1 <- c("<0.01", "<0.032", "0.10", "0.32", "1", ">4")
# to_compare <- list(c("SRF", "DCM"), c("DCM", "MES"), c("MES", "Malaspina"))
to_compare <- list(c("SRF", "SRF_trans"), c("DCM", "DCM_trans"), c("MES", "MES_trans"), c("Malaspina", "Malaspina_trans"))

get_pnps <- function(in_type){
  return(pn_ps_metagenome %>% filter(gene_type==in_type) %>% filter(pnps<4))$log_pnps}

for (depth in c("SRF", "DCM", "MES", "Malaspina")){
  print(t.test(get_pnps(paste(depth, "_trans", sep = "")), get_pnps("SRF"))$p.value)
}

t.test(get_pnps("SRF"), get_pnps("DCM"))
t.test(get_pnps("DCM"), get_pnps("MES"))
t.test(get_pnps("MES"), get_pnps("Malaspina"))

stat.test <- pn_ps_metagenome %>%
  group_by(gene_type) %>%
  t_test(len ~ supp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

p1 <- pn_ps_metagenome %>%
  filter(pnps < 4) %>%
  ggplot(aes(x= gene_type2, y = log_pnps, fill = fct_rev(depth))) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) +
  scale_y_continuous(breaks = scale1, labels = log_scale1, 
                     name = "pN/pS", limits = c(-1.5, 1.3)) +
  theme_classic() +
  xlab("transposase        normal ORFs") +
  stat_summary(fun.data = boxplot.give.n, geom = "text", 
               position=position_dodge(.74)) +
  scale_fill_manual(breaks=c("SRF", "DCM","MES","Malaspina"), 
                    values=c("#73C2FB", "#008081", "#0080FF", "#1034A6")) + 
  coord_flip() + 
  geom_signif(
    y_position = c(0.12, 0.4), 
    xmin = c(2.1, 1.28), 
    xmax = c(2.3, 2.3),
    annotation = c("****", "****"), 
    tip_length = 0.005) +
  guides(fill=guide_legend(title="Depth")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = c(0.9, 0.45),
        legend.background = element_rect(fill = "transparent", 
                                         color = "transparent"))

bin_taxon$log_trans_tograph <- ifelse(bin_taxon$percent_trans == 0, -2, log10(bin_taxon$percent_trans))
scale <- c(-2, -1.5, -1, -0.5, 0)
log_scale <- c("0%\nno hits", "0.03%", "0.10%", "0.33%", "1.0%")
p2 <- bin_taxon %>% 
  filter(!is.na(depth))%>%
  # mutate(is_MES = ifelse(depth == "MES", "MES", "SRF&DCM")) %>%
  ggplot(aes(x = median_bin_pnps, y = log_trans_tograph)) + 
  geom_point(aes(color=fct_rev(depth)), alpha = 0.8) +
  theme_classic() +
  scale_x_log10() +
  scale_y_continuous(breaks = scale, labels = log_scale) +
  xlab("Median MAG pN/pS (Tara Oceans only)") +
  ylab("% of transposase ORF in MAG") + 
  scale_color_manual(values=c("#0080FF", "#008081", "#73C2FB")) +
  theme(legend.position="none")

ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(0.55,0.45),labels = c("A", "B"))

ggsave("p2.pdf", plot = last_plot())




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
bin_taxon$log_bin_pnps <- log(bin_taxon$median_bin_pnps)
summary(lm(`transposase gene calls in genome (%)`~log_bin_pnps, data = bin_taxon))

# deep ocean bins dont have higher pnps
bin_taxon <- bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")) %>%
  mutate(is_MES = ifelse(depth == "MES", "MES", "SRF&DCM"))
t.test(log_bin_pnps ~ is_MES, data = bin_taxon %>% filter(depth %in% c("SRF", "DCM", "MES")))

