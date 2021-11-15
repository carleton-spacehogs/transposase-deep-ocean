setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(pn_ps_integron, pn_ps_bins, pn_ps_total, tara_integron_summary, malaspina_integron_summary, tara_integron_all, malaspina_integron_all) %=% init_integron()

make_summary_data_graphable <- function(integron_summary, reference_sample){
  summary_graph <- integron_summary %>%
    filter(COG_function != "total") %>%
    select(-c("integron_count", paste(reference_sample, "_count", sep = ""))) %>%
    melt(id.vars=c("COG_function", "ratio_prop"),
       variable.names = c("integron_prop", paste(reference_sample, "_prop", sep = "")))
  summary_graph$variable <- factor(summary_graph$variable, levels = rev(levels(summary_graph$variable)))
  return(summary_graph)
}

make_summary_data_graphable(tara_integron_summary, "normal") %>%
  filter(!COG_function %in% c("~Total", "Function unknown", 
                              "General function prediction only")) %>%
  ggplot(aes(fill=variable, y=value, x=reorder(COG_function, ratio_prop))) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() +
  #xlab("COG genes function categories") + 
  xlab("") +
  ylab("Proprotion to total ORF calls") + 
  theme_classic() +
  guides(fill = guide_legend(reverse = TRUE, title = "COG calls origin:")) +
  scale_fill_manual(values=c("dark gray", "orange"),
                    labels = c("all Tara Oceans\nmetagenomes (21 Millon)", 
                               "cassette ORFs (n = 1636)")) +
  theme(legend.position = c(0.75, 0.68),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

# same things above for deep ocean
make_summary_data_graphable(malaspina_integron_summary, "deep") %>%
  filter(!COG_function %in% c("~Total", "Function unknown", 
                              "General function prediction only")) %>%
  ggplot(aes(fill=variable, y=value, x=reorder(COG_function, ratio_prop))) + 
  geom_bar(position="dodge", stat="identity") + 
  coord_flip() +
  xlab("") + ylim(c(0, 0.15)) +
  theme_classic() +
  ylab("Proprotion to total gene calls") + 
  guides(fill = guide_legend(reverse = TRUE, title = "COG calls from:")) + 
  scale_fill_discrete(labels = c("the entire Malaspina\nmetagenome (2 Million)", "cassette genes (n = 338)")) +
  theme(legend.position = c(0.75, 0.45),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


tara_integron_all %>% filter(category=="Defense mechanisms") %>%
  count()
tara_integron_all %>% 
  filter(category=="Defense mechanisms") %>%
  filter(grepl("toxin",gene_function, ignore.case = TRUE)) %>%
  count()

pn_ps_integron %>% filter(category!="Defense mechanisms") %>%
  filter(pNpS < 10) %>%
  select("pNpS") %>% summary()

pn_ps_integron %>% 
  filter(pNpS < 10) %>%
  select("pNpS") %>% count()

t.test(pNpS~is_defense, data = pn_ps_integron %>% filter(pNpS < 10))

# from https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
t.test.from.summary.data <- function(mean1, sd1, n1, mean2, sd2, n2, ...) {
  data1 <- scale(1:n1)*sd1 + mean1
  data2 <- scale(1:n2)*sd2 + mean2
  t.test(data1, data2, ...)
}

defesen_pnps <- pn_ps_integron %>% filter(category=="Defense mechanisms") %>%
  filter(pNpS < 10) %>% select("pNpS")

m_x <- 0.18505400556663568
m_y <- mean(defesen_pnps$pNpS)
s_x <- 0.3826944115726565
s_y <- sd(defesen_pnps$pNpS)
t.test.from.summary.data(m_x, s_x, 3000, m_y, s_y, nrow(defesen_pnps))


pn_ps_integron %>% 
  filter(pNpS < 10) %>%
  filter(from_integron_type != "In0") %>%
  ggplot(aes(y=pNpS, x=is_defense)) +
  geom_boxplot() + 
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  stat_compare_means(comparisons = list(c("defense mechanisms", "non-defense"))) +
  coord_flip(ylim=c(0,7.5))
  
pn_ps_integron %>% 
  filter(pnps < 10) %>%
  ggplot(aes(y=log_pnps, x=category)) +
  geom_boxplot() + 
  stat_summary(fun.data = boxplot.give.n, geom = "text") +
  coord_flip()

pn_ps_integron <- pn_ps_integron %>% filter(!is.na(category)) %>%
  mutate(is_metabolism = ifelse(category %in% c("Lipid transport and metabolism","Nucleotide transport and metabolism"), "lipid", "non-lipid"))
  # mutate(is_metabolism = ifelse(category %like% "metabolism", "metabolism", "non-metabolism"))
  
t.test(log_pnps~ is_metabolism, data = pn_ps_integron %>% 
         filter(pnps < 10) %>% filter(!is.na(category)))

tmp <- pn_ps_integron%>%group_by(category) %>%
  summarize(count = n(), median=median(pnps), mean=mean(pnps))
tmp
