setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

graph_taxon <- as.factor(c("Verrucomicrobiae", "Flavobacteria", "Gammaproteobacteria"))
# graph_taxon <- factor(graph_taxon, levels= c("Acidimicrobidae", "Flavobacteria","Alphaproteobacteria", "Gammaproteobacteria"))

g_tax <- bin_taxon[,c("Class", "biofilm_count", "log_percent_trans", 
                      "depth", "percent_trans", "is_biofilm")] %>% 
  filter(Class %in% graph_taxon) %>% filter(!is.na(depth)) %>%
  mutate(depth = fct_rev(depth)) 
g_tax <- filter_outliers(g_tax, "percent_trans")

p_depth <- g_tax %>%
  ggplot(aes(x=percent_trans, y = depth)) +
  facet_wrap(~Class, ncol = 4) +
  geom_boxplot() + 
  # xlim(0, 0.8) +
  stat_summary(fun.data = boxplot.give.n, geom = "text",
               position=position_nudge(x = 0.1, y = 0)) + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
  # theme(strip.background =element_rect(fill=c("red", "blue", "blue"))) 

p_biofilm <- g_tax %>% 
  # filter(percent_trans < min(boxplot(bin_taxon_small$percent_trans)$out)) %>%
  # filter(biofilm_count < min(boxplot(bin_taxon_small$biofilm_count)$out)) %>%
  ggplot(aes(y = percent_trans, x = is_biofilm)) +
  facet_wrap(~Class, ncol = 4) +
  geom_boxplot() +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.1)) + 
  ylab("% of transposase ORFs in MAG") +
  # ylim(0, 0.8) +
  coord_flip() + 
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())

ggarrange(p_depth, p_biofilm, ncol = 1, nrow = 2)

# to_graph1 comes from bin_taxon.R
to_graph1  %>%
  ggplot(aes(y = mean_trans_percent, x = mean_biofilm_count)) +
  geom_point() + 
  # geom_smooth(method=lm, se = F)  +
  xlab("Mean biofilm-associated ORFs count") +
  ylab("Mean % of transposase ORFs") +
  xlim(-2,15) + ylim(-0.01,0.31) +
  theme_classic()+
  geom_text_repel(aes(label=`Class with more than 10 MAGs`), 
                  data = graph_df, box.padding = 0.6)



bin_taxon_small %>% 
  filter(Class == "Gammaproteobacteria") %>% 
  filter(is_biofilm == FALSE) %>% nrow()

bin_taxon_small %>% 
  filter(percent_trans < min(boxplot(bin_taxon_small$percent_trans)$out)) %>%
  filter(biofilm_count < min(boxplot(bin_taxon_small$biofilm_count)$out)) %>%
  ggplot(aes(y = percent_trans, x = percent_biofilm)) +
  facet_wrap(~Class) +
  geom_point(alpha = 0.15) +
  # geom_smooth(se = F) + 
  theme_classic()
  

summary(lm(percent_trans~biofilm_count, data = bin_taxon%>%
             filter(percent_trans < min(boxplot(bin_taxon_small$percent_trans)$out)) %>%
             filter(biofilm_count < min(boxplot(bin_taxon_small$biofilm_count)$out))))

table(bin_taxon$`Class with more than 10 MAGs`, bin_taxon$depth) %>% 
  as.data.frame.matrix() -> MAG_depth_count

table(bin_taxon$`Class with more than 10 MAGs`, bin_taxon$depth) %>%
  prop.table(margin = 1) %>%
  `*`(100) %>% round(2) %>%
  as.data.frame.matrix() -> temp
temp <- temp[order(temp$MES, decreasing = TRUE), ]

tran.lm <- lm(percent_trans~percent_biofilm, data = bin_taxon%>%filter(!is.na(Class)))
tran.lm1 <- update(tran.lm, .~.+Class)
tran.lm2 <- update(tran.lm, .~.*Class)
tran.lm3 <- update(tran.lm1, .~.+depth)
tran.lm4 <- update(tran.lm1, .~.*depth)

anova(tran.lm, tran.lm1)
anova(tran.lm1, tran.lm2)
anova(tran.lm1, tran.lm3)
anova(tran.lm3, tran.lm4)

anova(tran.lm3)

bin_taxon_class <- bin_taxon%>%filter(!is.na(Class))
tran.Classlast <- lm(percent_trans~biofilm_count+depth+Class, data = bin_taxon_class)
anova(tran.Classlast)

tran.biolast <- lm(percent_trans~depth+Class+percent_biofilm, data = bin_taxon_class)
anova(tran.biolast)

tran.depthlast <- lm(percent_trans~Class+biofilm_count+depth, data = bin_taxon_class)
anova(tran.depthlast)




