setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
summary = init_integron_category()

# total 1225 cassette "proteins" from Malaspina, 7294 from Tara Oceans metagenomes.
# 359 cassette from Malaspina receive COG function calls, compared to 1951 from Tara.

# Malaspina total ORFs: 2438506, with COG annotation: 1657370 -> annotation rate: 0.67967
# South Atlatnic total: 8263808, COG annota: 4405933 -> rate: 0.5331
# North Atlantic total: 4718433, COG annota: 2759276 -> rate: 0.5847
# Mediterrean total: 3818905, COG: 2270527 -> rate: 0.5945
# Red sea total: 1988319, COG: 1145209 -> rate: 0.57597

# source, total_cassettes, functional_cassettes
# total            25178   9929 -> 39.4%
# ARS               2417   840
# CPC               2194   747
# EAC               1680   677
# IN                 687   271
# MED                907   442
# NAT               1679   861
# NP                2701   1051
# RS                 509   267
# SAT               4205   1757
# SP                5174   1110
# deep              3025   1206 -> 40.0%

# COG_integrons[COG_integrons.function == "nan"] # 16177 rows, 16177/25178 = 0.6425054
# COG_integrons['function'] = COG_integrons.function.str.replace(pattern, "Function unknown", regex=True)
# COG_integrons[COG_integrons.function == "Function unknown"] # 17439 rows, 17439/25178 = 0.69262848518

graphing <- summary %>%
  filter(integron_prop > 0.001) %>%
  arrange(integron_prop) %>%
  mutate(integron_val_rank = integron_prop) %>%
  select(-c("integron_count", "normal_count")) %>%
  reshape2::melt(id.vars=c("COG_function","integron_val_rank"),
                 variable.names = c("integron_prop", "normal_prop"))
graphing$type <- factor(graphing$variable,levels = c("normal_prop", "integron_prop"))

num_cassette <- as.character(sum(summary$integron_count))
cas_legend <- sprintf("Cassette ORFs (n = %s)      ", num_cassette)
nor_legend <- expression(paste("all other ORFs in", italic("Tara"), " Oceans"))

out_category <- c("Chromatin structure and dynamics") # because there is only 1 in integrons

integron_order = rev(arrange(summary, integron_prop)$COG_function)

# p_cassette_category <- 
graphing %>% 
  # filter(!COG_function %in% out_category) %>%
  mutate(COG_function = factor(COG_function, levels = integron_order)) %>%
  ggplot(aes(fill=type, y=value, x=COG_function)) + # SORT
  geom_bar(position="dodge", stat="identity") + 
  coord_flip(ylim = c(0, 0.15), clip = "off") +
  ylab("Proprotion to total ORF calls") + 
  theme_classic() +
  guides(fill = guide_legend(reverse = TRUE, title = "COG gene-calls for:")) +
  scale_fill_manual(values=c("dark gray", "orange"), 
                    labels = c(nor_legend, cas_legend)) +
  annotate(geom="text", x=18.6, y=0.108, size = 3.3,
           label="& Malaspina metagenomes") +
  theme(legend.position = c(0.65, 0.93),
        axis.title.y=element_blank(),
        legend.background = element_rect(fill="transparent",color="transparent")) 
  # annotate(geom="text", x=0, y=-0.08, label="* TM: transport and metabolism", color="red")

# p_bd sees gen_figure5_biofilm_defense_corr.R
ggarrange(p_cassette_category, p_bd, labels = c("A", "B"), 
          ncol = 2, nrow = 1, widths = c(0.7, 0.35))

ggsave("F5_cassette-COG-cateogory_biofilm-vs-defense.png", 
       plot = last_plot(),
       height = 8,
       width = 8)
ggsave("F5_cassette-COG-cateogory_biofilm-vs-defense.pdf", 
       plot = last_plot(),
       height = 8,
       width = 8)

total = sum(summary$integron_count) + sum(summary$normal_count)
integron_total = sum(summary$integron_count)

options(scipen = 0)
COG_enrichment.stat = apply(summary, 1, function(row){
  group1 = integron_total
  int_count = as.numeric(row["integron_count"])
  nor_count = as.numeric(row["normal_count"])
  group2 = int_count + nor_count
  # test for over-representation (enrichment)
  p.val = phyper(int_count-1, group2, total-group2, group1,lower.tail= FALSE)
  return(p.val)
})


summary.stats = summary %>%
  mutate(p.val = COG_enrichment.stat) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() 


