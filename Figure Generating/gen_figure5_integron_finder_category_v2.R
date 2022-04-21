setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
summary <- init_integron_category()

# total 1225 cassette "proteins" from Malaspina, 7294 from Tara Oceans metagenomes.
# 359 cassette from Malaspina receive COG function calls, compared to 1951 from Tara.

# Malaspina total ORFs: 2438506, with COG annotation: 1657370 -> annotation rate: 0.67967
# South Atlatnic total: 8263808, COG annota: 4405933 -> rate: 0.5331
# North Atlantic total: 4718433, COG annota: 2759276 -> rate: 0.5847
# Mediterrean total: 3818905, COG: 2270527 -> rate: 0.5945
# Red sea total: 1988319, COG: 1145209 -> rate: 0.57597

graphing <- summary %>%
    mutate(integron_val_rank = integron_prop) %>%
    filter(COG_function != "total") %>%
    select(-c("integron_count", "normal_count")) %>%
    melt(id.vars=c("COG_function","integron_val_rank"), #  "ratio_prop", 
       variable.names = c("integron_prop", "normal_prop"))
graphing$type <- factor(graphing$variable,levels = c("normal_prop", "integron_prop"))

num_cassette <- as.character(summary[summary$COG_function == "~Total","integron_count"])
cas_legend <- sprintf("Cassette ORFs (n = %s)      ", num_cassette)
nor_legend <- expression(paste("all other ORFs in", italic("Tara"), " Oceans"))
out_category <- c("~Total")

p_cassette_category <- graphing %>% 
  filter(!COG_function %in% out_category) %>%
  ggplot(aes(fill=type, y=value, 
             x=reorder(COG_function, desc(integron_val_rank)))) + # SORT
  geom_bar(position="dodge", stat="identity") + 
  coord_flip(ylim = c(0, 0.15), clip = "off") +
  ylab("Proprotion to total ORF calls") + 
  # coord_cartesian(ylim = c(0, 0.15), clip = "off") +
  theme_classic() +
  guides(fill = guide_legend(reverse = TRUE, title = "COG gene-calls for:")) +
  scale_fill_manual(values=c("dark gray", "orange"), 
                    labels = c(nor_legend, cas_legend)) +
  annotate(geom="text", x=20.6, y=0.108, size = 3.3,
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
