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
cas_legend <- sprintf("Cassette ORFs (n = %s)\n", num_cassette)
nor_legend <- "The rest of Tara Oceans &\nMalaspina metagenomes"
lab_font <- c(rep(c("plain"),each=17), "bald", rep(c("plain"),each=3), "bald", "bald")
out_category <- c("~Total")

graphing %>% 
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
  theme(legend.position = c(0.7, 0.93),
        axis.title.y=element_blank(),
        # axis.text.y = element_text(face = lab_font),
        legend.background = element_rect(fill="transparent",color="transparent")) 
  # annotate(geom="text", x=0, y=-0.08, label="* TM: transport and metabolism", color="red")

ggsave("cassette_function_categories.png", plot = last_plot())
