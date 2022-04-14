setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()

bin_taxon$prop_metabolism <- bin_taxon$total_metabolic_ORFs/bin_taxon$`Total ORFs`

bin_taxon %>% 
  # mutate(size_fraction = factor(size_fraction, levels = ls_order)) %>% 
  ggplot(aes(x=size_fraction, y = prop_metabolism)) +
  # facet_wrap(~depth, ncol= 1) +
  #xlab("MAG lifestyle") +
  # coord_cartesian(ylim =c(0.7, 7.7)) +
#  stat_summary(fun.data = boxplot.give.n, geom = "text",
#               position=position_nudge(x = 0, y = 1.5)) +
  geom_boxplot(outlier.alpha = 0.1)+
  theme_classic() +
  theme()# strip.text.x = element_blank()

bin_taxon %>% 
  ggplot(aes(x=size_fraction, y = carbohydrate)) +
  geom_boxplot(outlier.alpha = 0.1,notch=TRUE)+
  theme_classic() +
  theme()# strip.text.x = element_blank()

bin_taxon %>% 
  ggplot(aes(x=size_fraction, y = shannon_base2)) +
  geom_boxplot(outlier.alpha = 0.1)+
  ylim(c(2.4, 2.7)) +
  theme_classic() +
  theme()# strip.text.x = element_blank()

t.test(simpson~size_fraction, 
       bin_taxon%>%filter(size_fraction %in% c("particle", "planktonic")))
t.test(shannon_base2~size_fraction, 
       bin_taxon%>%filter(size_fraction %in% c("particle", "planktonic")))



library(fmsb)

metabolic_names <- c("energy production and conversion",
                     "amino acid", "nucleotide", "carbohydrate",
                     "coenzyme", "lipid", "inorganic ion")
log_prop_names <- c("log_prop_energy","log_prop_amino", "log_prop_nucleotide", "log_prop_lipid",
                    "log_prop_coenzyme", "log_prop_carbohydrate", "log_prop_inorganicion")
c = 0.0001
metabolic_analy <- bin_taxon %>% select(c("size_fraction", metabolic_names, "Total ORFs"))

for (i in 1:7) {
  metabolic_analy <- cbind(metabolic_analy, 
    100* metabolic_analy[,metabolic_names[i]]/metabolic_analy[,"Total ORFs"] + c)
}

colnames(metabolic_analy) <- c("size_fraction", metabolic_names, "Total ORFs", log_prop_names)

metabolic_summary <- metabolic_analy %>%
  group_by(size_fraction) %>%
  summarise(# log_prop_energy = quantile(log_prop_energy, c(0.25, 0.5, 0.75), na.rm = TRUE),
    log_prop_nucleotide = quantile(log_prop_nucleotide, c(0.25, 0.5, 0.75), na.rm = TRUE),        
    log_prop_amino = quantile(log_prop_amino, c(0.25, 0.5, 0.75), na.rm = TRUE),
    log_prop_lipid = quantile(log_prop_lipid, c(0.25, 0.5, 0.75), na.rm = TRUE),  
    log_prop_coenzyme = quantile(log_prop_coenzyme, c(0.25, 0.5, 0.75), na.rm = TRUE),
    log_prop_carbohydrate = quantile(log_prop_carbohydrate, c(0.25, 0.5, 0.75), na.rm = TRUE),    
    log_prop_inorganicion= quantile(log_prop_inorganicion, c(0.25, 0.5, 0.75), na.rm = TRUE),
            prob = c(0.25, 0.5, 0.75)) %>%
  as.data.frame()


log_prop_names <- c("log_prop_nucleotide", "log_prop_amino", "log_prop_carbohydrate",
                    "log_prop_coenzyme", "log_prop_lipid", "log_prop_inorganicion")

particle_sum <- metabolic_summary%>%
  filter(size_fraction == "particle") %>%
  select(log_prop_names)
particle_sum <- rbind(rep(10,6) , rep(1,6), particle_sum)

plankton_sum <- metabolic_summary%>%
  filter(size_fraction == "planktonic") %>%
  select(log_prop_names)
plankton_sum <- rbind(rep(10,6) , rep(1,6), plankton_sum)

graph_name <- c("nucleotide","amino", "lipid", 
                "coenzyme", "carbohydrate", "inorganic\nion")

colnames(particle_sum) <- graph_name
colnames(plankton_sum) <- graph_name
par(mfrow=c(1,2))
create_beautiful_radarchart(particle_sum, title = "particle", caxislabels = c("2%", "4%", "6%", "8%", "10%"))
create_beautiful_radarchart(plankton_sum, title = "planktonic", caxislabels = c("2%", "4%", "6%", "8%", "10%"))










metabolic_names <- c("energy production and conversion",
                     "amino acid", "nucleotide", "carbohydrate",
                     "coenzyme", "lipid", "inorganic ion")

metabolic_summary <- metabolic_analy %>%
  group_by(size_fraction) %>%
  summarise(
    energy = quantile(`energy production and conversion`, c(0.25, 0.5, 0.75), na.rm = TRUE),
    nucleotide = quantile(nucleotide, c(0.25, 0.5, 0.75), na.rm = TRUE),        
    `amino acid` = quantile(`amino acid`, c(0.25, 0.5, 0.75), na.rm = TRUE),
    lipid = quantile(lipid, c(0.25, 0.5, 0.75), na.rm = TRUE),  
    coenzyme = quantile(coenzyme, c(0.25, 0.5, 0.75), na.rm = TRUE),
    carbohydrate = quantile(carbohydrate, c(0.25, 0.5, 0.75), na.rm = TRUE),    
    `inorganic ion`= quantile(`inorganic ion`, c(0.25, 0.5, 0.75), na.rm = TRUE),
    prob = c(0.25, 0.5, 0.75)) %>%
  as.data.frame()

graph_name <- c("energy",
                   "amino acid", "nucleotide", "carbohydrate",
                   "coenzyme", "lipid", "inorganic ion")

particle_sum <- metabolic_summary%>%
  filter(size_fraction == "particle") %>%
  select(graph_name)
particle_sum <- rbind(rep(300,7) , rep(1,7), particle_sum)

plankton_sum <- metabolic_summary%>%
  filter(size_fraction == "planktonic") %>%
  select(graph_name)
plankton_sum <- rbind(rep(300,7) , rep(1,7), plankton_sum)

par(mfrow=c(1,2))
create_beautiful_radarchart(particle_sum, title = "particle")
create_beautiful_radarchart(plankton_sum, title = "planktonic")










particle <- metabolic_analy%>%
  filter(size_fraction == "particle") %>%
  select(log_prop_names)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
particle <- rbind(rep(-1,10) , rep(-4,10), particle)

planktonic <- metabolic_analy%>%
  filter(size_fraction == "planktonic") %>%
  select(log_prop_names)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
planktonic <- rbind(rep(-1,10) , rep(-4,10), planktonic)

# The default radar chart 


colnames(planktonic) <- graph_name
colnames(particle) <- graph_name
par(mfrow=c(1,2))
radarchart(planktonic)
radarchart(particle)

create_beautiful_radarchart <- function(data, color1 = "#00AFBB", 
                                        color2 = "#FFFFFF", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color1, pfcol = scales::alpha(color1, 0.3), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
create_beautiful_radarchart(particle)
create_beautiful_radarchart(planktonic)



