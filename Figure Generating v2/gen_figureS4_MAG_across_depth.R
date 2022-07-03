setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins, low_trans, high_trans) %=% init_bins()

col_list <- c("size_fraction","depth", "Class", "class_trans")

all_tara <- bin_taxon[,col_list] %>%
  filter(!is.na(depth))

all_mala <- malaspina_bins[,col_list] %>%
  filter(size_fraction != "error")

all <- rbind(all_tara, all_mala)

# labels = c("Planktonic", "Mixed", "Particle-associated"),
# values = c('green','gray', "orange"

depth_size <- all %>% count(size_fraction, depth) 
depth_size$depth <- factor(depth_size$depth, levels=c("SRF","DCM","MES","BAT"))
depth_size$size_fraction <- factor(depth_size$size_fraction, 
                                   levels = c("particle", "mixed", "planktonic"))

MAG_l <- ggplot(depth_size, aes(fill=size_fraction, y=fct_rev(depth), x=n))  + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  xlab("Proportion of MAGs") +
  ylab("Depth") +
  scale_fill_manual(labels = rev(c("Planktonic", "Mixed", "Particle-associated")),
                    values = rev(c('green','gray', "orange")))+
  guides(fill=guide_legend(title="MAG Lifestyle", reverse = TRUE)) + #, reverse = TRUE
  theme_classic() +
  theme(legend.justification='left',
        legend.margin = margin(0),
        legend.position="top",
        legend.direction="vertical")

# ggsave("S4_MAG_lifestyle_depth.png", plot = last_plot())

depth_trans<- all %>% count(class_trans, depth) 
depth_trans$depth <- factor(depth_trans$depth, levels=c("SRF","DCM","MES","BAT"))

# g_tax_label <- c("High %-transposase classes (\u0251- \u03b2- \u03b3- proteobact. and Actinobact.)", 
#                  "Normal transposase abundance classes", 
#                  "Low %-transposase classes (Flavobact., Acidimicrobidae, etc.)")

g_tax_label <- c("High transposase abundance", 
                 "Average transposase abundance",
                 "Low transposase abundance")

# g_tax_col <- rev(c("light blue","sky blue","steelblue"))
color_code <- c("steelblue","cyan","yellow")

MAG_r <- ggplot(depth_trans, aes(fill=class_trans, y=fct_rev(depth), x=n)) + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  xlab("Proportion of MAGs") +
  ylab("Depth") +
  guides(fill=guide_legend(title="MAG taxon class category", reverse = TRUE)) +
  scale_fill_manual(guide = guide_legend(reverse = TRUE),
                    labels = g_tax_label, values = g_tax_col)+
  theme_classic() +
  theme(axis.title.y = element_blank(),
        # legend.title = element_blank(),
        legend.justification='left',
        legend.margin = margin(0),
        legend.position="top",
        legend.direction="vertical")

ggarrange(MAG_l, MAG_r, labels = c("A", "B"), 
          ncol = 2, nrow = 1, widths = c(0.5, 0.5))

# ggsave("S4_MAG_prop_depth_combine.svg", plot = last_plot())

ggsave("S4_MAG_prop_depth_combine.png", 
       plot = last_plot(),
       height = 3,
       width = 8)

