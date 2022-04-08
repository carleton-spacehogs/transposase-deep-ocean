setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()

# from bin_taxon_genomesize_statistics.R
low_trans <- c("Flavobacteria","Acidimicrobidae","novelClass_E",
               "Gemmatimonadetes","SAR202-2","Marinisomatia")
high_trans <- c("Alphaproteobacteria","Gammaproteobacteria",
                "Betaproteobacteria","Actinobacteria")

col_list <- c("size_fraction","depth", "Class", "class_trans")

all_tara <- bin_taxon[,col_list] %>%
  filter(!is.na(depth))

all_mala <- malaspina_bins[,col_list] %>%
  filter(size_fraction != "error")

all <- rbind(all_tara, all_mala)

depth_trans<- all %>% count(class_trans, depth) 
depth_trans$depth <- factor(depth_trans$depth, levels=c("SRF","DCM","MES","BAT"))

ggplot(depth_trans, aes(fill=class_trans, y=fct_rev(depth), x=n)) + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  xlab("Proportion of MAGs") +
  ylab("Depth") +
  guides(fill=guide_legend(title="MAG Class category", reverse = TRUE)) +
  scale_fill_manual(labels = c("High %-transposase", "Normal %-transposase", "Low %-transposase"),
                     values = c('orange','gray', "green"))+
  theme_classic()

depth_size <- all %>% count(size_fraction, depth) 
depth_size$depth <- factor(depth_size$depth, levels=c("SRF","DCM","MES","BAT"))

ggplot(depth_size, aes(fill=size_fraction, y=fct_rev(depth), x=n))  + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  xlab("Proportion of MAGs") +
  ylab("Depth") +
  guides(fill=guide_legend(title="MAG Lifestyle")) +
  scale_fill_manual(labels = c("Planktonic", "Mixed", "Particle-\nassociated"),
                    values = c('green','gray', "orange"))+
  theme_classic()
