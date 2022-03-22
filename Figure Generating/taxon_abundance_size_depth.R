mala_cov <- read_csv("data/Malaspina-transposase-coverage.csv")
mala_bin_cov <- read_csv("data/malaspina_bin_abundance_per_sample.csv")

mala_meta <- mala_cov[,c("sample1","lower_filter_size")]
mala_bin_cov <- merge(mala_meta, mala_bin_cov, by = "sample1")


setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()

bin_sel <- c("bin","depth","Class","prot","bact")
tara_bins <- bin_taxon[,bin_sel]
mala_bins <- malaspina_bins[,c("bin","depth","Class")]
mala_cov <- read_csv("data/Malaspina-Combined.RPKM.csv")[,c("bin","freeLiving_abundance","particle_abundance")]
mala_bins <- merge(mala_bins, mala_cov, by = "bin")

colnames(mala_bins) <- bin_sel

merged_bins <- rbind(tara_bins, mala_bins)

class_sel <- c("Actinobacteria", "Alphaproteobacteria","Betaproteobacteria","Gammaproteobacteria", # high transposase
            "Acidobacteria","Gemmatimonadetes", "Flavobacteria", "Marinisomatia") # low transposase abundance

# prot sum 3581, bact sum 7624
merged_bins%>% 
  filter(Class %in% class_sel) %>%
  mutate(bact = bact/7624) %>%
  mutate(prot = prot/3581) %>%
  group_by(depth) %>%
  summarise(across(where(is.numeric), list(sum = sum)))

merged_bins%>% 
  # filter(Class %in% class_sel) %>%
  # group_by(Class, depth) %>%
  summarise(across(where(is.numeric), list(sum = sum)))











