init_env <- function(){
  library(readr);
  library(tidyverse)
  library(ggpubr)
  library(purrr)
  library(ggpubr)
  library(ggrepel)
  library("ggplot2")
  library("ggpmisc")
  library("readxl")
  library("car")
  library("forcats")
  library(grid)
  library(gridExtra)
  library("GGally") 
  library(dplyr)
  library(data.table)
  library(gridExtra)
  library(png)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
  getwd()
}

boxplot.give.n <- function(x){ return(c(y = mean(x), label = length(x)))}


# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# code taken from https://stackoverflow.com/questions/7519790/assign-multiple-new-variables-on-lhs-in-a-single-line
# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

init_bins <- function(){
  init_env()
  RPKM <- read_excel("data/TOBG-Combined.RPKM.xlsx")
  virus <- read_csv("data/viral_contig_each_bin.csv", col_names = c("bin", "viral_contig_count"))
  taxon <- read_csv("data/bin_taxon2.csv")
  taxon <- merge(taxon, RPKM[c("bin", "Max_RPKM", "sum")])
  # HTG <- read_csv("data/HGT_genes_each_bin.csv")
  # taxon <- merge(taxon, HTG, by="bin", all.x = TRUE)
  # taxon$HTG_proportion <- taxon$HTG/taxon$`Total ORFs`
  taxon$depth <- factor(taxon$depth, levels = c("SRF", "DCM", "mixed", "epi", "MES"))
  taxon$`transposase gene calls in genome (%)` <- taxon$`Percent Transposases`*100
  taxon$`biofilm gene calls in genome (%)` <- taxon$`percent_biofilm`*100
  taxon$`defense mechanisms gene calls in genome (%)` <- taxon$`percent_defense`*100
  taxon$log_complete_bin_size <- log(taxon$`complete genome size (Mbp)`)
  taxon$complete_ORF_count <- taxon$`Total ORFs`/taxon$`Percent Complete`*100
  
  taxon_count_greater_than_10 <- taxon %>% 
    #filter(depth %in% c("SRF", "DCM", "MES")) %>%
    group_by(Class) %>%  filter(n() > 10) %>%
    select(Class) %>% unique()
  
  order_count_greater_than_5 <- taxon %>% 
    group_by(Order) %>% filter(n() > 5) %>%
    select(Order) %>% unique()
  
  taxon <- taxon %>% mutate(`Class with more than 10 MAGs` = 
                             ifelse(Class %in% c(taxon_count_greater_than_10)$Class, Class, "Others Or Unknown"))
  taxon <- taxon %>% mutate(`Class with more than 10 MAGs` = 
                              ifelse(!is.na(`Class with more than 10 MAGs`), `Class with more than 10 MAGs`, "Others Or Unknown"))
  taxon <- taxon %>% mutate(`Order>5` = 
                              ifelse(Order %in% c(order_count_greater_than_5)$Order, Order, "Others Or Unknown"))
  taxon <- taxon %>% mutate(`Order>5` = 
                              ifelse(!is.na(`Order>5`), `Order>5`, "Others Or Unknown"))
  
  rm(virus, RPKM)
  pn_ps_bins <- read_csv("data/bin_median_pnps.csv")
  taxon <- merge(taxon, pn_ps_bins, by="bin")
  
  depth_comparisons <- list(c("DCM", "MES"), c("SRF", "DCM"), c("SRF", "MES") )
  
  malaspina_bins <- read_csv("data/malaspina_bin_taxon.csv")
  pn_ps_malaspina_bins <- read_csv("data/malaspina_bin_median_pnps.csv")
  malaspina_bins$bin <- gsub('mp-deep_mag-', 'deep_MAG_', malaspina_bins$magId)
  malaspina_bins <- merge(malaspina_bins, pn_ps_malaspina_bins, by="bin")
  
  return(list(taxon, depth_comparisons, malaspina_bins))
}

init_integron <- function(){
  init_env()
  # pn_ps_integron <- read_csv("data/integron_pnps_exploded.csv")
  # pn_ps_integron <- pn_ps_integron %>%  mutate(gene_type = ifelse(category == "Defense mechanisms",  "defense mech", "non-defense cassette gene calls"))
  
  pn_ps_integron <- read_csv("data/all_integron_with_pnps_exploded.csv")
  pn_ps_bins <- read_csv("data/bin_median_pnps.csv")
  pn_ps_total <- read_csv("data/pNpS_total.csv")
  pn_ps_total <- pn_ps_total %>% mutate(gene_type = 
                              ifelse(integron == "Y", "all cassette genes", ifelse(transposase == "Y", "transposase", "normal") ))
  pn_ps_total$gene_type <- factor(pn_ps_total$gene_type, levels = c("normal", "transposase", "all cassette genes"))
  pn_ps_total <- pn_ps_total %>%  mutate(log_pnps = ifelse(pnps < 0.001, -3, log10(pnps)))
  
  # tara_integron_summary <- read_csv("data/all_integron_func_category_count.csv")
  tara_integron_summary <- read_csv("data/integron_func_category_count_known_function.csv")
  tara_integron_summary$ratio_prop = tara_integron_summary$integron_prop/tara_integron_summary$normal_prop
  tara_integron_summary <- tara_integron_summary[order(tara_integron_summary$COG_function),]
  tara_integron_all <- read_csv("data/all_ocean_merged_integrons.csv")
  
  deep_integron_summary <- read_csv("data/all_deep_integron_func_category_count.csv")
  deep_integron_summary$ratio_prop = deep_integron_summary$integron_prop/deep_integron_summary$deep_prop
  deep_integron_summary <- deep_integron_summary[order(deep_integron_summary$COG_function),]
  deep_integron_all <- read_csv("data/deep_integrons_merged_func_category.csv")
  
  return(list(pn_ps_integron, pn_ps_bins, pn_ps_total, tara_integron_summary, deep_integron_summary, tara_integron_all, deep_integron_all))
}

init_transposase_in_bins <- function(){
  init_env()
  tara_trans_in_bins <- read_csv("data/transposase_in_bins_pnps.csv")
  mala_trans_in_bins <- read_csv("data/malaspina_transposase_in_bins.csv")
  
  return(list(tara_trans_in_bins, mala_trans_in_bins))
}

init_tara <- function(){
  init_env()
  DNA_cov_file <- read_excel("data/DNA_Biofilm_Trans_ToxinAntitoxin_Coverage.xlsx")
  RNA_cov_file <- read_excel("data/RNA_Biofilm_Trans_ToxinAntitoxin_Coverage.xlsx")
  
  DNA_RNA_connector <- read_excel("data/DNA_RNA_connector.xlsx")
  DNA_Metadata <- read_excel("data/DNA_Location_Metadata.xlsx")
  DNA_Metadata$Layer_DNA <- factor(DNA_Metadata$Layer_DNA, levels = c("SRF", "DCM", "MIX", "MES"))
  DNA_Metadata$upper_size_dna <- factor(DNA_Metadata$upper_size_dna, levels = c("1.6", "3"))
  RNA_Metadata <- read_excel("data/RNA_Location_Metadata.xlsx")
  RNA_Metadata$Layer_RNA <- factor(RNA_Metadata$Layer_RNA, levels = c("SRF", "DCM", "MIX", "MES"))
  RNA_Metadata$upper_size_rna <- factor(RNA_Metadata$upper_size_rna, levels = c("1.6", "3"))
  
  DNA_tara <- merge(x=DNA_Metadata, y=DNA_cov_file, by ='connector_DNA', all = TRUE)
  RNA_tara <- merge(x=RNA_Metadata, y=RNA_cov_file, by ='connector_RNA', all = TRUE)
  
  DNA_Merged <- merge(x=DNA_RNA_connector, y=DNA_tara, by = 'connector_DNA', all = FALSE)
  DNA_RNA_tara <- merge(x=DNA_Merged, y=RNA_tara, by ='connector_RNA', all = FALSE)
  
  #DNA_RNA_tara$biofilm_exp_rate <- DNA_RNA_tara$log_rna_biofilm/DNA_RNA_tara$log_dna_biofilm
  #DNA_RNA_tara$trans_exp_rate <- DNA_RNA_tara$log_rna_trans/DNA_RNA_tara$log_dna_trans
  #DNA_RNA_tara$toxin_exp_rate <- DNA_RNA_tara$log_rna_toxin/DNA_RNA_tara$log_dna_toxin
  #DNA_RNA_tara$defense_exp_rate <- DNA_RNA_tara$log_rna_defense/DNA_RNA_tara$log_dna_defense
  
  #DNA_RNA_tara$biofilm_exp_rate <- log(DNA_RNA_tara$RNA_Biofilm/DNA_RNA_tara$DNA_Biofilm)
  #DNA_RNA_tara$trans_exp_rate <- log(DNA_RNA_tara$RNA_Trans/DNA_RNA_tara$DNA_Transposase)
  #DNA_RNA_tara$toxin_exp_rate <- log(DNA_RNA_tara$RNA_TA/DNA_RNA_tara$DNA_TA)
  #DNA_RNA_tara$defense_exp_rate <- log(DNA_RNA_tara$RNA_Defense/DNA_RNA_tara$DNA_Defense)
  
  DNA_RNA_tara$biofilm_exp_rate <- DNA_RNA_tara$RNA_Biofilm/DNA_RNA_tara$DNA_Biofilm
  DNA_RNA_tara$trans_exp_rate <- DNA_RNA_tara$RNA_Trans/DNA_RNA_tara$DNA_Transposase
  DNA_RNA_tara$toxin_exp_rate <- DNA_RNA_tara$RNA_TA/DNA_RNA_tara$DNA_TA
  DNA_RNA_tara$defense_exp_rate <- DNA_RNA_tara$RNA_Defense/DNA_RNA_tara$DNA_Defense
  
  DNA_RNA_tara <- DNA_RNA_tara%>%group_by(Layer_DNA)%>%
    mutate(mean_layer_trans_exp=mean(trans_exp_rate), 
           mean_layer_biofilm_exp=mean(biofilm_exp_rate),
           median_dna_trans = median(DNA_Transposase),
           median_dna_biofilm = median(DNA_Biofilm),
           median_dna_defense = median(DNA_Defense),
           mean_dna_trans = mean(DNA_Transposase),
           mean_dna_biofilm = mean(DNA_Biofilm),
           mean_dna_defense = mean(DNA_Defense),
           mean_dna_toxin = mean(DNA_TA),
           median_rna_trans = median(RNA_Trans),
           median_rna_biofilm = median(RNA_Biofilm),
           median_rna_defense = median(RNA_Defense)
           )
  
  rm(DNA_cov_file, RNA_cov_file, DNA_Merged, DNA_RNA_connector, DNA_Metadata, RNA_Metadata)
  
  depth_comparisons <- list( c("SRF", "DCM"), c("DCM", "MES"), c("SRF", "MES") )
  DNA_tara$Depth <- as.numeric(DNA_tara$Depth_DNA)
  RNA_tara$Depth <- as.numeric(RNA_tara$Depth_RNA)
  DNA_tara <- mutate(DNA_tara, is_MES = ifelse(Layer_DNA == "MES", "MES", "SRF, DCM"))
  
  return(list(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparisons))
}
