init_env <- function(){
  library(readr)
  library(scatterplot3d)
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
  library(stringr)
  library(png)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
  getwd()
}

boxplot.give.n <- function(x){ return(c(y = mean(x), label = length(x)))}
boxplot.give.nr <- function(x){ 
  if (length(x)==522366) {
    c(y = mean(x), label = 500000) # round up the 500,000 non transposase ORFs
  } else {
    return(c(y = mean(x), label = signif(length(x), 4)))
  }
}

filter_outliers <- function(df, colname){
  b <- colname
  out <- boxplot(df[,b])$out
  if (length(out) == 0){
    return(df) 
  } else {
    # need to use get() because this is a dynamic variable
    return(filter(df, get(b) < max(out) && get(b) > min(out)))
  }
}

filter_outliers_vct <- function(vector){
  out <- boxplot(vector)$out
  if (length(out) == 0){ return(vector) } else {
    return(vector[-out])
  }
}

less_than <- function(vct, upper) {
  out <- c()
  for (each in vct) {
    if (each[1] < upper) {
      out <- c(each, out)
    } else { a <- 1}
  }
  return (out)
}

get_r <- function(reglm) {return(summary(reglm)$r.squared)}

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

init_individual_metagenomes <- function(){
  init_env()
  SRF <- read_csv("data/trans_biofilm_defense_and_normal_pnps_SRF.csv")
  DCM <- read_csv("data/trans_biofilm_defense_and_normal_pnps_DCM.csv")
  MES <- read_csv("data/trans_biofilm_defense_and_normal_pnps_MES.csv")
  Malaspina <- read_csv("data/trans_biofilm_defense_and_normal_pnps_BAT.csv")
  all <- rbind(SRF, DCM, MES, Malaspina)
  all <- all%>%filter(!is.na(pnps))
  all <- all%>%filter(pnps<4)
  all$log_pnps <- ifelse(all$pnps < 0.01, -2, log10(all$pnps)) 
  all$gene_type <- factor(all$gene_type, levels = c("biofilm", "defense", "n", "transposase"))
  all$depth <- factor(all$depth, levels=c("SRF", "DCM", "MES", "BAT")) 
  return(list(all))
}

init_mala_cov <- function(){
  init_env()
  mala_cov <- read_csv("data/Malaspina-transposase-coverage.csv")
  mala_cov$log_dna_trans <- log10(mala_cov$trans_prop)
  mala_cov$log_dna_biofilm <- log10(mala_cov$biofilm_prop)
  mala_cov$log_dna_defense <- log10(mala_cov$defense_prop)
  mala_cov$log_dna_toxin <- log10(mala_cov$toxin_prop)
  mala_cov$DNA_Transposase <- mala_cov$trans_prop
  mala_cov$DNA_Biofilm <- mala_cov$biofilm_prop
  mala_cov$DNA_Defense <- mala_cov$defense_prop
  mala_cov$Layer_DNA <- "BAT"
  mala_cov$is_MES = "MES, BAT"
  mala_cov$Depth <- mala_cov$Depth * (-1)
  return(mala_cov)
}

init_bins <- function(){
  init_env()
  c <- 0.0000000001 # for log(0)
  # low_trans <- c("Flavobacteria","Acidimicrobidae","novelClass_E",
  #                "Gemmatimonadetes","SAR202-2","Marinisomatia")
  low_trans <- c("Flavobacteria","Acidimicrobidae",
                 "Gemmatimonadetes","Marinisomatia")
  high_trans <- c("Alphaproteobacteria","Gammaproteobacteria",
                  "Betaproteobacteria","Actinobacteria")
  taxon <- read_csv("data/bin_taxon.csv")
  COG_diversity <- read_csv("data/tara_bin_COG_diversity2.csv")
  taxon$depth <- factor(taxon$depth, levels = c("SRF", "DCM", "MES"))
  taxon$`transposase gene calls in genome (%)` <- taxon$percent_trans
  taxon$`biofilm gene calls in genome (%)` <- taxon$percent_biofilm
  taxon$`defense mechanisms gene calls in genome (%)` <- taxon$percent_defense
  taxon$log_complete_bin_size <- log10(taxon$`complete genome size (Mbp)` + c)
  taxon$log_percent_trans <- log10(taxon$percent_trans + c)
  taxon$log_percent_biofilm <- log10(taxon$percent_biofilm + c)
  taxon$complete_ORF_count <- taxon$`Total ORFs`/taxon$`Percent Complete`*100
  taxon <- taxon %>% mutate(is_deep_sea = ifelse(depth == "unsure", "unsure", ifelse(depth == "MES", "deep_sea", "SRF&DCM")))
  taxon$is_biofilm <- ifelse(taxon$biofilm_count > 0, "present", "absent")
  taxon$is_tara <- "Tara"
  taxon$graphing_log_trans <- ifelse(taxon$log_percent_trans < -9, -2.5, taxon$log_percent_trans)
  taxon$graphing_log_biofilm <- ifelse(taxon$log_percent_biofilm < -9, -2.5, taxon$log_percent_biofilm)
  taxon$class_trans <- ifelse(taxon$Class %in% low_trans, "low",
                            ifelse(taxon$Class %in% high_trans, "high", "normal"))
  taxon$class_trans <- factor(taxon$class_trans, levels = c("high", "normal", "low"))
  taxon$Class <- gsub('Verrucomicrobiae', 'Verrucomicrobia', taxon$Class)
  
  pn_ps_bins <- read_csv("data/bin_median_pnps.csv")
  pn_ps_bins <- pn_ps_bins %>% filter(count >= 100)
  taxon <- merge(taxon, pn_ps_bins, by="bin", all.x = TRUE)
  taxon <- merge(taxon, COG_diversity, by="bin", all.x = TRUE)
  taxon$log_median_bin_pnps <- log10(taxon$median_bin_pnps)
  taxon$size_fraction <- factor(taxon$size_fraction, levels = c("planktonic", "mixed", "particle"))
  
  depth_comparisons <- list(c("DCM", "MES"), c("SRF", "DCM"), c("SRF", "MES") )
  
  malaspina_bins <- read_csv("data/malaspina_bin_taxon_over70complete.csv")
  pn_ps_malaspina_bins <- read_csv("data/malaspina_bin_median_pnps.csv")
  malaspina_bins_trans_biofilm_TA <- read_csv("data/Malaspina_origin_biofilm_trans_TA_each_bin.csv")
  malaspina_bins$bin <- gsub('mp-deep_mag-', 'deep_MAG_', malaspina_bins$magId)
  malaspina_bins <- merge(malaspina_bins, pn_ps_malaspina_bins, by="bin", all.x = TRUE)
  malaspina_bins <- merge(malaspina_bins, malaspina_bins_trans_biofilm_TA, by="bin", all.x = TRUE)
  malaspina_bins$log_percent_trans <- log10(malaspina_bins$percent_trans + c)
  malaspina_bins$log_percent_biofilm <- log10(malaspina_bins$percent_biofilm + c)
  malaspina_bins$log_median_bin_pnps <- log10(malaspina_bins$median_bin_pnps)
  malaspina_bins$graphing_log_trans <- ifelse(malaspina_bins$log_percent_trans < -9, -2.5, malaspina_bins$log_percent_trans)
  malaspina_bins$graphing_log_biofilm <- ifelse(malaspina_bins$log_percent_biofilm < -9, -2.5, malaspina_bins$log_percent_biofilm)
  malaspina_bins$Class <- gsub('Acidimicrobiia', 'Acidimicrobidae', malaspina_bins$Class)
  malaspina_bins$Class <- gsub('Verrucomicrobiae', 'Verrucomicrobia', malaspina_bins$Class)
  malaspina_bins$class_trans <- ifelse(malaspina_bins$Class %in% low_trans, "low",
                              ifelse(malaspina_bins$Class %in% high_trans, "high", "normal"))
  malaspina_bins$class_trans <- factor(malaspina_bins$class_trans, levels = c("high", "normal", "low"))
  
  tmp <- c(taxon[, "Class"], malaspina_bins[, "Class"])
  tmp<-tmp[!is.na(tmp)]
  tmp2 <- as.data.frame(table(tmp)) %>% filter(Freq > 9)
  big_taxa <- tmp2[ , "tmp"]
  taxon$`Class with more than 10 MAGs` <- ifelse(is.na(taxon$Class), "Others Or Unknown",
                                          ifelse(taxon$Class %in% big_taxa, taxon$Class,"Others Or Unknown"))
  malaspina_bins$`Class with more than 10 MAGs` <- ifelse(is.na(malaspina_bins$Class), "Others Or Unknown",
                                          ifelse(malaspina_bins$Class %in% big_taxa, malaspina_bins$Class,"Others Or Unknown"))
  malaspina_bins$is_biofilm <- ifelse(malaspina_bins$biofilm_count > 0, "present", "absent")
  malaspina_bins$depth <- "BAT"
  malaspina_bins$g_depth <- "Deep\nMalaspina"
  malaspina_bins$is_deep_sea <- "deep_sea"
  malaspina_bins$size_fraction <- factor(malaspina_bins$size_fraction, levels = c("planktonic", "mixed", "particle"))
  return(list(taxon, depth_comparisons, malaspina_bins))
}

init_integron <- function(){
  init_env()
  pn_ps_integron <- read_csv("data/all_integron_with_pnps_exploded.csv")
  pn_ps_integron <- pn_ps_integron %>% mutate(log_pnps = ifelse(pnps < 0.01, -2, log10(pnps)))
  pn_ps_bins <- read_csv("data/bin_median_pnps.csv")
  pn_ps_total <- read_csv("data/pNpS_total.csv")
  # pn_ps_total <- pn_ps_total %>% mutate(gene_type = ifelse(integron == "Y", "all cassette genes", 
  #                             ifelse(transposase == "Y", "transposase", "normal") ))
  tara_non_cassette <- pn_ps_total %>% 
    filter(integron == "N" & transposase == "N") %>% 
    mutate(gene_type = "normal")
  tara_non_cassette <- tara_non_cassette %>% 
    filter(pnps < 4) %>% filter(!is.na(pnps)) %>%
    mutate(log_pnps = ifelse(pnps < 0.01, -2, log10(pnps)))
  
  malaspina_total <- read_csv("data/trans_biofilm_defense_and_normal_pnps_BAT.csv")
  
  deep_integron <- read_csv("data/malaspina_pNpS2_integron.csv")
  deep_integron <- deep_integron %>% 
    mutate(log_pnps = ifelse(pnps < 0.01, -2, log10(pnps))) %>%
    mutate(gene_type = ifelse(COG20_CATEGORY == "nan", "no_call", 
                       ifelse(COG20_CATEGORY == "V", "defense", "non_defense")))
  deep_non_cassette <- read_csv("data/malaspina_pNpS2_non_integrons_subsampled.csv")
  deep_non_cassette <- deep_non_cassette %>% 
    filter(pnps < 4) %>% filter(!is.na(pnps)) %>%
    mutate(log_pnps = ifelse(pnps < 0.01, -2, log10(pnps)))
  # deep_non_trans <- deep_non_cassette %>% filter(gene_type != "_T")
  deep_non_cassette$gene_type <- "normal"
  
  # 364,720 Tara Oceans ORFs and 135,280 Malaspina ORFs
  # following the ratio of 1,286 and 477 ORFs from the 
  # Tara Oceans and Malaspina metagenomes with pN/pS
  
  return(list(pn_ps_integron, pn_ps_bins, sample_n(tara_non_cassette, 364720), 
              sample_n(deep_non_cassette, 135280), deep_integron))
}

init_integron_category <- function(){
  summary <- read_csv("data/tara_malaspina_integron_func_category_count.csv")
  
  col1 <- c('Secondary metabolites biosynthesis, transport and catabolism',
            "Posttranslational modification, protein turnover, chaperones", 
            "Intracellular trafficking, secretion, and vesicular transport",
            "Cell cycle control, cell division, chromosome partitioning",
            "Translation, ribosomal structure and biogenesis",
            "transport and metabolism")
  col2 <- c("Secondary metabolites...", "Posttranslational modification...",
            "Intracellular trafficking...", "Cell cycle control...",
            "Translation, ribosomal structure...", "TM*")
  replace <- data.frame(col1, col2)
  
  for (r in 1:nrow(replace)){
    summary$COG_function <- gsub(replace[r,1],replace[r,2],summary$COG_function)
  }
  
  return(summary)
}

init_transposase_in_bins <- function(){
  init_env()
  tara_trans_in_bins <- read_csv("data/transposase_in_bins_pnps.csv")
  mala_trans_in_bins <- read_csv("data/malaspina_transposase_in_bins.csv")
  
  return(list(tara_trans_in_bins, mala_trans_in_bins))
}

init_tara <- function(){
  init_env()
  DNA_cov <- read_excel("data/DNA_Biofilm_Trans_Defense_Coverage.xlsx")
  DNA_cov$log_dna_biofilm <- log10(DNA_cov$DNA_Biofilm)
  DNA_cov$log_dna_trans <- log10(DNA_cov$DNA_Transposase)
  DNA_cov$log_dna_defense <- log10(DNA_cov$DNA_Defense)
  
  RNA_cov <- read_excel("data/RNA_Biofilm_Trans_Defense_Coverage.xlsx")
  RNA_cov$log_rna_biofilm <- log10(RNA_cov$RNA_Biofilm)
  RNA_cov$log_rna_trans <- log10(RNA_cov$RNA_Transposase)
  RNA_cov$log_rna_defense <- log10(RNA_cov$RNA_Defense)
  
  DNA_RNA_connector <- read_excel("data/DNA_RNA_connector.xlsx")
  DNA_Metadata <- read_excel("data/DNA_Location_Metadata.xlsx")
  DNA_Metadata$Depth <- as.numeric(DNA_Metadata$Depth_DNA)
  DNA_Metadata$Layer_DNA <- factor(DNA_Metadata$Layer_DNA, levels = c("SRF", "DCM", "MIX", "MES"))
  DNA_Metadata$Ocean_short<-factor(DNA_Metadata$Ocean_short)
  DNA_Metadata$upper_size_dna <- factor(DNA_Metadata$upper_size_dna, levels = c("1.6", "3"))
  RNA_Metadata <- read_excel("data/RNA_Location_Metadata.xlsx")
  RNA_Metadata <- filter(RNA_Metadata, Depth_RNA != "NA")
  RNA_Metadata$Depth <- as.numeric(RNA_Metadata$Depth_RNA)
  RNA_Metadata <- RNA_Metadata %>% mutate(Layer_RNA = ifelse(Layer_RNA != "MIX", Layer_RNA,
                                              ifelse(Depth_RNA <= 120, "DCM", "MES")))
  RNA_Metadata$Layer_RNA <- factor(RNA_Metadata$Layer_RNA, levels = c("SRF", "DCM", "MES"))
  RNA_Metadata$upper_size_rna <- factor(RNA_Metadata$upper_size_rna, levels = c("1.6", "3"))
  
  DNA_tara <- merge(x=DNA_Metadata, y=DNA_cov, by ='connector_DNA', all = TRUE)
  RNA_tara <- merge(x=RNA_Metadata, y=RNA_cov, by ='connector_RNA', all = TRUE)
  
  DNA_Merged <- merge(x=DNA_RNA_connector, y=DNA_tara, by = 'connector_DNA', all = FALSE)
  DNA_RNA_tara <- merge(x=DNA_Merged, y=RNA_tara, by ='connector_RNA', all = FALSE)
  
  DNA_RNA_tara$biofilm_exp_rate <- DNA_RNA_tara$RNA_Biofilm/DNA_RNA_tara$DNA_Biofilm
  DNA_RNA_tara$trans_exp_rate <- DNA_RNA_tara$RNA_Trans/DNA_RNA_tara$DNA_Transposase
  DNA_RNA_tara$defense_exp_rate <- DNA_RNA_tara$RNA_Defense/DNA_RNA_tara$DNA_Defense
  
  rm(DNA_cov, RNA_cov, DNA_Merged, DNA_RNA_connector, DNA_Metadata, RNA_Metadata)
  
  depth_comparisons <- list( c("SRF", "DCM"), c("DCM", "MES"), c("SRF", "MES") )
  DNA_tara <- mutate(DNA_tara, is_MES = ifelse(Layer_DNA == "MES", "MES, BAT", 
                                               ifelse(Layer_DNA == "MIX", "MIX", "SRF, DCM")))
  
  return(list(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparisons))
}
