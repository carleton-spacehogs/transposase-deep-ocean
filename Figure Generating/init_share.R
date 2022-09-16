init_env <- function(){
  library(readr)
  library(scatterplot3d)
  library(tidyverse)
  library(purrr)
  library(ggpubr)
  library(ggrepel)
  library(ggtext)
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
  library(rstatix)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
  getwd()
}

googleDrive.loc = "G:/My Drive/Rika SpaceHog Lab - Jimmy's Transposase Notebook (2020-2022)/Transposase-MS-v2"

boxplot.give.n <- function(x){ return(c(y = mean(x), label = length(x)))}
boxplot.give.nr <- function(x){ 
  if (length(x)==522366) {
    c(y = mean(x), label = 500000) # round up the 500,000 non transposase ORFs
  } else {
    return(c(y = mean(x), label = signif(length(x), 4)))
  }
}

get.fudge = function(df, var){
  vec = unlist(df[var])
  vec = vec[!is.na(vec)]
  return(min(vec[vec > 0]))
}

filter_outliers <- function(df, colname){
  b <- colname
  out <- boxplot.stats(df[,b])$out
  print(out)
  if (length(out) == 0){
    return(df) 
  } else {
    # need to use get() because this is a dynamic variable
    return(filter(df, !get(b) %in% out))
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

get_r <- function(reglm) {return(unlist(summary(reglm)['r.squared']))}

get_rp <- function(reglm) {
  r <- summary(reglm)$r.squared
  p <- rev(anova(reglm)$`Pr(>F)`)[2]
  out <- format(signif(p, digits = 3), scientific = TRUE)
  return(out)
}

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

`%--%` <- function(x, y) { # for string replacement
  do.call(sprintf, c(list(x), y))
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

get_logs = function(df, name_vct1, name_vct2, DNA_or_RNA){
  if (length(name_vct1) != length(name_vct2)) {
    print("error, 2 name vectors must have the same lenght")
    return("error")
  }
  for (i in 1:length(name_vct1)){
    old = paste(toupper(DNA_or_RNA), name_vct1[i], sep="_")
    new = paste("log", tolower(DNA_or_RNA), name_vct2[i], sep="_")
    df[new] = log2(df[old])
  }
  return(df)
}

get_CAZpep_percent = function(df, DNA_or_RNA){
  sect_CAZ = paste(toupper(DNA_or_RNA), "secretory_CAZyme", sep="_")
  sect_pep = paste(toupper(DNA_or_RNA), "secretory_peptidase", sep="_")
  total_CAZ = paste(toupper(DNA_or_RNA), "CAZyme", sep="_")
  total_pep = paste(toupper(DNA_or_RNA), "peptidase", sep="_")
  prop_CAZ = paste("percent", tolower(DNA_or_RNA), "sect_CAZ", sep="_")
  prop_pep = paste("percent", tolower(DNA_or_RNA), "sect_pep", sep="_")
  prop_both = paste("percent", tolower(DNA_or_RNA), "sect_CAZpep", sep="_")
  
  df[prop_CAZ] = df[sect_CAZ]/df[total_CAZ] * 100
  df[prop_pep] = df[sect_pep]/df[total_pep] * 100
  df[prop_both] = (df[sect_CAZ]+df[sect_pep])/(df[total_CAZ] + df[total_pep]) * 100
  return(df)
}

init_individual_metagenomes <- function(){
  init_env()
  SRF <- read_csv("data/trans_biofilm_defense_and_normal_pnps_SRF_all.csv")
  DCM <- read_csv("data/trans_biofilm_defense_and_normal_pnps_DCM_all.csv")
  MES <- read_csv("data/trans_biofilm_defense_and_normal_pnps_MES_all.csv")
  MAL <- read_csv("data/trans_biofilm_defense_and_normal_pnps_BAT_all.csv")
  all <- rbind(SRF, DCM, MES, MAL)
  all <- all%>%filter(!is.na(pnps))
  all <- all%>%filter(pnps<4)
  all$log_pnps <- ifelse(all$pnps < 0.01, -2, log10(all$pnps)) 
  all$gene_type <- factor(all$gene_type, levels = c("biofilm", "defense", "n", "transposase"))
  all$depth <- factor(all$depth, levels=c("SRF", "DCM", "MES", "BAT"))
  rm(SRF, DCM, MES, MAL)
  return(list(all))
}

# newvct = c("trans","defense","signalT","replication","sect_CAZ","sect_pep","aminoTM",
#            "lipidTM","coenzyme","energyPC","cellWall")

newvct = c("trans", "CAZyme", "peptidase", "sect_CAZ", "sect_pep",
"Lipid_TM", "Coenzyme_TM", "Signal_transduction_mechanisms", 
"defense", "Energy_production_and_conversion", "Replication_recombination_and_repair", "Cell_wall",
"Amino_acid_TM","Cell_motility", "Extracellular_struct", "Cytoskeleton",
"Carbohydrate_TM", "Inorganic_ion_TM", "Nucleotide_TM",
"Chromatin_structure_and_dynamics", "Cell_cycle_control", "Intracellular_trafficking", "Secondary_metabolites", 
"Translation_ribosomal_structure", "Transcription", "Posttranslational_modification", "Mobilome_prophages_transposons")

oldvct = c("transposase", "CAZyme", "peptidase", "secretory_CAZyme", "secretory_peptidase",
"Lipid_TM", "Coenzyme_TM", "Signal_transduction_mechanisms", 
"Defense_mechanisms", "Energy_production_and_conversion", "Replication_recombination_and_repair", "Cell_wall",
"Amino_acid_TM","Cell_motility", "Extracellular_struct", "Cytoskeleton",
"Carbohydrate_TM", "Inorganic_ion_TM", "Nucleotide_TM",
"Chromatin_structure_and_dynamics", "Cell_cycle_control", "Intracellular_trafficking", "Secondary_metabolites", 
"Translation_ribosomal_structure", "Transcription", "Posttranslational_modification", "Mobilome_prophages_transposons")

init_mala_cov <- function(){
  init_env()
  mala_cov = read_csv("../OM-RGC-and-abundance/Malaspina-genes-coverage.csv")
  # mala_cov = read_csv("../OM-RGC-and-abundance/Malaspina-genes-coverage-OM-RGC.csv")
  mala_cov = get_logs(mala_cov, oldvct, newvct, "DNA")
  mala_cov = get_CAZpep_percent(mala_cov, "DNA")
  
  mala_cov = mala_cov %>%
    mutate(Layer_DNA = "BAT",
      is_MES = "MES, BAT",
      Depth = Depth * (-1))
  return(mala_cov)
}

# assume the dataframe to be DNA_RNA_tara
get_expression = function(df, name_vct = oldvct){
  for (i in 1:length(name_vct)){
    DNA = paste0("DNA_", name_vct[i])
    RNA = paste0("RNA_", name_vct[i])
    expression = paste0(name_vct[i], "_exp")
    print(DNA)
    print(RNA)
    df[expression] = df[RNA]/df[DNA]
  }
  return(df)
}


init_MAGs_pnps_depths <- function(gene){
  init_env()
  defense="../toxin-db/all_oceans_defense_mech_depth_pnps.csv"
  toxin="../toxin-db/all_oceans_toxin_depth_pnps.csv"
  indi_toxin="../toxin-db/all_oceans_toxin_individual_depth_pnps.csv"
  COG_toxin="../toxin-db/all_oceans_COG_toxin_depth_pnps.csv"
  mobilome="../toxin-db/all_oceans_mobilome_depth_pnps.csv"
  trans="../toxin-db/all_oceans_transposase_individual_depth_pnps.csv"
  prep_pnps = function(filename){
    g=read_csv(filename)
    g$SRF_ratio = g$SRF_gene_median-g$SRF_bin_median
    g$DCM_ratio = g$DCM_gene_median-g$DCM_bin_median
    g$MES_ratio = g$MES_gene_median-g$MES_bin_median
    g$BAT_ratio = g$deep_gene_median-g$deep_bin_median
    return (g)
  }
  prep_pnps2 = function(filename){
    g=read_csv(filename)
    g$SRF_ratio = g$SRF_pnps-g$SRF_bin_median
    g$DCM_ratio = g$DCM_pnps-g$DCM_bin_median
    g$MES_ratio = g$MES_pnps-g$MES_bin_median
    g$BAT_ratio = g$deep_pnps-g$deep_bin_median
    return (g)  
  }
  return(list(prep_pnps(defense), prep_pnps(toxin), 
              prep_pnps2(indi_toxin), prep_pnps(COG_toxin),
              prep_pnps(mobilome), prep_pnps2(trans)))
}

merge_MAGs_eukaryote = function(bin_df){
  euk = read_csv("../identify_eukaryotes/MAG_composition_summary.csv")
  euk_short = euk[c("bin","eukarya_prop")]
  bin_df = merge(bin_df, euk_short, by = "bin")
  return(bin_df)
}

init_bins <- function(){
  init_env()
  sel_col = c("bin", "Class", "Order", "Family", "complete_genome_size", #"Genus", "Species", 
              "percent_trans", "percent_defense", "Total ORFs", "depth")

  taxon = read_csv("data/bin_taxon.csv")
  origin = read_csv("../MAG-related/Tara_bins_origin.csv")[,c("bin","depth","size_fraction")]
  taxon = merge(origin, taxon, by = "bin")
  
  malaspina_bins = read_csv("data/malaspina_bin_taxon_over70complete.csv")
  # pn_ps_malaspina_bins <- read_csv("data/malaspina_bin_median_pnps.csv")
  gene_abun = read_csv("data/Malaspina_origin_biofilm_trans_TA_each_bin.csv")
  malaspina_bins$bin = gsub('mp-deep_mag-', 'deep_MAG_', malaspina_bins$magId)
  malaspina_bins = merge(malaspina_bins, gene_abun, by="bin", all.x = TRUE)
  malaspina_bins$depth = "BAT"
  
  all_MAGs = rbind(taxon[sel_col], malaspina_bins[sel_col])
  all_MAGs$Class = gsub('Acidimicrobiia', 'Acidimicrobidae',
                   gsub('Verrucomicrobiae', 'Verrucomicrobia', all_MAGs$Class))
  
  Class_count = table(all_MAGs$Class[!is.na(all_MAGs$Class)])
  big_taxa = Class_count[Class_count >= 10]
  
  all_MAGs = all_MAGs %>% mutate(
    large_classes = ifelse(is.na(Class), "Others Or Unknown",
                    ifelse(Class %in% names(big_taxa), Class,"Others Or Unknown")),
    depth = factor(depth, # gsub("unsure", "MIXED", depth),
                   levels = c("SRF", "DCM", "MES", "BAT")),
    `Total ORFs` = as.numeric(`Total ORFs`)
  )
  all_MAGs = merge_MAGs_eukaryote(all_MAGs)
  return(all_MAGs)
}

merge_MAGs_CAZyme_peptide = function(bin_df){
  pa <- read_csv("../particle_association_redux/MAGs_diamond_secretory_CAZyme_peptidase.csv")
  pa[is.na(pa)] = 0
  pa <- merge(bin_df, pa, by = "bin")
  pa = pa %>% mutate( percent_sect_CAZpep =
    (signal_CAZ_count + signal_pep_count)/(total_CAZ_count + total_pep_count) * 100,                 
    percent_sect_CAZ = signal_CAZ_count/total_CAZ_count * 100,
    percent_sect_pep = signal_pep_count/total_pep_count * 100)
  pa$percent_sect_pep[is.na(pa$percent_sect_pep)] = 0 # no peptidase at all
  return(pa)
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
  
  malaspina_total <- read_csv("data/trans_biofilm_defense_and_normal_pnps_BAT_all.csv")
  
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
  # summary = read_csv('../integron_finder_v2/integron_known_COG_category_summary.csv')
  
  # col1 <- c('Secondary metabolites biosynthesis, transport and catabolism',
  #           "Mobilome: prophages, transposons",
  #           "Posttranslational modification, protein turnover, chaperones", 
  #           "Intracellular trafficking, secretion, and vesicular transport",
  #           "Cell cycle control, cell division, chromosome partitioning",
  #           "Replication, recombination and repair",
  #           "Translation, ribosomal structure and biogenesis")
  # col2 <- c("Secondary metabolites", "Mobilome prophages transposons",
  #           "Posttranslational modification",
  #           "Intracellular trafficking", "Cell cycle control",
  #           "Replication recombination and repair",
  #           "Translation, ribosomal structure")
  
  col1 <- c('Secondary metabolites biosynthesis, transport and catabolism',
            "Posttranslational modification, protein turnover, chaperones", 
            "Intracellular trafficking, secretion, and vesicular transport",
            "Cell cycle control, cell division, chromosome partitioning",
            "Translation, ribosomal structure and biogenesis",
            "transport and metabolism")
  col2 <- c("Secondary metabolites...", "Posttranslational modification...",
            "Intracellular trafficking...", "Cell cycle control...",
            "Translation, ribosomal structure...", "T.M.")

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

init_tara_md = function(factorize = TRUE){
  DNA_Metadata <- read_excel("data/DNA_Location_Metadata.xlsx")
  DNA_Metadata$Depth <- as.numeric(DNA_Metadata$Depth_DNA)
  
  RNA_Metadata <- read_excel("data/RNA_Location_Metadata.xlsx")
  RNA_Metadata = RNA_Metadata %>% mutate(
    Depth = as.numeric(Depth_RNA),
    Layer_RNA = ifelse(Layer_RNA != "MIX", Layer_RNA,
                ifelse(Depth_RNA <= 120, "DCM", "MES"))) %>% 
    filter(Depth_RNA != "NA") 
  
  if (factorize){
    DNA_Metadata = DNA_Metadata %>% mutate(
      Layer_DNA = factor(Layer_DNA, levels = c("SRF", "DCM", "MIX", "MES")),
      Ocean_short = factor(Ocean_short),
      upper_size_dna = factor(upper_size_dna, levels = c("1.6", "3")))
    RNA_Metadata$Layer_RNA <- factor(RNA_Metadata$Layer_RNA, levels = c("SRF", "DCM", "MES"))
    RNA_Metadata$upper_size_rna <- factor(RNA_Metadata$upper_size_rna, levels = c("1.6", "3")) 
  }
  return(list(DNA_Metadata, RNA_Metadata))
}

init_tara <- function(){
  init_env()
  DNA_cov = read_csv("../OM-RGC-and-abundance/tara_DNA_genes_abundance.csv")
  DNA_cov = get_logs(DNA_cov, oldvct, newvct, "DNA")
  DNA_cov = get_CAZpep_percent(DNA_cov, "DNA")
  
  RNA_cov = read_csv("../OM-RGC-and-abundance/tara_RNA_genes_abundance.csv")
  RNA_cov = get_logs(RNA_cov, oldvct, newvct, "RNA")
  RNA_cov = get_CAZpep_percent(RNA_cov, "RNA")

  g(DNA_Metadata, RNA_Metadata) %=% init_tara_md()
  DNA_tara <- merge(x=DNA_Metadata, y=DNA_cov, by ='connector_DNA', all = TRUE)
  RNA_tara <- merge(x=RNA_Metadata, y=RNA_cov, by ='connector_RNA', all = TRUE)
  
  DNA_tara = DNA_tara %>% mutate(
    is_MES = ifelse(Layer_DNA == "MES", "MES, BAT", 
             ifelse(Layer_DNA == "MIX", "MIX", "SRF, DCM")))
  
  DNA_RNA_connector = read_excel("data/DNA_RNA_connector.xlsx")
  tmp = merge(x=DNA_RNA_connector, y=DNA_tara, by = 'connector_DNA', all = FALSE)
  tmp = tmp[ ,!(colnames(tmp) %in% c("ENA_Run_ID", "Depth"))] # have it in RNA
  DNA_RNA_tara = merge(x=tmp, y=RNA_tara, by ='connector_RNA', all = FALSE)
  DNA_RNA_tara = DNA_RNA_tara[ ,!grepl("sum|read", colnames(DNA_RNA_tara))]
  
  return(list(DNA_tara, RNA_tara, DNA_RNA_tara))
}
