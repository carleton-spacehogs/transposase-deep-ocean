setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()

col_list <- c("complete genome size (Mbp)", "percent_trans", "log_percent_trans", "class_trans",
              "depth", "Class", "is_deep_sea", "is_biofilm", "percent_biofilm", "log_percent_biofilm",
              "Class with more than 10 MAGs", "size_fraction")
stat_all_tara <- bin_taxon[,col_list] # %>% filter(!is.na(depth))

stat_all_mala <- malaspina_bins[,col_list] %>%
  filter(size_fraction != "error")

stat_all <- rbind(stat_all_tara, stat_all_mala)

t.test(percent_trans~is_deep_sea, stat_all)
print(paste("95% CI:", 1/0.1522989, 1/0.1858968))

stat_all_deep <- stat_all %>% filter(percent_trans>0) %>% 
  filter(is_deep_sea == "deep_sea")

stat_all_s <- stat_all %>% filter(percent_trans>0) %>% 
  filter(is_deep_sea == "SRF&DCM")

summary(lm(log_percent_trans~Class, data=stat_all))
summary(lm(log_percent_trans~`complete genome size (Mbp)`, data=stat_all))
summary(lm(log_percent_trans~size_fraction, data=stat_all))


depth.lm <- lm(log_percent_trans~depth, data=stat_all)
gsize_depth.lm <- lm(log_percent_trans~depth + `complete genome size (Mbp)`, stat_all)
lifestyle_depth.lm <- lm(log_percent_trans~depth + size_fraction, stat_all)
add1.lm <- update(gsize_depth.lm, .~. +size_fraction)
add2.lm <- update(add1.lm, .~. +Class)
add3.lm <- update(add2.lm, .~. +Class*depth)
add4.lm <- update(add3.lm, .~. +I(`complete genome size (Mbp)`^2))

get_r(depth.lm)
get_r(gsize_depth.lm)
get_r(lifestyle_depth.lm)
get_r(add1.lm)
get_r(add2.lm)
get_r(add3.lm)
get_r(add4.lm)

anova(add1.lm)
anova(add3.lm)
anova(add4.lm)












taxon_depth.lm <- lm(log_percent_trans~depth+Class, data=stat_all)
anova(taxon_depth.lm)

taxon_depth2.lm <- lm(log_percent_trans~depth*Class, data=stat_all)
anova(taxon_depth2.lm)

biofilm_depth.lm <- lm(log_percent_trans~depth+log_percent_biofilm, stat_all)
anova(biofilm_depth.lm)

size_fraction_depth.lm <- lm(log_percent_trans~depth+size_fraction, stat_all)
get_r(size_fraction_depth.lm)

gsize_depth.lm <- lm(log_percent_trans~depth+`complete genome size (Mbp)`, stat_all)
anova(gsize_depth.lm)

# binary is_biofilm not as good as log_percent_biofilm
# log_percent_biofilm is not as good as size_fraction

stat_allq.lm <- update(taxon_depth2.lm, .~. +size_fraction)
stat_all1.lm <- update(stat_allq.lm, .~. +`complete genome size (Mbp)`)
stat_all2.lm <- update(stat_all1.lm, .~. + I(`complete genome size (Mbp)`^2), stat_all)
summary(stat_all2.lm)

anova(stat_all2.lm)
get_r(stat_all2.lm)

# get high/low abundance class
coe <- as.data.frame(summary(quadra_bin_size.lm)$coefficients)
colnames(coe)<- c('estimate', 'error', 'tvalue', "pvalue")
coe <- tibble::rownames_to_column(coe, "varname")
coe_class <- coe %>%
  slice(grep("^Clas", varname)) %>% 
  filter(pvalue < 0.05)

library("car")
car::vif(stat_ally.lm, singular.ok = TRUE)
stat_allx.lm <- update(taxon_depth.lm, .~. +`complete genome size (Mbp)`)
stat_ally.lm <- update(stat_allx.lm, .~. +log_percent_biofilm)
anova(stat_ally.lm)


stat_all$class_trans1 <- ifelse(stat_all$class_trans == "low", "low", "not_low")
stat_all$class_trans2 <- ifelse(stat_all$class_trans == "high", "high", "not_high")

biofilm_stat_all <- filter_outliers(stat_all, "percent_biofilm")
t.test(percent_biofilm ~ class_trans2, stat_all)
boxplot(percent_biofilm ~ class_trans, biofilm_stat_all)
t.test(percent_biofilm ~ class_trans, biofilm_stat_all %>% filter(class_trans != "high"))
boxplot(`complete genome size (Mbp)`~is_biofilm, stat_all, ylim = c(0, 9))


stat_all %>% 
  with(table(class_trans1, is_deep_sea)) %>% 
  chisq.test()

stat_all %>% 
  with(table(class_trans2, is_deep_sea)) %>% 
  chisq.test()

stat_all %>% 
  with(table(class_trans, is_deep_sea)) %>% 
  chisq.test()



big_taxa <- unique(stat_all$`Class with more than 10 MAGs`)
big_taxa <- big_taxa[!big_taxa == "Others Or Unknown"]
bi_var <- stat_all %>% select(c("Class with more than 10 MAGs", "percent_trans"))

more_taxa <- c()
less_taxa <- c()
for (taxa in big_taxa){
  this_taxa <- bi_var %>% filter(`Class with more than 10 MAGs`==taxa)
  others <- bi_var %>% filter(`Class with more than 10 MAGs`!=taxa)
  res <- t.test(this_taxa$percent_trans, others$percent_trans)
  if (res$p.value < 0.05){
    if (res$estimate[1] < res$estimate[2]){
      less_taxa <- c(less_taxa, taxa)
      print(paste(taxa, "->", "less transposases than the rest"))
    } else {
      more_taxa <- c(more_taxa, taxa)
      print(paste(taxa, "->", "more transposases than the rest"))
    }
  }
}

# > more_taxa
# [1] "Alphaproteobacteria" "Gammaproteobacteria" "Betaproteobacteria"  "Actinobacteria"     
# > less_taxa
# [1] "Flavobacteria"    "Acidimicrobidae"  "novelClass_E"     "Gemmatimonadetes"
# [5] "SAR202-2"         "Marinisomatia" 

class_trans <-bi_var %>% group_by(`Class with more than 10 MAGs`) %>%
  summarise(across(where(is.numeric), list(mean = mean)))

class_trans <- class_trans[order(-class_trans$percent_trans_mean),]

depth_stats_taxon <- stat_all %>% 
  with(table(`Class with more than 10 MAGs`, depth)) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column(var = "Class with more than 10 MAGs")
write.csv(depth_stats_taxon,"taxon_depth_summary.csv", row.names = FALSE)

SRF <- stat_all %>% filter(depth == "SRF")
DCM <- stat_all %>% filter(depth == "DCM")
MES <- stat_all %>% filter(depth == "MES")
Mala <- stat_all %>% filter(depth == "Deep Malaspina")
t.test(MES$percent_trans, Mala$percent_trans)

t.test(`complete genome size (Mbp)`~is_deep_sea, stat_all)

SRF[,"percent_trans"]

for (each in list(SRF, DCM, MES, Mala)){
  print(mean(each[,"percent_trans"]))
}

