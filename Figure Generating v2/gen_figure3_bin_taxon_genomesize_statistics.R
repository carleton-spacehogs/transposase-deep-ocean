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

sum <- summary(add4.lm)

get_rp(depth.lm)
get_rp(gsize_depth.lm)
get_rp(lifestyle_depth.lm)
get_rp(add1.lm)
get_rp(add2.lm)
get_rp(add3.lm)
get_rp(add4.lm)


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


# which taxa class is transposase enriched?
big_taxa <- stat_all %>% 
  group_by(Class) %>%
  summarise(count = n()) %>%
  filter(count > 10) %>%
  drop_na() %>%
  as.data.frame()

big_taxa <- big_taxa$Class

bi_var <- stat_all %>% select(c("Class", "log_percent_trans"))

more_taxa <- c()
less_taxa <- c()
for (taxa in big_taxa){
  this_taxa <- bi_var %>% filter(Class==taxa)
  others <- bi_var %>% filter(Class!=taxa)
  res <- t.test(this_taxa$log_percent_trans, others$log_percent_trans)
  if (res$p.value < 0.05){
    if (res$estimate[1] < res$estimate[2]){
      less_taxa <- c(less_taxa, taxa)
      print(paste(taxa, "->", "less transposases")) # than the rest 
    } else {
      more_taxa <- c(more_taxa, taxa)
      print(paste(taxa, "->", "more transposases")) # than the rest
    }
  }
}

# > more_taxa
# [1] "Alphaproteobacteria" "Gammaproteobacteria" "Betaproteobacteria"  "Deltaproteobacteria"
# > less_taxa
# [1] "Verrucomicrobia" "Flavobacteria"   "Acidimicrobidae" "novelClass_E"    "Opitutae"       
# [6] "SAR202-2" 

stat_all %>% group_by(Class) %>%
  summarise(mean.log.trans = mean(log_percent_trans),
            count = n()) %>%
  filter(Class %in% big_taxa) %>% 
  arrange(mean.log.trans)

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

