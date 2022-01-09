setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()

all_tara <- bin_taxon %>% 
  select("complete genome size (Mbp)", "percent_trans", "log_percent_trans", "class_trans",
         "depth", "Class", "is_deep_sea", "is_biofilm", "percent_biofilm", "log_percent_biofilm",
         "Class with more than 10 MAGs") %>%
  filter(!is.na(depth))

all_mala <- malaspina_bins %>% 
  select(c("complete genome size (Mbp)", "percent_trans", "log_percent_trans", "class_trans",
           "depth", "Class", "is_deep_sea", "is_biofilm","percent_biofilm", "log_percent_biofilm",
           "Class with more than 10 MAGs"))

all  <- filter_outliers(rbind(all_tara, all_mala),"percent_trans")
# all <- rbind(all_tara, all_mala)

depth.lm <- lm(log_percent_trans~depth, data=all)
summary(depth.lm)

taxon_depth.lm <- lm(log_percent_trans~depth+Class, data=all)
anova(taxon_depth.lm)

taxon_depth2.lm <- lm(log_percent_trans~depth*Class, data=all)
anova(taxon_depth2.lm)

biofilm_depth.lm <- lm(log_percent_trans~depth+log_percent_biofilm, all)
anova(biofilm_depth.lm)

gsize_depth.lm <- lm(log_percent_trans~depth+`complete genome size (Mbp)`, all)
anova(gsize_depth.lm)

# binary is_biofilm not as good as log_percent_biofilm
allq.lm <- update(taxon_depth2.lm, .~. +log_percent_biofilm)
anova(allq.lm)
get_r(allq.lm)

all1.lm <- update(allq.lm, .~. +`complete genome size (Mbp)`)
anova(all1.lm)
get_r(all1.lm)

quadra_all.lm <- lm(log_percent_trans~depth*Class +log_percent_biofilm 
                    +`complete genome size (Mbp)` + I(`complete genome size (Mbp)`^2), all)
summary(quadra_bin_size.lm)

anova(quadra_bin_size.lm)
get_r(linear_bin_size.lm)
get_r(quadra_bin_size.lm)

# get high/low abundance class
coe <- as.data.frame(summary(quadra_bin_size.lm)$coefficients)
colnames(coe)<- c('estimate', 'error', 'tvalue', "pvalue")
coe <- tibble::rownames_to_column(coe, "varname")
coe_class <- coe %>%
  slice(grep("^Clas", varname)) %>% 
  filter(pvalue < 0.05)


anova(lm(log_percent_trans~depth*Class + log_percent_biofilm +`complete genome size (Mbp)`, all))
summary(lm(`complete genome size (Mbp)`~percent_biofilm, all))

library("car")
car::vif(ally.lm, singular.ok = TRUE)
allx.lm <- update(taxon_depth.lm, .~. +`complete genome size (Mbp)`)
ally.lm <- update(allx.lm, .~. +log_percent_biofilm)
anova(ally.lm)


all$class_trans1 <- ifelse(all$Class %in% low_trans, "low", "not_low")
all$class_trans2 <- ifelse(all$Class %in% high_trans, "high", "not_high")

all %>% 
  with(table(class_trans1, is_deep_sea)) %>% 
  chisq.test()

all %>% 
  with(table(class_trans2, is_deep_sea)) %>% 
  chisq.test()

all %>% 
  with(table(class_trans, is_deep_sea)) %>% 
  chisq.test()

big_taxa <- unique(all$`Class with more than 10 MAGs`)
big_taxa <- big_taxa[!big_taxa == "Others Or Unknown"]
bi_var <- all %>% select(c("Class with more than 10 MAGs", "percent_trans"))

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

depth_stats_taxon <- all %>% 
  with(table(`Class with more than 10 MAGs`, depth)) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column(var = "Class with more than 10 MAGs")
write.csv(depth_stats_taxon,"taxon_depth_summary.csv", row.names = FALSE)

MES <- all %>% filter(depth == "MES")
Mala <- all %>% filter(depth == "Deep Malaspina")
t.test(MES$percent_trans, Mala$percent_trans)

t.test(`complete genome size (Mbp)`~is_deep_sea, all)

