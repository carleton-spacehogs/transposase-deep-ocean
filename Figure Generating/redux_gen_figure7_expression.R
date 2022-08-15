setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
options(scipen=10000)

DNA_RNA_tara <- DNA_RNA_tara %>%filter(Layer_DNA %in% c("SRF", "DCM", "MES"))
# cols <- c("biofilm_exp_rate", "trans_exp_rate", "defense_exp_rate", "Layer_RNA",
#           "log_dna_biofilm", "log_dna_defense", "log_dna_trans",
#           "log_rna_biofilm", "log_rna_defense", "log_rna_trans", "upper_size_rna")
DNA_RNA = DNA_RNA_tara

gen_plot <- function(gene_name, log_scale, percent_scale, upscale, 
                     DNA_color, RNA_color, txt_top, txt_bot){
  log_scale <- log_scale + upscale
  DNA_gene <- paste("log_dna",gene_name, sep="_")
  RNA_gene <- paste("log_rna",gene_name, sep="_")
  tmp1 <- DNA_RNA[,c(DNA_gene, "Layer_RNA")]
  tmp2 <- DNA_RNA[,c(RNA_gene, "Layer_RNA")]
  tmp1$type <- "DNA"
  tmp2$type <- "RNA"
  
  tmp1 <- filter_outliers(tmp1, DNA_gene)
  tmp2 <- filter_outliers(tmp2, RNA_gene)
  
  g_cols <- c("log_abundance", "depth", "type")
  DNA_RNA_bar <- rbind(`colnames<-`(tmp1, g_cols),
                       `colnames<-`(tmp2, g_cols))
  
  sum_bar <- DNA_RNA_bar %>% 
    group_by(type, depth) %>% 
    summarise(median_abun = median(log_abundance) + upscale, 
    mean_abun = mean(log_abundance),
    lower = quantile(log_abundance, probs = c(.25)) + upscale, 
    upper = quantile(log_abundance, probs = c(.75)) + upscale)
  
  x_max <- max(max(DNA_RNA[,c(DNA_gene)]), max(DNA_RNA[,c(RNA_gene)]))+ upscale
  labD = "DNA"
  labR = "RNA"
  
  p <- ggplot(sum_bar, aes(y=fct_rev(depth), x=median_abun, fill=fct_rev(type))) + 
    geom_bar(stat="identity", position=position_dodge()) +
    # scale_x_continuous(breaks = log_scale, labels = percent_scale) +# 
    geom_errorbar(aes(xmin=lower, xmax=upper), width=.2,
                  position=position_dodge(.9)) + 
    scale_fill_manual(values = c("DNA" = DNA_color, "RNA" = RNA_color))+
    annotate(geom="text", x=x_max*txt_top, y=3.25, label=labD) + 
    annotate(geom="text", x=x_max*txt_bot, y=2.8, label=labR) + 
    theme_classic() +
    
    theme(axis.title.y = element_blank(),
          legend.position = "none")
  
  if (gene_name == "trans") {
    return(p + xlab("% reads mapped to transposases"))
  } else if (gene_name == "defense") {
    return(p + xlab("% reads mapped to defense mechanism ORFs"))
  } else if (gene_name == "biofilm"){
    p <- p + xlab("% reads mapped to biofilm ORFs")
    p <- p + scale_y_discrete(labels=c("MES  \n(n = 22)",
                                       "DCM  \n(n = 36)",
                                       "SRF   \n(n = 64)"))
  }
  return(p)
}

exp_biofilm <- gen_plot("biofilm", c(-4.301, -4, -3.698, -3.3979),
                        c("0.005%", "0.01%", "0.02%", "0.04%"), 4.3,
                        "gray", "green", 0.54, 0.44)

exp_trans <- gen_plot("trans", c(-4.523, -4, -3.523, -3, -2.523), 
                      c("0.003%", "0.01%", "0.03%", "0.10%","0.3%"),
                      4.52, "gray", "light blue", 0.54, 0.46)

exp_defense <- gen_plot("defense", c(-3, -2.699, -2.398, -2.097),
                        c("0.1%", "0.2%", "0.4%","0.8%"),
                        3, "gray", "red", 0.1, 0.1)

gen_plot("sect_CAZpep", c(-3, -2.699, -2.398, -2.097),
         c("0.1%", "0.2%", "0.4%","0.8%"),
         3, "gray", "red", 0.1, 0.1)


sel_cols = c("biofilm_exp_rate", "defense_exp_rate", "trans_exp_rate", "Layer_RNA")
melted <- DNA_RNA[, sel_cols] %>%
  `colnames<-`(c("Biofilm", "Defense mech.", "Transposase", "Depth")) %>%
  reshape2::melt(id.vars=c('Depth'),var='Gene')

library(ggpubr)
library(rstatix)
melted$Depth <- fct_rev(melted$Depth)
melted$log_exp <- log10(melted$value)

# three gene's expression looks very normal with log transform
# hist((filter(melted, Gene == "Biofilm"))$log_exp) 

stat.test <- melted %>%
  group_by(Gene) %>%
  t_test(log_exp~Depth) %>%
  add_significance("p") %>% 
  add_xy_position()
stat.test2 <- filter(stat.test, group1 =="MES")
stat.test2$y.position <- log10(c(3.2, 3.8, 1.9, 1.2, 5.6, 25)) # bio, def, trans

exp_summary <- melted %>% # 
  ggplot(aes(x=Depth, y=value, fill = Gene)) +
  scale_fill_manual(values = c("green", "red","light blue"))+
  facet_wrap(~ Gene, ncol = 3, scales = "free_x") +
  geom_boxplot(outlier.colour = "gray") +
  stat_pvalue_manual(stat.test2, label = "p.signif", tip.length = 0.01) +
  ylab("RNA/DNA ratio (%RNA / %DNA per sample)") +
  scale_y_log10() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  coord_flip() 

ggarrange(nrow = 2, ncol = 1,
  ggarrange(exp_biofilm, exp_defense, widths = c(0.5,0.46), labels = c("A","B"), ncol = 2),
  ggarrange(exp_trans, exp_summary, widths = c(0.4,0.6), labels = c("C","D"), ncol = 2))

ggsave("F7_DNA-RNA-ratio_target-genes.png", 
       plot = last_plot(),
       height = 4.3,
       width = 7.5)
ggsave("F7_DNA-RNA-ratio_target-genes.pdf", 
       plot = last_plot(),
       height = 4.3,
       width = 7.5)

for (g in c("biofilm", "defense", "trans")){
  DNA_exp <- paste(g,"exp_rate", sep="_")
  print(DNA_exp)
  dcm<-DNA_RNA%>%filter(Layer_RNA == "DCM")
  srf<-DNA_RNA%>%filter(Layer_RNA == "SRF")
  mes<-DNA_RNA%>%filter(Layer_RNA == "MES")
  srf_data <- less_than(srf[,c(DNA_exp)],out_exp)
  dcm_data <- less_than(dcm[,c(DNA_exp)],out_exp)
  mes_data <- less_than(mes[,c(DNA_exp)],out_exp)
  print(paste("SRF vs MES, P:", t.test(srf_data, mes_data)$p.value))
  print(paste("SRF vs DCM, P:", t.test(srf_data, dcm_data)$p.value))
}














# trashed 

gen_plot1 <- function(gene_name, log_scale, percent_scale, upscale, 
                      DNA_color, RNA_color, txt_top, txt_bot){
  log_scale <- log_scale + upscale
  DNA_gene <- paste("log_dna",gene_name, sep="_")
  RNA_gene <- paste("log_rna",gene_name, sep="_")
  tmp1 <- DNA_RNA[,c(DNA_gene, "upper_size_rna")]
  tmp2 <- DNA_RNA[,c(RNA_gene, "upper_size_rna")]
  tmp1$type <- "DNA"
  tmp2$type <- "RNA"
  
  tmp1 <- filter_outliers(tmp1, DNA_gene)
  tmp2 <- filter_outliers(tmp2, RNA_gene)
  
  g_cols <- c("log_abundance", "size_fraction", "type")
  DNA_RNA_bar <- rbind(`colnames<-`(tmp1, g_cols),
                       `colnames<-`(tmp2, g_cols))
  
  sum_bar <- DNA_RNA_bar %>% 
    group_by(type, size_fraction) %>% 
    summarise(median_abun = median(log_abundance) + upscale, 
              mean_abun = mean(log_abundance),
              lower = quantile(log_abundance, probs = c(.25)) + upscale, 
              upper = quantile(log_abundance, probs = c(.75)) + upscale)
  
  x_max <- max(max(DNA_RNA[,c(DNA_gene)]), max(DNA_RNA[,c(RNA_gene)]))+ upscale
  labD = "DNA"
  labR = "RNA"
  
  p <- ggplot(sum_bar, aes(y=fct_rev(size_fraction), x=median_abun, fill=fct_rev(type))) + 
    geom_bar(stat="identity", position=position_dodge()) +
    scale_x_continuous(breaks = log_scale, labels = percent_scale) +# 
    geom_errorbar(aes(xmin=lower, xmax=upper), width=.2,
                  position=position_dodge(.9)) + 
    # scale_fill_manual(values = c("DNA" = DNA_color, "RNA" = RNA_color))+
    # annotate(geom="text", x=x_max*txt_top, y=1.25, label=labD) + 
    # annotate(geom="text", x=x_max*txt_bot, y=1.8, label=labR) + 
    theme_classic() +
    
    theme(axis.title.y = element_blank()) # ,legend.position = "none"
  
  if (gene_name == "biofilm") {
    p <- p + xlab("% reads mapped to biofilm ORFs")
    return(p)
  } else if (gene_name == "trans") {
    return(p + xlab("% reads mapped to transposases"))
  } else {
    return(p + xlab("% reads mapped to defense mechanism ORFs"))
  }
}

exp_biofilm1 <- gen_plot1("biofilm", c(-4.301, -4, -3.698, -3.3979),
                          c("0.005%", "0.01%", "0.02%", "0.04%"), 4.3,
                          "gray", "green", 0.54, 0.44)

exp_trans1 <- gen_plot1("trans", c(-4.523, -4, -3.523, -3, -2.523), 
                        c("0.003%", "0.01%", "0.03%", "0.10%","0.3%"),
                        4.52, "gray", "light blue", 0.54, 0.46)

exp_defense1 <- gen_plot1("defense", c(-3, -2.699, -2.398, -2.097),
                          c("0.1%", "0.2%", "0.4%","0.8%"),
                          3, "gray", "red", 0.1, 0.1)

ggarrange(exp_biofilm1,exp_defense1,exp_trans1)

out_exp <- min(boxplot(c(DNA_RNA_tara$biofilm_exp_rate,DNA_RNA_tara$trans_exp_rate,DNA_RNA_tara$defense_exp_rate))$out)

out_b <- min(boxplot(DNA_RNA_tara$trans_exp_rate)$out)
out_d <- min(boxplot(DNA_RNA_tara$defense_exp_rate)$out)
out_t <- min(boxplot(DNA_RNA_tara$biofilm_exp_rate)$out)

cal1 <- function(vec) { 
  out_exp <- min(c(boxplot(vec)$out, 1000)) # in case no outliers
  return(signif(mean(less_than(vec, out_exp)),3))
}
cal2 <- function(vec) { 
  out_exp <- min(c(boxplot(vec)$out, 1000)) # in case no outliers
  return(signif(sd(less_than(vec, out_exp)),3))
}

DNA_RNA_tara[,c("biofilm_exp_rate", "trans_exp_rate", "defense_exp_rate", "Layer_RNA")] %>%
  group_by(Layer_RNA) %>%
  summarise(Biofilm = cal1(biofilm_exp_rate),
            Defense = cal1(defense_exp_rate),
            Transposase = cal1(trans_exp_rate),
            B_sd = cal2(biofilm_exp_rate),
            D_sd = cal2(defense_exp_rate),
            T_sd = cal2(trans_exp_rate),
            B_median = median(biofilm_exp_rate),
            D_median = median(defense_exp_rate),
            T_median = median(trans_exp_rate)) %>%
  as.data.frame(row.names = FALSE)

library(png)
img <- readPNG("./F7_exp_ratio_table.png")
exp_table_img <- ggplot() + background_image(img) + 
  theme(plot.margin = margin(t=0, l=0.05, r=0.05, b=0, unit = "cm"))


colnames(exp_table)[1] <- "Depth"

exp_table <- DNA_RNA_tara[,c("biofilm_exp_rate", "trans_exp_rate", "defense_exp_rate", "Layer_RNA")] %>%
  group_by(Layer_RNA) %>%
  summarise(`   Biofilm   ` = signif(median(biofilm_exp_rate, out_exp),3),
            `  Defense  ` = signif(median(defense_exp_rate, out_exp),3),
            Transposase = signif(median(trans_exp_rate, out_exp),3)) %>%
  as.data.frame(row.names = FALSE)

library(condformat)
exp_table_img <- condformat(exp_table) %>%
  rule_fill_bar(`   Biofilm   `, limits = c(0, NA), low = "white", high = "green") %>% 
  rule_fill_bar(`  Defense  `, limits = c(0, NA), low = "white", high = "red") %>% 
  rule_fill_bar(Transposase, limits = c(0, NA), low = "white", high = "light blue") %>%
  theme_grob(rows = NULL) %>%
  theme_caption(caption = "My Caption") %>%
  condformat2grob()

B_text <- text_grob("          Mean expression ratio (%RNA / %DNA)\n  of target genes in each sample")
ggarrange(
  ggarrange(pCassette, annotate_figure(exp_table_img, top = B_text), 
            ncol = 1, heights = c(0.7, 0.3), labels = c("A", "B")),
  ggarrange(exp_biofilm, exp_defense, exp_trans, 
            ncol =1, labels = c("C", "D", "E")),
  widths = c(0.5, 0.4),
  ncol = 2
)
