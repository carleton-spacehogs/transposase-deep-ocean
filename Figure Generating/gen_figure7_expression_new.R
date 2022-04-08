setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()
options(scipen=10000)

DNA_RNA_tara <- DNA_RNA_tara %>%filter(Layer_DNA %in% c("SRF", "DCM", "MES"))
scale <- c(-4.5, -4, -3.5, -3, -2.602, -2)
percent_scale <- c("0.003%", "0.01%", "0.03%", "0.10%","0.25%","1.00%")
def_scale <- c(-2.8239, -2.673389, -2.522879, -2.221849, -2.372364, -2.09691)
def_percent_scale <- c("0.15%", "0.21%", "0.30%", "0.60%", "0.42%", "0.80%")

gen_plot <- function(gene_name, y_text, log_scale, percent_scale, upscale, 
                     DNA_color, RNA_color, txt_top, txt_bot){
  log_scale <- log_scale + upscale
  DNA_gene <- paste("log_dna",gene_name, sep="_")
  RNA_gene <- paste("log_rna",gene_name, sep="_")
  tmp1 <- DNA_RNA_tara[,c(DNA_gene, "Layer_RNA")]
  tmp2 <- DNA_RNA_tara[,c(RNA_gene, "Layer_RNA")]
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
    lower = quantile(log_abundance, probs = c(.25)) + upscale, 
    upper = quantile(log_abundance, probs = c(.75)) + upscale)
    # summarise(sd = sd(log_abundance),
    #  mean_abun = mean(log_abundance) + upscale, 
    #  upper = mean_abun + sd,
    #  lower = mean_abun - sd)
  
  x_max <- max(max(DNA_RNA_tara[,c(DNA_gene)]), max(DNA_RNA_tara[,c(RNA_gene)]))+ upscale
  labD = "DNA"
  labR = "RNA"
  # if (gene_name == "biofilm"){
  #   labD = "DNA (metagenomes)"
  #   labR = "RNA (metatranscriptomes)"
  # }
  
  ggplot(sum_bar, aes(y=fct_rev(depth), x=median_abun, fill=fct_rev(type))) + 
    geom_bar(stat="identity", position=position_dodge()) +
    scale_x_continuous(breaks = log_scale, labels = percent_scale) + 
    geom_errorbar(aes(xmin=lower, xmax=upper), width=.2,
                  position=position_dodge(.9)) + 
    scale_fill_manual(values = c("DNA" = DNA_color, "RNA" = RNA_color))+
    xlab(y_text) +
    annotate(geom="text", x=x_max*txt_top, y=3.25, label=labD) + 
    annotate(geom="text", x=x_max*txt_bot, y=2.8, label=labR) + 
    theme_classic() +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
}

y1 <- "% reads mapped to biofilm ORFs"
y2 <- paste("% reads mapped to defense mech. ORFs",strrep(" ",4))
# y3 <- paste("% DNA/RNA reads mapped to transposases", strrep(" ",4)) # , strrep(" ",8)
y3 <- "% reads mapped to transposases"

exp_biofilm <- gen_plot("biofilm", y1, scale, percent_scale, 
                        4.5, "gray", "green", 0.5, 0.48)
exp_defense <- gen_plot("defense", y2, def_scale, def_percent_scale, 
                        2.85, "gray", "red", 0.7, 0.4)
exp_trans <- gen_plot("trans", y3, scale, percent_scale, 
                      4.5, "gray", "light blue", 0.5, 0.38)

out_exp <- min(boxplot(c(DNA_RNA_tara$biofilm_exp_rate,DNA_RNA_tara$trans_exp_rate,DNA_RNA_tara$defense_exp_rate))$out)

exp_table <- DNA_RNA_tara[,c("biofilm_exp_rate", "trans_exp_rate", "defense_exp_rate", "Layer_RNA")] %>%
  group_by(Layer_RNA) %>%
  summarise(`   Biofilm   ` = signif(mean(less_than(biofilm_exp_rate, out_exp)),3),
            `  Defense  ` = signif(mean(less_than(defense_exp_rate, out_exp)),3),
            Transposase = signif(mean(less_than(trans_exp_rate, out_exp)),3)) %>%
  as.data.frame(row.names = FALSE)

library(png)
img <- readPNG("./F7_exp_ratio_table.png")
exp_table_img <- ggplot() + background_image(img) + 
  theme(plot.margin = margin(t=0, l=0.05, r=0.05, b=0, unit = "cm"))

ggarrange(exp_biofilm, exp_defense,exp_trans, exp_table_img,
          labels = c("A","B","C","D"), ncol = 1, heights = c(1, 1, 1, 0.63))

ggsave("F7_expression_bars.png", plot = last_plot())

for (g in c("biofilm", "defense", "trans")){
  DNA_exp <- paste(g,"exp_rate", sep="_")
  print(DNA_exp)
  dcm<-DNA_RNA_tara%>%filter(Layer_RNA == "DCM")
  srf<-DNA_RNA_tara%>%filter(Layer_RNA == "SRF")
  mes<-DNA_RNA_tara%>%filter(Layer_RNA == "MES")
  srf_data <- less_than(srf[,c(DNA_exp)],out_exp)
  dcm_data <- less_than(dcm[,c(DNA_exp)],out_exp)
  mes_data <- less_than(mes[,c(DNA_exp)],out_exp)
  print(paste("SRF vs MES, P:", t.test(srf_data, mes_data)$p.value))
  print(paste("SRF vs DCM, P:", t.test(srf_data, dcm_data)$p.value))
}


# trashed 

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
