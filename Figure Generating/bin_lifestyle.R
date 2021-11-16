setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparisons, malaspina_bins) %=% init_bins()

#pB from metagenome_pnps.R. pA is deprecated

biofilm<-select(bin_taxon, "lifestyle", "bin", `biofilm gene calls in genome (%)`)
trans<-select(bin_taxon, "lifestyle", "bin", `transposase gene calls in genome (%)`)
defense<-select(bin_taxon, "lifestyle", "bin", `defense mechanisms gene calls in genome (%)`)

set_colnames=c("lifestyle", "bin", "prop", "gene")
defense_scale <- 8
biofilm$gene <- "biofilm"
colnames(biofilm) <- set_colnames
trans$gene <- "transposase"
colnames(trans) <- set_colnames
defense$gene <- "defense"
colnames(defense) <- set_colnames
defense$prop <- defense$prop/defense_scale

to_graph <- rbind(trans, biofilm, defense)
to_graph$gene <- factor(to_graph$gene, levels = c("transposase", "biofilm", "defense"))
options(scipen=10000)
pC <- ggplot(to_graph, aes(x=`lifestyle`, y=prop, fill=gene)) + 
  geom_boxplot(outlier.color = "gray") + 
  #geom_violin(adjust = 5) +
  theme_classic() +
  scale_y_continuous(name = "% of transposases, biofilm in total ORFs", # first axis
    sec.axis = sec_axis(~.*defense_scale, name="% of defense mechanisms in total ORFs") # second axis
  ) +
  xlab("Lifestyles of MAGs") +
  scale_x_discrete(labels=c("free-living\n(n = 96)", "mix\n(n = 1081)", "particle-assoc\n(n = 130)")) +
  coord_flip(ylim= c(0,0.65)) +
  scale_fill_manual(breaks=c("defense","biofilm","transposase"), values=c('red','darkolivegreen1','cadetblue1'))

ggarrange(pB, pC, ncol = 2, labels = c("A", "B"), widths = c(0.5, 0.5), heights = c(0.9, 1))


particle_biofilm_prop <- bin_taxon[bin_taxon$lifestyle == "particle", "biofilm gene calls in genome (%)"]
particle_trans_prop <- bin_taxon[bin_taxon$lifestyle == "particle", "transposase gene calls in genome (%)"]
particle_defense_prop <- bin_taxon[bin_taxon$lifestyle == "particle", "defense mechanisms gene calls in genome (%)"]

free_living_biofilm_prop <- bin_taxon[bin_taxon$lifestyle == "free_living", "biofilm gene calls in genome (%)"]
free_living_trans_prop <- bin_taxon[bin_taxon$lifestyle == "free_living", "transposase gene calls in genome (%)"]
free_living_defense_prop <- bin_taxon[bin_taxon$lifestyle == "free_living", "defense mechanisms gene calls in genome (%)"]

t.test(particle_biofilm_prop, free_living_biofilm_prop)
t.test(particle_trans_prop, free_living_trans_prop)
t.test(particle_defense_prop, free_living_defense_prop)


bin_taxon%>%
  ggplot(aes(x=lifestyle,y=HTG_proportion)) +
  geom_boxplot() +
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0))

bin_taxon$log_HTG_prop <- log(bin_taxon$HTG_proportion)
compare_means(log_HTG_prop~lifestyle, data = bin_taxon, method = "t.test")

# trashed
ggarrange(ggarrange(pB, pC, ncol = 2, labels = c("A", "B"), widths = c(0.45, 0.55)), # Second row with box and dot plots
          pA,
          nrow = 2, 
          heights = c(0.5, 0.5),
          labels = c("", "    C"))
