setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
g(bin_taxon, depth_comparison, malaspina_bins) %=% init_bins()
transposase_in_bins <- init_transposase_in_bins()

col_list <- c("complete genome size (Mbp)", "percent_trans", "log_percent_trans", 
              "percent_biofilm", "is_biofilm", "size_fraction",
              "depth", "Class", "class_trans", "graphing_log_trans", "graphing_log_biofilm")
all_tara <- bin_taxon[,col_list] %>%
  filter(!is.na(depth))

all_mala <- malaspina_bins[,col_list] %>%
  filter(size_fraction != "error")

all <- rbind(all_tara, all_mala) %>%
  mutate(biofilm_bin=cut(percent_biofilm, breaks = c(-1,0.001,0.2,5)))%>%
  mutate(size_bin=cut(`complete genome size (Mbp)`, breaks = c(0,2,3,4,5,20)))

fact <- c("None", '<0-0.2%', '>0.2%')
all$biofilm_bin <- as.character(all$biofilm_bin)
all$biofilm_bin[all$biofilm_bin == '(-1,0.001]'] <- fact[1]
all$biofilm_bin[all$biofilm_bin == '(0.001,0.2]'] <- fact[2]
all$biofilm_bin[all$biofilm_bin == '(0.2,5]'] <- fact[3]
all$biofilm_bin <- factor(all$biofilm_bin, levels = fact)

all$size_bin <- as.character(all$size_bin)
all$size_bin[all$size_bin == '(0,2]'] <- '0-2'
all$size_bin[all$size_bin == '(2,3]'] <- '2-3'
all$size_bin[all$size_bin == '(3,4]'] <- '3-4'
all$size_bin[all$size_bin == '(4,5]'] <- '4-5'
all$size_bin[all$size_bin == '(5,20]'] <- '>5'
all$size_bin <- factor(all$size_bin, levels = c(">5","4-5","3-4","2-3","0-2"))


all %>% 
  ggplot(aes(x=biofilm_bin, y = percent_trans)) +
  ylim(c(0,0.8))+
  ## facet_wrap(~depth, nrow= 1) +
  geom_boxplot(outlier.alpha = 0.1) +
  xlab("%-biofilm ORF in MAG") +
  ylab("%-transposase") +
  stat_summary(fun.data = boxplot.give.n, geom = "text",
               position=position_nudge(x = 0, y = 0.2))+
  theme_classic() # +

p_top <- all %>% 
  mutate(size_fraction = factor(size_fraction, levels = ls_order)) %>% 
  ggplot(aes(x=size_fraction, y = percent_trans)) +
  ylim(c(0,0.8))+
  ## facet_wrap(~depth, nrow= 1) +
  geom_boxplot(outlier.alpha = 0.1) +
  xlab("MAG lifestyle") +
  ylab("%-transposase") +
  stat_summary(fun.data = boxplot.give.n, geom = "text",
               position=position_nudge(x = 0, y = 0.2))+
  theme_classic() # +
  # theme(axis.title.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.line.x = element_blank(),
  #      axis.ticks.x = element_blank())

all_dummy <- all
all_dummy$size_fraction <- "all"
ls_order <- c("planktonic", "mixed", "particle", "all")
p_mid <- all %>% 
  mutate(size_fraction = factor(size_fraction, levels = ls_order)) %>% 
  ggplot(aes(x=size_fraction, y = `complete genome size (Mbp)`)) +
  facet_wrap(~depth, nrow= 1) +
  scale_y_continuous(breaks = c(0, 2.5, 5, 7.5), limits = c(0, 8.2))+
  xlab("MAG lifestyle") +
  ylab("complete genome size")+
  stat_summary(fun.data = boxplot.give.n, geom = "text",
               position=position_nudge(x = 0, y = 1.5)) +
  geom_boxplot(outlier.alpha = 0.1)+
  theme_classic() +
  theme(strip.text.x = element_blank())

p_bot<- all %>% 
  filter(!is.na(depth)) %>%
  ggplot(aes(x=fct_rev(size_bin), y=percent_trans)) + 
  geom_boxplot(outlier.alpha = 0.1) +
  ylim(c(0, 0.8))+
  stat_summary(fun.data = boxplot.give.n, geom = "text", position=position_nudge(x = 0, y = 0.03)) + 
  facet_wrap(~depth, nrow= 1) +
  ylab("%-transposase") +
  xlab("complete genome size (Mbp)") +
  theme_classic() +
  theme(strip.text.x = element_blank())

ggarrange(p_top, p_mid, p_bot, ncol = 1, heights = c(3, 3, 3),
          labels = c("A", "B", "C"))



