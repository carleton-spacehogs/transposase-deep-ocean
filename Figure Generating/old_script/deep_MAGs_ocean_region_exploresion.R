

mala_bin_abun = read_csv("data/Malaspina-Combined.RPKM.csv")
mala_metadata = read_csv("data/Malaspina-genes-coverage.csv")[c("sample","Ocean")]
names(mala_bin_abun) <- sub("_0.2", "", sub("_0.8", "", names(mala_bin_abun)))

tmp=cbind(mala_bin_abun["bin"], mala_bin_abun[grep("SRR", names(mala_bin_abun))]) %>%
  remove_rownames %>% column_to_rownames(var="bin")
md=mala_metadata %>% remove_rownames %>% column_to_rownames(var="sample")

mala_bin_abun = merge(t(tmp), md, by = 0) %>% remove_rownames %>%
  column_to_rownames(var="Row.names")

out = as.data.frame(t(aggregate(. ~ Ocean, mala_bin_abun, sum)))
names(out) <- out[1,]
out <- out[-1,]
out$Atlantic = as.numeric(out$Atlantic)
out$Indian = as.numeric(out$Indian)
out$Pacific = as.numeric(out$Pacific)

out = out %>%
  # filter(Atlantic+Indian+Pacific > 5) %>%
  mutate(origin = ifelse(Atlantic > (Indian + Pacific)*2, "Atlantic",
                  ifelse(Indian > (Atlantic + Pacific)*2, "Indian",
                  ifelse(Pacific > (Atlantic + Indian)*2, "Pacific", "undecided"))))
