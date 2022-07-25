setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source("./init_share.R")
# g(DNA_tara, RNA_tara, DNA_RNA_tara, depth_comparison) %=% init_tara()
mala_cov <- init_mala_cov()
g(DNA_Metadata, RNA_Metadata) %=% init_tara_md(factorize = FALSE)

sum_reads_into_row = function(hit_table_f, base_path){
  sum_reads_per_file = function(singel_hit_table, base_path){
    sample_name = str_replace(singel_hit_table, paste(base_path,"/",sep = ""), "")
    sample_name = str_replace(sample_name, "_all_hits_table.txt", "")
    if (file.info(singel_hit_table)$size < 60) { # no reads mapped
      print(paste("empty table:",singel_hit_table))
      return(c(sample_name, 0))
    } else {
      ht = suppressMessages(read_delim(singel_hit_table, delim = "\t",
                                       escape_double = FALSE, trim_ws = TRUE))
      return(c(sample_name, sum(ht$Read_counts)))
    }
    return(df)
  }
  df = mapply(sum_reads_per_file, hit_table_f, base_path)
  df = t(as.data.frame(df))
  rownames(df) <- 1:nrow(df)
  colnames(df) = c("sample","euk_read_count")
  df = as.data.frame(df)
  df = transform(df, euk_read_count = as.numeric(euk_read_count))
  return(df)
}

hit_table = "*_all_hits_table.txt"
base02="../identify_eukaryotes/EukDetect-deep/02_size_fraction/filtering"
base08="../identify_eukaryotes/EukDetect-deep/08_size_fraction/filtering"
deep_frac02 = Sys.glob(file.path(base02, hit_table))
deep_frac08 = Sys.glob(file.path(base08, hit_table))

euk_read_count02 = sum_reads_into_row(deep_frac02, base02)
euk_read_count02$filter_size <- "size02-08"
euk_read_count08 = sum_reads_into_row(deep_frac08, base08)
euk_read_count08$filter_size <- "size08-05"
deep_euk_count = rbind(euk_read_count02, euk_read_count08)

mala_cov = merge(mala_cov, deep_euk_count, by = "sample")
mala_cov$euk_prop = (mala_cov$euk_read_count*150)/(mala_cov$Sequencing_depth_Gbp*1e+9)
fudge_factor = min(mala_cov$euk_prop[mala_cov$euk_prop > 0])/2
mala_cov$log_euk_prop = log10(mala_cov$euk_prop + fudge_factor)

base_tara="../identify_eukaryotes/EukDetect-normal/filtering"
base_arct="../identify_eukaryotes/EukDetect-arctic/filtering"
tara_normal = Sys.glob(file.path(base_tara, hit_table))
tara_arctic = Sys.glob(file.path(base_arct, hit_table))
euk_tara = rbind(sum_reads_into_row(tara_normal, base_tara),
                 sum_reads_into_row(tara_arctic, base_arct))

sel_col = c("Layer_DNA","upper_size_dna","Raw reads","Ocean_DNA","ENA_Run_ID")
match_md_row = function(x){
  print(x)
  md_row = DNA_Metadata[DNA_Metadata$ENA_Run_ID %like% x["sample"], sel_col]
  multi_factor = as.numeric(str_count(md_row$ENA_Run_ID, "ERR"))
  df = cbind(x["sample"], as.numeric(x["euk_read_count"])*multi_factor, md_row)
  return(df)
}

apply(euk_tara, 1, match_md_row)


tara_merged = as.data.frame(t(apply(euk_tara, 1, match_md_row)))
tara_merged = transform(tara_merged, euk_count = as.numeric(euk_count), 
                        total_reads = as.numeric(`Raw reads`)) %>%
  mutate(filter_size = ifelse(upper_size_dna == "1.6", "size0.22-1.6", "size0.22-3.0"))
tara_merged$euk_prop = tara_merged$euk_count/tara_merged$total_reads
tara_merged$log_euk_prop = log10(tara_merged$euk_prop)

sel_col2 = c("filter_size", "log_euk_prop","euk_prop","Ocean_DNA", "Layer_DNA")
depth_cat = rbind(tara_merged[sel_col2], mala_cov[sel_col2])
depth_cat$Ocean = unlist(depth_cat$Ocean_DNA)
depth_cat$Layer_DNA = unlist(depth_cat$Layer_DNA)
boxplot(log_euk_prop ~ filter_size, depth_cat)
boxplot(log_euk_prop ~ Ocean, depth_cat, las = 3)

depth_cat %>%
  filter(Layer_DNA != "MIX") %>%
  mutate(size2 = ifelse(filter_size %in% c("size0.22-1.6","size02-08"), 
                        "planktonic", "particle"),
         Depth = factor(Layer_DNA, levels = c("BAT", "MES", "DCM", "SRF"))) %>%
  # group_by(Layer_DNA, size2) %>% count()
  ggplot(aes(x = log_euk_prop, y = Depth)) +
  geom_boxplot(aes(fill = filter_size)) +
  # facet_wrap(~size2, ncol = 1)+
  geom_point(aes(color = Ocean)) +
  xlab("Proportion of eukaryotic reads (log10)")

mala_cov %>%
  filter(filter_size == "size02-08") %>%
  ggplot(aes(y=log_dna_defense, x=percent_sect_CAZ)) +
  geom_point(aes(color = filter_size)) +
  geom_smooth(method = "lm", se = F, color = "orange")



