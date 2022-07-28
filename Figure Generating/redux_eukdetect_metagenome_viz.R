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
euk_read_count02$filter_size <- "size0.2-0.8"
euk_read_count08 = sum_reads_into_row(deep_frac08, base08)
euk_read_count08$filter_size <- "size0.8-5.0"
deep_euk_count = rbind(euk_read_count02, euk_read_count08)
deep_euk_count$depth = "BAT"

base_big_SRF = "../identify_eukaryotes/EukDetect-SRF_5-20/filtering"
base_big_DCM = "../identify_eukaryotes/EukDetect-DCM_5-20/filtering"
big_size_SRF_f = Sys.glob(file.path(base_big_SRF, hit_table))
big_size_DCM_f = Sys.glob(file.path(base_big_DCM, hit_table))

big_size_SRF = sum_reads_into_row(big_size_SRF_f, base_big_SRF)
big_size_SRF$depth = "SRF"
big_size_DCM = sum_reads_into_row(big_size_DCM_f, base_big_DCM)
big_size_DCM$depth = "DCM"
big_size_count = rbind(big_size_SRF, big_size_DCM)
big_size_count$filter_size <- "size5-20"

# mala_cov = merge(mala_cov, deep_euk_count, by = "sample")
# mala_cov$euk_prop = (mala_cov$euk_read_count*150)/(mala_cov$Sequencing_depth_Gbp*1e+9)
# fudge_factor = min(mala_cov$euk_prop[mala_cov$euk_prop > 0])/2
# mala_cov$log_euk_prop = log10(mala_cov$euk_prop + fudge_factor)

base_tara="../identify_eukaryotes/EukDetect-normal/filtering"
base_arct="../identify_eukaryotes/EukDetect-arctic/filtering"
tara_normal = Sys.glob(file.path(base_tara, hit_table))
tara_arctic = Sys.glob(file.path(base_arct, hit_table))
euk_tara = rbind(sum_reads_into_row(tara_normal, base_tara),
                 sum_reads_into_row(tara_arctic, base_arct))
euk_tara$euk_read_count = as.numeric(euk_tara$euk_read_count)

sel_col = c("Layer_DNA","upper_size_dna")
match_md_row = function(x){
  md_row = DNA_Metadata[DNA_Metadata$ENA_Run_ID %like% x["sample"], sel_col]
  if (nrow(md_row) > 0) {
    df = c(x, unlist(md_row))
    return(df)
  } else {
    print(paste("cannot find", x["sample"]))
  }
}

merged_tara = as.data.frame(apply(euk_tara, 1, match_md_row))
merged_tara = as.data.frame(t(merged_tara)) %>%
  filter(sample != "cannot find ERR594346") %>%
  mutate(depth = Layer_DNA,
         filter_size = ifelse(upper_size_dna == "1.6", "size0.22-1.6", "size0.22-3.0"),
         euk_read_count = as.numeric(as.character(euk_read_count)))


sel_col2 = c("sample", "filter_size", "euk_read_count", "depth")
depth_cat = rbind(merged_tara[sel_col2], deep_euk_count[sel_col2], big_size_count[sel_col2])
total_reads = read_csv("../identify_eukaryotes/count-reads.out")
all_cat = merge(depth_cat, total_reads, by = "sample")
all_cat$euk_prop = all_cat$euk_read_count/all_cat$num_reads
fudge_all = min(all_cat$euk_prop[all_cat$euk_prop > 0])/2
all_cat$log_euk_prop = log10(all_cat$euk_prop + fudge_all)

all_cat %>%
  filter(depth != "MIX") %>%
  mutate(size2 = ifelse(filter_size %in% c("size0.22-1.6","size02-08"), 
                        "planktonic", "particle"),
         Depth = factor(depth, levels = c("BAT", "MES", "DCM", "SRF"))) %>%
  ggplot(aes(x = log_euk_prop, y = Depth, fill = filter_size)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0), alpha = 0.5) +
  xlab("Proportion of eukaryotic marker reads (log10)") +
  theme_classic()


