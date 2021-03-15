#!/usr/bin/Rscript

suppressMessages({
  source("sleuth_methods.R")
})

num_cores <- 4
host <- "useast.ensembl.org"
gene_column <- 'ens_gene'

transcript_gene_mapping <- get_mouse_gene_names(host = host)

sample_file <- '../metadata/sample_list.txt'
exp_info_file <- '../metadata/exp_info.txt'
in_dir <- '../results/kallisto'
out_dir <- '../results/sleuth_kal'
out_prefix <- 'L2KD'
control_name <- "control"

sample <- read.table(sample_file, sep = "\t", stringsAsFactors = F,
                       col.names = c("sample"))
s2c <- read.table(file.path(exp_info_file), header = T, sep="\t",
                  stringsAsFactors = F)
s2c$condition <- relevel(as.factor(s2c$condition), control_name)

kal_dirs <- file.path(in_dir, sample[, 1])
kal_dirs <- kal_dirs[match(s2c$sample, sample[, 1])]
s2c$path <- kal_dirs

beta <- "conditionL2KD"


message("Running sleuth")
sleuth_res <- run_sleuth(s2c, max_bootstrap = 100, beta = beta,
               gene_column = NULL,
               num_cores = num_cores, which_var = "obs_tpm",
               read_bootstrap_tpm = TRUE)
output_sleuth(out_dir, out_prefix, sleuth_res)

message("Running sleuth-CN on the data using Actb as reference gene")
denom <- "ENSRNOT00000080216.1"
real_prefix <- paste(out_prefix, "cn_actb", sep = "_")
cn_res <- run_cn(s2c, max_bootstrap = 100, delta = 0.01, beta = beta,
                   denom = denom, gene_column = NULL,
                   num_cores = num_cores, which_var = "obs_tpm")
output_sleuth(out_dir, real_prefix, cn_res)

message("Running sleuth-CN on the data using Hprt as reference gene")
denom <- "ENSRNOT00000045153.3"
real_prefix <- paste(out_prefix, "cn_hprt1", sep = "_")
cn_res <- run_cn(s2c, max_bootstrap = 100, delta = 0.01, beta = beta,
                   denom = denom, gene_column = NULL,
                   num_cores = num_cores, which_var = "obs_tpm")
output_sleuth(out_dir, real_prefix, cn_res)

message("Running both")
denom <- c("ENSRNOT00000080216.1", "ENSRNOT00000045153.3")
real_prefix <- paste(out_prefix, "cn_both", sep = "_")
cn_res <- run_cn(s2c, max_bootstrap = 100, delta = 0.01, beta = beta,
                   denom = denom, gene_column = NULL,
                   num_cores = num_cores, which_var = "obs_tpm")
output_sleuth(out_dir, real_prefix, cn_res)
