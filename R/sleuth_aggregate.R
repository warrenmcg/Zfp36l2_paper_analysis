.libPaths(c('~/R_library', .libPaths()))

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
so_file <- args[1]
out_prefix <- args[2]

beta <- 'conditionL2KD'

kal_so <- readRDS(so_file)
kal_so$gene_column <- "ens_gene"
kal_so$target_mapping <- dplyr::select(kal_so$target_mapping,
                                       -length)
kal_res <- sleuth::sleuth_results(kal_so, beta, pval_aggregate = T, weight_func = exp)
kal_res$sum_mean_obs_counts <- log2(exp(kal_res$sum_mean_obs_counts))
kal_res <- dplyr::rename(kal_res, log2_mean_observed_tpms = sum_mean_obs_counts)

# This identifies whether a significant gene can be described as overall
# up-regulated, down-regulated, or "mixed"
# If it does not have any significantly changing transcripts
# (i.e. "num_sig_trans" = 0), then "dir" = NA
trans_res <- sleuth::sleuth_results(kal_so, beta)
gene_res <- trans_res %>% group_by(ens_gene) %>%
  filter(qval<=0.05, .preserve = TRUE) %>%
  summarise(num_sig_trans = n(),
            dir = ifelse(n() == 0, NA,
                         ifelse(all(b > 0), "up",
                                ifelse(all(b < 0), "down", "mixed"))),
            .groups = "keep")

# Merge the number of sig. transcripts & direction into the main
# results table
kal_res <- merge(kal_res, gene_res, by = "ens_gene", all.x = TRUE)
kal_res$num_sig_trans[is.na(kal_res$num_sig_trans)] <- 0

out_file <- paste0(out_prefix, "_gene_aggregate_allResults.txt")
write.table(kal_res, out_file,
            sep = "\t", quote = F, row.names = F)
sig_kal_res <- kal_res[which(kal_res$qval <= 0.05), ]
out_file <- paste0(out_prefix, "_gene_aggregate_sigResults.txt")
write.table(sig_kal_res, out_file,
            sep = "\t", quote = F, row.names = F)
