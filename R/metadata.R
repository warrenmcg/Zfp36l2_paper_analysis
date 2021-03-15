#!/usr/bin/Rscript

sample_names <- list.files('../data', pattern = 'fastq.gz')
sample_names <- gsub('.fastq.gz', '', sample_names)
small_names <- simplify2array(strsplit(sample_names, '_'))[1,]

dir.create('../metadata', showWarnings = FALSE)
write.table(sample_names, '../metadata/sample_list.txt', row.names = F, col.names = F, quote = F, sep = "\t")

exp_df <- data.frame(sample = sample_names, alias = small_names)
exp_df$condition <- ifelse(grepl('Cont', sample_names), 'control', 'L2KD')
write.table(exp_df, '../metadata/exp_info.txt', row.names = F, quote = F, sep = "\t")
