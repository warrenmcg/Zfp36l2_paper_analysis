.libPaths(c('~/R_library', .libPaths()))
library("sleuth")#, lib.loc = "~/test_R_library")
library("sleuthALR")

get_rat_gene_names <- function(host = "dec2016.archive.ensembl.org") {
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "rnorvegicus_gene_ensembl",
    host = host)
  ttg <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "transcript_version",
    "ensembl_gene_id", "external_gene_name", "version", "transcript_length"),
    mart = mart)
  ttg <- dplyr::rename(ttg, ext_gene = external_gene_name, length = transcript_length)
  ttg <- dplyr::mutate(ttg, ens_gene = paste(ensembl_gene_id, version, sep="."),
    target_id = paste(ensembl_transcript_id, transcript_version, sep="."))
  ttg <- dplyr::select(ttg, target_id, ens_gene, ext_gene, length)
  ttg
}

output_sleuth <- function(out_dir, out_prefix, sleuth_res) {
  so <- sleuth_res$so
  s2c <- so$sample_to_covariates

  message("saving the sleuth object")
  so_file <- file.path(out_dir, paste0(out_prefix, '.rds'))
  sleuth_save(so, so_file)

  message("printing out results and plots")
  wt <- sleuth_results(so, beta, pval_aggregate = FALSE)
  wt_file <- file.path(out_dir, paste0(out_prefix, '_wtAllResults.txt'))
  wt_sig_file <- file.path(out_dir, paste0(out_prefix, '_wtSigResults.txt'))
  write.table(wt, file = wt_file, row.names=F, quote=F, sep="\t")
  write.table(wt[which(wt$qval <= 0.05), ], file = wt_sig_file, row.names=F, quote=F, sep="\t")

  lrt <- sleuth_results(so, "reduced:full", test_type = "lrt", pval_aggregate = FALSE)
  lrt_file <- file.path(out_dir, paste0(out_prefix, '_lrtAllResults.txt'))
  lrt_sig_file <- file.path(out_dir, paste0(out_prefix, '_lrtSigResults.txt'))
  write.table(lrt, file = lrt_file, row.names=F, quote=F, sep="\t")
  write.table(lrt[which(lrt$qval <= 0.05), ], file = lrt_sig_file, row.names=F, quote=F, sep="\t")

  grDevices::pdf(paste(out_dir,"/", out_prefix, "_pcaClusterPlot.pdf", sep=""))
  print(plot_pca(so, text_labels=T, color_by="condition", units = "tpm"))
  grDevices::dev.off()

  grDevices::pdf(paste(out_dir,"/", out_prefix, "_wtPvalHist.pdf", sep=""))
  graphics::hist(wt$pval, xlab="p-values",
       main = paste(levels(s2c$condition)[1],"vs",levels(s2c$condition)[2]))
  grDevices::dev.off()

  grDevices::pdf(paste(out_dir,"/", out_prefix, "_lrtPvalHist.pdf", sep=""))
  graphics::hist(lrt$pval, xlab="p-values",
       main = paste(levels(s2c$condition)[1],"vs",levels(s2c$condition)[2]))
  grDevices::dev.off()

  invisible(NULL)
}

run_sleuth_prep <- function(sample_info, max_bootstrap = 30, gene_column = NULL,
  filter_target_id = NULL, ...) {
  so <- sleuth_prep(sample_info, ~ condition, max_bootstrap = max_bootstrap,
    target_mapping = transcript_gene_mapping, filter_target_id = filter_target_id,
    ...)
  so
}

run_cn <- function(sample_info,
  max_bootstrap = 30,
  gene_column = NULL,
  denom = NULL,
  filter_target_id = NULL,
  delta = NULL,
  impute_proportion = 0.65,
  which_var = 'obs_tpm',
  method = "multiplicative",
  shrink_fun = sleuth::basic_shrink_fun,
  run_ash = FALSE,
  beta = "conditionB",
  ...) {

  if (length(denom) > 1) {
    lr_type <- 'alr'
  } else if (denom == 'best') {
    denom_var <- ifelse(which_var == "obs_tpm", "tpm", "est_counts")
    denom <- sleuthALR::choose_denom(sample_info = sample_info, num_cores = 1, which_var = denom_var,
                                     filter_length = TRUE, target_mapping = transcript_gene_mapping)
    lr_type <- 'alr'
  } else if (denom == 'iqlr') {
    lr_type <- 'iqlr'
  } else if (denom == 'all' | denom == 'clr') {
    lr_type <- 'clr'
  } else {
    lr_type <- 'alr'
  }

  so <- sleuthALR::make_lr_sleuth_object(sample_info,
    target_mapping = transcript_gene_mapping,
    beta = beta,
    denom_name = denom, aggregate_column = gene_column,
    max_bootstrap = max_bootstrap,
    filter_target_id = filter_target_id,
    lr_type = lr_type,
    impute_method = method,
    delta = delta,
    impute_proportion = impute_proportion,
    run_models = FALSE,
    ...)

  so <- sleuth_fit(so, ~ condition, 'full', which_var = which_var, shrink_fun = shrink_fun)
  so <- sleuth_wt(so, beta)
  so <- sleuth_fit(so, ~ 1, 'reduced', which_var = which_var, shrink_fun = shrink_fun)
  so <- sleuth_lrt(so, 'reduced', 'full')

  denom_names <- sleuthALR::get_denom_names(so)
  lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt',
    show_all = FALSE)[, c('target_id', 'pval', 'qval', 'test_stat')]
  wt <- sleuth_results(so, beta,
    show_all = FALSE)[, c('target_id', 'pval', 'qval', 'b', 'se_b')]

  if (run_ash) {
    message('running ash to get posterior beta estimates')
    # sink is necessary to hide a whole bunch of cat output from ash
    sink(tempfile())
    df <- so$fits$full$models$df.residual
    ash_res <- ashr::ash(wt$b, wt$se_b, df = df, method = "shrink")
    res <- ash_res$result
    new_wald <- res$PosteriorMean / res$PosteriorSD
    new_pval <- 2*pnorm(abs(new_wald), lower.tail = FALSE)
    new_qval <- p.adjust(new_pval, method = "BH")
    wt$pval <- new_pval
    wt$qval <- new_qval
    wt$b <- res$PosteriorMean
    sink()
  }

  wt$se_b <- NULL
  res <- list(sleuthALR.lrt = lrt, sleuthALR.wt = wt, denoms = denom_names, so = so)
  res
}

run_sleuth <- function(sample_info,
  max_bootstrap = 30,
  gene_column = NULL,
  filter_target_id = NULL,
  shrink_fun = sleuth::basic_shrink_fun,
  run_ash = FALSE,
  beta = "conditionB",
  ...) {


  so <- NULL
  if (!is.null(gene_column)) {
    so <- run_sleuth_prep(sample_info, max_bootstrap = max_bootstrap,
      aggregation_column = gene_column, filter_target_id = filter_target_id,
      ...)
  } else {
    so <- run_sleuth_prep(sample_info, max_bootstrap = max_bootstrap,
      filter_target_id = filter_target_id, ...)
  }

  so <- sleuth_fit(so, so$full_formula, 'full', shrink_fun = shrink_fun)
  so <- sleuth_wt(so, beta)
  so <- sleuth_fit(so, ~ 1, 'reduced', shrink_fun = shrink_fun)
  so <- sleuth_lrt(so, 'reduced', 'full')

  res <- NULL
  lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt',
    show_all = FALSE)[, c('target_id', 'pval', 'qval', 'test_stat')]
  wt <- sleuth_results(so, beta,
    show_all = FALSE)[, c('target_id', 'pval', 'qval', 'b', 'se_b')]

  if (run_ash) {
    message('running ash to get posterior beta estimates')
    # sink is necessary to hide a whole bunch of cat output from ash
    sink(tempfile())
    df <- so$fits$full$models$df.residual
    ash_res <- ashr::ash(wt$b, wt$se_b, df = df, method = "shrink")
    res <- ash_res$res
    new_wald <- res$PosteriorMean / res$PosteriorSD
    new_pval <- 2*pnorm(abs(new_wald), lower.tail = FALSE)
    new_qval <- p.adjust(new_pval, method = "BH")
    wt$pval <- new_pval
    wt$qval <- new_qval
    wt$b <- res$PosteriorMean
    sink()
  }

  wt$se_b <- NULL
  res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  res$so <- so

  res
}
