if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    stop('Usage: Rscript r_pkg_install.R R_DIR')
  }
  r_dir <- args[1]
} else {
  stopifnot(!is.null(r_dir))
}

if (!file.exists('r_pkg_install_success.txt')) {
  if(compareVersion(as.character(getRversion()), '3.5.0') < 0) {
    source("https://bioconductor.org/biocLite.R")
    ## Install Bioconductor packages not available on bioconda
    biocLite(c('SRAdb'))
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")    
    BiocManager::install(c('SRAdb'), type = 'source', checkBuilt = TRUE)
  }

  ## place non-standard packages into a separate directory
  dir.create(r_dir, showWarnings = F)
  .libPaths(c(r_dir, .libPaths()))

  ## Install modified repos
  devtools::install_github('warrenmcg/sleuth', ref = 'speedy_fit')

  ## Install other github packages
  devtools::install_github('warrenmcg/sleuth-CN')

  message <- "all packages have been installed correctly"
  sink('r_pkg_install_success.txt')
  cat(message)
  sink()
} else {
  message('It appears that all the R packages have been installed')
  message('There is nothing to do')
}
