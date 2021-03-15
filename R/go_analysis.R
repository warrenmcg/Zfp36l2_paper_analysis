.libPaths(c('~/R_library', .libPaths()))

suppressMessages(library(topGO))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_isoform.R INPUT OUT_PREFIX')
}

species <- 'Mouse'
gene_file <- args[1]
out_prefix <- args[2]

gene_data <- read.table(gene_file, sep = "\t", quote = "", header = T, stringsAsFactors = F)
gene_data <- gene_data[!is.na(gene_data$pval), ]
all_genes <- factor(as.integer(gene_data$qval <= 0.05))
names(all_genes) <- substr(gene_data$target_id, 1, 18)
GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes,
               annot = annFUN.org, mapping = "org.Rn.eg.db",
               ID = "ensembl", nodeSize = 5)
weight01.fisher <- runTest(GOdata, statistic = "fisher")
fisher_scores <- score(weight01.fisher)
goterms <- Term(GOTERM)
go_res <- data.frame(go_id = names(fisher_scores), go_term = goterms[names(fisher_scores)],
                     pval = fisher_scores)

mart <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
  host = 'useast.ensembl.org',
  dataset = 'rnorvegicus_gene_ensembl')
annos <- biomaRt::getBM(mart = mart,
  attributes = c('ensembl_gene_id', 'external_gene_name'))
gene_names <- annos$external_gene_name
names(gene_names) <- annos$ensembl_gene_id
rm(annos)

nodes <- names(GOdata@graph@nodeData@data)
scores <- GOdata@allScores
names(scores) <- GOdata@allGenes
gene_list <- lapply(nodes, function(id) {
  genes <- ls(GOdata@graph@nodeData@data[[id]]$genes)
  scores <- all_genes[genes]
  sig_genes <- names(scores[which(scores == 1)])
  names <- gene_names[genes]
  sig_names <- gene_names[sig_genes]
  string <- paste(names, collapse = ", ")
  sig_string <- paste(sig_names, collapse = ", ")
  res <- data.frame(go_id = id,
    size = length(names),
    significant = length(sig_names),
    sig_genes = sig_string,
    all_genes = string,
    stringsAsFactors = F)
  res
})
gene_df <- dplyr::bind_rows(gene_list)
go_res <- merge(go_res, gene_df, by = "go_id")
sig_go_res <- go_res[which(go_res$pval <= 0.05), ]

write.table(go_res, file = paste(out_prefix, "all_go_res.txt", sep = "_"), sep = "\t", row.names = F, quote = F)
write.table(sig_go_res, file = paste(out_prefix, 'sig_go_res.txt', sep = "_"), sep = "\t", row.names = F, quote = F)

gene_go_list <- lapply(nodes, function(id) {
  genes <- ls(GOdata@graph@nodeData@data[[id]]$genes)
  go_name <- goterms[id]
  names(go_name) <- NULL
  scores <- all_genes[genes]
  symbols <- gene_names[genes]
  names(symbols) <- NULL
  sig_genes <- names(scores[which(scores == 1)])
  data.frame(go_id = id, go_name = go_name,
    ens_gene_id = genes,
    gene_symbol = symbols,
    significant = genes %in% sig_genes,
    stringsAsFactors = F)
})
gene_go_df <- dplyr::bind_rows(gene_go_list)
write.table(gene_go_df, file = paste(out_prefix, "long_go_res.txt", sep = "_"),
  sep = "\t", row.names = F, quote = F)
