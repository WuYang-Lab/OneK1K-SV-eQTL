#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: add_qvalues.R <cis_qtl.txt.gz> <output.txt.gz>")
}

input <- args[[1]]
output <- args[[2]]
if (!file.exists(input)) {
  stop("Input file not found: ", input)
}

dt <- fread(input)
if (!"pval_beta" %in% names(dt)) {
  stop("Input does not contain pval_beta: ", input)
}

dt[, pval_beta := as.numeric(pval_beta)]
valid <- which(is.finite(dt$pval_beta) & dt$pval_beta >= 0 & dt$pval_beta <= 1)
dt[, qval := NA_real_]

method <- "BH"
if (length(valid) > 0L && requireNamespace("qvalue", quietly = TRUE)) {
  qvalues <- tryCatch(
    qvalue::qvalue(dt$pval_beta[valid])$qvalues,
    error = function(e) NULL
  )
  if (!is.null(qvalues)) {
    dt$qval[valid] <- qvalues
    method <- "Storey qvalue"
  } else {
    dt$qval[valid] <- p.adjust(dt$pval_beta[valid], method = "BH")
  }
} else if (length(valid) > 0L) {
  dt$qval[valid] <- p.adjust(dt$pval_beta[valid], method = "BH")
}

fwrite(dt, output, sep = "\t")
message("Added qval using ", method, "; output: ", output)
