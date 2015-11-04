# loadMothurTaxSummary: load Mothur tax.summary file into data frame with columns taxa, rows samples
loadMothurTaxSummary = function(taxSummaryFile, excludeSamples = c()) {
  counts = read.table(taxSummaryFile, header = TRUE)
  rownames(counts) = make.names(counts$taxon, unique = TRUE)
  # get rid of excess columns in mothur tax.summary file
  counts = counts[,-grep("taxon|taxlevel|rankID|daughterlevels|total", colnames(counts))]
  counts = counts[,!colnames(counts) %in% excludeSamples]
  counts = counts[!rownames(counts) %in% c("Root"),]
  as.data.frame(t(counts))
}
