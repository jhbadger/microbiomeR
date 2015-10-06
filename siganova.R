
# siganova: return data frame of taxa significantly abundant in respect to a group.
# inputs: tax.summary file from Mothur, numeric taxonomic level (1=kingdom,2=phylum, etc), group to check
# optional arrays of samples to exclude (e.g. Mock) and taxa to exclude (e.g. unclassified)
# optional normalization to 100000 reads/sample (default)
# optional correction method (default "fdr")

siganova = function(taxSummaryFile, taxLevel, group, excludeSamples=c(), excludeTaxa=c(), normalize=TRUE, correction = "fdr") {
  counts = read.table(taxSummaryFile, header = TRUE)
  counts = counts[counts[1]==taxLevel,]
  rownames(counts) = make.names(counts$taxon, unique=TRUE)
  counts = counts[,-grep("taxon|taxlevel|rankID|daughterlevels|total", colnames(counts))]
  counts = counts[,!colnames(counts) %in% excludeSamples]
  counts = counts[!rownames(counts) %in% excludeTaxa,]
  if (normalize) {
    counts = round((counts/colSums(counts))*100000)
  }
  if (length(group) == 0) { # no group provided
    group = colnames(counts)
  }
  aof=function(x){anova(aov(x~group))}
  anovaresults = apply(t(counts), 2, aof)
  pvalues=t(data.frame(lapply(anovaresults, function(x) {x["Pr(>F)"][1,]})))
  pvalues = as.matrix(pvalues[complete.cases(pvalues),])
  pvalues=data.frame(pvalues[order(pvalues),])
  pvalues[,"fdr"]=p.adjust(p = pvalues[,1], method=correction)
  sig_fdr = pvalues[pvalues$fdr<0.05,]
  colnames(sig_fdr)=c("raw pvalue", correction)
  if (length(row.names(sig_fdr) > 0)) {
    sig_fdr
  }
  else {
    "No significantly differentially abundant taxa" 
  }
}