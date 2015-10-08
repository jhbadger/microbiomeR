
# siganova: return data frame of taxa significantly abundant in respect to a group.
# inputs: tax.summary file from Mothur, numeric taxonomic level (1=kingdom,2=phylum, etc), group to check
# optional arrays of samples to exclude (e.g. Mock) and taxa to exclude (e.g. unclassified)
# optional normalization to 100000 reads/sample (default)
# optional correction method (default "fdr")

siganova = function(taxSummaryFile, taxLevel = 2, group, excludeSamples=c(), excludeTaxa=c(), 
                    normalize=TRUE, correction = "fdr") {
  counts = read.table(taxSummaryFile, header = TRUE)
  counts = counts[counts[1]==taxLevel,]
  rownames(counts) = make.names(counts$taxon, unique=TRUE)
  counts = counts[,-grep("taxon|taxlevel|rankID|daughterlevels|total", colnames(counts))]
  counts = counts[,!colnames(counts) %in% excludeSamples]
  if (length(excludeTaxa) > 0) {
    counts = counts[!rownames(counts) %in% 
                      rownames(counts)[grep(paste(excludeTaxa,collapse="|"),rownames(counts))],]
  }
  if (length(group) == 0 ) {
    group = colnames(counts)
  }
  if (normalize) {
    counts = round((counts/colSums(counts))*100000)
  }
  if (length(group) == 0) { # no group provided
    group = colnames(counts)
  }
  fractions = round(1000*rowSums(counts)/sum(colSums(counts)))/1000
  aof=function(x){anova(aov(x~group))}
  anovaresults = apply(t(counts), 2, aof)
  pvalues=t(data.frame(lapply(anovaresults, function(x) {x["Pr(>F)"][1,]})))
  pvalues = as.matrix(pvalues[complete.cases(pvalues),])
  pvalues=data.frame(pvalues[order(pvalues),])
  pvalues[,"fdr"]=p.adjust(p = pvalues[,1], method=correction)
  pvalues[,"fraction"] = fractions[rownames(pvalues)]
  sig_fdr = pvalues[pvalues$fdr<0.05,]
  colnames(sig_fdr)=c("raw pvalue", correction, "fraction")
  if (length(row.names(sig_fdr) > 0)) {
    sig_fdr
  }
  else {
    "No significantly differentially abundant taxa" 
  }
}