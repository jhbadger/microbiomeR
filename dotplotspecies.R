# siganova: plot dot plot of provided taxon in respect to a group (or by sample if no group provided).
# inputs: tax.summary file from Mothur, taxon string, group to group by
# optional arrays of samples to exclude (e.g. Mock)

dotplotspecies = function(taxSummaryFile, taxon, group=c(), excludeSamples=c(),normalize=TRUE) {
  library(ggplot2)
  library(reshape2)
  counts = read.table(taxSummaryFile, header = TRUE)
  counts = counts[counts$taxon==taxon,]
  rownames(counts) = counts$taxon
  # get rid of excess columns in mothur tax.summary file
  counts = counts[,-grep("taxon|taxlevel|rankID|daughterlevels|total", colnames(counts))]
  counts = counts[,!colnames(counts) %in% excludeSamples]
  if (normalize) {
    counts = round((counts/rowSums(counts))*100000)
  }
  if (length(group) == 0 ) {
    group = colnames(counts)
  }
  counts=melt(counts)
  counts$group = group
  colnames(counts) = c("Sample", "Value","Group")
  ggplot(counts,aes(x=Group,y=Value)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
    theme(axis.text.x=element_text(angle=-90,hjust=0)) + theme(legend.position="none") +
    scale_x_discrete(name="") + scale_y_continuous("Taxon Counts") + ggtitle(taxon)
}