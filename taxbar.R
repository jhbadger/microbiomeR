# taxbar: make bar graph of taxonomic distribution of groups of samples
# inputs: tax.summary file from Mothur, numeric taxonomic level (1=kingdom,2=phylum, etc), group linking samples to groups
# optional minimum fraction to include (default 0.05)
# optional arrays of samples to exclude (e.g. Mock) and taxa to exclude (e.g. unclassified)

taxbar = function(taxSummaryFile, taxLevel, group=c(), minFraction=0.05, excludeSamples=c(), excludeTaxa=c()) {
  library(ggplot2)
  library(reshape2)
  counts = read.table(taxSummaryFile, header = TRUE)
  counts = counts[counts[1]==taxLevel,]
  rownames(counts) = make.names(counts$taxon, unique = TRUE)
  # get rid of excess columns in mothur tax.summary file
  counts = counts[,-grep("taxon|taxlevel|rankID|daughterlevels|total", colnames(counts))]
  counts = counts[,!colnames(counts) %in% excludeSamples]
  if (length(excludeTaxa) > 0) {
    counts = counts[!rownames(counts) %in% 
                      rownames(counts)[grep(paste(excludeTaxa,collapse="|"),rownames(counts))],]
  }
  if (length(group) == 0 ) {
    group = colnames(counts)
  }
  # summarize by group
  counts = t(sapply(by(t(counts),group,colSums),identity))
  # eliminate taxa below minFraction
  counts = counts[,colSums(counts)/sum(counts)>=minFraction]
  counts = counts/rowSums(counts)
  mcounts=melt(counts)
  colnames(mcounts)=c("Group","Taxon", "Fraction")
  mcounts$Sample=droplevels(mcounts$Group)
  mcounts$Taxon=droplevels(mcounts$Taxon)
  names=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  ggplot(mcounts, aes(x = Group, y = Fraction, fill = Taxon)) + 
    geom_bar(stat = "identity") + theme(axis.text.x=element_text(size=8)) +
    theme(axis.text.x=element_text(angle=-90,hjust=0)) +
    scale_x_discrete(name="") + scale_y_continuous("Relative Abdundance") +
    ggtitle(paste(names[taxLevel],"abundance")) +
    guides(fill=guide_legend(title="Legend"))
}
