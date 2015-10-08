# unifracpca: return scatter plot graph of unifrac distances, colored and/or shaped by a factor
# inputs: square Unifrac distance matrixfrom Mothur, optional group factor and colors, shapes,
# and sizes to use for the various factor levels, plus an optional title string if default
# not appropriate dist flag allows returning PCA distance matrix for inspection rather than plotting

unifracpca = function(unifracDistSquare, group=c(), excludeSamples=c(), colors=c(), shapes=c(),
                      sizes=c(), exstring="PCA on unweighted unifrac distances", dist=FALSE) {
  library(ggplot2)
  unifrac = read.table(unifracDistSquare, skip=1, row.names = 1)
  colnames(unifrac) = rownames(unifrac)
  unifrac = unifrac[!colnames(unifrac) %in% excludeSamples,!colnames(unifrac) %in% excludeSamples]
  obj = prcomp(unifrac, scale. = T)
  percents = as.integer((obj$sdev)^2 / sum(obj$sdev^2)*1000)/10.0
  scores = obj$x[,1:2]
  lambda = obj$sdev[1:2]
  d	= data.frame(id=1:nrow(scores), t(t(scores)/lambda))
  if (length(group) > 0) {
    d$group = group
  }
  if (dist) {
    return(d)
  }
  plot = ggplot(d,	aes(x=PC1,	y=PC2))
  if (length(group) > 0) {
    plot = plot + geom_point(aes(color=group, shape=group))
  }
  else {
    plot = plot+geom_point()
  }
  if (length(colors) > 0) {
    plot = plot + scale_color_manual(values=colors, labels=levels(group))
  }
  if (length(shapes) > 0) {
    plot = plot + scale_shape_manual(values=shapes, labels=levels(group)) 
  }
  if (length(sizes) > 0) {
    plot = plot + scale_size_manual(values=sizes, labels=levels(group)) 
  }
  plot = plot + scale_x_continuous(name=paste("PC1 ",percents[1], "%", sep="")) 
  plot = plot + scale_y_continuous(name=paste("PC2 ",percents[2], "%", sep=""))
  plot = plot + ggtitle(exstring)
  plot
}