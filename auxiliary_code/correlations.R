# Okay, a little bit of explanation before I forget.
# The first function, correlation_matrix(), is straightforward: it makes a correlation matrix
# of all samples in a QDNAseq-object. All the way on the bottom of the code is an example of how
# you might want to visualize this data. The correlations are based on the correlation of segment
# values of all bins. But this analysis may be hampered by poor segmentation! Instead, you could
# just take the raw copynumber values, but there is too much noise that is relatively equal between
# samples. You would, for instance, see pretty good correlation between two samples with absolutely 
# no copy numer aberrations. That is where the correlation_matrix_adjusted() comes in. It either
# gives everything equal segments, or it performs pair-wise matching of segment break points for 
# all pairs. It therefore calls the TemplateFromEqualSegments() and two_sample_compare() functions, 
# respectively. 

correlation_matrix <- function(object, trncname=FALSE) {
  samples <- sampleNames(object)
  if(trncname==TRUE) {samples <- gsub("_.*","",samples)}
  if(trncname!=FALSE&&trncname!=TRUE) {samples <- gsub(trncname,"",samples)}
  cormat <- cor(na.exclude(object@assayData$segmented))
  rownames(cormat) <- samples
  colnames(cormat) <- samples
  return(cormat)
}

correlation_matrix_adjusted <- function(object, trncname=FALSE, equalsegments=FALSE, funtype = 'mean') {
  samples <- sampleNames(object)
  n <- length(samples)
  if(trncname==TRUE) {samples <- gsub("_.*","",samples)}
  if(trncname!=FALSE&&trncname!=TRUE) {samples <- gsub(trncname,"",samples)}
  if (equalsegments==FALSE) {
    cormat <- matrix(nrow=n,ncol=n)
    colnames(cormat) <- samples
    rownames(cormat) <- samples
    cormat[n,n] <- 1
    for (i in 1:(n-1)) {
      cormat[i,i] <- 1
      for (j in (i+1):n) {
        cormat[i,j] <- two_sample_compare(object,i,index2=j)$correlation
        cormat[j,i] <- cormat[i,j]
      }
    }
  } else {
    size <- length(object@featureData@data$chromosome)
    tempmat <- matrix(nrow = size, ncol = n)
    for (i in 1:length(samples)) {
      tempmat[,i] <- TemplateFromEqualSegments(object,i,equalsegments=equalsegments,funtype=funtype)$segments
    }
    segmentdf <- na.exclude(as.data.frame(tempmat))
    colnames(segmentdf) <- samples
    cormat <- cor(segmentdf)
  }
  
  return(cormat)
}

# library(corrplot)
# corrplot(cormat, method = "circle")
