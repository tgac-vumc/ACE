TemplateFromEqualSegments <- function(object, index = 1, equalsegments = 20, funtype = 'mean') {
  fun <- get(funtype)
  library(Biobase)
  library(QDNAseq)
  fd <- fData(object)
  pd <- pData(object)
  copynumbers <- as.vector(object@assayData$copynumber[,index])
  segments <- copynumbers
  chr <- as.vector(fd$chromosome)
  start <- as.vector(fd$start)
  end <- as.vector(fd$end)
  bin <- 1:length(chr)
  
  for (c in 1:22) {
    bincounter <- 1 
    nos <- floor(sum(!is.na(copynumbers)&chr==c)/equalsegments)
    leftovers <- sum(!is.na(copynumbers)&chr==c)%%equalsegments
    nobs <- rep(equalsegments,nos)
    if (nos == 0) {
      nos <- 1
      nobs <- leftovers
      leftovers <- 0
    }
    if (leftovers > 0) {
      for (a in 0:(leftovers-1)) {
        divvy <- a%%nos+1
        nobs[divvy] <- nobs[divvy] + 1
      }
    }
    
    for (i in 1:nos) {
      segments[which(!is.na(copynumbers)&chr==c)[bincounter:(bincounter+nobs[i]-1)]] <- rep(fun(segments[which(!is.na(copynumbers)&chr==c)[bincounter:(bincounter+nobs[i]-1)]]),nobs[i])
      bincounter <- bincounter + nobs[i]
    }
    
  }
  
  template <- data.frame(bin,chr,start,end,copynumbers,segments)
return(template)
}