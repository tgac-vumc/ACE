# this code can be used to more systematically use the squaremodel function
# squaremodel_summary returns a list of plots starting with the matrix plots and followed by up to 7 CNPs of the best fits
# besides the squaremodel, it also requires the sample data, either as template dataframe or as QDNAseq-object with the appropriate sample index
# you can save the plots to a variable or directly print them (which currently is the default!)
# I figured a 2x4 summary page would be nice. You can't currently specify this with parameters, 
# so if you want to adjust this you have to manually chance the code below.

squaremodel_summary <- function(template,squaremodel,QDNAseqobjectsample=FALSE,samplename,printplots=TRUE,outputdir,imagetype='pdf',trncname=FALSE) {
  library(ggplot2)
  library(Biobase)
  binsize <- 0
  if (QDNAseqobjectsample) {
    fd <- fData(template)
    pd <- pData(template)
    binsize <- fd$end[1]/1000
  } else {binsize <- template$end[1]/1000}
  if (missing(outputdir)) {outputdir<-"."}
  if (missing(samplename)) {
    if (QDNAseqobjectsample) {
      samplename<-pd$name[QDNAseqobjectsample]
      if(trncname==TRUE) {samplename <- gsub("_.*","",samplename)}
      if(trncname!=FALSE&&trncname!=TRUE) {samplename <- gsub(trncname,"",samplename)}
    } else {samplename<-"Sample"}
  }
  imagefunction <- get(imagetype)
  
  cnplots <- min(7,length(unique(squaremodel$minimadf$error)))
  plots <- vector(mode = 'list', length = 8)
  plots[[1]] <- squaremodel$matrixplot + ggtitle(paste0(samplename,", bin size ",binsize,
                                                        " kbp, method = ",squaremodel$method,
                                                        ", penalty = ",squaremodel$penalty,
                                                        ", penploidy = ",squaremodel$penploidy))
  if (cnplots!=0) {
    for (i in 1:cnplots) {
      index <- tail(which(squaremodel$minimadf$error==unique(sort(squaremodel$minimadf$error))[i]),1)
      plots[[1+i]] <- singleplot(template=template,QDNAseqobjectsample=QDNAseqobjectsample,standard=1,
                                 cellularity=squaremodel$minimadf$cellularity[index],
                                 error=squaremodel$minimadf$error[index],
                                 ploidy=squaremodel$minimadf$ploidy[index],
                                 title=paste0("fit ",i, ", ploidy = ",squaremodel$minimadf$ploidy[index]))
    }
  }
  if (printplots == TRUE) {
    if (imagetype == 'pdf') {
      imagefunction(file.path(outputdir,paste0(samplename,".pdf")),width=10.5)
      print(plots)
      dev.off()
    } else {
      if (length(plots)==1) {
        imagefunction(file.path(outputdir,paste0(samplename,".",imagetype)),width = 720, height = 480)
        print(plots)
        dev.off()
      } else {
        imagefunction(file.path(outputdir,paste0(samplename,".",imagetype)),width = 1440, height = 1920)
        print(multiplot(plotlist = plots, cols=2))
        dev.off()
      }
    }
  }
  return(plots)
  
}

# loop through all samples of your object and make squaremodel summaries for those
# since this only works on QDNAseq-objects, I use some of the info in those objects for the filenames and graphs
loop_squaremodel <- function(object,ptop=5,pbottom=1,prows=100,method='RMSE',penalty=0,penploidy=0,
                             outputdir,imagetype='pdf',trncname=FALSE, printsummaries=TRUE) {
  
  imagefunction <- get(imagetype)
  if (missing(outputdir)) {outputdir<-"."}
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  
  library(ggplot2)
  library(Biobase)
  
  
  fd <- fData(object)
  pd <- pData(object)
  if(trncname==TRUE) {pd$name <- gsub("_.*","",pd$name)}
  if(trncname!=FALSE&&trncname!=TRUE) {pd$name <- gsub(trncname,"",pd$name)}
  binsize <- fd$end[1]/1000
  
  matrixplots <- vector(mode = 'list', length = length(pd$name))
  
  for (i in 1:length(pd$name)) {
    model <- squaremodel(object,QDNAseqobjectsample=i,ptop=ptop,pbottom=pbottom,prows=prows,
                         penalty=penalty,penploidy=penploidy,method=method)
    squaremodel_summary(object,model,QDNAseqobjectsample=i,outputdir=outputdir,imagetype=imagetype,trncname=trncname)
    matrixplots[[i]] <- model$matrixplot + ggtitle(paste0(pd$name[i]))
  }
  
  if(printsummaries == TRUE) {
    if(imagetype == 'pdf') {
      imagefunction(file.path(outputdir,paste0("matrixplots.",imagetype)),width=10.5)
      print(matrixplots)
      dev.off()
    } else {
      imagefunction(file.path(outputdir,paste0("matrixplots.",imagetype)),width=1440,height = 480*ceiling(length(matrixplots)/2))
      print(multiplot(plotlist = matrixplots, cols=2))
      dev.off()
    }
  }
}
