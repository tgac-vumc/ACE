singlemodel_exclude <- function(template,ploidy = 2, standard, QDNAseqobjectsample = FALSE, method = 'RMSE', penalty = 0, exclude=c(), highlightminima = TRUE) {
  library(ggplot2)
  if(QDNAseqobjectsample) {template <- ObjectsampleToTemplate(template, QDNAseqobjectsample)}
  template <- template[!template$chr %in% exclude,]
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  if(missing(standard)) { standard <- median(rep(segmentdata$values,segmentdata$lengths)) }
  
  fraction <- c()
  expected <- c()
  temp <- c()
  errorlist <- c()
  for (i in 5:100) {
    fraction[i-4] <- i/100
    for (p in 1:12) {
      expected[p] <- standard*(1+(p-ploidy)*fraction[i-4]/(fraction[i-4]*(ploidy-2)+2))
    }
    for (j in 1:length(segmentdata$values)) {
      if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))^2}
      else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))}
      else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty)}
      else {print("Not a valid method")}
    }
    if(method=='RMSE') {errorlist[i-4] <- sqrt(sum(rep(temp,segmentdata$lengths))/sum(segmentdata$lengths))}
    else if(method=='SMRE') {errorlist[i-4] <- sum(rep(temp,segmentdata$lengths))/sum(segmentdata$lengths)^2}
    else if(method=='MAE') {errorlist[i-4] <- sum(rep(temp,segmentdata$lengths))/sum(segmentdata$lengths)}
    
  }
  
  minima <- c()
  rerror <- c()
  
  if (round(errorlist[1], digits = 10) < round(errorlist[2], digits = 10)) {
    lastminimum <- fraction[1]
    minima[1] <- fraction[1]
    rerror[1] <- errorlist[1]/max(errorlist)
  }
  for (l in 6:99) {
    if (round(errorlist[l-4], digits = 10) < round(errorlist[l-5], digits = 10) & round(errorlist[l-4], digits = 10) < round(errorlist[l-3], digits =10)) { 
      lastminimum <- fraction[l-4]
      minima <- append(minima,fraction[l-4])
      rerror <- append(rerror,(errorlist[l-4]/max(errorlist)))
    }
  }
  if (errorlist[100-4] <= errorlist[100-5]) {
    lastminimum <- fraction[100-4]
    minima <- append(minima, fraction[100-4])
    rerror <- append(rerror, errorlist[100-4]/max(errorlist))
  }
  
  cellularity <- 5:100
  tempdf <- data.frame(cellularity,errorlist=errorlist/max(errorlist))
  minimadf <- data.frame(minima=minima*100,rerror)
  if(highlightminima==TRUE) {
    tempplot <- ggplot() +
      scale_y_continuous(name = "relative error", limits = c(0,1.05), expand=c(0,0)) +
      scale_x_continuous(name = "cellularity (%)") +
      geom_vline(xintercept = seq(from = 10, to = 100, by = 10), color = "#666666", linetype = "dashed") +
      geom_point(aes(y=errorlist, x=cellularity), data=tempdf) +
      geom_point(aes(y=rerror, x=minima), data=minimadf, color = 'red') +
      theme_classic() + theme(
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
      ggtitle("errorlist") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    tempplot <- ggplot() +
      scale_y_continuous(name = "relative error", limits = c(0,1.05), expand=c(0,0)) +
      scale_x_continuous(name = "cellularity (%)") +
      geom_vline(xintercept = seq(from = 10, to = 100, by = 10), color = "#666666", linetype = "dashed") +
      geom_point(aes(y=errorlist, x=cellularity), data=tempdf) +
      theme_classic() + theme(
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
      ggtitle("errorlist") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(list(ploidy=ploidy,standard=standard,method=method,penalty=penalty,minima=minima,rerror=rerror,errorlist=errorlist,errorplot=tempplot))
  
}

squaremodel_exclude <- function(template, QDNAseqobjectsample = FALSE, prows=100, ptop=5, pbottom=1, method = 'RMSE', penalty = 0, penploidy = 0, exclude=c(), highlightminima = TRUE) {
  library(ggplot2)
  library(Biobase)
  if(QDNAseqobjectsample) {template <- ObjectsampleToTemplate(template, QDNAseqobjectsample)}
  template <- template[!template$chr %in% exclude,]
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  
  fraction <- c()
  errormatrix <- matrix(nrow=(prows+1),ncol=96)
  listofploidy <- c()
  listofcellularity <- c()
  listoferrors <- c()
  for (t in 0:prows) {
    ploidy <- ptop-((ptop-pbottom)/prows)*t
    listofploidy <- append(listofploidy, rep(ploidy,96))
    expected <- c()
    temp <- c()
    errorlist <- c()
    for (i in 5:100) {
      fraction[i-4] <- i/100
      for (p in 1:12) {
        expected[p] <- 1+(p-ploidy)*fraction[i-4]/(fraction[i-4]*(ploidy-2)+2)
      }
      for (j in 1:length(segmentdata$values)) {
        if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i-4]^penalty))^2}
        else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i-4]^penalty))}
        else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i-4]^penalty)}
        else {print("Not a valid method")}
      }
      if(method=='RMSE') {errorlist[i-4] <- sqrt(sum(rep(temp,segmentdata$lengths))/sum(segmentdata$lengths))}
      else if(method=='SMRE') {errorlist[i-4] <- sum(rep(temp,segmentdata$lengths))/sum(segmentdata$lengths)^2}
      else if(method=='MAE') {errorlist[i-4] <- sum(rep(temp,segmentdata$lengths))/sum(segmentdata$lengths)}
    }
    listofcellularity <- append(listofcellularity, fraction)
    listoferrors <- append(listoferrors, errorlist)
    errormatrix[t+1,] <- errorlist
    
  }
  minimat <- matrix(nrow=(prows+1),ncol=96)
  for (i in 1:(prows+1)) {
    for (j in 1:96) {
      if (i==1|i==(prows+1)) {
        minimat[i,j] <- FALSE
      } else if (j==96) {
        minimat[i,j] <- errormatrix[i,j]==min(errormatrix[(i-1):(i+1),(j-1):j])
      } else {
        minimat[i,j] <- errormatrix[i,j]==min(errormatrix[(i-1):(i+1),(j-1):(j+1)])
      }
    }
  }
  round(errormatrix, digits = 10)
  round(listoferrors, digits = 10)
  errordf <- data.frame(ploidy=listofploidy,
                        cellularity=listofcellularity,
                        error=listoferrors/max(listoferrors),
                        minimum=as.vector(t(minimat)))
  minimadf <- errordf[which(errordf$minimum==TRUE),]
  if(highlightminima==TRUE){
    tempplot <- ggplot() +
      geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
      geom_point(data=minimadf, aes(x=cellularity, y=ploidy, alpha=min(error)/error), shape=16) +
      scale_fill_gradient(low="green", high="red") +
      ggtitle("Matrix of errors") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    tempplot <- ggplot() +
      geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
      scale_fill_gradient(low="green", high="red") +
      ggtitle("Matrix of errors") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(list(method=method, 
              penalty=penalty, 
              penploidy=penploidy,
              errormatrix=errormatrix, 
              minimatrix = minimat, 
              errordf=errordf, 
              minimadf=minimadf, 
              matrixplot=tempplot))
  
}