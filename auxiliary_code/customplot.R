# This is basically a copy of the singleplot function, with some minor adjustments to line sizes et cetera
# Customize your plots to your liking with this function. You can call it singleplot and source it to overwrite the ACNE singleplot function
# Note that this will not change the output of acne, because it does not use singleplot, but it will affect postanalysisloop!
# Note that this code still depends on ACNE function ObjectsampleToTemplate if you are using a QDNAseq-object as template


customplot <- function(template,cellularity = 1, error, ploidy = 2, standard, title = "Plot",QDNAseqobjectsample = FALSE, cap = 12, chrsubset) {
  library(ggplot2)
  if(QDNAseqobjectsample) {template <- ObjectsampleToTemplate(template, QDNAseqobjectsample)}
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  if(missing(standard)) { standard <- median(rep(segmentdata$values,segmentdata$lengths)) }
  adjustedcopynumbers <- ploidy + (template$copynumbers-standard)/((1/ploidy)*cellularity*standard)
  adjustedsegments <- ploidy + (template$segments-standard)/((1/ploidy)*cellularity*standard)
  df <- data.frame(bin=template$bin,adjustedcopynumbers,adjustedsegments)
  
  rlechr <- rle(as.vector(template$chr))
  binchrend <- c()
  currentbin <- 0
  binchrmdl <- c()
  for (i in 1:length(rlechr$values)) {
    currentmiddle <- currentbin+rlechr$lengths[i]/2
    currentbin <- currentbin+rlechr$lengths[i]
    binchrend <- append(binchrend, currentbin)
    binchrmdl <- append(binchrmdl, currentmiddle)
  } 
  colnames(df)[2] <- "copynumbers"
  colnames(df)[3] <- "segments"
  
  dfna <- na.exclude(df)
  cappedcopynumbers <- dfna[which(dfna$copynumbers > cap),]
  if(length(cappedcopynumbers$copynumbers)>0) {cappedcopynumbers$copynumbers <- cap-0.1}
  cappedsegments <- dfna[which(dfna$segments > cap),]
  if(length(cappedsegments$segments)>0) {cappedsegments$segments <- cap-0.1}
  toppedcopynumbers <- dfna[which(dfna$copynumbers <= 0),]
  if(length(toppedcopynumbers$copynumbers)>0) {toppedcopynumbers$copynumbers <- 0.1}
  toppedsegments <- dfna[which(dfna$segments <= 0),]
  if(length(toppedsegments$segments)>0) {toppedsegments$segments <- 0.1}
  
  line1 <- paste0("Cellularity: ", cellularity)
  line2 <- paste0("")
  if(!missing(error)) {
    line2 <- paste0("Relative error: ", round(error, digits = 3))
  }
  if(missing(chrsubset)){
    ggplot() +
      scale_y_continuous(name = "copies", limits = c(0,cap), breaks = c(0:cap), expand=c(0,0)) +
      scale_x_continuous(name = "chromosome", limits = c(0,binchrend[22]), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
      geom_hline(yintercept = c(0:4), color = '#333333', size = 1) +
      geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 1) +
      geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed", size = 1) +
      geom_point(aes(x = bin,y = copynumbers),data=df, size = 0.5, color = 'black') +
      geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
      geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
      geom_point(aes(x = bin,y = segments),data=df, size = 2, color = 'darkorange') +
      geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 2, color = 'darkorange', shape = 24) +
      geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 2, color = 'darkorange', shape = 25) +
      theme_classic() + theme(
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-0.9, ymax = cap), fill = 'white') +
      annotate("text", x = binchrmdl[2], y = cap-0.3, label = line1) +
      annotate("text", x = binchrmdl[2], y = cap-0.7, label = line2)
  } else {
    firstchr <- range(chrsubset)[1]
    lastchr <- range(chrsubset)[2]
    if(firstchr==1){firstbin<-0
    } else {firstbin<-binchrend[firstchr-1]+1}
    ggplot() +
      scale_y_continuous(name = "copies", limits = c(0,cap), breaks = c(0:cap), expand=c(0,0)) +
      scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), breaks = binchrmdl[firstchr:lastchr], labels = firstchr:lastchr, expand = c(0,0)) +
      geom_hline(yintercept = c(0:4), color = '#333333', size = 1) +
      geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 1) +
      geom_vline(xintercept = binchrend[firstchr:lastchr], color = "#666666", linetype = "dashed", size = 1) +
      geom_point(aes(x = bin,y = copynumbers),data=df, size = 0.5, color = 'black') +
      geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
      geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
      geom_point(aes(x = bin,y = segments),data=df, size = 2, color = 'darkorange') +
      geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 2, color = 'darkorange', shape = 24) +
      geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 2, color = 'darkorange', shape = 25) +
      theme_classic() + theme(
        axis.line = element_line(color='black', size=1), axis.ticks = element_line(color='black', size=1), axis.text = element_text(color='black')) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
}
