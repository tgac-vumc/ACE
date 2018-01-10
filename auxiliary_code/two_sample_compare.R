two_sample_compare <- function(template1,index1=FALSE,ploidy1=2,cellularity1=1,standard1,name1,
                               template2,index2=FALSE,ploidy2=2,cellularity2=1,standard2,name2,
                               equalsegments=FALSE,plot=TRUE,cap=12,qcap=12,trncname=FALSE,legend=TRUE,chrsubset) {
  library(ggplot2)
  library(Biobase)
  if (missing(template2)) {template2<-template1}
  if(index1) {
    pd1<-pData(template1)
    if(trncname==TRUE) {pd1$name <- gsub("_.*","",pd1$name)}
    if(trncname!=FALSE&&trncname!=TRUE) {pd1$name <- gsub(trncname,"",pd1$name)}
    template1 <- ObjectsampleToTemplate(template1, index1)
    if(missing(name1)) {name1<-pd1$name[index1]}
  } else if(missing(name1)) {name1<-"sample1"}
  if(index2) {
    pd2<-pData(template2)
    if(trncname==TRUE) {pd2$name <- gsub("_.*","",pd2$name)}
    if(trncname!=FALSE&&trncname!=TRUE) {pd2$name <- gsub(trncname,"",pd2$name)}
    template2 <- ObjectsampleToTemplate(template2, index2)
    if(missing(name2)) {name2<-pd2$name[index2]}
  } else if(missing(name2)) {name2<-"sample2"}
  segmentdata1 <- rle(as.vector(na.exclude(template1$segments)))
  segmentdata2 <- rle(as.vector(na.exclude(template2$segments)))
  if(missing(standard1)) { standard1 <- median(rep(segmentdata1$values,segmentdata1$lengths)) }
  adjustedcopynumbers1 <- ploidy1 + ((template1$copynumbers-standard1)*(cellularity1*(ploidy1-2)+2))/(cellularity1*standard1)
  adjcnna1 <- as.vector(na.exclude(adjustedcopynumbers1))
  if(missing(standard2)) { standard2 <- median(rep(segmentdata2$values,segmentdata2$lengths)) }
  adjustedcopynumbers2 <- ploidy2 + ((template2$copynumbers-standard2)*(cellularity2*(ploidy2-2)+2))/(cellularity2*standard2)
  adjcnna2 <- as.vector(na.exclude(adjustedcopynumbers2))
  template1na <- na.exclude(template1)
  template1na <- template1na[,1:4]
  template1na$copynumbers <- adjcnna1
  template2na <- na.exclude(template2)
  template2na <- template2na[,1:4]
  template2na$copynumbers <- adjcnna2
  
  if(length(template1na$chr)!=length(template2na$chr)){print("bins not matching")}
  
  Chromosome <- c()
  Start <- c()
  End <- c()
  Num_Bins <- c()
  Mean1 <- c()
  Mean2 <- c()
  SE1 <- c()
  SE2 <- c()
  p_value <- c()
  q_value <- c()
  
  if(equalsegments==FALSE){
    segments1 <- c()
    for (s in 1:length(segmentdata1$lengths)) {
      segments1[s] <- sum(segmentdata1$lengths[1:s]) 
    }
    segments2 <- c()
    for (t in 1:length(segmentdata2$lengths)) {
      segments2[t] <- sum(segmentdata2$lengths[1:t]) 
    }
    allsegments <- sort(unique(c(segments1,segments2)))
    primer <- 1
    for (i in 1:(length(allsegments))) {
      Chromosome[i] <- as.vector(template1na$chr)[primer]
      Start[i] <- as.vector(template1na$start)[primer]
      End[i] <- as.vector(template1na$end)[allsegments[i]]
      Num_Bins[i] <- allsegments[i]-primer+1
      Mean1[i] <- mean(template1na$copynumbers[primer:allsegments[i]])
      SE1[i] <- sd(template1na$copynumbers[primer:allsegments[i]])/sqrt(Num_Bins[i])
      Mean2[i] <- mean(template2na$copynumbers[primer:allsegments[i]])
      SE2[i] <- sd(template2na$copynumbers[primer:allsegments[i]])/sqrt(Num_Bins[i])
      if (Num_Bins[i]>1) {p_value[i] <- t.test(template1na$copynumbers[primer:allsegments[i]],
                           template2na$copynumbers[primer:allsegments[i]])$p.value
      } else {p_value[i]<-1}
      primer <- allsegments[i]+1
    }
    
  } else {
    bincounter <- 1
    segmentcounter <- 0
    if(class(equalsegments)!="numeric"){equalsegments<-20}
    for (c in 1:22) {
      nos <- floor(length(template1na$chr[which(template1na$chr==c)])/equalsegments)
      leftovers <- length(template1na$chr[which(template1na$chr==c)])%%equalsegments
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
        Chromosome[segmentcounter + i] <- as.vector(template1na$chr)[bincounter]
        Start[segmentcounter + i] <- as.vector(template1na$start)[bincounter]
        End[segmentcounter + i] <- as.vector(template1na$end)[bincounter+nobs[i]-1]
        Num_Bins[segmentcounter + i] <- nobs[i]
        Mean1[segmentcounter + i] <- mean(template1na$copynumbers[bincounter:(bincounter+nobs[i]-1)])
        SE1[segmentcounter + i] <- sd(template1na$copynumbers[bincounter:(bincounter+nobs[i]-1)])/sqrt(nobs[i])
        Mean2[segmentcounter + i] <- mean(template2na$copynumbers[bincounter:(bincounter+nobs[i]-1)])
        SE2[segmentcounter + i] <- sd(template2na$copynumbers[bincounter:(bincounter+nobs[i]-1)])/sqrt(nobs[i])
        p_value[segmentcounter + i] <- t.test(template1na$copynumbers[bincounter:(bincounter+nobs[i]-1)],
                                              template2na$copynumbers[bincounter:(bincounter+nobs[i]-1)])$p.value
        bincounter <- bincounter + nobs[i]
        }
      segmentcounter <- segmentcounter+nos
    }
    
  }
  q_value <- p.adjust(p_value,method="BH")
  combinedsegmentsdf <- data.frame(Chromosome,Start,End,Num_Bins,Mean1,SE1,Mean2,SE2,p_value,q_value)
  if (plot==TRUE) {
    q_capped <- sapply(q_value, function(x) max(10^-qcap,x))
    q_corrected <- -log(q_capped,10)*cap/qcap
    df <- data.frame(bin=template1na$bin,Mean1=rep(Mean1,Num_Bins),Mean2=rep(Mean2,Num_Bins),q_value=rep(q_corrected,Num_Bins))
    rlechr <- rle(as.vector(template1$chr))
    binchrend <- c()
    currentbin <- 0
    binchrmdl <- c()
    for (i in 1:length(rlechr$values)) {
      currentmiddle <- currentbin+rlechr$lengths[i]/2
      currentbin <- currentbin+rlechr$lengths[i]
      binchrend <- append(binchrend, currentbin)
      binchrmdl <- append(binchrmdl, currentmiddle)
    }
    if(missing(chrsubset)){
      compareplot <- ggplot() +
        scale_y_continuous(name = "copies", limits = c(0,cap), breaks = c(0:cap), expand=c(0,0), sec.axis = sec_axis(~.*qcap/cap, name = "-log10(q-value)")) +
        scale_x_continuous(name = "chromosome", limits = c(0,binchrend[22]), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
        geom_bar(aes(x=bin, y = q_value), data=df, fill='green', stat='identity') +
        geom_hline(yintercept = c(0:4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = Mean1),data=df, size = 1, color = 'red') +
        geom_point(aes(x = bin,y = Mean2),data=df, size = 1, color = 'blue') +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black'))
        
      if (legend==TRUE){
        compareplot <- compareplot +
          geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-0.9, ymax = cap), fill = 'white') +
          annotate("text", x = binchrend[2], y = cap-0.2, label = name1, color = 'red') +
          annotate("text", x = binchrend[2], y = cap-0.6, label = name2, color = 'blue')
      }
    } else {
      firstchr <- range(chrsubset)[1]
      lastchr <- range(chrsubset)[2]
      if(firstchr==1){firstbin<-0
      } else {firstbin<-binchrend[firstchr-1]+1}
      compareplot <- ggplot() +
        scale_y_continuous(name = "copies", limits = c(0,cap), breaks = c(0:cap), expand=c(0,0), sec.axis = sec_axis(~.*qcap/cap, name = "-log10(q-value)")) +
        scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), breaks = binchrmdl[firstchr:lastchr], labels = firstchr:lastchr, expand = c(0,0)) +
        geom_bar(aes(x=bin, y = q_value), data=df, fill='green', stat='identity') +
        geom_hline(yintercept = c(0:4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend[firstchr:lastchr], color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = Mean1),data=df, size = 1, color = 'red') +
        geom_point(aes(x = bin,y = Mean2),data=df, size = 1, color = 'blue') +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black'))
      
      if (legend==TRUE){
        compareplot <- compareplot +
          geom_rect(aes(xmin=firstbin+1, xmax = firstbin+(binchrend[lastchr]-firstbin)/3.5, ymin = cap-0.9, ymax = cap), fill = 'white') +
          annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-0.2, label = name1, color = 'red') +
          annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-0.6, label = name2, color = 'blue')
      }
    }  
    
    return(list(two_sample_df=combinedsegmentsdf,compareplot=compareplot))
  } else {return(combinedsegmentsdf)}
}


ObjectsampleToTemplate <- function(copyNumbersSegmented, index = 1) {
  library(Biobase)
  library(QDNAseq)
  fd <- fData(copyNumbersSegmented)
  segments <- as.vector(copyNumbersSegmented@assayData$segmented[,index])
  copynumbers <- as.vector(copyNumbersSegmented@assayData$copynumber[,index])
  chr <- as.vector(fd$chromosome)
  start <- as.vector(fd$start)
  end <- as.vector(fd$end)
  bin <- 1:length(chr)
  template <- data.frame(bin,chr,start,end,copynumbers,segments)
  return(template)
}