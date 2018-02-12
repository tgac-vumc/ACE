# This function allows the user to compare copy number profiles of two samples. You can enter the parameters of each 
# sample, but depending on the intended use, this may not be necessary. Alternatively, you can use the cellularity-agnostic
# method "SD", in which case bin and segment values are converted to a Z-score, i.e. how many standard deviations they are 
# removed from the mean ("MAD" is also available but I think this approach is flawed). Using altmethod="SD" basically obviates 
# sample parameters (i.e. ploidy, cellularity, and standard). The crux of two_sample_compare() is that it gives the 
# two samples in question the same segment breakpoints. It can do this by taking all breakpoints of both samples, or by 
# splitting up the chromosomes in artificial segments of a specified number of bins (by means of the parameter equalsegments).
# This way, two samples can be compared regardless of segmentation. All resulting segments come with a mean value and a 
# standard error of the mean. Corresponding segments of the two samples are then tested for having significantly different 
# means as determined with a two-sided t-test. In addition to the reported p-value, the multiple-testing adjusted q-value is
# returned. This value is calculated using p.adjust(p-value, method="BH"). A dataframe with all segment features and values
# is returned, similar to getadjustedsegments(). Additionally a correlation of segments is given (but only if plot=TRUE).
# Finally, the combined copy number plot is given. Segment values are converted to absolute copies as in other ACE functions.
# To keep the plots from becoming too busy, individual bins are left out. Negative log10 q-values are given on a secondary 
# axis. When using altmethod="SD", keep in mind the y-axis now displays the Z-score. You will have to adjust your y-axis limits,
# to include negative numbers, for instance ymin=-3.

two_sample_compare <- function(template1,index1=FALSE,ploidy1=2,cellularity1=1,standard1,name1,
                               template2,index2=FALSE,ploidy2=2,cellularity2=1,standard2,name2,
                               equalsegments=FALSE,altmethod=FALSE,ymin=0,
                               plot=TRUE,cap=12,qcap=12,trncname=FALSE,legend=TRUE,chrsubset) {
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
  if(altmethod==FALSE){
    template1na <- template1na[,1:4]
    template1na$copynumbers <- adjcnna1
  }
  template2na <- na.exclude(template2)
  if(altmethod==FALSE){
    template2na <- template2na[,1:4]
    template2na$copynumbers <- adjcnna2
  }
  
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
      if(altmethod=="MAD"){Mean1[i] <- median(template1na$copynumbers[primer:allsegments[i]])}
      SE1[i] <- sd(template1na$copynumbers[primer:allsegments[i]])/sqrt(Num_Bins[i])
      if(altmethod=="MAD"){SE1[i] <- mad(template1na$copynumbers[primer:allsegments[i]])/sqrt(Num_Bins[i])}
      Mean2[i] <- mean(template2na$copynumbers[primer:allsegments[i]])
      if(altmethod=="MAD"){Mean2[i] <- median(template2na$copynumbers[primer:allsegments[i]])}
      SE2[i] <- sd(template2na$copynumbers[primer:allsegments[i]])/sqrt(Num_Bins[i])
      if(altmethod=="MAD"){SE2[i] <- mad(template2na$copynumbers[primer:allsegments[i]])/sqrt(Num_Bins[i])}
      if (Num_Bins[i]>1) {p_value[i] <- t.test(template1na$copynumbers[primer:allsegments[i]],
                           template2na$copynumbers[primer:allsegments[i]])$p.value
      } else {p_value[i]<-1}
      primer <- allsegments[i]+1
    }
    if(altmethod=="SD"){
      primer <- 1
      ns1 <- mean(rep(Mean1,Num_Bins))
      sd1 <- sd(rep(Mean1,Num_Bins))
      Mean1 <- (Mean1-ns1)/sd1
      SE1 <- SE1/sd1
      ns2 <- mean(rep(Mean2,Num_Bins))
      sd2 <- sd(rep(Mean2,Num_Bins))
      Mean2 <- (Mean2-ns2)/sd2
      SE2 <- SE2/sd2
      for (i in 1:(length(allsegments))) {
        cn1 <- (template1na$copynumbers[primer:allsegments[i]]-ns1)/sd1
        cn2 <- (template2na$copynumbers[primer:allsegments[i]]-ns2)/sd2
        if (Num_Bins[i]>1) {p_value[i] <- t.test(cn1,cn2)$p.value
        } else {p_value[i]<-1}
        primer <- allsegments[i]+1
      }
    }
    if(altmethod=="MAD"){
      primer <- 1
      ns1 <- median(rep(Mean1,Num_Bins))
      mad1 <- mad(rep(Mean1,Num_Bins))
      Mean1 <- (Mean1-ns1)/mad1
      SE1 <- SE1/mad1
      ns2 <- median(rep(Mean2,Num_Bins))
      mad2 <- mad(rep(Mean2,Num_Bins))
      Mean2 <- (Mean2-ns2)/mad2
      SE2 <- SE2/mad2
      for (i in 1:(length(allsegments))) {
        cn1 <- (template1na$copynumbers[primer:allsegments[i]]-ns1)/mad1
        cn2 <- (template2na$copynumbers[primer:allsegments[i]]-ns2)/mad2
        if (Num_Bins[i]>1) {p_value[i] <- t.test(cn1,cn2)$p.value
        } else {p_value[i]<-1}
        primer <- allsegments[i]+1
      }
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
    if(altmethod=="SD"){
      ns1 <- mean(rep(Mean1,Num_Bins))
      sd1 <- sd(rep(Mean1,Num_Bins))
      Mean1 <- (Mean1-ns1)/sd1
      SE1 <- SE1/sd1
      ns2 <- mean(rep(Mean2,Num_Bins))
      sd2 <- sd(rep(Mean2,Num_Bins))
      Mean2 <- (Mean2-ns2)/sd2
      SE2 <- SE2/sd2
      for (i in 1:(length(Mean1))) {
        firstbin <- which(template1na$chr==Chromosome[i] & template1na$start==Start[i])
        lastbin <- which(template1na$chr==Chromosome[i] & template1na$end==End[i])
        cn1 <- (template1na$copynumbers[firstbin:lastbin]-ns1)/sd1
        cn2 <- (template2na$copynumbers[firstbin:lastbin]-ns2)/sd2
        if (Num_Bins[i]>1) {p_value[i] <- t.test(cn1,cn2)$p.value
        } else {p_value[i]<-1}
      }
    }
  }
  q_value <- p.adjust(p_value,method="BH")
  combinedsegmentsdf <- data.frame(Chromosome,Start,End,Num_Bins,Mean1,SE1,Mean2,SE2,p_value,q_value)
  if (plot==TRUE) {
    q_capped <- sapply(q_value, function(x) max(10^-qcap,x))
    q_corrected <- -log(q_capped,10)*cap/qcap
    df <- data.frame(bin=template1na$bin,Mean1=rep(Mean1,Num_Bins),Mean2=rep(Mean2,Num_Bins),q_value=rep(q_corrected,Num_Bins))
    correlation <- cor(Mean1,Mean2)
    rlechr <- rle(as.vector(template1$chr))
    binchrend <- c()
    currentbin <- 0
    binchrmdl <- c()
    if(altmethod==FALSE){yname<-"copies"} else {yname <- paste0(altmethod," units")}
    for (i in 1:length(rlechr$values)) {
      currentmiddle <- currentbin+rlechr$lengths[i]/2
      currentbin <- currentbin+rlechr$lengths[i]
      binchrend <- append(binchrend, currentbin)
      binchrmdl <- append(binchrmdl, currentmiddle)
    }
    if(missing(chrsubset)){
      compareplot <- ggplot() +
        scale_y_continuous(name = yname, limits = c(ymin,cap), breaks = c(ymin:cap), expand=c(0,0), sec.axis = sec_axis(~.*qcap/cap, name = "-log10(q-value)")) +
        scale_x_continuous(name = "chromosome", limits = c(0,binchrend[22]), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
        geom_bar(aes(x=bin, y = q_value), data=df, fill='green', stat='identity') +
        geom_hline(yintercept = c(ymin:4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = Mean1),data=df, size = 1, color = 'red') +
        geom_point(aes(x = bin,y = Mean2),data=df, size = 1, color = 'blue') +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black'))
        
      if (legend==TRUE){
        compareplot <- compareplot +
          geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-1.2, ymax = cap), fill = 'white') +
          annotate("text", x = binchrend[2], y = cap-0.2, label = name1, color = 'red') +
          annotate("text", x = binchrend[2], y = cap-0.6, label = name2, color = 'blue') +
          annotate("text", x = binchrend[2], y = cap-1, label = paste0("r = ",round(correlation,digits=3)), color = 'black')
      }
    } else {
      firstchr <- range(chrsubset)[1]
      lastchr <- range(chrsubset)[2]
      if(firstchr==1){firstbin<-0
      } else {firstbin<-binchrend[firstchr-1]+1}
      compareplot <- ggplot() +
        scale_y_continuous(name = yname, limits = c(ymin,cap), breaks = c(ymin:cap), expand=c(0,0), sec.axis = sec_axis(~.*qcap/cap, name = "-log10(q-value)")) +
        scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), breaks = binchrmdl[firstchr:lastchr], labels = firstchr:lastchr, expand = c(0,0)) +
        geom_bar(aes(x=bin, y = q_value), data=df, fill='green', stat='identity') +
        geom_hline(yintercept = c(ymin:4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend[firstchr:lastchr], color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = Mean1),data=df, size = 1, color = 'red') +
        geom_point(aes(x = bin,y = Mean2),data=df, size = 1, color = 'blue') +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black'))
      
      if (legend==TRUE){
        compareplot <- compareplot +
          geom_rect(aes(xmin=firstbin+1, xmax = firstbin+(binchrend[lastchr]-firstbin)/3.5, ymin = cap-1.2, ymax = cap), fill = 'white') +
          annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-0.2, label = name1, color = 'red') +
          annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-0.6, label = name2, color = 'blue') +
          annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-1, label = paste0("r = ",round(correlation,digits=3)), color = 'black')
      }
    }  
    
    return(list(two_sample_df=combinedsegmentsdf,correlation=correlation,compareplot=compareplot))
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