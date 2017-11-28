two_sample_compare <- function(template1,index1=FALSE,ploidy1=2,cellularity1=1,standard1,
                               template2,index2=FALSE,ploidy2=2,cellularity2=1,standard2,
                               equalsegments=FALSE,plot=TRUE,cap=12,qcap=12) {
  library(ggplot2)
  library(Biobase)
  if(index1) {template1 <- ObjectsampleToTemplate(template1, index1)}
  if(index2) {template2 <- ObjectsampleToTemplate(template2, index2)}
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
  template2na <- na.exclude(template2na[,1:4])
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
          index <- a%%nos+1
          nobs[index] <- nobs[index] + 1
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
    q_capped <- sapply(q_value, function(x) max(10^-cap,x))
    df <- data.frame(bin=template1na$bin,Mean1=rep(Mean1,Num_Bins),Mean2=rep(Mean2,Num_Bins),q_value=rep(q_capped,Num_Bins))
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
    compareplot <- ggplot() +
      scale_y_continuous(name = "copies", limits = c(0,cap), breaks = c(0:cap), expand=c(0,0), sec.axis = sec_axis(~.*qcap/cap, name = "-log10(q-value)")) +
      scale_x_continuous(name = "chromosome", limits = c(0,binchrend[22]), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
      geom_bar(aes(x=bin, y = -log(q_value,10)), data=df, fill='green', stat='identity') +
      geom_hline(yintercept = c(0:4), color = '#333333', size = 0.5) +
      geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 0.5) +
      geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
      geom_point(aes(x = bin,y = Mean1),data=df, size = 1, color = 'red') +
      geom_point(aes(x = bin,y = Mean2),data=df, size = 1, color = 'blue') +
      theme_classic() + theme(
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black'))
    
    return(list(two_sample_df=combinedsegmentsdf,compareplot=compareplot))
  } else {return(combinedsegmentsdf)}
}