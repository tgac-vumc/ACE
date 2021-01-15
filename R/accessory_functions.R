# What does ACEcall contribute to currently available callers for CNAs? Well ...
# you can set a cutoff for absolute copy number gains or losses. A gain of 0.1 
# can be "irrelevant" in a very pure tumor sample, but represents a full gain in
# a sample with only 10% tumor material! This function will output a new 
# "template" data frame that contains the adjusted copynumber and segment 
# information, a call, and a p-value and q-value assessing the probability that 
# if the true segment mean equals the ploidy, the resulting or a more extreme
# segment mean would be found. The plot color codes the segments according to
# the calls.

ACEcall <- function(template, QDNAseqobjectsample=FALSE, cellularity=1, ploidy=2, standard, plot=TRUE, title, 
                    pcutoff, qcutoff=-3, subcutoff=0.2, trncname=FALSE, cap=12, bottom=0, chrsubset, 
                    onlyautosomes = TRUE, sgc = c()) {
  if(QDNAseqobjectsample) {
    if(missing(title)) {
      pd<-Biobase::pData(template)
      if(trncname==TRUE) {pd$name <- gsub("_.*","",pd$name)}
      if(trncname!=FALSE&&trncname!=TRUE) {pd$name <- gsub(trncname,"",pd$name)}
      title <- pd$name[QDNAseqobjectsample]
    }
    template <- objectsampletotemplate(template, QDNAseqobjectsample)
  }
  gc <- rep(2, nrow(template))
  gc[template$chr %in% sgc] <- 1
  if(missing(title)) {title <- "Plot"}
  #segmentdata <- rle(as.vector(na.exclude(template$segments)))
  if(missing(standard) || !is(standard, "numeric")) { standard <- median(template$segments, na.rm = T) }
  adjustedcopynumbers <- template$copynumbers*(ploidy+2/cellularity-2)/standard - gc/cellularity + gc
  adjustedsegments <- template$segments*(ploidy+2/cellularity-2)/standard - gc/cellularity + gc
  #adjustedcopynumbers <- ploidy + ((template$copynumbers-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
  #adjcopynumberdata <- as.vector(na.exclude(adjustedcopynumbers))
  #adjustedsegments <- ploidy + ((template$segments-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
  template$chr <- gsub("chr","",template$chr,ignore.case = TRUE)
  df <- data.frame(bin=template$bin,chr=template$chr,adjustedcopynumbers,adjustedsegments,gc)
  colnames(df)[3] <- "copynumbers"
  colnames(df)[4] <- "segments"
  
  dfna <- na.exclude(df)
  bin <- dfna$bin
  copynumbers <- dfna$copynumbers
  segments <- dfna$segments
  if(bottom=="min"){bottom<-floor(min(dfna$copynumbers))}
  if(cap=="max"){cap<-ceiling(max(dfna$copynumbers))}
  if(!missing(chrsubset)){
    firstchr <- range(chrsubset)[1]
    lastchr <- range(chrsubset)[2]
    dfna <- dfna[dfna$chr %in% rlechr$values[seq(firstchr,lastchr)],]
  }
  segmentdata <- rle(as.vector(paste0(dfna$chr, "_", dfna$segments)))
  num_bins <- segmentdata$lengths
  segment_mean <- c()
  segment_SE <- c()
  calls <- c()
  probnorm <- c()
  qprobnorm <- c()
  probloss <- c()
  probgain <- c()
  probdloss <- c()
  probamp <- c()
  color <- c()
  p <- c()
  dl <- 1-subcutoff
  l <- round(ploidy-1,digits=0)+subcutoff
  sl <- round(ploidy,digits=0)-subcutoff
  sg <- round(ploidy,digits=0)+subcutoff
  g <- round(ploidy+1,digits=0)-subcutoff
  dg <- round(ploidy+2,digits=0)-subcutoff
  amp <- round(ploidy+3,digits=0)-subcutoff
  
  counter <- 1
  for (i in seq_along(segmentdata$lengths)) {
    segment_mean[i] <- mean(dfna$copynumbers[seq(counter, counter+segmentdata$lengths[i]-1)])
    segment_SE[i] <- sd(dfna$copynumbers[seq(counter, counter+segmentdata$lengths[i]-1)])/sqrt(segmentdata$lengths[i])
    p[i] <- 2*(1-pnorm(abs(round(ploidy,digits = 0)-segment_mean[i])/segment_SE[i]))
    probnorm[i] <- round(log(p[i],10),digits=3)
    counter <- counter+segmentdata$lengths[i]
  }
  qprobnorm <- round(log(p.adjust(p, method='BH'),10),digits=3)
  
  segment_m <- rep(segment_mean, num_bins)+2-dfna$gc
  calls[(segment_m<sl)] <- -0.5
  color[(segment_m<sl)] <- 'turquoise'
  calls[(segment_m<l)] <- -1
  color[(segment_m<l)] <- 'blue'
  calls[(segment_m<dl)] <- -2
  color[(segment_m<dl)] <- 'darkblue'
  calls[(segment_m>sg)] <- 0.5
  color[(segment_m>sg)] <- 'gold'
  calls[(segment_m>g)] <- 1
  color[(segment_m>g)] <- 'darkorange'
  calls[(segment_m>dg)] <- 2
  color[(segment_m>dg)] <- 'red'
  calls[(segment_m>amp)] <- 3
  color[(segment_m>amp)] <- 'purple'
  calls[(segment_m>=sl&segment_m<=sg)] <- 0
  color[(segment_m>=sl&segment_m<=sg)] <- 'black'
  if(missing(pcutoff)){
    calls[rep(qprobnorm,num_bins)>=qcutoff] <- 0
    color[rep(qprobnorm,num_bins)>=qcutoff] <- 'black'
  } else {
    calls[rep(probnorm,num_bins)>=pcutoff] <- 0
    color[rep(probnorm,num_bins)>=pcutoff] <- 'black'
  }
  
  # templatena <- na.exclude(template)
  # templatena$segment_mean <- rep(segment_mean,num_bins)
  # templatena$segment_SE <- rep(segment_SE,num_bins)
  # templatena$probnorm <- rep(probnorm,num_bins)
  # templatena$qprobnorm <- rep(qprobnorm,num_bins)
  # templatena$calls <- rep(calls,num_bins)
  # templatena$color <- rep(color,num_bins)
  df$segment_mean <- rep(NA,length(df$bin))
  df$segment_mean[dfna$bin] <- rep(segment_mean,num_bins) # templatena$segment_mean
  df$segment_SE <- rep(NA,length(df$bin))
  df$segment_SE[dfna$bin] <- rep(segment_SE,num_bins) # templatena$segment_SE
  df$pnorm_log10 <- rep(NA,length(df$bin))
  df$pnorm_log10[dfna$bin] <- rep(probnorm,num_bins) # templatena$probnorm
  df$qnorm_log10 <- rep(NA,length(df$bin))
  df$qnorm_log10[dfna$bin] <- rep(qprobnorm,num_bins) # templatena$qprobnorm
  df$calls <- rep(NA,length(df$bin))
  df$calls[dfna$bin] <- calls # templatena$calls
  df$color <- rep('black',length(df$bin))
  df$color[dfna$bin] <- color # templatena$color
  dfna <- na.exclude(df)
  
  if(plot==TRUE){
    if(onlyautosomes==TRUE) {
      rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1, 22)]))
    }	else if (onlyautosomes==FALSE) {rlechr <- rle(as.vector(template$chr))
    } else {rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1, onlyautosomes)]))}
    binchrend <- c()
    currentbin <- 0
    binchrmdl <- c()
    for (i in seq_along(rlechr$values)) {
      currentmiddle <- currentbin+rlechr$lengths[i]/2
      currentbin <- currentbin+rlechr$lengths[i]
      binchrend <- append(binchrend, currentbin)
      binchrmdl <- append(binchrmdl, currentmiddle)
    }
    
    if(bottom=="min"){bottom<-floor(min(dfna$adjustedcopynumbers))}
    if(cap=="max"){cap<-ceiling(max(dfna$adjustedcopynumbers))}
    cappedcopynumbers <- dfna[(dfna$adjustedcopynumbers > cap),]
    if(length(cappedcopynumbers$adjustedcopynumbers)>0) {cappedcopynumbers$adjustedcopynumbers <- cap-0.1}
    cappedsegments <- dfna[(dfna$adjustedsegments > cap),]
    if(length(cappedsegments$segments)>0) {cappedsegments$adjustedsegments <- cap-0.1}
    toppedcopynumbers <- dfna[(dfna$adjustedcopynumbers <= bottom),]
    if(length(toppedcopynumbers$adjustedcopynumbers)>0) {toppedcopynumbers$adjustedcopynumbers <- bottom+0.1}
    toppedsegments <- dfna[(dfna$adjustedsegments <= bottom),]
    if(length(toppedsegments$adjustedsegments)>0) {toppedsegments$adjustedsegments <- bottom+0.1}
    line1 <- paste0("Cellularity: ", cellularity)
    if(missing(chrsubset)){
      calledplot <- ggplot2::ggplot() +
        scale_y_continuous(name = "copies", limits = c(bottom,cap), breaks = seq(bottom, cap), expand=c(0,0)) +
        scale_x_continuous(name = "chromosome", limits = c(0,tail(binchrend,1)), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
        geom_hline(yintercept = seq(0, 4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = seq(5, cap-1), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = copynumbers),data=dfna[(dfna$copynumbers>bottom&dfna$copynumbers<cap),], size = 0.1, color = 'gray') +
        geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'gray', shape = 24) +
        geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'gray', shape = 25) +
        geom_point(aes(x = bin,y = segments),data=dfna[(dfna$segments>bottom&dfna$segments<cap),], size = 1, color = dfna$color[(dfna$segments>bottom&dfna$segments<cap)]) +
        geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'purple', shape = 24) +
        geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkblue', shape = 25) +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-0.9, ymax = cap), fill = 'white') +
        annotate("text", x = binchrmdl[2], y = cap-0.3, label = line1)
    } else {
      firstchr <- range(chrsubset)[1]
      lastchr <- range(chrsubset)[2]
      if(firstchr==1){firstbin<-0
      } else {firstbin<-binchrend[firstchr-1]+1}
      dfna <- dfna[dfna$chr %in% rlechr$values[seq(firstchr, lastchr)],]
      calledplot <- ggplot2::ggplot() +
        scale_y_continuous(name = "copies", limits = c(bottom,cap), breaks = seq(bottom, cap), expand=c(0,0)) +
        scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), breaks = binchrmdl[seq(firstchr,lastchr)], labels = rlechr$values[seq(firstchr,lastchr)], expand = c(0,0)) +
        geom_hline(yintercept = seq(0, 4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = seq(5, cap-1), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend[seq(firstchr,lastchr)], color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = copynumbers),data=dfna[(dfna$copynumbers>bottom&dfna$copynumbers<cap),], size = 0.1, color = 'gray') +
        geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'gray', shape = 24) +
        geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'gray', shape = 25) +
        geom_point(aes(x = bin,y = segments),data=dfna[(dfna$segments>bottom&dfna$segments<cap),], size = 1, color = dfna$color[(dfna$segments>bottom&dfna$segments<cap)]) +
        geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'purple', shape = 24) +
        geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkblue', shape = 25) +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    return(list(calledtemplate=df,calledplot=calledplot))
  } else {return(calledtemplate=df)}
}

# Finally, the combined copy number plot is given. Segment values are converted to absolute copies as in other ACE functions.
# To keep the plots from becoming too busy, individual bins are left out. Negative log10 q-values are given on a secondary 
# axis. When using altmethod="SD", keep in mind the y-axis now displays the Z-score. You will have to adjust your y-axis limits,
# to include negative numbers, for instance bottom=-3.

twosamplecompare <- function(template1, index1=FALSE, ploidy1=2, cellularity1=1, standard1, name1,
                             template2, index2=FALSE, ploidy2=2, cellularity2=1, standard2, name2,
                             equalsegments=FALSE, altmethod=FALSE, cap=12, qcap=12, bottom=0, 
                             plot=TRUE, trncname=FALSE, legend=TRUE, chrsubset, onlyautosomes=TRUE,
                             showcorrelation=TRUE) {
  if (missing(template2)) {template2<-template1}
  if(index1) {
    pd1<-Biobase::pData(template1)
    if(trncname==TRUE) {pd1$name <- gsub("_.*","",pd1$name)}
    if(trncname!=FALSE&&trncname!=TRUE) {pd1$name <- gsub(trncname,"",pd1$name)}
    template1 <- objectsampletotemplate(template1, index1)
    if(missing(name1)) {name1<-pd1$name[index1]}
  } else if(missing(name1)) {name1<-"sample1"}
  if(index2) {
    pd2<-Biobase::pData(template2)
    if(trncname==TRUE) {pd2$name <- gsub("_.*","",pd2$name)}
    if(trncname!=FALSE&&trncname!=TRUE) {pd2$name <- gsub(trncname,"",pd2$name)}
    template2 <- objectsampletotemplate(template2, index2)
    if(missing(name2)) {name2<-pd2$name[index2]}
  } else if(missing(name2)) {name2<-"sample2"}
  if(nrow(template1)!=nrow(template2)){print("bins not matching")}
  na1 <- apply(template1, 1, function(x) {any(is.na(x))})
  na2 <- apply(template2, 1, function(x) {any(is.na(x))})
  template1na <- template1[!(na1 | na2), ]
  template2na <- template2[!(na1 | na2), ]
  segmentdata1 <- rle(as.vector(na.exclude(template1na$segments)))
  segmentdata2 <- rle(as.vector(na.exclude(template2na$segments)))
  if(missing(standard1) || !is(standard1, "numeric")) { standard1 <- median(rep(segmentdata1$values,segmentdata1$lengths)) }
  adjustedcopynumbers1 <- ploidy1 + ((template1na$copynumbers-standard1)*(cellularity1*(ploidy1-2)+2))/(cellularity1*standard1)
  adjcnna1 <- as.vector(na.exclude(adjustedcopynumbers1))
  if(missing(standard2) || !is(standard2, "numeric")) { standard2 <- median(rep(segmentdata2$values,segmentdata2$lengths)) }
  adjustedcopynumbers2 <- ploidy2 + ((template2na$copynumbers-standard2)*(cellularity2*(ploidy2-2)+2))/(cellularity2*standard2)
  adjcnna2 <- as.vector(na.exclude(adjustedcopynumbers2))

  bin <- template1na$bin
  if(altmethod==FALSE){
    template1na <- template1na[, seq(1,4)]
    template1na$copynumbers <- adjcnna1
  }
  
  if(altmethod==FALSE){
    template2na <- template2na[, seq(1,4)]
    template2na$copynumbers <- adjcnna2
  }
  
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
    for (s in seq_along(segmentdata1$lengths)) {
      segments1[s] <- sum(segmentdata1$lengths[seq(1,s)]) 
    }
    segments2 <- c()
    for (t in seq_along(segmentdata2$lengths)) {
      segments2[t] <- sum(segmentdata2$lengths[seq(1,t)]) 
    }
    allsegments <- sort(unique(c(segments1,segments2)))
    primer <- 1
    for (i in seq_along(allsegments)) {
      Chromosome[i] <- as.vector(template1na$chr)[primer]
      Start[i] <- as.vector(template1na$start)[primer]
      End[i] <- as.vector(template1na$end)[allsegments[i]]
      Num_Bins[i] <- allsegments[i]-primer+1
      Mean1[i] <- mean(template1na$copynumbers[seq(primer,allsegments[i])])
      if(altmethod=="MAD"){Mean1[i] <- median(template1na$copynumbers[seq(primer,allsegments[i])])}
      SE1[i] <- sd(template1na$copynumbers[seq(primer,allsegments[i])])/sqrt(Num_Bins[i])
      if(altmethod=="MAD"){SE1[i] <- mad(template1na$copynumbers[seq(primer,allsegments[i])])/sqrt(Num_Bins[i])}
      Mean2[i] <- mean(template2na$copynumbers[seq(primer,allsegments[i])])
      if(altmethod=="MAD"){Mean2[i] <- median(template2na$copynumbers[seq(primer,allsegments[i])])}
      SE2[i] <- sd(template2na$copynumbers[seq(primer,allsegments[i])])/sqrt(Num_Bins[i])
      if(altmethod=="MAD"){SE2[i] <- mad(template2na$copynumbers[seq(primer,allsegments[i])])/sqrt(Num_Bins[i])}
      if (Num_Bins[i]>1) {p_value[i] <- t.test(template1na$copynumbers[seq(primer,allsegments[i])],
                                               template2na$copynumbers[seq(primer,allsegments[i])])$p.value
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
      for (i in seq_along(allsegments)) {
        cn1 <- (template1na$copynumbers[seq(primer,allsegments[i])]-ns1)/sd1
        cn2 <- (template2na$copynumbers[seq(primer,allsegments[i])]-ns2)/sd2
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
      for (i in seq_along(allsegments)) {
        cn1 <- (template1na$copynumbers[seq(primer,allsegments[i])]-ns1)/mad1
        cn2 <- (template2na$copynumbers[seq(primer,allsegments[i])]-ns2)/mad2
        if (Num_Bins[i]>1) {p_value[i] <- wilcox.test(cn1,cn2)$p.value
        } else {p_value[i]<-1}
        primer <- allsegments[i]+1
      }
    }
  } else {
    bincounter <- 1
    segmentcounter <- 0
    if(!is(equalsegments, "numeric")){equalsegments<-20}
    if(onlyautosomes==TRUE) {
      rlechr <- rle(as.vector(template1$chr[template1$chr %in% seq(1,22)]))
    }	else if (onlyautosomes==FALSE) {rlechr <- rle(as.vector(template1$chr))
    } else {rlechr <- rle(as.vector(template1$chr[template1$chr %in% seq(1,onlyautosomes)]))}
    
    for (c in rlechr$values) {
      nos <- floor(length(template1na$chr[(template1na$chr==c)])/equalsegments)
      leftovers <- length(template1na$chr[(template1na$chr==c)])%%equalsegments
      nobs <- rep(equalsegments,nos)
      if (nos == 0) {
        nos <- 1
        nobs <- leftovers
        leftovers <- 0
      }
      if (leftovers > 0) {
        for (a in seq(0, leftovers-1)) {
          divvy <- a%%nos+1
          nobs[divvy] <- nobs[divvy] + 1
        }
      }
      for (i in seq(1, nos)) {
        if(nobs[i]) {
          Chromosome[segmentcounter + i] <- as.vector(template1na$chr)[bincounter]
          Start[segmentcounter + i] <- as.vector(template1na$start)[bincounter]
          End[segmentcounter + i] <- as.vector(template1na$end)[bincounter+nobs[i]-1]
          Num_Bins[segmentcounter + i] <- nobs[i]
          Mean1[segmentcounter + i] <- mean(template1na$copynumbers[seq(bincounter, bincounter+nobs[i]-1)])
          SE1[segmentcounter + i] <- sd(template1na$copynumbers[seq(bincounter, bincounter+nobs[i]-1)])/sqrt(nobs[i])
          Mean2[segmentcounter + i] <- mean(template2na$copynumbers[seq(bincounter, bincounter+nobs[i]-1)])
          SE2[segmentcounter + i] <- sd(template2na$copynumbers[seq(bincounter, bincounter+nobs[i]-1)])/sqrt(nobs[i])
          p_value[segmentcounter + i] <- t.test(template1na$copynumbers[seq(bincounter, bincounter+nobs[i]-1)],
                                                template2na$copynumbers[seq(bincounter, bincounter+nobs[i]-1)])$p.value
          bincounter <- bincounter + nobs[i]
        }
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
      for (i in seq_along(Mean1)) {
        firstbin <- which(template1na$chr==Chromosome[i] & template1na$start==Start[i])
        lastbin <- which(template1na$chr==Chromosome[i] & template1na$end==End[i])
        cn1 <- (template1na$copynumbers[seq(firstbin, lastbin)]-ns1)/sd1
        cn2 <- (template2na$copynumbers[seq(firstbin, lastbin)]-ns2)/sd2
        if (Num_Bins[i]>1) {p_value[i] <- t.test(cn1,cn2)$p.value
        } else {p_value[i]<-1}
      }
    }
    if(altmethod=="MAD"){
      ns1 <- median(rep(Mean1,Num_Bins))
      mad1 <- mad(rep(Mean1,Num_Bins))
      Mean1 <- (Mean1-ns1)/mad1
      SE1 <- SE1/mad1
      ns2 <- median(rep(Mean2,Num_Bins))
      mad2 <- mad(rep(Mean2,Num_Bins))
      Mean2 <- (Mean2-ns2)/mad2
      SE2 <- SE2/mad2
      for (i in seq_along(Mean1)) {
        firstbin <- which(template1na$chr==Chromosome[i] & template1na$start==Start[i])
        lastbin <- which(template1na$chr==Chromosome[i] & template1na$end==End[i])
        cn1 <- (template1na$copynumbers[seq(firstbin,lastbin)]-ns1)/mad1
        cn2 <- (template2na$copynumbers[seq(firstbin,lastbin)]-ns2)/mad2
        if (Num_Bins[i]>1) {p_value[i] <- wilcox.test(cn1,cn2)$p.value
        } else {p_value[i]<-1}
      }
    }
  }
  q_value <- p.adjust(p_value,method="BH")
  combinedsegmentsdf <- data.frame(Chromosome,Start,End,Num_Bins,Mean1,SE1,Mean2,SE2,p_value,q_value)
  correlation <- cor(rep(Mean1,Num_Bins),rep(Mean2,Num_Bins))
  if (plot==TRUE) {
    if(bottom=="min"){bottom<-floor(min(Mean1, Mean2))}
    if(cap=="max"){cap<-ceiling(max(Mean1, Mean2))}
    q_capped <- sapply(q_value, function(x) max(10^-qcap,x))
    q_corrected <- -log(q_capped,10)*cap/qcap
    df <- data.frame(bin=template1na$bin,chr=template1na$chr,Mean1=rep(Mean1,Num_Bins),Mean2=rep(Mean2,Num_Bins),q_value=rep(q_corrected,Num_Bins))
    if(!missing(chrsubset)) {
      firstchr <- range(chrsubset)[1]
      lastchr <- range(chrsubset)[2]
      rlechr <- rle(as.vector(template1$chr))
      df <- df[df$chr %in% rlechr$values[seq(firstchr,lastchr)],]
    }
    subsetcorrelation <- cor(df$Mean1,df$Mean2)
    if(onlyautosomes==TRUE) {
      rlechr <- rle(as.vector(template1$chr[template1$chr %in% seq(1,22)]))
    }	else if (onlyautosomes==FALSE) {rlechr <- rle(as.vector(template1$chr))
    } else {rlechr <- rle(as.vector(template1$chr[template1$chr %in% seq(1,onlyautosomes)]))}
    binchrend <- c()
    currentbin <- 0
    binchrmdl <- c()
    if(altmethod==FALSE){yname<-"copies"} else {yname <- paste0(altmethod," units")}
    for (i in seq_along(rlechr$values)) {
      currentmiddle <- currentbin+rlechr$lengths[i]/2
      currentbin <- currentbin+rlechr$lengths[i]
      binchrend <- append(binchrend, currentbin)
      binchrmdl <- append(binchrmdl, currentmiddle)
    }
    if(missing(chrsubset)){
      compareplot <- ggplot2::ggplot() +
        scale_y_continuous(name = yname, limits = c(bottom,cap), breaks = seq(bottom,cap), expand=c(0,0), 
                           sec.axis = sec_axis(~.*qcap/cap, name = "-log10(q-value)")) +
        scale_x_continuous(name = "chromosome", limits = c(0,tail(binchrend,1)), 
                           breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
        geom_bar(aes(x=bin, y = q_value), data=df, fill='green', stat='identity') +
        geom_hline(yintercept = seq(0,4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = seq(5,cap-1), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = Mean1),data=df, size = 1, color = 'red') +
        geom_point(aes(x = bin,y = Mean2),data=df, size = 1, color = 'blue') +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
          axis.text = element_text(color='black'))
      
      if (legend==TRUE){
        compareplot <- compareplot +
          geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-1.2, ymax = cap), fill = 'white') +
          annotate("text", x = binchrend[2], y = cap-0.2, label = name1, color = 'red') +
          annotate("text", x = binchrend[2], y = cap-0.6, label = name2, color = 'blue')
        if (showcorrelation==TRUE) {
          compareplot <- compareplot +
            annotate("text", x = binchrend[2], y = cap-1, label = paste0("r = ",round(correlation,digits=3)), color = 'black')
        }
      }
    } else {
      firstchr <- range(chrsubset)[1]
      lastchr <- range(chrsubset)[2]
      if(firstchr==1){firstbin<-0
      } else {firstbin<-binchrend[firstchr-1]+1}
      compareplot <- ggplot2::ggplot() +
        scale_y_continuous(name = yname, limits = c(bottom,cap), breaks = seq(bottom,cap), expand=c(0,0), 
                           sec.axis = sec_axis(~.*qcap/cap, name = "-log10(q-value)")) +
        scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), 
                           breaks = binchrmdl[seq(firstchr,lastchr)], 
                           labels = rlechr$values[seq(firstchr,lastchr)], expand = c(0,0)) +
        geom_bar(aes(x=bin, y = q_value), data=df, fill='green', stat='identity') +
        geom_hline(yintercept = seq(0,4), color = '#333333', size = 0.5) +
        geom_hline(yintercept = seq(5,cap-1), color = 'lightgray', size = 0.5) +
        geom_vline(xintercept = binchrend[seq(firstchr,lastchr)], color = "#666666", linetype = "dashed") +
        geom_point(aes(x = bin,y = Mean1),data=df, size = 1, color = 'red') +
        geom_point(aes(x = bin,y = Mean2),data=df, size = 1, color = 'blue') +
        theme_classic() + theme(
          axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black'))
      
      if (legend==TRUE){
        compareplot <- compareplot +
          geom_rect(aes(xmin=firstbin+1, xmax = firstbin+(binchrend[lastchr]-firstbin)/3.5, ymin = cap-1.2, ymax = cap), fill = 'white') +
          annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-0.2, label = name1, color = 'red') +
          annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-0.6, label = name2, color = 'blue')
        if(showcorrelation==TRUE) {
          compareplot <- compareplot +
            annotate("text", x = firstbin+(binchrend[lastchr]-firstbin)/7, y = cap-1, label = paste0("r = ",round(subsetcorrelation,digits=3)), color = 'black')
        }
      }
    }  
    
    if (!missing(chrsubset)) {
      return(list(twosampledf=combinedsegmentsdf,
                  correlation=correlation,
                  subsetcorrelation=subsetcorrelation,
                  compareplot=compareplot))
    } else {return(list(twosampledf=combinedsegmentsdf,
                        correlation=correlation,
                        compareplot=compareplot))}
  } else if (showcorrelation==TRUE) {return(list(twosampledf=combinedsegmentsdf,
                                                 correlation=correlation))
  } else {return(combinedsegmentsdf)}
}

# This code can be used to more systematically use the squaremodel function.
# squaremodelsummary() returns a list of plots starting with the matrix plots and followed by up to 7 CNPs of the best fits.
# Beside the squaremodel, it also requires the sample data, either as template dataframe or as QDNAseq-object with the appropriate sample index
# You can save the plots to a variable or directly print them (which currently is the default!).
# I figured a 2x4 summary page would be nice. You can't currently specify this with parameters, 
# so if you want to adjust this you have to manually change the code below.

squaremodelsummary <- function(template,QDNAseqobjectsample=FALSE,squaremodel,samplename,printplots=TRUE,outputdir,imagetype='pdf',trncname=FALSE) {
  binsize <- 0
  if (QDNAseqobjectsample) {
    fd <- Biobase::fData(template)
    pd <- Biobase::pData(template)
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
    for (i in seq(1,cnplots)) {
      index <- tail(which(squaremodel$minimadf$error==unique(sort(squaremodel$minimadf$error))[i]),1)
      plots[[1+i]] <- singleplot(template=template,QDNAseqobjectsample=QDNAseqobjectsample,
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
loopsquaremodel <- function(object,ptop=5,pbottom=1,prows=100,method='RMSE',penalty=0,penploidy=0,
                            outputdir,imagetype='pdf',trncname=FALSE, returnmodels = FALSE, 
                            printplots = TRUE, printobjectsummary=TRUE) {
  
  imagefunction <- get(imagetype)
  if (missing(outputdir)) {outputdir<-"."}
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  
  fd <- Biobase::fData(object)
  pd <- Biobase::pData(object)
  if(trncname==TRUE) {pd$name <- gsub("_.*","",pd$name)}
  if(trncname!=FALSE&&trncname!=TRUE) {pd$name <- gsub(trncname,"",pd$name)}
  binsize <- fd$end[1]/1000
  
  if (returnmodels[1] != FALSE) {sqmlist <- vector(mode = 'list', length = length(pd$name))}
  matrixplots <- vector(mode = 'list', length = length(pd$name))
  
  for (i in seq_along(pd$name)) {
    model <- squaremodel(object,QDNAseqobjectsample=i,ptop=ptop,pbottom=pbottom,prows=prows,
                         penalty=penalty,penploidy=penploidy,method=method)
    sms <- squaremodelsummary(object,model,QDNAseqobjectsample=i,outputdir=outputdir,imagetype=imagetype,trncname=trncname,printplots=printplots)
    matrixplots[[2*(i-1)+1]] <- model$matrixplot + ggtitle(paste0(pd$name[i]))
    matrixplots[[2*(i-1)+2]] <- sms[[2]]
    if (returnmodels[1] !=  FALSE) {
      sqmlist[[i]]$samplename <- pd$name[i]
      if (returnmodels[1] != TRUE) {model <- model[name = returnmodels]}
      sqmlist[[i]] <- append(sqmlist[[i]], model)
    }
    
  }
  
  if(printobjectsummary == TRUE) {
    if(imagetype == 'pdf') {
      imagefunction(file.path(outputdir,paste0("objectsummary.",imagetype)),width=10.5)
      print(matrixplots)
      dev.off()
    } else {
      imagefunction(file.path(outputdir,paste0("objectsummary.",imagetype)),width=1440,height = 480*ceiling(length(matrixplots)/2))
      print(multiplot(plotlist = matrixplots, cols=2))
      dev.off()
    }
  }
  if (returnmodels[1] != FALSE) {return(sqmlist)}
}

# The first function, correlationmatrix(), is straightforward: it makes a correlation matrix
# of all samples in a QDNAseq-object. The correlations are based on the correlation of segment
# values of all bins. But this analysis may be hampered by poor segmentation! Instead, you could
# just take the raw copynumber values, but there is too much noise that is relatively equal between
# samples. You would, for instance, see pretty good correlation between two samples with absolutely 
# no copy numer aberrations. That is where the correlationmatrixadjusted() comes in. It either
# gives everything equal segments, or it performs pair-wise matching of segment break points for 
# all pairs. It therefore calls the templatefromequalsegments() and twosamplecompare() functions, 
# respectively. 

correlationmatrix <- function(object, trncname=FALSE) {
  samples <- Biobase::sampleNames(object)
  if (trncname==TRUE) {samples <- gsub("_.*","",samples)}
  if (trncname!=FALSE&&trncname!=TRUE) {samples <- gsub(trncname,"",samples)}
  cormat <- cor(na.exclude(assayData(object)$segmented))
  rownames(cormat) <- samples
  colnames(cormat) <- samples
  return(cormat)
}

correlationmatrixadjusted <- function(object, trncname=FALSE, equalsegments=FALSE, funtype = 'mean') {
  samples <- Biobase::sampleNames(object)
  n <- length(samples)
  if (trncname==TRUE) {samples <- gsub("_.*","",samples)}
  if (trncname!=FALSE&&trncname!=TRUE) {samples <- gsub(trncname,"",samples)}
  if (equalsegments==FALSE) {
    cormat <- matrix(nrow=n,ncol=n)
    colnames(cormat) <- samples
    rownames(cormat) <- samples
    cormat[n,n] <- 1
    for (i in seq(1, n-1)) {
      cormat[i,i] <- 1
      for (j in seq(i+1, n)) {
        cormat[i,j] <- twosamplecompare(object,i,index2=j)$correlation
        cormat[j,i] <- cormat[i,j]
      }
    }
  } else {
    size <- length(fData(object)$chromosome)
    tempmat <- matrix(nrow = size, ncol = n)
    for (i in seq_along(samples)) {
      tempmat[,i] <- templatefromequalsegments(object,i,equalsegments=equalsegments,funtype=funtype)$segments
    }
    segmentdf <- na.exclude(as.data.frame(tempmat))
    colnames(segmentdf) <- samples
    cormat <- cor(segmentdf)
  }
  
  return(cormat)
}

# This function creates a template (used in many ACE functions) from a segment
# file (created by many other programs, and available for download from the TCGA
# amongst others!!!). The argument segmentdf can either be a data frame with
# segmented data or the location of a tab-delimited text-file. It is important
# that chromosomes are either integers or chr followed by the number of the
# chromosome (e.g. chr1) It doesn't hurt if the file contains segmented data on
# X and Y chromosome, as long as they are listed after the autosomes. The
# function expects the following order of columns: chromosome, start, end,
# number of bins, segment value. You can use the arguments of the function to
# specify the index of the relevant columns, if necessary The argument log
# actually transforms log-values back to linear scaling; if set to TRUE, it
# converts back from natural logarithm. Alternatively you can specify the base
# of the log-conversion (commonly 2). The standard error and standard deviation
# columns are "cheat" arguments. Segment files do not specify the copynumber
# value per bin, which is one of the columns that the template needs. If
# omitted, copynumbers will take the same values as segments, and will disappear
# behind the segment values in the copy number plot. It is possible that
# standard error or standard deviation data is available per segment. In that
# case, you can create "simulated" copynumber values that are taken from a
# normal distribution corresponding to the segment mean and the given standard
# deviation or standard error. Obviously, this is merely meant as a
# visualization aide!!!

segmentstotemplate <- function(segmentdf, chrci=1, startci=2, endci=3, binsci=4, meanci=5, seci, sdci, log=FALSE) {
  if(is(segmentdf, "character")) {segmentdf <- try(read.table(segmentdf, header = TRUE, comment.char = "", sep = "\t"))}
  if (inherits(segmentdf, "try-error")) {print(paste0("could not find ", segmentdf))
  } else {
    if(log==TRUE){log<-exp(1)}
    if(log) {
      segmentdf[,meanci] <- log^segmentdf[,meanci]
    }
    segments <- rep(segmentdf[,meanci],segmentdf[,binsci])
    chr <- rep(segmentdf[,chrci],segmentdf[,binsci])
    chr <- gsub("chr","",chr,ignore.case = TRUE)
    if(!missing(sdci)) {
      copynumbers <- c()
      for (i in seq_along(segmentdf[,1])) {copynumbers <- append(copynumbers,rnorm(n=segmentdf[i,binsci],
                                                                                  mean=segmentdf[i,meanci],
                                                                                  sd=segmentdf[i,sdci]))}
    } else if(!missing(seci)) {
      copynumbers <- c()
      for (i in seq_along(segmentdf[,1])) {copynumbers <- append(copynumbers,rnorm(n=segmentdf[i,binsci],
                                                                                  mean=segmentdf[i,meanci],
                                                                                  sd=(segmentdf[i,seci]*sqrt(segmentdf[i,binsci]))))}
    } else {
      copynumbers <- segments
    }
    bin <- seq_along(segments)
    start <- c()
    end <- c()
    for (i in seq_along(segmentdf[,1])) {
      start <- append(start,round(seq(from=segmentdf[i,startci],by=(segmentdf[i,endci]+1-segmentdf[i,startci])/segmentdf[i,binsci],
                                length.out=segmentdf[i,binsci])))
      end <- append(end,round(seq(to=segmentdf[i,endci],by=(segmentdf[i,endci]+1-segmentdf[i,startci])/segmentdf[i,binsci],
                            length.out=segmentdf[i,binsci])))
    }
    template <- data.frame(bin,chr,start,end,copynumbers,segments)
    return(template)
  }
}

compresstemplate <- function(template, factor = 20, funtype = 'median'){
  fun <- get(funtype)
  chr <- c()
  start <- c()
  end <- c()
  bin <- c()
  copynumbers <- c()
  segments <- c()
  bincounter <- 1
  nbcounter <- 1
  
  rlechr <- rle(as.vector(template$chr))
  
  for (c in seq_along(rlechr$values)) {
    nob <- rlechr$lengths[c]
    nnb <- floor(nob/factor)
    vob <- rep(factor, nnb)
    leftovers <- nob%%factor
    if (leftovers > 0) {
      for (a in seq(0, leftovers-1)) {
        divvy <- a%%nnb+1
        vob[divvy] <- vob[divvy] + 1
      }
    }
    for (i in (vob)) {
      bin[nbcounter] <- nbcounter
      chr[nbcounter] <- rlechr$values[c]
      start[nbcounter] <- template$start[bincounter]
      end[nbcounter] <- template$end[bincounter+i-1]
      copynumbers[nbcounter] <- fun(template$copynumbers[seq(bincounter, bincounter+i-1)], na.rm = TRUE)
      segments[nbcounter] <- fun(template$segments[seq(bincounter, bincounter+i-1)], na.rm = TRUE)
      nbcounter <- nbcounter + 1
      bincounter <- bincounter + i
    }
  }
  compressedtemplate <- data.frame(bin,chr,start,end,copynumbers,segments)
  return(compressedtemplate)
}

# This function was inspired by twosamplecompare, which also has the option to
# "resegment" copy number data using equally sized segments. I figured the
# option to do this was useful enough to make it its own function, since the
# output can be used as input for other functions (singlemodel, squaremodel,
# singleplot, ACEcall) You obviously want equal segments if you use this
# function ... the argument equalsegments specifies how many bins (or probes)
# should be in this segment. Each chromosome is divided into segments of roughly
# this number of bins. The function expects an ACE template, but obviously the
# value of segments is not necessary. You can also feed it a data frame which
# contains the columns "bin", "chr", "start", "end", "copynumbers", with
# copynumbers being the value of a bin or probe. The columns "bin", "start", and
# "end" are not used, but they are integral to templates for other ACE
# functions.

templatefromequalsegments <- function(template, QDNAseqobjectsample = FALSE, equalsegments = 20, funtype = 'mean', chrsubset, onlyautosomes = TRUE) {
  fun <- get(funtype)
  if(QDNAseqobjectsample) {
    template <- objectsampletotemplate(template, QDNAseqobjectsample)
  } else { template$chr <- gsub("chr","",template$chr,ignore.case = TRUE)} 
  copynumbers <- as.vector(template$copynumbers)
  if("segments" %in% colnames(template)) {
    segments <- as.vector(template$segments)
    } else {segments <- copynumbers}
  chr <- as.vector(template$chr)
  start <- as.vector(template$start)
  end <- as.vector(template$end)
  bin <- as.vector(template$bin)
  
  if(onlyautosomes==TRUE) {
    rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1,22)]))
  }	else if (onlyautosomes==FALSE) {rlechr <- rle(as.vector(template$chr))
  } else {rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1, onlyautosomes)]))}
  
  if(!missing(chrsubset)){
    chromosomes <- rlechr$values[chrsubset]
  } else {chromosomes <- rlechr$values}

  for (c in chromosomes) {
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
      for (a in seq(0, leftovers-1)) {
        divvy <- a%%nos+1
        nobs[divvy] <- nobs[divvy] + 1
      }
    }
    
    for (i in seq(1, nos)) {
      if(nobs[i]) {
        segments[which(!is.na(copynumbers)&chr==c)[seq(bincounter, bincounter+nobs[i]-1)]] <- rep(fun(copynumbers[which(!is.na(copynumbers)&chr==c)[seq(bincounter, bincounter+nobs[i]-1)]]),nobs[i])
        bincounter <- bincounter + nobs[i]
      }
    }
    
  }
  
  template <- data.frame(bin,chr,start,end,copynumbers,segments)
  return(template)
  
}

# Another function inspired by twosamplecompare! The idea of this one is
# actually more basic. The name of the function says it all. This time there are
# two types of input for the segments and an important option. First off: you
# can give your segments as something like a GRanges object or a data frame
# obtained from getadjustedsegments, or you can use an ACE template and use that
# (similar as in twosamplecompare). Second, you now have the option to either
# only use the segments from the segmentinput, or you combine segments as in
# twosamplecompare. Another nice feature of this function when compared to
# twosamplecompare is that you do not have to have equal number of bins, since
# bins of the segmentinput are ignored: it will just use segment info with
# "Chromosome", "Start", and "End" specified in aptly named columns. Since you
# are forcing segments on a template, the template does not need to have segment
# values itself if combinesegments = FALSE. This functions returns the
# resegmented template.

forcesegmentsontemplate <- function(segmentinput, template, QDNAseqobjectsample = FALSE, 
                                    combinesegments = FALSE, funtype = 'mean') {
    
    fun <- get(funtype)
    if (QDNAseqobjectsample) {template <- objectsampletotemplate(template, QDNAseqobjectsample)}
    if (is(segmentinput, "GRanges")) {
        segmentdf <- data.frame(Chromosome = GenomicRanges::seqnames(segmentinput),
                                Start = GenomicRanges::start(segmentinput),
                                End = GenomicRanges::end(segmentinput))
    } else if (colnames(segmentinput)[1] == "bin") {
        segmentdf <- getadjustedsegments(segmentinput)
    } else {
        Chromosome <- segmentinput[,grep("chr", colnames(segmentinput), ignore.case = TRUE)[1]]
        Start <- segmentinput[,grep("start", colnames(segmentinput), ignore.case = TRUE)[1]]
        End <- segmentinput[,grep("end", colnames(segmentinput), ignore.case = TRUE)[1]]
        segmentdf <- data.frame(Chromosome, Start, End)
    }
    
    template$chr <- as.vector(gsub("chr", "", template$chr,ignore.case = TRUE))
    segmentdf$Chromosome <- as.vector(gsub("chr", "", segmentdf$Chromosome,ignore.case = TRUE))
    
    if (combinesegments == TRUE) {
        ownsegments <- getadjustedsegments(template)
        allsegments <- rbind(segmentdf[, seq(1, 3)], ownsegments[, seq(1, 3)])
        allsegments <- allsegments[order(allsegments$Chromosome, allsegments$Start, allsegments$End),]
        Chromosome <- c()
        Start <- c()
        End <- c()
        for (c in rle(as.vector(allsegments$Chromosome))$values) {
            starts <- sort(unique(allsegments$Start[allsegments$Chromosome == c]))
            ends <- sort(unique(allsegments$End[allsegments$Chromosome == c]))
            if (length(starts != length(ends))) {
                starts <- append(starts, ends[-length(ends)]+1)
                ends <- append(ends, starts[-1]-1)
                starts <- sort(unique(starts))
                ends <- sort(unique(ends))
            }
            Chromosome <- append(Chromosome, rep(c, length(starts)))
            Start <- append(Start, starts)
            End <- append(End, ends)
        }
        segmentdf <- data.frame(Chromosome, Start, End)
    }
    
    template$segments <- NA
    
    for (s in seq_along(segmentdf$Chromosome)) {
        indices <- which(template$chr == segmentdf$Chromosome[s] &
                             template$start >= segmentdf$Start[s] &
                             template$end <= segmentdf$End[s])
        template$segments[indices] <- fun(na.exclude(template$copynumbers[indices]))
    }
    template$segments[is.na(template$copynumbers)] <- NA
    return(template)
}