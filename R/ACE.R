##### ACE has arrived! #####

# There are a number of functions: ACE, ploidyplotloop, objectsampletotemplate,
# singlemodel, squaremodel, singleplot, getadjustedsegments, linkvariants,
# postanalysisloop, analyzegenomiclocations. I have also added the "multiplot"
# function for making summary sheets.

# ACE takes a folder with rds-files (default, make sure all are segmented) or
# bam-files and returns plots for the most likely errors in convenient
# subfolders. In case of bam-files it will run binsizes 30, 100, 500 and 1000
# kbp as default, but a vector of desired binsizes can be used as input. Bin
# sizes available are 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp. Note that
# output will be larger and programs will run longer with the smaller binsizes
# (1 kbp is pretty hardcore). ploidies specifies for which ploidies fits will be
# made: default is 2 but you can input a vector with all desired ploidies, e.g.
# c(2,4). Beggers be choosers! ACE likes pdf-files, but for large datasets or
# small binsizes you might want to go with 'png'. A standard method for error
# calculation is the root mean squared error (RMSE), but you can argue you
# should actually punish more the long segments that are a little bit off
# compared to the short segments that are way off. 'MAE' weighs every error
# equally, whereas 'SMRE' does the latter; these will generally generate more
# fits (I see that as a downside). but it is more likely that the segments of
# which you are sure are sticking tighter to the integer ploidy in the best fit.
# The parameter "penalty" sets a penalty for the error calculated at lower
# cellularity. Errors are divided by the cellularity to the power of the
# penalty: default = 0 (no penalty). The parameters cap and bottom determine the
# limits of the y-axis on the resulting copy number plots.

# New functionality for large data sets (many samples!): Fits should be chosen
# manually, but ... you can now open the file with all errorlists to pick the
# most likely fit without looking at the profile. Use the generated
# tab-delimited file (fitpicker.tsv) in the ploidy subdirectory to fill in the
# column with the likely fit. If you set the argument autopick to TRUE, ACE will
# fill in the likely fit column ... thanks for the vote of confidence! Added
# argument trncname (truncate name) which truncates the name to everything
# before the first instance of _ - set TRUE if applicable, or specify regular
# expression, e.g. "-.*" to REMOVE everything after and including the first dash
# found in a sample name. Added argument printsummaries: superbig summary files
# may crash the program, so you can set this argument to FALSE or 2 if you still
# want the error plots. ACE will create a file (parameters.tsv) in the input
# directory, specifying the supplied (or default when omitted) parameters.

# ploidyplotloop is the meat of ACE, and it can be run as a separate function as
# well; it takes a QDNAseq-object as input and also needs a folder to write the
# files to. This is particularly handy if you want to analyze a whole
# QDNAseq-object that is loaded in your R environment.

# objectsampletotemplate is pretty much there to parse QDNAseqobjects into the
# dataframe structure used by the singlemodel and singleplot functions. These
# latter functions call objectsampletotemplate itself when necessary, but it can
# be handy to make a template if you expect some repeated use of the functions
# or if you want to make the template and then do your own obscure manipulations
# to it (I know you want to!).

# singlemodel: like the name implies it runs the fitting algorithm on a single
# sample and spits out the info you want! Not strictly necessary to save it to a
# variable, but that might still be handy. singlemodel comes with the
# awesomeness of manual input: you can restrain the model to the ploidy you
# expect (default 2, but hey, ploidy 5 happens, right?) and you can tell the
# program if you think there is a different standard (this should only be
# necessary when the median segment happens to be a subclonal variant; very
# rare).

# squaremodel: runs the algorithm for each tested ploidy being the standard. If
# you feel you are doing too much tinkering with variables using the singlemodel
# function, then try this! You can choose the range of ploidies it should test;
# the default being all ploidies between 5 and 1 in 100 decrements of 0.04. Like
# singlemodel, the function returns a list, which you can save to a variable.

# singleplot: again, you already know what it does! But did you know you can
# change the chart title to anything you like? Well, I'm here to tell you it
# can. And there is more! Don't believe the cellularities the program calculated
# for you? You don't have to! Just fill in the cellularity you believe to be
# true, and singleplot will take your word for it!

# getadjustedsegments: get info of the actual segments in a handy dataframe!
# Contains start and end location, "number of probes" (number of bins supporting
# the segment value), value of the segment (Segment_Mean is the value from
# QDNAseq, Segment_Mean2 is calculated from adjustedcopynumbers within the
# segment), and the nearest ploidy with the chance (log10 p-value) that the
# segment has this ploidy. A very low p-value indicates a high chance of
# subclonality, although this value should be approached with extreme caution
# (and is not yet corrected for multiple testing)

# linkvariants: now we're talking! Give a tab-delimited file with variant data,
# and this function will tell you what the copy number is at that genomic
# location, and it will guess how many copies are mutant! Read more about the
# options for this function in the code. Now includes an option to analyze
# germline heterozygous SNPs! Also make sure to indicate if X and Y are single
# germline copy (sgc = c("X","Y")) if you are linking variants on these
# chromosomes in male subjects.

# postanalysisloop: you've run ACE, you've picked your models, you have your
# mutation data ... Now if you could only ... Say no more! This function
# combines the power of all above functions and your own brain (you still choose
# the models). For the brave of heart.

# analyzegenomiclocations: a sort of simplified version of linkvariants, you can
# manually input a single genomic location as specified by chromosome and
# position you can also input multiple locations using vectors of the same
# length. The output will be a data frame. You can also enter frequencies to
# quickly calculate mutant copies, but then you also need to enter cellularity!
# Note: the function requires the data frame with adjusted segments (output of
# getadjustedsegments)

# Recently a lot of functions were complemented with the ability to take into
# account single germline copies of chromosomes. This is of course particularly
# relevant for the sex chromosomes in males. Firstly this was necessary to get
# correct values for these chromosomes, but secondly it can even be helpful in
# deciding the correct ploidy of a tumor when it is used in the fitting
# functions! For now it is too involving to include this option in analyses for
# entire objects, such as runACE and loopsquaremodel, so it is only available in
# single sample analyses. postanalysisloop does have an option for it, because
# this is post model fitting. Please consult the documentation.

# That's pretty much it for now, let me know if you run into some errors or
# oddities j.poell@amsterdamumc.nl

runACE <- function(inputdir = "./", outputdir, filetype = 'rds', genome = "hg19", binsizes, ploidies = 2, 
                   imagetype = 'pdf', method = 'RMSE', penalty = 0, cap = 12, bottom = 0, 
                   trncname = FALSE, printsummaries = TRUE, savereadcounts = FALSE, autopick = FALSE) { 
	imagefunction <- get(imagetype)
	if(substr(inputdir,nchar(inputdir),nchar(inputdir)) != "/") {inputdir <- paste0(inputdir,"/")}
	if(missing(outputdir)) { outputdir <- substr(inputdir,0,nchar(inputdir)-1) }
	if(!dir.exists(outputdir)) {dir.create(outputdir)}
	if(filetype=='bam'){
		if(missing(binsizes)) { binsizes <- c(100,500,1000) }
	  parameters <- data.frame(options = c("inputdir","outputdir","filetype","binsizes","ploidies","imagetype","method","penalty","cap","bottom","trncname","printsummaries","autopick"), 
	                           values = c(inputdir,outputdir,filetype,paste0(binsizes,collapse=", "),paste0(ploidies,collapse=", "),imagetype,method,penalty,cap,bottom,trncname,printsummaries,autopick))
	  for (b in binsizes) {
		  currentdir <- file.path(outputdir,paste0(b,"kbp"))
		  dir.create(currentdir)
		  bins <- QDNAseq::getBinAnnotations(binSize = b, genome = genome)
		  readCounts <- QDNAseq::binReadCounts(bins, path = inputdir)
		  if (savereadcounts == TRUE) {
		    saveRDS(readCounts, file = file.path(outputdir, paste0(b, "kbp-raw.rds")))
		  }
		  readCountsFiltered <- QDNAseq::applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
		  readCountsFiltered <- QDNAseq::estimateCorrection(readCountsFiltered)
		  # the default correctBins will output a ratio; using method = 'median' will return a corrected readcount
		  copyNumbers <- QDNAseq::correctBins(readCountsFiltered)
		  copyNumbers <- QDNAseq::normalizeBins(copyNumbers)
		  copyNumbers <- QDNAseq::smoothOutlierBins(copyNumbers)
		  copyNumbersSegmented <- QDNAseq::segmentBins(copyNumbers, transformFun = 'sqrt') # the transformFun is not available in older versions of QDNAseq!
		  copyNumbersSegmented <- QDNAseq::normalizeSegmentedBins(copyNumbersSegmented)
		  saveRDS(copyNumbersSegmented, file = file.path(outputdir,paste0(b,"kbp.rds")))
		  
		  ploidyplotloop(copyNumbersSegmented,currentdir,ploidies,imagetype,method,penalty,cap,bottom,trncname,printsummaries,autopick)
		  
		}
	}
	else if(filetype=='rds'){
	  parameters <- data.frame(options = c("inputdir","outputdir","filetype","ploidies","imagetype","method","penalty","cap","bottom","trncname","printsummaries","autopick"), 
	                           values = c(inputdir,outputdir,filetype,paste0(ploidies,collapse=", "),imagetype,method,penalty,cap,bottom,trncname,printsummaries,autopick))
	  files <- list.files(inputdir, pattern = "\\.rds$")
	  for (f in seq_along(files)) {
			currentdir <- file.path(outputdir,paste0(substr(files[f],0,nchar(files[f])-4)))
			dir.create(currentdir)
			copyNumbersSegmented <- readRDS(file.path(inputdir,files[f]))
			ploidyplotloop(copyNumbersSegmented,currentdir,ploidies,imagetype,method,penalty,cap,bottom,trncname,printsummaries,autopick)
		}
	}
	else {
	print("not a valid filetype")
	}
	write.table(parameters, file=file.path(outputdir,"parameters.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
}

ploidyplotloop <- function(copyNumbersSegmented,currentdir,ploidies=2,imagetype='pdf',method='RMSE',penalty=0,
                           cap=12,bottom=0,trncname=FALSE,printsummaries=TRUE,autopick=FALSE) {
  imagefunction <- get(imagetype)
	fd <- Biobase::fData(copyNumbersSegmented)
	pd <- Biobase::pData(copyNumbersSegmented)
	if(trncname==TRUE) {pd$name <- gsub("_.*","",pd$name)}
	if(trncname!=FALSE&&trncname!=TRUE) {pd$name <- gsub(trncname,"",pd$name)}
	binsize <- fd$end[1]/1000
# I have commented out the best fit and last minimum plots, they are all in likely fits anyway	
#	bfplots <- vector(mode = 'list', length = 4*length(pd$name))
#	lmplots <- vector(mode = 'list', length = 4*length(pd$name))
	
	for (q in ploidies) {
	  qdir <- file.path(currentdir,paste0(q,"N"))
	  dir.create(qdir)
	  
  	likelyplots <- vector(mode = 'list', length = 3*length(pd$name))
  	listerrorplots <- vector(mode = 'list', length = length(pd$name))
  	fitpicker <- matrix(ncol = 16, nrow = length(pd$name))
  	colnames(fitpicker) <- c("sample","likely_fit","ploidy","standard","fit_1","fit_2","fit_3","fit_4","fit_5","fit_6","fit_7","fit_8","fit_9","fit_10","fit_11","fit_12")
  	dir.create(file.path(qdir,"likelyfits"))  
  	
  	for (a in seq_along(pd$name)) {
  		segmentdata <- rle(as.vector(na.exclude(assayData(copyNumbersSegmented)$segmented[,a])))
  		standard <- median(rep(segmentdata$values,segmentdata$lengths))
  			
  		fraction <- c()
  		expected <- c()
  		temp <- c()
  		errorlist <- c()
  		
  		for (i in seq(5,100)) {
  		  fraction[i-4] <- i/100
  		  for (p in seq(1,12)) {
  		    expected[p] <- standard*(p*fraction[i-4] + 2*(1-fraction[i-4]))/(fraction[i-4]*q + 2*(1-fraction[i-4]))
  		  }
  		  # the maximum error 0.5 was added to make sure hyperamplifications (p>12) don't get ridiculous errors
  		  for (j in seq_along(segmentdata$values)) {
  		    if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))^2}
  		    else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))}
  		    else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty)}
  		    else {print("Not a valid method")}
  		  }
  		  if(method=='RMSE') {errorlist[i-4] <- sqrt(sum(temp*segmentdata$lengths)/sum(segmentdata$lengths))}
  		  else if(method=='SMRE') {errorlist[i-4] <- sum(temp*segmentdata$lengths)/sum(segmentdata$lengths)^2}
  		  else if(method=='MAE') {errorlist[i-4] <- sum(temp*segmentdata$lengths)/sum(segmentdata$lengths)}
  		  
  		}
  			
  		minima <- c()
  		rerror <- c()
  		
  		if (round(errorlist[1], digits = 10) < round(errorlist[2], digits = 10)) {
  		  lastminimum <- fraction[1]
  		  minima[1] <- fraction[1]
  		  rerror[1] <- errorlist[1]/max(errorlist)
  		}
  		
  		for (l in seq(6,99)) {
  		  if (round(errorlist[l-4], digits = 10) < round(errorlist[l-5], digits = 10) & round(errorlist[l-4], digits = 10) < round(errorlist[l-3], digits = 10)) { 
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
  		 			
  		expsd <- sqrt(pd$expected.variance[a])
  		obssd <- QDNAseq:::sdDiffTrim(assayData(copyNumbersSegmented)$copynumber[,a], na.rm=TRUE)
  		chr <- fd$chromosome[fd$chromosome %in% seq(1,22)]
  		bin <- seq_along(chr)
  		rlechr <- rle(chr)
  		lastchr <- length(rlechr$values)
  		binchrend <- c()
  		currentbin <- 0
  		binchrmdl <- c()
  		for (i in seq_along(rlechr$values)) {
  		  currentmiddle <- currentbin+rlechr$lengths[i]/2
  		  currentbin <- currentbin+rlechr$lengths[i]
  		  binchrend <- append(binchrend, currentbin)
  		  binchrmdl <- append(binchrmdl, currentmiddle)
  		}
  		
  		# create a subdirectory for your sample
  		fp <- file.path(qdir,pd$name[a])
  		if(!dir.exists(fp)) {
  		  dir.create(fp)
  		}
  			
  		dir.create(file.path(fp,"graphs"))
  		
  		cellularity <- seq(5,100)
  		tempdf <- data.frame(cellularity,errorlist=errorlist/max(errorlist))
  		minimadf <- data.frame(minima=minima*100,rerror)
  		tempplot <- ggplot2::ggplot() +
  		  scale_y_continuous(name = "relative error", limits = c(0,1.05), expand=c(0,0)) +
  		  scale_x_continuous(name = "cellularity (%)") +
  		  geom_vline(xintercept = seq(from = 10, to = 100, by = 10), color = "#666666", linetype = "dashed") +
  		  geom_point(aes(y=errorlist, x=cellularity), data=tempdf) +
  		  geom_point(aes(y=rerror, x=minima), data=minimadf, color = 'red') +
  		  theme_classic() + theme(
  		    axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
  		  ggtitle(paste0(pd$name[a], " - errorlist")) +
  		  theme(plot.title = element_text(hjust = 0.5))
  		
  #		bfplots[[(4*(a-1)+2)]] <- tempplot
  #		lmplots[[(4*(a-1)+2)]] <- tempplot
  		likelyplots[[(3*(a-1)+3)]] <- tempplot
  		listerrorplots[[a]] <- tempplot
  		
  		imagefunction(file.path(fp,paste0(pd$name[a],"_errorlist.",imagetype)))
  		print(tempplot)
  		dev.off()

  		plots <- vector(mode = 'list', length = (length(minima)))
  		
  		bfi <- tail(which(rerror==min(rerror)), n = 1)
  		
  		fitpicker[a,1] <- pd$name[a]
  		if (autopick==TRUE) { fitpicker[a,2] <- minima[bfi] }	else { fitpicker[a,2] <- "" }
  		fitpicker[a,3] <- q
  		fitpicker[a,4] <- standard
  		

  		for (m in seq_along(minima)) {
  		  fitpicker[a,m+4] <- minima[m]
  		  adjustedcopynumbers <- q + ((assayData(copyNumbersSegmented)$copynumber[,a]-standard)*(minima[m]*(q-2)+2))/(minima[m]*standard)
  		  adjustedsegments <- q + ((assayData(copyNumbersSegmented)$segmented[,a]-standard)*(minima[m]*(q-2)+2))/(minima[m]*standard)
  		  df <- as.data.frame(cbind(bin,adjustedcopynumbers[seq_along(bin)],adjustedsegments[seq_along(bin)]))
  		  colnames(df)[2] <- "copynumbers"
  		  colnames(df)[3] <- "segments"
  		  dfna <- na.exclude(df)
  		  copynumbers <- dfna$copynumbers
  		  segments <- dfna$segments
  		  cappedcopynumbers <- dfna[dfna$copynumbers > cap,]
  		  if(length(cappedcopynumbers$copynumbers)>0) {cappedcopynumbers$copynumbers <- cap-0.1}
  		  cappedsegments <- dfna[dfna$segments > cap,]
  		  if(length(cappedsegments$segments)>0) {cappedsegments$segments <- cap-0.1}
  		  toppedcopynumbers <- dfna[dfna$copynumbers <= bottom,]
  		  if(length(toppedcopynumbers$copynumbers)>0) {toppedcopynumbers$copynumbers <- bottom+0.1}
  		  toppedsegments <- dfna[dfna$segments <= bottom,]
  		  if(length(toppedsegments$segments)>0) {toppedsegments$segments <- bottom+0.1}
  		  line1 <- paste0("Cellularity: ", minima[m])
  		  line2 <- paste0("Relative error: ", round(rerror[m], digits = 3))
  		  line3 <- paste0("Expected Noise: ", round(expsd,digits = 4))
  		  line4 <- paste0("Observed Noise: ", round(obssd,digits = 4))
  		  fn <- file.path(fp,"graphs",paste0(pd$name[a], " - ",q,"N fit ", m, ".",imagetype))
  		  
  		  tempplot <- ggplot2::ggplot() +
    			scale_y_continuous(name = "copies", limits = c(bottom,cap), breaks = seq(bottom,cap), expand=c(0,0)) +
    			scale_x_continuous(name = "chromosome", limits = c(0,binchrend[lastchr]), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
    			geom_hline(yintercept = seq(0,4), color = '#333333', linewidth = 0.5) +
    			geom_hline(yintercept = seq(5,cap-1), color = 'lightgray', linewidth = 0.5) +
    			geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
    			geom_point(aes(x = bin,y = copynumbers),data=dfna[dfna$copynumbers>bottom&dfna$copynumbers<cap,], size = 0.1, color = 'black') +
  		    geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
  		    geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
    			geom_point(aes(x = bin,y = segments),data=dfna[dfna$segments>bottom&dfna$segments<cap,], size = 1, color = 'darkorange') +
  		    geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'darkorange', shape = 24) +
  		    geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkorange', shape = 25) +
  		    theme_classic() + theme(
  		      axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
  		    ggtitle(paste0(pd$name[a], " - binsize ", binsize, " kbp - ", pd$used.reads[a], " reads - ",q,"N fit ", m)) +
    			geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-0.9, ymax = cap), fill = 'white') +
    			annotate("text", x = binchrmdl[2], y = cap-0.3, label = line1) +
    			annotate("text", x = binchrmdl[2], y = cap-0.7, label = line2) +
    			geom_rect(aes(xmin=binchrmdl[12], xmax = binchrmdl[lastchr], ymin = cap-0.9, ymax = cap), fill = 'white') +
    			annotate("text", x = binchrmdl[16], y = cap-0.3, label = line3) +
    			annotate("text", x = binchrmdl[16], y = cap-0.7, label = line4) +
  		    theme(plot.title = element_text(hjust = 0.5))
  		  plots[[m]] <- tempplot
  		  if(bfi==m) {
  		    likelyplots[[(3*(a-1)+1)]] <- tempplot
  		    if(imagetype %in% c('pdf','svg')) {
  		      imagefunction(file.path(qdir,"likelyfits",paste0(pd$name[a],"_bestfit_",q,"N.",imagetype)),width=10.5)
  		    } else {
  		      imagefunction(file.path(qdir,"likelyfits",paste0(pd$name[a],"_bestfit_",q,"N.",imagetype)),width=720)
  		    }  
  		    print(tempplot)
  		    dev.off()  
  		  }
  		  if(m==length(minima)) {
  		    likelyplots[[(3*(a-1)+2)]] <- tempplot
  		    if(imagetype %in% c('pdf','svg')) {
  		      imagefunction(file.path(qdir,"likelyfits",paste0(pd$name[a],"_lastminimum_",q,"N.",imagetype)),width=10.5)
  		    } else {
  		      imagefunction(file.path(qdir,"likelyfits",paste0(pd$name[a],"_lastminimum_",q,"N.",imagetype)),width=720)
  		    }
  		    print(tempplot)
  		    dev.off() 
  		  }
  		  if(imagetype %in% c('pdf','svg')) {
  		    imagefunction(fn,width=10.5)
  		  } else {
  		    imagefunction(fn,width=720)
  		  }
  		  print(tempplot)
  		  dev.off()
  		  

  		}
  			
  		if(imagetype %in% c('pdf','svg')) {
  		  pdf(file.path(fp,paste0("summary_",pd$name[a],".pdf")),width=10.5)
  		  print(plots)
  		  dev.off()
  		} else {
  		  if (length(plots)==1) {
  		    imagefunction(file.path(fp,paste0("summary_",pd$name[a],".",imagetype)), width = 720)
  		    print(plots)
  		    dev.off()
  		  } else {
    		  imagefunction(file.path(fp,paste0("summary_",pd$name[a],".",imagetype)), width = 2160, height = 480*ceiling(length(plots)/3))
    		  print(multiplot(plotlist = plots, cols=3)) 
    		  dev.off()
    		}
    		  
  		}
  		  
  	}
	
	
  	
  	if(printsummaries == TRUE) {
    	if(imagetype %in% c('pdf','svg')) {
    	  pdf(file.path(qdir,"summary_likelyfits.pdf"),width=10.5)
      	  print(likelyplots)
      	  dev.off()
      	pdf(file.path(qdir,"summary_errors.pdf"))
      	  print(listerrorplots)
      	  dev.off()
      	} else { 
    	  imagefunction(file.path(qdir,paste0("summary_likelyfits.",imagetype)), width = 2160, height = 480*length(pd$name))
      	  print(multiplot(plotlist = likelyplots, cols=3))
      	  dev.off()
      	imagefunction(file.path(qdir,paste0("summary_errors.",imagetype)), width = 1920, height = 480*ceiling(length(pd$name)/4))
      	if (length(listerrorplots)==1){print(listerrorplots)} else {print(multiplot(plotlist = listerrorplots, cols=4))}
      	  dev.off()
      	}
  	} else if(printsummaries == 2) {
  	  if(imagetype %in% c('pdf','svg')) {
  	    pdf(file.path(qdir,"summary_errors.pdf"))
  	    print(listerrorplots)
  	    dev.off()
  	    } else { 
  	    imagefunction(file.path(qdir,paste0("summary_errors.",imagetype)), width = 1920, height = 480*ceiling(length(pd$name)/4))
  	    if (length(listerrorplots)==1){print(listerrorplots)} else {print(multiplot(plotlist = listerrorplots, cols=4))}
  	      dev.off()
  	   }
  	} 
	
  	write.table(fitpicker, file=file.path(qdir,paste0("fitpicker_",q,"N.tsv")), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
  	
	}
}


# this code makes the generic input template dataframe for singlemodel and singleplot
# those function can take QDNAseq objects when specified, so it is not necessary to run this separately
objectsampletotemplate <- function(copyNumbersSegmented, index = 1) {
  fd <- Biobase::fData(copyNumbersSegmented)
	try(segments <- as.vector(assayData(copyNumbersSegmented)$segmented[,index]))
	copynumbers <- as.vector(assayData(copyNumbersSegmented)$copynumber[,index])
	chr <- as.vector(fd$chromosome)
	start <- as.vector(fd$start)
	end <- as.vector(fd$end)
	bin <- seq_along(chr)
	if(inherits(segments,"NULL")){
	  template <- data.frame(bin,chr,start,end,copynumbers)
	} else {	template <- data.frame(bin,chr,start,end,copynumbers,segments)}
	return(template)
}


# this function takes a template dataframe as specified in the above function
# you can use a QDNAseq object as "template", then specify QDNAseqobjectsample by the sample number!
# e.g. model <- singlemodel(copyNumbersSegmented, QDNAseqobjectsample = 3)
singlemodel <- function(template,QDNAseqobjectsample = FALSE, ploidy = 2, standard, method = 'RMSE', 
                        exclude = c("X","Y"), sgc = c(), penalty = 0, highlightminima = TRUE) {
  if(QDNAseqobjectsample) {template <- objectsampletotemplate(template, QDNAseqobjectsample)}
  template <- template[!template$chr %in% exclude,]
  if(missing(standard) || !is(standard, "numeric")) { standard <- median(template$segments, na.rm = T) }
  templateh <- template[template$chr %in% sgc,]
  template <- template[!template$chr %in% sgc,]
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  segmentdatah <- rle(as.vector(na.exclude(templateh$segments)))
  fraction <- c()
  expected <- c()
  hexpected <- c()
  temp <- c()
  errorlist <- c()
  for (i in seq(5,100)) {
    fraction[i-4] <- i/100
    for (p in seq(1,12)) {
      expected[p] <- standard*(p*fraction[i-4] + 2*(1-fraction[i-4]))/(fraction[i-4]*ploidy + 2*(1-fraction[i-4]))
      hexpected[p] <- standard*(p*fraction[i-4] + 1*(1-fraction[i-4]))/(fraction[i-4]*ploidy + 2*(1-fraction[i-4]))
      #		  expected[p] <- standard*(1+(p-ploidy)*fraction[i-4]/(fraction[i-4]*(ploidy-2)+2))
    }
    for (j in seq_along(segmentdata$values)) {
      if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))^2}
      else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))}
      else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty)}
      else {print("Not a valid method")}
    }
    for (k in seq_along(segmentdatah$values)) {
      if(method=='RMSE') {temp[j+k] <- (min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i-4]^penalty))^2}
      else if(method=='SMRE') {temp[j+k] <- sqrt(min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i-4]^penalty))}
      else if(method=='MAE') {temp[j+k] <- min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i-4]^penalty)}
      else {print("Not a valid method")}
    }
    lengths <- c(segmentdata$lengths, segmentdatah$lengths)
    if(method=='RMSE') {errorlist[i-4] <- sqrt(sum(temp*lengths)/sum(lengths))}
    else if(method=='SMRE') {errorlist[i-4] <- sum(temp*lengths)/sum(lengths)^2}
    else if(method=='MAE') {errorlist[i-4] <- sum(temp*lengths)/sum(lengths)}
    
  }
  
  minima <- c()
  rerror <- c()
  
  if (round(errorlist[1], digits = 10) < round(errorlist[2], digits = 10)) {
    lastminimum <- fraction[1]
    minima[1] <- fraction[1]
    rerror[1] <- errorlist[1]/max(errorlist)
  }
  for (l in seq(6,99)) {
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
  
  cellularity <- seq(5,100)
  tempdf <- data.frame(cellularity,errorlist=errorlist/max(errorlist))
  minimadf <- data.frame(minima=minima*100,rerror)
  if(highlightminima==TRUE) {
    tempplot <- ggplot2::ggplot() +
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
    tempplot <- ggplot2::ggplot() +
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


# Turns out, tumors are often quite messy, genomically speaking
# Occasionally, you might want to be able to compare the fits at 2N directly with the fits at other ploidies
# Also, with pooring quality or lots of subclonality, it is more common for the standard to be messed up
# In line with graphs presented by ASCAT and CELLULOID, squaremodel will do fitting and plot a graph for 
#   fits using two variables: ploidy and cellularity
# On top of the penalty for low cellularity, you can also add a penalty for ploidies
# It is implemented as follows: error*(1+abs(ploidy-2))^penploidy
# The plot has the two variables as axis and the color code indicates the relative error
# To make the minima pop out, the color code is the inverse of the relative error
# Minima are found by checking each value for neighboring values, and will only return true if its the lowest error
squaremodel <- function(template, QDNAseqobjectsample = FALSE, prows=100, ptop=5, pbottom=1, method = 'RMSE',
                        exclude = c("X","Y"), sgc = c(), penalty = 0, penploidy = 0, 
                        cellularities = seq(5,100), highlightminima = TRUE, standard) {
  if(QDNAseqobjectsample) {template <- objectsampletotemplate(template, QDNAseqobjectsample)}
  template <- template[!template$chr %in% exclude,]
  if(missing(standard) || !is(standard, "numeric")) { standard <- median(template$segments, na.rm = T) }
  templateh <- template[template$chr %in% sgc,]
  template <- template[!template$chr %in% sgc,]
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  segmentdatah <- rle(as.vector(na.exclude(templateh$segments)))
  fraction <- sort(cellularities)/100
  error <- c()
  errormatrix <- matrix(nrow=(prows+1),ncol=length(fraction))
  listofploidy <- c()
  listofcellularity <- c()
  listoferrors <- c()
  for (t in seq(0,prows)) {
    ploidy <- ptop-((ptop-pbottom)/prows)*t
    listofploidy <- append(listofploidy, rep(ploidy,length(fraction)))
    expected <- c()
    hexpected <- c()
    temp <- c()
    errorlist <- c()
    for (i in seq_along(fraction)) {
      for (p in seq(1,12)) {
        expected[p] <- standard*(p*fraction[i] + 2*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
        hexpected[p] <- standard*(p*fraction[i] + 1*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
      }
      for (j in seq_along(segmentdata$values)) {
        if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))^2}
        else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))}
        else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty)}
        else {print("Not a valid method")}
      }
      for (k in seq_along(segmentdatah$values)) {
        if(method=='RMSE') {temp[j+k] <- (min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i]^penalty))^2}
        else if(method=='SMRE') {temp[j+k] <- sqrt(min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i]^penalty))}
        else if(method=='MAE') {temp[j+k] <- min(abs(segmentdatah$values[k]-hexpected),0.5)/(fraction[i]^penalty)}
        else {print("Not a valid method")}
      }
      lengths <- c(segmentdata$lengths, segmentdatah$lengths)
      if(method=='RMSE') {errorlist[i] <- sqrt(sum(temp*lengths)/sum(lengths))}
      else if(method=='SMRE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)^2}
      else if(method=='MAE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)}
    }
    listofcellularity <- append(listofcellularity, fraction)
    listoferrors <- append(listoferrors, errorlist)
    errormatrix[t+1,] <- errorlist
    
  }
  minimat <- matrix(nrow=(prows+1),ncol=length(fraction))
  for (i in seq(1,prows+1)) {
    for (j in seq(1,length(fraction))) {
      if (i==1|i==(prows+1)) {
        minimat[i,j] <- FALSE
      } else if (j==length(fraction)) {
        minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j)])
      } else {
        minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j+1)])
      }
    }
  }
  round(errormatrix, digits = 10)
  round(listoferrors, digits = 10)
  errordf <- data.frame(ploidy=listofploidy,
                        cellularity=listofcellularity,
                        error=listoferrors/max(listoferrors),
                        minimum=as.vector(t(minimat)))
  minimadf <- errordf[errordf$minimum==TRUE,]
  minimadf <- minimadf[order(minimadf$error,-minimadf$cellularity),]
  if(highlightminima==TRUE){
    tempplot <- ggplot2::ggplot() +
      geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
      geom_point(data=minimadf, aes(x=cellularity, y=ploidy, alpha=min(error)/error), shape=16) +
      scale_fill_gradient(low="green", high="red") +
      ggtitle("Matrix of errors") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    tempplot <- ggplot2::ggplot() +
      geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
      scale_fill_gradient(low="green", high="red") +
      ggtitle("Matrix of errors") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(list(method=method, 
              penalty=penalty, 
              penploidy=penploidy,
              standard=standard,
              errormatrix=errormatrix, 
              minimatrix = minimat, 
              errordf=errordf, 
              minimadf=minimadf, 
              matrixplot=tempplot))
  
}


# this function takes a template dataframe as specified in objectsampletotemplate
# you can use a QDNAseq object as "template", then specify QDNAseqobjectsample by the sample number!
# don't forget to specify cellularity, error, ploidy, and standard; you can find these in the model output
# chrsubset restricts the plot to the selected CONTIGUOUS chromosomes: the function takes the lowest and highest 
#  integer from the entered vector and plots the range of chromosomes
#  e.g. chrsubset = c(3,6) and chrsubset = 3:6 result in the same plot
#  subsetting chromosomes messes up the position for the legend, 
#  using chrsubset=1:22 is therefore a "hack" to remove the legend but still get the same plot
# you can scale your upper and lower y-axis limit to the data using cap="max" and bottom="min" respectively
singleplot <- function(template, QDNAseqobjectsample = FALSE, cellularity = 1, error, ploidy = 2, standard, title,
                       trncname = FALSE, cap = 12, bottom = 0, chrsubset, onlyautosomes = TRUE, sgc = c()) {
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
  if(missing(standard) || !is(standard, "numeric")) { standard <- median(template$segments[!template$chr %in% c("X", "Y")], na.rm = T) }
  # formula below is a simplification of the one that is commented out
  adjustedcopynumbers <- template$copynumbers*(ploidy+2/cellularity-2)/standard - gc/cellularity + gc
  adjustedsegments <- template$segments*(ploidy+2/cellularity-2)/standard - gc/cellularity + gc
  # adjustedcopynumbers <- (((template$copynumbers/standard)*(cellularity*ploidy+2*(1-cellularity)))-gc*(1-cellularity))/cellularity
  # adjustedsegments <- (((template$segments/standard)*(cellularity*ploidy+2*(1-cellularity)))-gc*(1-cellularity))/cellularity
  template$chr <- gsub("chr","",template$chr,ignore.case = TRUE)
  df <- data.frame(bin=template$bin,chr=template$chr,adjustedcopynumbers,adjustedsegments)
  
  if(onlyautosomes==TRUE) {
    rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1,22)]))
  }	else if (onlyautosomes==FALSE) {rlechr <- rle(as.vector(template$chr))
  } else {rlechr <- rle(as.vector(template$chr[template$chr %in% seq(1,onlyautosomes)]))}
  binchrend <- c()
  currentbin <- 0
  binchrmdl <- c()
  for (i in seq_along(rlechr$values)) {
    currentmiddle <- currentbin+rlechr$lengths[i]/2
    currentbin <- currentbin+rlechr$lengths[i]
    binchrend <- append(binchrend, currentbin)
    binchrmdl <- append(binchrmdl, currentmiddle)
  } 
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
  cappedcopynumbers <- dfna[dfna$copynumbers > cap,]
  if(length(cappedcopynumbers$copynumbers)>0) {cappedcopynumbers$copynumbers <- cap-0.1}
  cappedsegments <- dfna[dfna$segments > cap,]
  if(length(cappedsegments$segments)>0) {cappedsegments$segments <- cap-0.1}
  toppedcopynumbers <- dfna[dfna$copynumbers <= bottom,]
  if(length(toppedcopynumbers$copynumbers)>0) {toppedcopynumbers$copynumbers <- bottom+0.1}
  toppedsegments <- dfna[dfna$segments <= bottom,]
  if(length(toppedsegments$segments)>0) {toppedsegments$segments <- bottom+0.1}
  
  
  
  line1 <- paste0("Cellularity: ", cellularity)
  line2 <- paste0("Ploidy: ", ploidy)
  if(!missing(error)) {
    line2 <- paste0("Relative error: ", round(error, digits = 3))
  }
  if(missing(chrsubset)){
    ggplot2::ggplot() +
      scale_y_continuous(name = "copies", limits = c(bottom,cap), breaks = seq(bottom,cap), expand=c(0,0)) +
      scale_x_continuous(name = "chromosome", limits = c(0,tail(binchrend,1)), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
      geom_hline(yintercept = seq(0,4), color = '#333333', linewidth = 0.5) +
      geom_hline(yintercept = seq(5,(cap-1)), color = 'lightgray', linewidth = 0.5) +
      geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
      geom_point(aes(x = bin,y = copynumbers),data=dfna[dfna$copynumbers>bottom&dfna$copynumbers<cap,], size = 0.1, color = 'black') +
      geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
      geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
      geom_point(aes(x = bin,y = segments),data=dfna[dfna$segments>bottom&dfna$segments<cap,], size = 1, color = 'darkorange') +
      geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'darkorange', shape = 24) +
      geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkorange', shape = 25) +
      theme_classic() + theme(
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-0.9, ymax = cap), fill = 'white') +
      annotate("text", x = binchrmdl[2], y = cap-0.3, label = line1) +
      annotate("text", x = binchrmdl[2], y = cap-0.7, label = line2)
  } else {
    
    if(firstchr==1){firstbin<-0
    } else {firstbin<-binchrend[firstchr-1]+1}
    ggplot2::ggplot() +
      scale_y_continuous(name = "copies", limits = c(bottom,cap), breaks = seq(bottom,cap), expand=c(0,0)) +
      scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), 
                         breaks = binchrmdl[seq(firstchr,lastchr)], 
                         labels = rlechr$values[seq(firstchr,lastchr)], expand = c(0,0)) +
      geom_hline(yintercept = seq(0,4), color = '#333333', linewidth = 0.5) +
      geom_hline(yintercept = seq(5,(cap-1)), color = 'lightgray', linewidth = 0.5) +
      geom_vline(xintercept = binchrend[seq(firstchr,lastchr)], color = "#666666", linetype = "dashed") +
      geom_point(aes(x = bin,y = copynumbers),data=dfna[dfna$copynumbers>bottom&dfna$copynumbers<cap,], size = 0.1, color = 'black') +
      geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
      geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
      geom_point(aes(x = bin,y = segments),data=dfna[dfna$segments>bottom&dfna$segments<cap,], size = 1, color = 'darkorange') +
      geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'darkorange', shape = 24) +
      geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkorange', shape = 25) +
      theme_classic() + theme(
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
}


# similarly to singleplot, this function takes a model and its template dataframe
# it returns a new dataframe with information about the individual segments: genomic location (chromosome, start, end),
# the segment value and the corresponding standard error and the p-value (in log10) associated with the nearest absolute copy number
getadjustedsegments <- function(template, QDNAseqobjectsample = FALSE, cellularity = 1, 
                                ploidy = 2, sgc = c(), standard, log=FALSE) {
  if(QDNAseqobjectsample) {template <- objectsampletotemplate(template, QDNAseqobjectsample)}
  template.na <- na.exclude(template)
  gc <- rep(2, nrow(template.na))
  gc[template.na$chr %in% sgc] <- 1
  if(log!=FALSE) {
    segmentdata <- rle(as.vector(paste0(template.na$chr, "_", template.na$segments)))
    segmentdata$values <- as.numeric(gsub(".*_", "", segmentdata$values))
  }
  # if(missing(standard) || !is(standard, "numeric")) { standard <- median(rep(segmentdata$values,segmentdata$lengths)) }
  if(missing(standard) || !is(standard, "numeric")) { standard <- median(template.na$segments) }
  # adjustedcopynumbers <- ploidy + ((template.na$copynumbers-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
  # adjustedsegments <- ploidy + ((template.na$segments-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
  adjustedcopynumbers <- template.na$copynumbers*(ploidy+2/cellularity-2)/standard - gc/cellularity + gc
  adjustedsegments <- template.na$segments*(ploidy+2/cellularity-2)/standard - gc/cellularity + gc
  adjsegmentdata <- rle(as.vector(paste0(template.na$chr, "_", adjustedsegments)))
  adjsegmentdata$values <- as.numeric(gsub(".*_", "", adjsegmentdata$values))
  adjcopynumberdata <- as.vector(adjustedcopynumbers)
  Chromosome <- c()
  Start <- c()
  End <- c()
  Num_Probes <- c()
  Segment_Mean <- c()
  Segment_Mean2 <- c()
  Segment_SE <- c()
  Copies <- c()
  P_log10 <- c()
  counter <- 1
  for (i in seq_along(adjsegmentdata$lengths)) {
    Chromosome[i] <- as.vector(template.na$chr)[counter]
    Start[i] <- as.vector(template.na$start)[counter]
    End[i] <- as.vector(template.na$end)[counter+adjsegmentdata$lengths[i]-1]
    Num_Probes[i] <- adjsegmentdata$lengths[i]
    if(log!=FALSE) {
      if(log==TRUE) {log <- 2}
      Segment_Mean[i] <- log(segmentdata$values[i]/standard, log)}
    if(log==FALSE) {
      Segment_Mean[i] <- adjsegmentdata$values[i]
      Segment_Mean2[i] <- mean(adjcopynumberdata[seq(counter, counter+adjsegmentdata$lengths[i]-1)])
      Segment_SE[i] <- sd(adjcopynumberdata[seq(counter, counter+adjsegmentdata$lengths[i]-1)])/sqrt(Num_Probes[i])
      Copies[i] <- round(Segment_Mean2[i])
      P_log10[i] <- round(log(2*(1-pnorm(abs(Copies[i]-Segment_Mean2[i])/Segment_SE[i])),10),digits=3)
    }
    counter <- counter+adjsegmentdata$lengths[i]
  }
  if(log!=FALSE) {df <- data.frame(Chromosome,Start,End,Num_Probes,Segment_Mean)}
  if(log==FALSE) {df <- data.frame(Chromosome,Start,End,Num_Bins=Num_Probes,Segment_Mean,Segment_Mean2,Segment_SE,Copies,P_log10)}
  return(df)
}


# If you have mutation data, you can use this function to find the number of
# copies of the corresponding genomic location. It expects a data frame with
# columns specifying Chromosome, Position, and Frequency (variant allele
# frequency in percentage). It can also be a file path to a tab-delimited
# file.linkvariants uses a data frame with segment information, such as created
# by getadjustedsegments. Again, this can be either a data frame or a file path
# of a tab-delimited text file. The function is also able to calculate variant
# copies in case it concerns heterozygous germline variants. Finally, the
# function can provide a confidence interval, but it needs an indication of read
# depth for this. When available, it will use the binom.test function to
# calculate upper and lower bounds of the chosen confidence level.

linkvariants <- function(variantdf, segmentdf, cellularity = 1, hetSNPs = FALSE,
                         chrindex=1,posindex=2,freqindex,altreadsindex,
                         totalreadsindex,refreadsindex,confidencelevel=FALSE,
                         append=TRUE, outputdir, sgc = c()){
  if(is(variantdf, "character")) {
    filename <- variantdf
    variantdf <- try(read.table(filename, header = TRUE, comment.char = "", sep = "\t"))
    writefiles <- TRUE
  } else {writefiles <- FALSE}
  if(is(segmentdf, "character")) {segmentdf <- try(read.table(segmentdf, header = TRUE, comment.char = "", sep = "\t"))}
  if(!inherits(variantdf, "try-error")) {
    Chromosome <- as.vector(variantdf[,chrindex])
    Chromosome <- gsub("chr","",Chromosome,ignore.case = TRUE)
    Position <- as.vector(variantdf[,posindex])
    if(!missing(altreadsindex)){
      altreads <- round(as.vector(variantdf[,altreadsindex]))
      if(!missing(freqindex)) {
        Frequency <- as.numeric(gsub("%","",as.vector(variantdf[,freqindex])))
        if(!missing(totalreadsindex)) {
          totalreads <- round(as.vector(variantdf[,totalreadsindex]))
        } else {
          totalreads <- round(altreads*100/Frequency)
        }
        
      } else if(!missing(totalreadsindex)) {
        totalreads <- round(as.vector(variantdf[,totalreadsindex]))
        Frequency <- 100*altreads/totalreads
      } else if(!missing(refreadsindex)){
        refreads <- round(as.vector(variantdf[,refreadsindex]))
        totalreads <- round(refreads+altreads)
        Frequency <- 100*altreads/totalreads
      } else {
        stop("Insufficient input data: supply variant frequency, total reads, or reference reads")
      }
    } else if(!missing(freqindex)&&missing(totalreadsindex)) {
      if (confidencelevel != FALSE) {
        warning("Not enough information for estimating confidence intervals")
      }
      confidencelevel <- FALSE
      Frequency <- as.numeric(gsub("%","",as.vector(variantdf[,freqindex])))
    } else if(!missing(freqindex)&&!missing(totalreadsindex)) {
      Frequency <- as.numeric(gsub("%","",as.vector(variantdf[,freqindex])))
      totalreads <- round(as.vector(variantdf[,totalreadsindex]))
      altreads <- round((Frequency/100)*totalreads)
    } else {
      stop("This function requires SOME kind of input data")
    }
    
    Copynumbers <- c()
    Variant_copies <- c()
    Variant_copies_low <- c()
    Variant_copies_high <- c()
    if (confidencelevel != FALSE) {
      ci <- apply(cbind(altreads,totalreads), 1, function(x) {
        binom.test(x[1], x[2], conf.level = confidencelevel)$conf.int })
      Frequency_low <- 100*ci[1,]
      Frequency_high <- 100*ci[2,]
    }
    
    if(length(Chromosome!=0)) {
      for (s in seq_along(Chromosome)) {
        cn <- segmentdf$Segment_Mean2[Chromosome[s] == segmentdf$Chromosome & Position[s] >= segmentdf$Start & Position[s] <= segmentdf$End]
        if(length(cn)==1){
          Copynumbers[s] <- cn
          if (hetSNPs == TRUE) {
            Variant_copies[s] <- (Frequency[s]/100)*(cn-2+2/cellularity) - (1/cellularity - 1)
            if (confidencelevel != FALSE) {
              Variant_copies_low[s] <- (Frequency_low[s]/100)*(cn-2+2/cellularity) - (1/cellularity - 1)
              Variant_copies_high[s] <- (Frequency_high[s]/100)*(cn-2+2/cellularity) - (1/cellularity - 1)
            }
          } else {
            if (Chromosome[s] %in% sgc) {
              Variant_copies[s] <- (Frequency[s]/100)*(cn-1+1/cellularity)
              if (confidencelevel != FALSE) {
                Variant_copies_low[s] <- (Frequency_low[s]/100)*(cn-1+1/cellularity)
                Variant_copies_high[s] <- (Frequency_high[s]/100)*(cn-1+1/cellularity)
              }
            } else {
              Variant_copies[s] <- (Frequency[s]/100)*(cn-2+2/cellularity)
              if (confidencelevel != FALSE) {
                Variant_copies_low[s] <- (Frequency_low[s]/100)*(cn-2+2/cellularity)
                Variant_copies_high[s] <- (Frequency_high[s]/100)*(cn-2+2/cellularity)
              }
            }
          }
        } else {
          Copynumbers[s] <- NA
          Variant_copies[s] <- NA
          if (confidencelevel != FALSE) {
            Variant_copies_low[s] <- NA
            Variant_copies_high[s] <- NA
          }
        }
      }
    }
    
    if(append==TRUE){
      output <- variantdf
      output$Copynumbers <- Copynumbers
      if (confidencelevel != FALSE) {
        output$Frequency_low <- Frequency_low
        output$Frequency_high <- Frequency_high
      }
      output$Variant_copies <- Variant_copies
      if (confidencelevel != FALSE) {
        output$Variant_copies_low <- Variant_copies_low
        output$Variant_copies_high <- Variant_copies_high
      }
    } else if (confidencelevel == FALSE){output <- data.frame(Chromosome,Position,Frequency,
                                                              Copynumbers,Variant_copies)
    } else {output <- data.frame(Chromosome,Position,Frequency,Frequency_low,Frequency_high,
                                 Copynumbers,Variant_copies,Variant_copies_low,Variant_copies_high)}
    if(writefiles) {
      if (missing(outputdir)) {
        fn <- gsub(".csv","_ACE.csv",filename)
        fn <- gsub(".tsv","_ACE.tsv",fn)
        fn <- gsub(".txt","_ACE.txt",fn)
        fn <- gsub(".xls","_ACE.xls",fn)
        write.table(output, file = fn ,sep="\t",row.names = FALSE, quote=FALSE)
      } else {
        if (dir.exists(outputdir)==FALSE) {dir.create(outputdir)}
        if (dirname(filename)==".") {newpath <- file.path(outputdir,filename)}
        else {newpath <- sub(dirname(filename),outputdir,filename)}
        fn <- gsub(".csv","_ACE.csv",newpath)
        fn <- gsub(".tsv","_ACE.tsv",fn)
        fn <- gsub(".txt","_ACE.txt",fn)
        fn <- gsub(".xls","_ACE.xls",fn)
        write.table(output, file = fn ,sep="\t",row.names = FALSE, quote=FALSE)
      }
    } 
    return(output)
  } else {print(paste0("failed to link variant data for ", filename))}
}

# YES! This monstruous function had to be updated due to the new and improved
# linkvariants function, a more-is-better version of linkmutationdata that also
# provides confidence intervals and the option to analyze heterozygous germline
# variants. The changes to postanalysisloop are quite subtle, but a lot of new
# arguments were necessary to accomodate linkvariants. These arguments are
# generally only used in the linkvariants call, so consult linkvariants
# documentation for those arguments! Note that the default for confidencelevel
# is FALSE in this function, so you will not get the confidence intervals, even
# if you provide read counts, unless you specify the confidence level (e.g. 0.95
# for 95% CI). One other notable change is that I have changed the argument
# mutationdata to variantdata, and it now thinks it should look for the
# directory "variantdata" instead. Of course, as before, you can specify the
# name of the directory! The implementation of single germline copy chromosomes
# is a bit different compared to the other functions. In this case, the
# modelsfile (or models dataframe) needs to have an additional column for
# gender. The index of this column is to be indicated using the genderci option.
# If the gender starts with the letter m, X and Y are assumed single copy. The
# option to plot sex chromosomes has also been added using the onlyautosomes
# argument (see singleplot).

postanalysisloop <- function(copyNumbersSegmented,modelsfile,variantdata,prefix="",postfix="",trncname=FALSE,inputdir=FALSE,
                             hetSNPs=FALSE,chrindex=1,posindex=2,freqindex,altreadsindex,totalreadsindex,refreadsindex,
                             confidencelevel=FALSE,append=TRUE,dontmatchnames=FALSE,printsegmentfiles=TRUE,printnewplots=TRUE,
                             imagetype='pdf',onlyautosomes=TRUE,outputdir="./",log=FALSE, segext='tsv', genderci) {
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  if(inputdir!=FALSE){
    if(missing(copyNumbersSegmented)) {
      files <- list.files(inputdir, pattern = "\\.rds$")
      if(length(files)==0){print("no rds-file detected")}
      if(length(files)>1){print(paste0("multiple rds-files detected, using first file: ",files[1]))}
      copyNumbersSegmented <- readRDS(file.path(inputdir,files[1]))
    }
    if(missing(modelsfile)){models <- try(read.table(file.path(inputdir,"models.tsv"), header = TRUE, comment.char = "", sep = "\t"))
    } else {models <- read.table(modelsfile, header = TRUE, comment.char = "", sep = "\t")}
    if(inherits(models, "try-error")) {print("failed to read modelsfile")}
    if(missing(variantdata)) {
      if (dir.exists(file.path(inputdir,"variantdata"))) {
        variantdata <- file.path(inputdir,"variantdata")
      } else { variantdata <- inputdir}
    }
  }
  if(missing(copyNumbersSegmented)){print("this function requires a QDNAseq-object")}
  if(is(copyNumbersSegmented, "character")) {copyNumbersSegmented <- try(readRDS(copyNumbersSegmented))}
  if(inherits(copyNumbersSegmented, "try-error")) {print("failed to read RDS-file")}
  if(missing(modelsfile)&&inputdir==FALSE){print("this function requires a tab-delimited file with model information per sample")}
  if(!missing(modelsfile)&&is(modelsfile, "character")) {models <- try(read.table(modelsfile, header = TRUE, comment.char = "", sep = "\t"))
  } else {models <- modelsfile}
  if(inherits(models, "try-error")) {print("failed to read modelsfile")}
  if(missing(variantdata)){print("not linking variant data")}
  if(!missing(variantdata)) {
    if (length(list.files(variantdata, pattern = "\\.csv$"))>0) {varext<-".csv"
    } else if (length(list.files(variantdata, pattern = "\\.txt$"))>0) {varext<-".txt"
    } else if (length(list.files(variantdata, pattern = "\\.tsv$"))>0) {varext<-".tsv"
    } else if (length(list.files(variantdata, pattern = "\\.xls$"))>0) {varext<-".xls"
    } else {print("file extension of variant files not supported: use .csv, .txt, .tsv, or .xls")}
  }
  fd <- Biobase::fData(copyNumbersSegmented)
  pd <- Biobase::pData(copyNumbersSegmented)
  if(trncname==TRUE) {pd$name <- gsub("_.*","",pd$name)}
  if(trncname!=TRUE&&trncname!=FALSE) {pd$name <- gsub(trncname,"",pd$name)}
  newplots <- vector(mode = 'list', length = (length(pd$name)))
  for (a in seq_along(pd$name)) {
    if(dontmatchnames==TRUE) { currentindex <- a
    } else { currentindex <- which(models[,1]==pd$name[a]) }
    if(length(currentindex)==1){
      cellularity <- as.numeric(models[currentindex,2])
      if(is.na(cellularity)){
        print(paste0("no valid cellularity given for ", pd$name[a]))
        newplots[[a]] <- ggplot2::ggplot() + 
          annotate("text", x=1, y=1, label = paste0("no plot available for ",pd$name[a]," - missing cellularity")) + 
          theme_classic() + 
          theme(line = element_blank(), axis.title =  element_blank(), axis.text = element_blank())
      } else {
        ploidy <- as.numeric(models[currentindex,3])
        if(is.na(ploidy)||length(ploidy)==0) {
          print(paste0("ploidy assumed 2 for ", pd$name[a]))
          ploidy <- 2
        }
        standard <- as.numeric(models[currentindex,4])
        if(is.na(standard)||length(standard)==0) {
          print(paste0("standard calculated from data for ", pd$name[a]))
          standard <- "standard"
        }
        gender <- "F"
        if (!missing(genderci)) {
          gender <- as.vector(models[currentindex,genderci])
        }
        if (onlyautosomes==TRUE) {
          newplots[[a]] <- singleplot(copyNumbersSegmented,cellularity=cellularity,ploidy=ploidy,standard=standard,QDNAseqobjectsample=a,title=pd$name[a])
        } else if (grepl("^m", gender, ignore.case = TRUE)) {
          newplots[[a]] <- singleplot(copyNumbersSegmented,cellularity=cellularity,ploidy=ploidy,standard=standard,QDNAseqobjectsample=a,title=pd$name[a],
                                      onlyautosomes=onlyautosomes, sgc = c("X", "Y"))
        } else {
          newplots[[a]] <- singleplot(copyNumbersSegmented,cellularity=cellularity,ploidy=ploidy,standard=standard,QDNAseqobjectsample=a,title=pd$name[a],
                                      onlyautosomes=onlyautosomes)
        }
        
        imagefunction <- get(imagetype)
        if(printnewplots==TRUE) {
          if (!dir.exists(file.path(outputdir,"newplots"))) {dir.create(file.path(outputdir,"newplots"))}
          if(imagetype %in% c('pdf','svg')) {
            imagefunction(file.path(outputdir,"newplots",paste0(pd$name[a],".",imagetype)),width=10.5)
            print(newplots[[a]])
            dev.off()
          } else {
            imagefunction(file.path(outputdir,"newplots",paste0(pd$name[a],".",imagetype)), width=720)
            print(newplots[[a]])
            dev.off()
          }
        }
        if (grepl("^m", gender, ignore.case = TRUE)) {
          segmentdf <- getadjustedsegments(copyNumbersSegmented,cellularity=cellularity,ploidy=ploidy,
                                           standard=standard,QDNAseqobjectsample=a,log=log,sgc=c("X","Y"))
        } else {
          segmentdf <- getadjustedsegments(copyNumbersSegmented,cellularity=cellularity,ploidy=ploidy,
                                           standard=standard,QDNAseqobjectsample=a,log=log)
        }
        if(!missing(variantdata)) {
          variantfile <- file.path(variantdata,paste0(prefix,pd$name[a],postfix,varext))
          folder <- file.path(outputdir,"variantdata")
          if(file.exists(variantfile)) {
            if (grepl("^m", gender, ignore.case = TRUE)) {
              linkvariants(variantfile,segmentdf,cellularity,hetSNPs,chrindex,posindex,freqindex,altreadsindex,
                           totalreadsindex,refreadsindex,confidencelevel,append,outputdir=folder,sgc=c("X","Y"))              
            } else {
              linkvariants(variantfile,segmentdf,cellularity,hetSNPs,chrindex,posindex,freqindex,altreadsindex,
                           totalreadsindex,refreadsindex,confidencelevel,append,outputdir=folder)
            }
          } else {warning(paste0("Cannot find file ", variantfile))}
        }
        if(printsegmentfiles==TRUE){
          if (!dir.exists(file.path(outputdir,"segmentfiles"))) {dir.create(file.path(outputdir,"segmentfiles"))}
          fn <- file.path(outputdir,"segmentfiles",paste0(pd$name[a],"_segments.",segext))
          write.table(segmentdf,fn,sep="\t",row.names = FALSE, quote=FALSE)
        }
      }
    } else {print(paste0("none or multiple models for ",pd$name[a]))}
  }
  return(newplots)
}

# Function to quickly look up adjusted copynumer and (optionally) mutant copies of specified genomic locations
# Functionality is very similar to linkmutationdata, but does not require "external" files
# Requires dataframe with adjusted segments
# Cellularity is only required to calculate mutant copies, use same number as used to generate segmentdf
# Chromosome, Position, and (optionally) Frequency can be single numbers or numeric vectors of the equal length
analyzegenomiclocations <- function(segmentdf, Chromosome, Position, Frequency, cellularity, sgc = c()){
  Copynumbers <- c()
  if(!missing(Frequency)){Mutant_copies <- c()}
  for (s in seq_along(Chromosome)) {
    cn <- segmentdf$Segment_Mean2[Chromosome[s] == segmentdf$Chromosome & Position[s] >= segmentdf$Start & Position[s] <= segmentdf$End]
    if(length(cn)==1){
      Copynumbers[s] <- cn
      if(!missing(Frequency)) {
        if (Chromosome[s] %in% sgc) {
          Mutant_copies[s] <- (Frequency[s]/100)*(cn-1+1/cellularity) 
        } else {
          Mutant_copies[s] <- (Frequency[s]/100)*(cn-2+2/cellularity) 
        }
      }
    } else {
      Copynumbers[s] <- NA
      if(!missing(Frequency)){Mutant_copies[s] <- NA}
    }
  }
  if(missing(Frequency)){return(data.frame(Chromosome, Position, Copynumbers))
  } else {return(data.frame(Chromosome, Position, Frequency, Copynumbers, Mutant_copies))}
}


# the code below for the multiplot function has been copy-pasted from www.cookbook-r.com which was made available under a CC0 license, courtesy of Winston Chang

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


