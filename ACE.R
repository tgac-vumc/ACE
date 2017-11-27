##### ACE has arrived! #####
# There are a number of functions: ACE, ploidyplotloop, ObjectsampleToTemplate, singlemodel, squaremodel, singleplot, 
# getadjustedsegments, linkmutationdata, postanalysisloop, analyzegenomiclocations.
# I have also added the "multiplot" function for making summary sheets.

# ACE takes a folder with rds-files (default, make sure all are segmented) or bam-files and returns plots for the most likely errors in convenient subfolders
# in case of bam-files it will run binsizes 30, 100, 500 and 1000 kbp as default, but a vector of desired binsizes can be used as input
# binsizes available are 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp
# note that output will be larger and programs will run longer with the smaller binsizes (1 kbp is pretty hardcore)
# ploidies specifies for which ploidies fits will be made: default is 2 but you can input a vector with all desired ploidies, e.g. c(2,4)
# beggers be choosers! ACE likes pdf-files, but for large datasets or small binsizes you might want to go with 'png'
# a standard method for error calculation is the root mean squared error (RMSE), but you can argue you should actually punish more the 
#   long segments that are a little bit off compared to the short segments that are way off
# 'MAE' weighs every error equally, whereas 'SMRE' does the latter; these will generally generate more fits (I see that as a downside)
#   but it is more likely that the segments of which you are sure are sticking tighter to the integer ploidy in the best fit
# the parameter "penalty" sets a penalty for the error calculated at lower cellularity
#   errors are divided by the cellularity to the power of the penalty: default = 0 (no penalty)
# new functionality for large data sets (many samples!):
# fits will always have to be chosen manually, but you can now open the file with all errorlists to pick the most likely fit without looking at the profile
# use the generated tab-delimited file to fill in the column with the likely fit
# added argument trncname (truncate name) which truncates the name to everything before the first instance of _ - set TRUE if applicable, 
#   or specify regular expression, e.g. "-.*" to REMOVE everything after and including the first dash found in a sample name
# added argument printsummaries: superbig summary files may crash the program, so you can set this argument to FALSE or 2 if you still want the error plots

# ploidyplotloop is the meat of ACE, and it can be run as a separate function as well; 
# it takes a QDNAseq-object as input and also needs a folder to write the files to
# this is particularly handy if you want to analyze a whole QDNAseq-object that is loaded in your R environment

# ObjectsampleToTemplate is pretty much there to parse QDNAseqobjects into the dataframe structure used by the singlemodel and singleplot functions
# these latter functions call ObjectsampleToTemplate itself when necessary, but it can be handy to make a template if you expect some repeated use of the functions
# or if you want to make the template and then do your own obscure manipulations to it (I know you want to!)

# singlemodel: like the name implies it runs the fitting algorithm on a single sample and spits out the info you want!
# not strictly necessary to save it to a variable, but that might still be handy
# singlemodel comes with the awesomeness of manual input: you can restrain the model to the ploidy you expect (default 2, but hey, ploidy 5 happens, right?)
# and you can tell the program if you think there is a different standard (this should only be necessary when the median segment happens to be a subclonal variant; very rare)

# squaremodel: this function drops the assumption of a standard, and instead runs the algorithm for each ploidy being the standard.
# if you feel you are doing too much tinkering with variables using the singlemodel function, then try this!
# you can choose the range of ploidies it should test; the default being all ploidies between 5 and 1 in 100 decrements of 0.04
# like singlemodel, the function returns a list, which you can save to a variable

# singleplot: again, you already know what it does! But did you know you can change the chart title to anything you like? Well, I'm here to tell you it can
# and there is more! Don't believe the cellularities the program calculated for you? You don't have to! Just fill in the cellularity you believe to be true, and singleplot will take your word for it!

# getadjustedsegments: get info of the actual segments in a handy dataframe! Contains start and end location, "number of probes" (number of bins supporting the segment value),
# value of the segment (Segment_Mean is the value from QDNAseq, Segment_Mean2 is calculated from adjustedcopynumbers within the segment),
# and the nearest ploidy with the chance (log10 p-value) that the segment has this ploidy. A very low p-value indicates a high chance of subclonality, 
# although this value should be approached with extreme caution (and is at the moment not corrected for multiple testing)

# linkmutationdata: now we're talking! Give a tab-delimited file with mutation data, and this function will tell you what the copy number is
# at that genomic location, and it will guess how many copies are mutant! Read more about the options for this function in the code

# postanalysisloop: you've run ACE, you've picked your models, you have your mutation data ... Now if you could only ...
# say no more! this function combines the power of all above functions and your own brain (you still choose the models)
# only for the brave of heart

# analyzegenomiclocations: a sort of simplified version of linkmutationdata, you can manually input a single genomic location as specified by chromosome and position
# you can also input multiple locations using vectors of the same length. The output will be a dataframe. You can also enter frequencies 
# to quickly calculate mutant copies, but then you also need to enter cellularity! Note: the function requires the dataframe with
# adjusted segments (output of getadjustedsegments)

# That's pretty much it for now, let me know if you run into some errors or oddities j.poell@vumc.nl

ACE <- function(inputdir = "./", outputdir, filetype = 'rds', binsizes, ploidies = 2, imagetype = 'pdf', method = 'RMSE', penalty = 0, cap = 12, trncname = FALSE, printsummaries = TRUE) { 
	imagefunction <- get(imagetype)
  library(QDNAseq)
	if(substr(inputdir,nchar(inputdir),nchar(inputdir)) != "/") {inputdir <- paste0(inputdir,"/")}
	if(missing(outputdir)) { outputdir <- substr(inputdir,0,nchar(inputdir)-1) }
	if(!dir.exists(outputdir)) {dir.create(outputdir)}
	if(filetype=='bam'){
		if(missing(binsizes)) { binsizes <- c(30,100,500,1000) }
	  parameters <- data.frame(options = c("inputdir","outputdir","filetype","binsizes","ploidies","imagetype","method","penalty","cap","trncname","printsummaries"), 
	                           values = c(inputdir,outputdir,filetype,paste0(binsizes,collapse=", "),paste0(ploidies,collapse=", "),imagetype,method,penalty,cap,trncname,printsummaries))
	  for (b in binsizes) {
		  currentdir <- file.path(outputdir,paste0(b,"kbp"))
		  dir.create(currentdir)
		  bins <- getBinAnnotations(binSize = b)
		  readCounts <- binReadCounts(bins, path = inputdir)
		  readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
		  readCountsFiltered <- estimateCorrection(readCountsFiltered)
		  # the default correctBins will output a ratio; using method = 'median' will return a corrected readcount
		  copyNumbers <- correctBins(readCountsFiltered)
		  copyNumbers <- normalizeBins(copyNumbers)
		  copyNumbers <- smoothOutlierBins(copyNumbers)
		  copyNumbersSegmented <- segmentBins(copyNumbers, transformFun = 'sqrt') # the transformFun is not available in older versions of QDNAseq!
		  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
		  saveRDS(copyNumbersSegmented, file = file.path(outputdir,paste0(b,"kbp.rds")))
		  
		  ploidyplotloop(copyNumbersSegmented,currentdir,ploidies,imagetype,method,penalty,cap,trncname,printsummaries)
		  
		}
	}
	else if(filetype=='rds'){
	  parameters <- data.frame(options = c("inputdir","outputdir","filetype","ploidies","imagetype","method","penalty","cap","trncname","printsummaries"), 
	                           values = c(inputdir,outputdir,filetype,paste0(ploidies,collapse=", "),imagetype,method,penalty,cap,trncname,printsummaries))
	  files <- list.files(inputdir, pattern = "\\.rds$")
	  for (f in 1:length(files)) {
			currentdir <- file.path(outputdir,paste0(substr(files[f],0,nchar(files[f])-4)))
			dir.create(currentdir)
			copyNumbersSegmented <- readRDS(file.path(inputdir,files[f]))
			ploidyplotloop(copyNumbersSegmented,currentdir,ploidies,imagetype,method,penalty,cap,trncname,printsummaries)
		}
	}
	else {
	print("not a valid filetype")
	}
	write.table(parameters, file=file.path(outputdir,"log.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
}

ploidyplotloop <- function(copyNumbersSegmented,currentdir,ploidies=2,imagetype='pdf',method='RMSE',penalty=0,cap=12,trncname=FALSE,printsummaries=TRUE) {
  imagefunction <- get(imagetype)
  library(ggplot2)
  library(Biobase)
	fd <- fData(copyNumbersSegmented)
	pd <- pData(copyNumbersSegmented)
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
  	
  	for (a in 1:length(pd$name)) {
  		segmentdata <- rle(as.vector(na.exclude(copyNumbersSegmented@assayData$segmented[,a])))
  		standard <- median(rep(segmentdata$values,segmentdata$lengths))
  			
  		fraction <- c()
  		expected <- c()
  		temp <- c()
  		errorlist <- c()
  		
  		for (i in 5:100) {
  		  fraction[i-4] <- i/100
  		  for (p in 1:12) {
  		    expected[p] <- standard*(1+(p-q)*fraction[i-4]/(fraction[i-4]*(q-2)+2))
  		  }
  		  # the maximum error 0.5 was added to make sure hyperamplifications (p>12) don't get ridiculous errors
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
  		obssd <- QDNAseq:::sdDiffTrim(copyNumbersSegmented@assayData$copynumber[,a], na.rm=T)
  		bin <- 1:length(fd$chromosome)
  		chr <- fd$chromosome
  		rlechr <- rle(chr)
  		binchrend <- c()
  		currentbin <- 0
  		binchrmdl <- c()
  		for (i in 1:length(rlechr$values)) {
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
  		
  		cellularity <- 5:100
  		tempdf <- data.frame(cellularity,errorlist=errorlist/max(errorlist))
  		tempplot <- ggplot() +
  		  scale_y_continuous(name = "relative error", limits = c(0,1.05), expand=c(0,0)) +
  		  scale_x_continuous(name = "cellularity (%)") +
  		  geom_vline(xintercept = seq(from = 10, to = 100, by = 10), color = "#666666", linetype = "dashed") +
  		  geom_point(aes(y=errorlist, x=cellularity), data=tempdf) +
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
  		
  		fitpicker[a,1] <- pd$name[a]
  		fitpicker[a,2] <- ""
  		fitpicker[a,3] <- q
  		fitpicker[a,4] <- standard
  		
  		bfi <- tail(which(rerror==min(rerror)))

  		for (m in 1:length(minima)) {
  		  fitpicker[a,m+4] <- minima[m]
  		  adjustedcopynumbers <- q + ((copyNumbersSegmented@assayData$copynumber[,a]-standard)*(minima[m]*(q-2)+2))/(minima[m]*standard)
  		  adjustedsegments <- q + ((copyNumbersSegmented@assayData$segmented[,a]-standard)*(minima[m]*(q-2)+2))/(minima[m]*standard)
  		  df <- as.data.frame(cbind(bin,adjustedcopynumbers,adjustedsegments))
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
  		  line1 <- paste0("Cellularity: ", minima[m])
  		  line2 <- paste0("Relative error: ", round(rerror[m], digits = 3))
  		  line3 <- paste0("Expected Noise: ", round(expsd,digits = 4))
  		  line4 <- paste0("Observed Noise: ", round(obssd,digits = 4))
  		  fn <- file.path(fp,"graphs",paste0(pd$name[a], " - ",q,"N fit ", m, ".",imagetype))
  		  
  		  tempplot <- ggplot() +
    			scale_y_continuous(name = "copies", limits = c(0,12), breaks = c(0:12), expand=c(0,0)) +
    			scale_x_continuous(name = "chromosome", limits = c(0,binchrend[22]), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
    			geom_hline(yintercept = c(0:4), color = '#333333', size = 0.5) +
    			geom_hline(yintercept = c(5:cap-1), color = 'lightgray', size = 0.5) +
    			geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
    			geom_point(aes(x = bin,y = copynumbers),data=df, size = 0.1, color = 'black') +
  		    geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
  		    geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
    			geom_point(aes(x = bin,y = segments),data=df, size = 1, color = 'darkorange') +
  		    geom_point(aes(x = bin,y = segments),data=cappedsegments, size = 1, color = 'darkorange', shape = 24) +
  		    geom_point(aes(x = bin,y = segments),data=toppedsegments, size = 1, color = 'darkorange', shape = 25) +
  		    theme_classic() + theme(
  		      axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
  		    ggtitle(paste0(pd$name[a], " - binsize ", binsize, " kbp - ", pd$used.reads[a], " reads - ",q,"N fit ", m)) +
    			geom_rect(aes(xmin=binchrmdl[1]/10, xmax = binchrmdl[4], ymin = cap-0.9, ymax = cap), fill = 'white') +
    			annotate("text", x = binchrmdl[2], y = cap-0.3, label = line1) +
    			annotate("text", x = binchrmdl[2], y = cap-0.7, label = line2) +
    			geom_rect(aes(xmin=binchrmdl[12], xmax = binchrmdl[22], ymin = cap-0.9, ymax = cap), fill = 'white') +
    			annotate("text", x = binchrmdl[16], y = cap-0.3, label = line3) +
    			annotate("text", x = binchrmdl[16], y = cap-0.7, label = line4)
  		  plots[[m]] <- tempplot
  		  if(bfi==m) {
  #		    bfplots[[(4*(a-1)+1)]] <- tempplot
  		    likelyplots[[(3*(a-1)+1)]] <- tempplot
  		    imagefunction(file.path(qdir,"likelyfits",paste0(pd$name[a],"_bestfit_",q,"N.",imagetype)))
  		    print(tempplot)
  		    dev.off()  
  		  }
  		  if(m==length(minima)) {
  #		    lmplots[[(4*(a-1)+1)]] <- tempplot
  		    likelyplots[[(3*(a-1)+2)]] <- tempplot
  		    imagefunction(file.path(qdir,"likelyfits",paste0(pd$name[a],"_lastminimum_",q,"N.",imagetype)))
  		    print(tempplot)
  		    dev.off() 
  		  }
  		  imagefunction(fn)
  		  # the print command is necessary when ggplotting in loops!
  		  print(tempplot)
  		  dev.off()
  		}
  			
  		if(imagetype == 'pdf') {
  		  imagefunction(file.path(fp,paste0("summary_",pd$name[a],".",imagetype)))
  		  print(plots)
  		  dev.off()
  		} else { 
  		  imagefunction(file.path(fp,paste0("summary_",pd$name[a],".",imagetype)), width = 1920, height = 480*ceiling(length(plots)/4))
  		  if (length(plots)==1) {print(plots)} else { print(multiplot(plotlist = plots, cols=4)) }
  		  dev.off()
  		}
  		  
  	}
	
	
  	
  	if(printsummaries == TRUE) {
    	if(imagetype == 'pdf') {
  #    	  imagefunction(file.path(qdir,paste0("summary_bestfits.",imagetype)))
  #    	  print(bfplots)
  #    	  dev.off()
  #  	  imagefunction(file.path(qdir,paste0("summary_lastminima.",imagetype)))
  #    	  print(lmplots)
  #    	  dev.off()
    	  imagefunction(file.path(qdir,paste0("summary_likelyfits.",imagetype)))
      	  print(likelyplots)
      	  dev.off()
      	imagefunction(file.path(qdir,paste0("summary_errors.",imagetype)))
      	  print(listerrorplots)
      	  dev.off()
      	} else { 
  #    	  imagefunction(file.path(qdir,paste0("summary_bestfits.",imagetype)), width = 1920, height = 480*length(pd$name))
  #    	  print(multiplot(plotlist = bfplots, cols=4))
  #    	  dev.off()
  #  	  imagefunction(file.path(qdir,paste0("summary_lastminima.",imagetype)), width = 1920, height = 480*length(pd$name))
  #    	  print(multiplot(plotlist = lmplots, cols=4))
  #    	  dev.off()
    	  imagefunction(file.path(qdir,paste0("summary_likelyfits.",imagetype)), width = 1440, height = 480*length(pd$name))
      	  print(multiplot(plotlist = likelyplots, cols=3))
      	  dev.off()
      	imagefunction(file.path(qdir,paste0("summary_errors.",imagetype)), width = 1920, height = 480*ceiling(length(pd$name)/4))
      	if (length(listerrorplots)==1){print(listerrorplots)} else {print(multiplot(plotlist = listerrorplots, cols=4))}
      	  dev.off()
      	}
  	} else if(printsummaries == 2) {
  	  if(imagetype == 'pdf') {
  	    imagefunction(file.path(qdir,paste0("summary_errors.",imagetype)))
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


# this function takes a template dataframe as specified in the above function
# you can use a QDNAseq object as "template", then specify QDNAseqobjectsample by the sample number!
# e.g. model <- singlemodel(copyNumbersSegmented, QDNAseqobjectsample = 3)
singlemodel <- function(template,ploidy = 2, standard, QDNAseqobjectsample = FALSE, method = 'RMSE', penalty = 0, highlightminima = TRUE) {
  library(ggplot2)
  if(QDNAseqobjectsample) {template <- ObjectsampleToTemplate(template, QDNAseqobjectsample)}
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


# Turns out, tumors are often quite messy, genomically speaking
# Occasionally, you might want to be able to compare the fits at 2N directly with the fits at other ploidies
# Also, with pooring quality or lots of subclonality, it is more common for the standard to be messed up
# In line with graphs presented by ASCAT and CELLULOID, squaremodel will do fitting and plot a graph for 
#   fits using two variables: ploidy and cellularity
# The use of a standard is no longer necessary: the fits are all using standard = 1
# This is relevant to keep in mind, because you have to specify standard = 1 if you use the fit in singleplot
# On top of the penalty for low cellularity, you can also add a penalty for ploidies
# It is implemented as follows: error*(1+abs(ploidy-2))^penploidy
# The plot has the two variables as axis and the color code indicates the relative error
# To make the minima pop out, the color code is the inverse of the relative error
# Minima are found by checking each value for neighboring values, and will only return true if its the lowest error
squaremodel <- function(template, QDNAseqobjectsample = FALSE, prows=100, ptop=5, pbottom=1, method = 'RMSE', penalty = 0, penploidy = 0, highlightminima = TRUE) {
  library(ggplot2)
  library(Biobase)
  if(QDNAseqobjectsample) {template <- ObjectsampleToTemplate(template, QDNAseqobjectsample)}
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
        expected[p] <- standard*(1+(p-ploidy)*fraction[i-4]/(fraction[i-4]*(ploidy-2)+2))
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

# this function takes a template dataframe as specified in ObjectsampleToTemplate
# you can use a QDNAseq object as "template", then specify QDNAseqobjectsample by the sample number!
# don't forget to specify cellularity, error, ploidy, and standard; you can find these in the model output
singleplot <- function(template,cellularity = 1, error, ploidy = 2, standard, title = "Plot",QDNAseqobjectsample = FALSE, cap = 12, chrsubset) {
  library(ggplot2)
  if(QDNAseqobjectsample) {template <- ObjectsampleToTemplate(template, QDNAseqobjectsample)}
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  if(missing(standard)) { standard <- median(rep(segmentdata$values,segmentdata$lengths)) }
  adjustedcopynumbers <- ploidy + ((template$copynumbers-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
	adjustedsegments <- ploidy + ((template$segments-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
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
  		geom_hline(yintercept = c(0:4), color = '#333333', size = 0.5) +
  		geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 0.5) +
  		geom_vline(xintercept = binchrend, color = "#666666", linetype = "dashed") +
  		geom_point(aes(x = bin,y = copynumbers),data=df, size = 0.1, color = 'black') +
  	  geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
  	  geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
  		geom_point(aes(x = bin,y = segments),data=df, size = 1, color = 'darkorange') +
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
	  firstchr <- range(chrsubset)[1]
	  lastchr <- range(chrsubset)[2]
	  if(firstchr==1){firstbin<-0
	  } else {firstbin<-binchrend[firstchr-1]+1}
	  ggplot() +
	    scale_y_continuous(name = "copies", limits = c(0,cap), breaks = c(0:cap), expand=c(0,0)) +
	    scale_x_continuous(name = "chromosome", limits = c(firstbin,binchrend[lastchr]), breaks = binchrmdl[firstchr:lastchr], labels = firstchr:lastchr, expand = c(0,0)) +
	    geom_hline(yintercept = c(0:4), color = '#333333', size = 0.5) +
	    geom_hline(yintercept = c(5:(cap-1)), color = 'lightgray', size = 0.5) +
	    geom_vline(xintercept = binchrend[firstchr:lastchr], color = "#666666", linetype = "dashed") +
	    geom_point(aes(x = bin,y = copynumbers),data=df, size = 0.1, color = 'black') +
	    geom_point(aes(x = bin,y = copynumbers),data=cappedcopynumbers, size = 0.5, color = 'black', shape = 24) +
	    geom_point(aes(x = bin,y = copynumbers),data=toppedcopynumbers, size = 0.5, color = 'black', shape = 25) +
	    geom_point(aes(x = bin,y = segments),data=df, size = 1, color = 'darkorange') +
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
getadjustedsegments <- function(template,cellularity = 1, ploidy = 2, standard, QDNAseqobjectsample = FALSE, log=FALSE) {
  if(QDNAseqobjectsample) {template <- ObjectsampleToTemplate(template, QDNAseqobjectsample)}
  segmentdata <- rle(as.vector(na.exclude(template$segments)))
  if(missing(standard)) { standard <- median(rep(segmentdata$values,segmentdata$lengths)) }
  adjustedcopynumbers <- ploidy + ((template$copynumbers-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
  adjustedsegments <- ploidy + ((template$segments-standard)*(cellularity*(ploidy-2)+2))/(cellularity*standard)
  template.na <- na.exclude(template)
  adjsegmentdata <- rle(as.vector(na.exclude(adjustedsegments)))
  adjcopynumberdata <- as.vector(na.exclude(adjustedcopynumbers))
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
  for (i in 1:length(adjsegmentdata$lengths)) {
    Chromosome[i] <- as.vector(template.na$chr)[counter]
    Start[i] <- as.vector(template.na$start)[counter]
    End[i] <- as.vector(template.na$end)[counter+adjsegmentdata$lengths[i]-1]
    Num_Probes[i] <- adjsegmentdata$lengths[i]
    if(log==TRUE) {Segment_Mean[i] <- log2(adjsegmentdata$values[i])}
    if(log==FALSE) {
      Segment_Mean[i] <- adjsegmentdata$values[i]
      Segment_Mean2[i] <- mean(adjcopynumberdata[counter:(counter+adjsegmentdata$lengths[i]-1)])
      Segment_SE[i] <- sd(adjcopynumberdata[counter:(counter+adjsegmentdata$lengths[i]-1)])/sqrt(Num_Probes[i])
      Copies[i] <- round(Segment_Mean2[i])
      P_log10[i] <- round(log(2*(1-pnorm(abs(Copies[i]-Segment_Mean2[i])/Segment_SE[i])),10),digits=3)
    }
    counter <- counter+adjsegmentdata$lengths[i]
  }
  if(log==TRUE) {df <- data.frame(Chromosome,Start,End,Num_Probes,Segment_Mean)}
  if(log==FALSE) {df <- data.frame(Chromosome,Start,End,Num_Bins=Num_Probes,Segment_Mean,Segment_Mean2,Segment_SE,Copies,P_log10)}
  return(df)
}


# If you have mutation data, you can use this function to find the ploidy of the corresponding genomic location
# Currently it expects a tab-delimited file with Chromosome, Position, and Frequency in columns 1, 2, and 3 respectively
# You can use the arguments to specify which columns if this is more convenient
# Frequency is expected to be a percentage! You can multiply Mutant_copies with 100 if it is given as a fraction
# The file should end with .txt, .csv or .tsv and not have any other dots beside the one before the extension!
# It expects the file to have a header row, which may be commented with a #
# You also need your segmentvalues using the getadjustedsegments function, either as a dataframe or file
# You can append the new columns to the original file (it will still save it under a new name) 
# or only output the used columns and new columns (filename is the same)
# In case of "append": note that R messes up column names if they have "special" characters or start with numbers
# outputdir can be specified, so you don't have to save in the same directory (especially handy in loops)
linkmutationdata <- function(filename, segmentdf, cellularity = 1,chrindex=1,posindex=2,freqindex=3, append=TRUE, outputdir){
  inputvcf <- try(read.table(filename, header = TRUE, comment.char = "", sep = "\t"))
  if(class(segmentdf)=="character") {segmentdf <- try(read.table(segmentdf, header = TRUE, comment.char = "", sep = "\t"))}
  if (!inherits(inputvcf, "try-error")) {
    Chromosome <- as.vector(inputvcf[,chrindex])
    Chromosome <- gsub("chr","",Chromosome,ignore.case = TRUE)
    Position <- as.vector(inputvcf[,posindex])
    Frequency <- as.numeric(gsub("%","",as.vector(inputvcf[,freqindex])))
    Copynumbers <- c()
    Mutant_copies <- c()
    if(length(Chromosome!=0)) {
      for (s in 1:length(Chromosome)) {
        cn <- segmentdf$Segment_Mean2[which(Chromosome[s] == segmentdf$Chromosome & Position[s] >= segmentdf$Start & Position[s] <= segmentdf$End)]
        if(length(cn)==1){
          Copynumbers[s] <- cn
          Mutant_copies[s] <- (Frequency[s]/100)*(cn-2+2/cellularity)
        } else {
          Copynumbers[s] <- NA
          Mutant_copies[s] <- NA
        }
      }
    }
    #Mutant_copies <- (Frequency/100)*(Copynumbers-2)+0.02*Frequency/cellularity
    if(append==TRUE){
      output <- inputvcf
      output$Copynumbers <- Copynumbers
      output$Mutant_copies <- Mutant_copies
    } else {output <- data.frame(Chromosome,Position,Frequency,Copynumbers,Mutant_copies)}
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
  } else {print(paste0("failed to link mutation data for ", filename))}
}

# This is the most involving function from the user perspective: it process estimated segment copynumbers and links mutationdata for an entire QDNAseq-object
# Additionally, you can print the plots of all selected models: the function will save the plots in a list and return the list
# The function needs models for all individual samples in the object you want analyzed
# This should be supplied by the user as a tab-delimited text with sample name, cellularity, ploidy, and standard
# If ploidy and standard are not given, they are assumed to be 2 and 1 respectively. Makes sure this holds true!
# Linking mutation data is optional if you only want to output the segment files or get new plots
# You can also set printsegmentfiles and printnewplots to FALSE to reduce output
# The argument mutationdata expects a directory path, e.g. "./" or "/vcf/mutationdata/"
# The directory should contain files that have the sample name in them, and end with .csv, .txt, .tsv or .xls (but all the same!)
# If the files have a prefix or postfix (such as mutations_sample.txt or sample_SNVs.csv), you specify this with the respective arguments
# copyNumbersSegmented can either be a QDNAseq-object already loaded into R or a file path to an RDS-file of the object
# trncname: get rid of everything from and including the underscore in the samplename; make sure it still matches the sample names in the mutation files!
# dontmatchnames: if you set this to TRUE, the function will just use the sample number in the QDNAseq-object as the row number in the modelsfile. Risky!
# inputdir is a bit of a convenience function. You can give a file path, and the function will look for the appropriate files there. 
#   It will still check for the first arguments, which will take priority over anything it finds in the directory. 
#   If not specified, it will look for a subdirectory mutationdata, if it doesn't find it, it will assume the files are in the inputdir
#   If it doesn't find the mutationdata, it should give the appropriate error
#   The models need to be in "models.txt" or otherwise specified
#   If not specified, it will open the first rds-file it finds in the inputdir: it will warn if it finds more than one rds-file, but still use the first!
# outputdir is there to make sure all your output ends in this directory, it will give mutationdata and segmentfiles in separate subdirectories
#   Specify the entire path for outputdir! If you don't start with /, it will make folders in your active working directory, not inputdir!

postanalysisloop <- function(copyNumbersSegmented,modelsfile,mutationdata,prefix="",postfix="",trncname=FALSE,inputdir=FALSE,
                             chrindex=1,posindex=2,freqindex=3,append=TRUE,dontmatchnames=FALSE,
                             printsegmentfiles=TRUE,printnewplots=TRUE,imagetype='pdf',outputdir="./"){
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  if(inputdir!=FALSE){
    if(missing(copyNumbersSegmented)) {
      files <- list.files(inputdir, pattern = "\\.rds$")
      if(length(files)==0){print("no rds-file detected")}
      if(length(files)>1){print(paste0("multiple rds-files detected, using first file: ",files[1]))}
      copyNumbersSegmented <- readRDS(file.path(inputdir,files[1]))
    }
    if(missing(modelsfile)){models <- try(read.table(file.path(inputdir,"models.txt"), header = TRUE, comment.char = "", sep = "\t"))
    } else {models <- read.table(modelsfile, header = TRUE, comment.char = "", sep = "\t")}
    if(inherits(models, "try-error")) {print("failed to read modelsfile")}
    if(missing(mutationdata)) {
      if (dir.exists(file.path(inputdir,"mutationdata"))) {
        mutationdata <- file.path(inputdir,"mutationdata")
      } else { mutationdata <- inputdir}
    }
  }
  if(missing(copyNumbersSegmented)){print("this function requires a QDNAseq-object")}
  if(class(copyNumbersSegmented)=="character") {copyNumbersSegmented <- try(readRDS(copyNumbersSegmented))}
  if(inherits(copyNumbersSegmented, "try-error")) {print("failed to read RDS-file")}
  if(missing(modelsfile)&&inputdir==FALSE){print("this function requires a tab-delimited file with model information per sample")}
  if(!missing(modelsfile)&&class(modelsfile)=="character") {models <- try(read.table(modelsfile, header = TRUE, comment.char = "", sep = "\t"))}
  if(inherits(models, "try-error")) {print("failed to read modelsfile")}
  if(missing(mutationdata)){print("not linking mutation data")}
  if(!missing(mutationdata)) {
    if (length(list.files(mutationdata, pattern = "\\.csv$"))>0) {mutext<-".csv"
  } else if (length(list.files(mutationdata, pattern = "\\.txt$"))>0) {mutext<-".txt"
  } else if (length(list.files(mutationdata, pattern = "\\.tsv$"))>0) {mutext<-".tsv"
  } else if (length(list.files(mutationdata, pattern = "\\.xls$"))>0) {mutext<-".xls"
  } else {print("file extension of mutation files not supported: use .csv, .txt, .tsv, or .xls")}
  }
  library(Biobase)
  fd <- fData(copyNumbersSegmented)
  pd <- pData(copyNumbersSegmented)
  if(trncname==TRUE) {pd$name <- gsub("_.*","",pd$name)}
  if(trncname!=TRUE&&trncname!=FALSE) {pd$name <- gsub(trncname,"",pd$name)}
  newplots <- vector(mode = 'list', length = (length(pd$name)))
  for (a in 1:length(pd$name)) {
    if(dontmatchnames==TRUE) { currentindex <- a
    } else { currentindex <- which(models[,1]==pd$name[a]) }
    if(length(currentindex)==1){
      cellularity <- as.numeric(models[currentindex,2])
      if(is.na(cellularity)){
        print(paste0("no valid cellularity given for ", pd$name[a]))
        newplots[[a]] <- ggplot() + 
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
          print(paste0("standard assumed 1 for ", pd$name[a]))
          standard <- 1
        }
        newplots[[a]] <- singleplot(copyNumbersSegmented,cellularity=cellularity,ploidy=ploidy,standard=standard,QDNAseqobjectsample=a,title=pd$name[a])
        imagefunction <- get(imagetype)
        if(printnewplots==TRUE) {
          if (!dir.exists(file.path(outputdir,"newplots"))) {dir.create(file.path(outputdir,"newplots"))}
          if(imagetype=='pdf'){
            imagefunction(file.path(outputdir,"newplots",paste0(pd$name[a],".",imagetype)),width=10.5)
            print(newplots[[a]])
            dev.off()
          } else {
            imagefunction(file.path(outputdir,"newplots",paste0(pd$name[a],".",imagetype)), width=720)
            print(newplots[[a]])
            dev.off()
          }
        }
        segmentdf <- getadjustedsegments(copyNumbersSegmented,cellularity=cellularity,ploidy=ploidy,standard=standard,QDNAseqobjectsample=a)
        if(!missing(mutationdata)) {
          mutationfile <- file.path(mutationdata,paste0(prefix,pd$name[a],postfix,mutext))
          folder <- file.path(outputdir,"mutationdata")
          linkmutationdata(mutationfile,segmentdf,cellularity,chrindex,posindex,freqindex,append,outputdir=folder)
        }
        if(printsegmentfiles==TRUE){
          if (!dir.exists(file.path(outputdir,"segmentfiles"))) {dir.create(file.path(outputdir,"segmentfiles"))}
          fn <- file.path(outputdir,"segmentfiles",paste0(pd$name[a],"_segments.tsv"))
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
analyzegenomiclocations <- function(segmentdf, Chromosome, Position, Frequency, cellularity){
  Copynumbers <- c()
  if(!missing(Frequency)){Mutant_copies <- c()}
  for (s in 1:length(Chromosome)) {
    cn <- segmentdf$Segment_Mean2[which(Chromosome[s] == segmentdf$Chromosome & Position[s] >= segmentdf$Start & Position[s] <= segmentdf$End)]
    if(length(cn)==1){
      Copynumbers[s] <- cn
      if(!missing(Frequency)) {Mutant_copies[s] <- (Frequency[s]/100)*(cn-2+2/cellularity) }
    } else {
      Copynumbers[s] <- NA
      if(!missing(Frequency)){Mutant_copies[s] <- NA}
    }
  }
  if(missing(Frequency)){return(data.frame(Chromosome, Position, Copynumbers))
  } else {return(data.frame(Chromosome, Position, Frequency, Copynumbers, Mutant_copies))}
}


# the code below for the multiplot functino has been copy-pasted www.cookbook-r.com which was made available under a CC0 license, courtesy of Winston Chang

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
  library(grid)
  
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
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
