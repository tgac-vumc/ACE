\name{squaremodelsummary}
\alias{squaremodelsummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a graphical summary of the result of squaremodel fitting
}
\description{
\code{squaremodelsummary} creates a graphical summary of a squaremodel by plotting the matrixplot and absolute copy number profiles corresponding with the 7 best fits. The list of plots can be saved to a variable and/or a file.
}
\usage{
squaremodelsummary(template, QDNAseqobjectsample = FALSE, 
  squaremodel, samplename, printplots = TRUE, outputdir, 
  imagetype = 'pdf', trncname = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{Object. Either a data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to analyze from the QDNAseqobject. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{squaremodel}{List of objects returned by \code{\link{squaremodel}}}
  \item{samplename}{Character string. Use this sample name in the title of the matrixplot. If the sample comes from a QDNAseq-object and this argument is omitted, the sample name is taken from the QDNAseq-object}
  \item{printplots}{Logical. Print the plots to file. Default = TRUE}
  \item{outputdir}{Character string. Print the plots to this directory}
  \item{imagetype}{Character string. Plots are printed to file using this graphics device. "pdf" will result in 8 pages with individual plots, while the other devices print a 2x4 mosaic. Default = "pdf"}
  \item{trncname}{Logical. If set to TRUE, \code{trncname} truncates the sample name from the first instance of "_" in the name. You can also specify the regular expression here, e.g. trncname = "-.*" truncates the name from the first dash. Default = FALSE}
}

\value{
Returns a list with eight plots (ggplot2-objects): the matrixplot and the copy number plots corresponding to the seven best fits. If printplots is TRUE, the plots will be printed to file. 
}

\author{
Jos B. Poell
}

\seealso{
\code{\link{squaremodel}}
}
\examples{
## simulated data assuming each chromosome comprises 100 bins
s <- jitter(c(1, 1, 0.8, 1.2, rep(1, 5), 1.4, rep(1, 13)), amount = 0)
n <- c(100, 100, 40, 60, rep(100, 5), 100, rep(100, 13))
bin <- 1:2200
chr <- rep(1:22, each = 100)
start <- rep(0:99*1000000+1, 22)
end <- rep(1:100*1000000, 22)
copynumbers <- jitter(rep(s,n), amount = 0.05)
segments <- rep(s, n)
template <- data.frame(bin = bin, chr = chr, start = start, end = end, 
  copynumbers = copynumbers, segments = segments)
sm <- squaremodel(template, method = 'MAE', penalty = 0.5, 
  penploidy = 0.5)
sms <- squaremodelsummary(template, squaremodel = sm, 
  samplename = "sim", printplots = FALSE)
sms[[1]]
sms[[2]]
## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
sqm <- squaremodel(copyNumbersSegmented, QDNAseqobjectsample = 2, 
  penalty = 0.5, penploidy = 0.5, 
  ptop = 4.3, pbottom = 1.8, prows = 250)
sqms <- squaremodelsummary(copyNumbersSegmented, 2, 
  squaremodel = sqm, printplots = FALSE)
sqms[[1]]
sqms[[2]] + ggplot2::ggtitle("Top fit for sample2")
}

