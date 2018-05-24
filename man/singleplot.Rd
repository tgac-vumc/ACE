\name{singleplot}
\alias{singleplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot an absolute copy number profile for a single sample
}
\description{
\code{singleplot} is the core plotting function of ACE. Input can be either a template or a QDNAseq-object with the index of the sample specified. Several of the arguments are parameters obtained from model fitting. Returns a \code{ggplot2} graph with absolute copies on the y-axis and genomic position on the x-axis. 
}
\usage{
singleplot(template, QDNAseqobjectsample = FALSE, 
  cellularity = 1, error, ploidy = 2, standard, title, 
  trncname = FALSE, cap = 12, bottom = 0, chrsubset, 
  onlyautosomes = TRUE)
}
\arguments{
  \item{template}{Object. Either a data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to analyze from the QDNAseq-object. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{cellularity}{Numeric. Used for rescaling bin and segment values. Printed on graph. Default = 1}
  \item{error}{Numeric. When given, it is printed on the graph below cellularity.}
  \item{ploidy}{Integer. Assume the median of segments has this absolute copy number. Default = 2}
  \item{standard}{Numeric. Force the given ploidy to represent this raw value. When omitted, the standard will be calculated from the data. When using parameters obtained from \code{squaremodel}, specify standard = 1}
  \item{title}{Character string. Overwrites the automatically generated title}
  \item{trncname}{Logical. In case of a QDNAseq object, the name of the sample is retrieved from the object and used as title. If set to TRUE, \code{trncname} truncates the sample name from the first instance of "_" in the name. You can also specify the regular expression here, e.g. trncname = "-.*" truncates the name from the first dash. Default = FALSE}
  \item{cap}{Integer. Influences your output copy number graphs. The upper limit of the y-axis is set at this number. When set to "max", it sets the cap to the maximum absolute copynumber value, rounded up. Bins and segments that exceed the cap are represented by a special mark. Recommended use between 8 and 16. Default = 12}
  \item{bottom}{Integer. Similar to cap, but for the lower limit of the y-axis. When set to "min", it sets the bottom to the minimum absolute copynumber value, rounded down. Bins and segments that subceed the bottom are represented by a special mark. Default = 0}
  \item{chrsubset}{Integer vector. Specify the chromosomes you want to plot. It will always take the full range of chromosomes in your subset, so specifying chrsubset = c(4, 8) will give the same plot as chrsubset = 4:8. When using a subset, singleplot will not plot the cellularity and error on the plot. Therefore, you can use this to make a copy number plot without this information by specifying chrsubset = 1:22}
  \item{onlyautosomes}{Logical or integer. You can fill in an integer to specify how many autosomes your species has. When TRUE, \code{singleplot} defaults to 22 (human) autosomes. When FALSE, \code{singleplot} will also plot whichever other chromosomes are specified in the template, e.g. "X", "Y", "MT". You can combine this argument with chrsubset, for instance chrsubset[1:23] to only include chromosome X (provided this is the 23rd chromosome)}
}

\value{
Returns a graph generated through the \code{ggplot2} package.
}

\author{
Jos B. Poell
}

\note{
  \code{singleplot} expects chromosome names, as specified in the \code{chr} column of the template, to be either just the integer chromosome number, or "chr" followed by the chromosome number. This is strictly required when onlyautosomes = TRUE.
}

\seealso{
  \code{\link{objectsampletotemplate}}, \code{\link{squaremodel}}, \code{\link{singlemodel}}
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
model <- singlemodel(template)
bestfit <- model$minima[model$rerror==min(model$rerror)]
singleplot(template, cellularity = tail(bestfit, 1), title = "sim")

## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
singlemodel(copyNumbersSegmented, QDNAseqobjectsample = 1)
singleplot(copyNumbersSegmented, QDNAseqobjectsample = 1, 
  cellularity = 0.79)
## QDNAseq 'blacklists' sex chromosomes, but singleplot can plot them
singleplot(copyNumbersSegmented, QDNAseqobjectsample = 1, 
  cellularity = 0.79, chrsubset = 12:24, onlyautosomes = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
