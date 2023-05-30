\name{twosamplecompare}
\alias{twosamplecompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Overlay copy number data of two samples and compare segment values
}
\description{
\code{twosamplecompare} "resegments" two samples to have the same breakpoints. Both samples' means of the resulting segments are tested for equality using the two-sided Welch two sample t-test. \code{twosamplecompare} returns a data frame with the comparisons per segment, it returns the correlation of segments, and a copy number plot with an overlay of (scaled) segment values of both samples and the associated -log10-transformed q-values. 
}
\usage{
twosamplecompare(template1, index1 = FALSE, ploidy1 = 2,
  cellularity1 = 1, standard1, name1, template2, index2 = FALSE, 
  ploidy2 = 2, cellularity2 = 1, standard2, name2,
  equalsegments = FALSE, altmethod = FALSE, cap = 12, qcap = 12, 
  bottom = 0, plot = TRUE, trncname = FALSE, legend = TRUE, 
  chrsubset, onlyautosomes = TRUE, sgc = c(),
  showcorrelation = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template1}{Object. Either a data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{index1}{Integer. If template1 is a QDNAseqobject, this specifies the index of the first sample. Default = FALSE}
  \item{ploidy1}{Integer. Assume the median of segments of the first sample has this absolute copy number. Default = 2}
  \item{cellularity1}{Numeric. Used for rescaling bin and segment values of the first sample. Default = 1}
  \item{standard1}{Numeric. Forces ploidy1 to represent this raw value. When omitted, the standard will be calculated from the data. When using parameters obtained from \code{squaremodel}, specify standard1 = 1}
  \item{name1}{Character string. Name of the first sample. Printed on graph}
  \item{template2}{Object. Either a data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object. When omitted, template1 will be used}
  \item{index2}{Integer. Specifies the index of the second sample in template2 or, when template2 is omitted, in template1. Default = FALSE}
  \item{ploidy2}{Integer. Assume the median of segments of the second sample has this absolute copy number. Default = 2}
  \item{cellularity2}{Numeric. Used for rescaling bin and segment values of the second sample. Default = 1}
  \item{standard2}{Numeric. Forces ploidy2 to represent this raw value. When omitted, the standard will be calculated from the data. When using parameters obtained from \code{squaremodel}, specify standard2 = 1}
  \item{name2}{Character string. Name of the second sample. Printed on graph}
  \item{equalsegments}{Logical or integer. If TRUE, \code{twosamplecompare} "resegments" both samples simply by creating segments containing roughly 20 bins, or as many bins as specified in this argument. When FALSE, both samples are resegmented by combining the break points and applying them to both samples. Default = FALSE}
  \item{altmethod}{Logical or character string. Instead of scaling the sample segments to absolute copies, scale them to standard units. There are two options: "SD" and "MAD". In the first case, the mean of segments is set to 0 and for each segment the distance (in standard deviations or "SD units" from the segment mean to the mean of segments is calculated in standard deviations. In case of "MAD", instead the median of segments, segment median, and median absolute deviation is used. Adjust the y-axis with the cap and bottom arguments for better visualization. Default = FALSE}
  \item{cap}{Integer. Influences your output copy number graphs. The upper limit of the y-axis is set at this number. When set to "max", it sets the cap to the maximum absolute copynumber value, rounded up. Bins and segments that exceed the cap are represented by a special mark. Recommended use between 8 and 16. Default = 12}
  \item{qcap}{Integer. Sets the upper limit of the secondary y-axis. Default = 12}
  \item{bottom}{Integer. Similar to cap, but for the lower limit of the y-axis. When set to "min", it sets the bottom to the minimum absolute copynumber value, rounded down. Bins and segments that subceed the bottom are represented by a special mark. Default = 0}
  \item{plot}{Logical. Produce a two-sample copy number plot. Default = TRUE}
  \item{trncname}{Logical. In case of a QDNAseq object, the name of the sample is retrieved from the object and used as title. If set to TRUE, \code{trncname} truncates the sample name from the first instance of "_" in the name. You can also specify the regular expression here, e.g. trncname = "-.*" truncates the name from the first dash. Default = FALSE}
  \item{legend}{Logical. Add the legend to the two-sample plot. Default = TRUE}
  \item{chrsubset}{Integer vector. Specify the chromosomes you want to plot. It will always take the full range of chromosomes in your subset, so specifying chrsubset = c(4, 8) will give the same plot as chrsubset = 4:8. When using a subset, \code{twosamplecompare} will not plot the cellularity and error on the plot.}
  \item{onlyautosomes}{Logical or integer. You can fill in an integer to specify how many autosomes your species has. When TRUE, \code{twosamplecompare} defaults to 22 (human) autosomes. When FALSE, \code{twosamplecompare} will also plot whichever other chromosomes are specified in the template, e.g. "X", "Y", "MT"}
  \item{sgc}{Integer or character vector. Specify which chromosomes occur with only a single copy in the germline. Note that this is assumed for both samples.}
  \item{showcorrelation}{Logical. Add the correlation to the two-sample plot. Default = TRUE}
  \item{quieterrors}{Logical. Mute warnings about errors that the statistical tests produce. Default = TRUE}
}
\details{
This function can be used for different types of comparisons. It can be used to compare a tumor sample with a healthy (preferably matched) control. In this case, it may not be necessary to fill in the cellularity, because it will not make a difference for the statistical tests. In this ability the function helps to determine which (if any) segments are significantly different from normal. The other major use is to compare two tumors from potentially the same origin, but that were separated in space or time. You can then assess if changes have occurred, or even whether the two samples are from different clonal origin. In this case it is important to achieve maximum similarity in segments. Now the argument altmethod may come in handy, because it does not require model fitting and optimization. The q-values that are obtained with this function should be interpreted with caution. The two-sample statistical tests will easily reach significance when the sample sizes, in this case bins per segment, are large. By creating equal segment sizes with the argument equalsegments, these biases disappear. 
}
\value{
\itemize{
\item twosampledf - data frame with the newly created segments and the information and comparison of both samples
\item correlation - Pearson correlation of the segment values of all bins between both samples
\item subsetcorrelation - same as correlation, but only applying to subset of chromosomes specified by the argument chrsubset
\item compareplot - ggplot2-graph of both samples with segment values in red (first sample) and blue (second sample). Green bars indicate q-values of the segments, scaled on the secondary axis
}
}

\author{
Jos B. Poell
}
\note{
The data frame, plot, and subsetcorrelation all use the same selection of chromosomes. The correlation in the plot corresponds to the displayed chromosomes. Note that the returned value \code{correlation} uses all segments in the data, also from the sex chromosomes when available. However, if there is no useful data for an entire chromosome, it will not constitute a segment and thus be excluded from the data frame, even though the chromsome may be included in the plot. 

If you want to get rid of the green significance bars in the plot, you can set qcap = 0. If you insist on getting rid of the entire secondary axis, save the plot to an object, then run: plotobject + scale_y_continuous(name = "copies", sec.axis = sec_axis(~., labels = NULL), expand = c(0,0))
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{templatefromequalsegments}}
}
\examples{
## simulated data assuming each chromosome comprises 100 bins
s1 <- jitter(c(1, 1, 0.8, 1.2, rep(1, 5), 1.4, rep(1, 13)), amount = 0)
s2 <- jitter(c(1, 1, 1.25, rep(1, 5), 1.5, rep(1, 13)), amount = 0)
n1 <- c(100, 100, 40, 60, rep(100, 5), 100, rep(100, 13))
n2 <- c(rep(100, 22))
bin <- 1:2200
chr <- rep(1:22, each = 100)
start <- rep(0:99*1000000+1, 22)
end <- rep(1:100*1000000, 22)
copynumbers1 <- jitter(rep(s1,n1), amount = 0.05)
copynumbers2 <- jitter(rep(s2,n2), amount = 0.05)
segments1 <- rep(s1, n1)
segments2 <- rep(s2, n2)
template1 <- data.frame(bin = bin, chr = chr, start = start, end = end, 
  copynumbers = copynumbers1, segments = segments1)
template2 <- data.frame(bin = bin, chr = chr, start = start, end = end, 
  copynumbers = copynumbers2, segments = segments2)
twosamplecompare(template1 = template1, template2 = template2, 
  cellularity1 = 0.4, cellularity2 = 0.5)
twosamplecompare(template1 = template1, template2 = template2, 
  cellularity1 = 0.4, cellularity2 = 0.5, equalsegments = 20)
## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
## for derivations of the parameters for this fit, see squaremodel
twosamplecompare(copyNumbersSegmented, index1 = 1, cellularity1 = 0.4, 
  standard1 = 1, index2 = 2, cellularity2 = 0.41, ploidy2 = 2.08, 
  standard2 = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
