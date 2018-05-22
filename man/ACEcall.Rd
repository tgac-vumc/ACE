\name{ACEcall}
\alias{ACEcall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Categorize and plot subclonal, single, and double gains / losses
}
\description{
ACEcall is the calling algorithm for ACE that utilizes the absolute copy number scaling to distinguish subclonal, single, and double losses and gains. For a segment to be "called", its segment mean needs to deviate from the general tumor ploidy both based on statistical significance and a large enough difference (both parameters can be set by the user). 
}
\usage{
ACEcall(template, QDNAseqobjectsample = FALSE, cellularity = 1, ploidy = 2, standard, 
        plot = TRUE, title, pcutoff = -3, qcutoff = -3, subcutoff = 0.2, 
        trncname = FALSE, cap = 12, bottom = 0, chrsubset, onlyautosomes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{Object. Either a template data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to analyze from the QDNAseq-object. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{cellularity}{Numeric. Used for rescaling bin and segment values. Printed on graph. Default = 1}
  \item{ploidy}{Integer. Assume the median of segments has this absolute copy number. Default = 2}
  \item{standard}{Numeric. Forces the given ploidy to represent this raw value. When omitted, the standard will be calculated from the data. When using parameters obtained from \code{squaremodel}, specify standard = 1}
  \item{plot}{Logical. Specifies whether the plot is created. When FALSE, only the data frame will be returned. Default = TRUE}
  \item{title}{Character string. Overwrites the automatically generated title}
  \item{pcutoff}{Numeric. Specifies the cutoff for statistical significance, without multiple testing correction, required to call a segment. Provide the desired log10-transformed p-value cutoff. When omitted (default), the q-value cutoff will be used instead.}
  \item{qcutoff}{Numeric. Specifies the cutoff for statistical significance, with multiple testing correction, required to call a segment. Provide the desired log10-transformed q-value cutoff. Default = -3}
  \item{subcutoff}{Numeric. Specifies the difference from the general tumor ploidy required to call a segment. Also used as margin for the other calls. Default = 0.2}
  \item{trncname}{Logical. In case of a QDNAseq object, the name of the sample is retrieved from the object and used as title. If set to TRUE, \code{trncname} truncates the sample name from the first instance of "_" in the name. You can also specify the regular expression here, e.g. trncname = "-.*" truncates the name from the first dash. Default = FALSE}
  \item{cap}{Integer. Influences your output copy number graphs. The upper limit of the y-axis is set at this number. When set to "max", it sets the cap to the maximum absolute copynumber value, rounded up. Bins and segments that exceed the cap are represented by a special mark. Recommended use between 8 and 16. Default = 12}
  \item{bottom}{Integer. Similar to cap, but for the lower limit of the y-axis. When set to "min", it sets the bottom to the minimum absolute copynumber value, rounded down. Bins and segments that subceed the bottom are represented by a special mark. Default = 0}
  \item{chrsubset}{Integer vector. Specify the chromosomes you want to plot. It will always take the full range of chromosomes in your subset, so specifying chrsubset = c(4, 8) will give the same plot as chrsubset = 4:8. When using a subset, \code{ACEcall} will not plot the cellularity and error on the plot.}
  \item{onlyautosomes}{Logical or integer. You can fill in an integer to specify how many autosomes your species has. When TRUE, \code{ACEcall} defaults to 22 (human) autosomes. When FALSE, \code{ACEcall} will also plot whichever other chromosomes are specified in the template, e.g. "X", "Y", "MT"}
}
\details{
The color code in the plot is as follows (in parentheses the call in the data frame):
  \itemize{
    \item{\strong{black}       - not called (0)}
    \item{\strong{gold}        - subclonal gain (0.5)}
    \item{\strong{turquoise}   - subclonal loss (-0.5)}
    \item{\strong{dark orange} - single copy gain (1)}
    \item{\strong{blue}        - single copy loss (-1)}
    \item{\strong{red}         - double copy gain (2)}
    \item{\strong{dark blue}   - double copy or full loss (-2)}
    \item{\strong{purple}      - amplification (3)}
  }
}
\value{
ACEcall returns a data frame that is similar to the input template, but it is supplemented with adjusted copy numbers, adjusted segment values, the segment mean and associated standard error, and the log10 p-value and q-value, which signifies the probability that if the true segment mean equals the general tumor ploidy (rounded to an integer), the resulting or a more extreme segment mean would be found. The q-value is the result of a Benjamini-Hochberg correction of the p-value, taking into account all tested segments. If the argument plot = TRUE, ACEcall will also return an absolute copy number plot (a ggplot2-object) with the segments color-coded to specify the calls. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jos B. Poell
}
\note{
For biological reasons, a segment is called a double loss when it is highly likely that some cells do not harbor a single copy of the segment of interest. Therefore, a segment with exactly one copy in a triploid tumor will be called as a single loss, not a double loss. Note that \code{ACEcall} compares segments to integer copy numbers. It therefore needs to round the given ploidy to an integer as well. It uses the \code{round} function for this. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
ACEcall(template, cellularity = 0.4, title = "sim")
## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
## for derivations of the parameters for this fit, see documentation on squaremodel
ACEcall(copyNumbersSegmented, QDNAseqobjectsample = 2, cellularity = 0.41, ploidy = 2.08, standard = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
