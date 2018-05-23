\name{squaremodel}
\alias{squaremodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate potential fits for a single sample using ploidy as a variable
}
\description{
\code{squaremodel} performs a "two-dimensional" fitting algorithm on a single sample. It calculates the error of the fit at each cellularity over a range of "ploidies".  Input can be either a template or a QDNAseq-object with the index of the sample specified. Returns a list with input parameters (method, penalty, and penploidy) and model characteristics (an error matrix, a logical matrix specifying minima, a data frame with all information, a data frame with only minima, and a graphical representation of the error matrix).
}
\usage{
squaremodel(template, QDNAseqobjectsample = FALSE, prows=100, ptop=5, pbottom=1, method = 'RMSE', 
            exclude = c(), penalty = 0, penploidy = 0, highlightminima = TRUE)
}

\arguments{
  \item{template}{Object. Either a template data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to analyze from the QDNAseqobject. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{prows}{Integer. Sets the resolution of the ploidy-axis. Determines how many decrements are used to get from ptop to pbottom (see below). Therefore, the actual number of rows is actually prows + 1. Default = 100}
  \item{ptop}{Numeric. Sets the highest ploidy at which to start testing fits. Default = 5}
  \item{pbottom}{Numeric. Sets the lowest ploidy to be tested. Default = 1}
  \item{method}{Character string specifying which error method to use. For more documentation, consult the vignette. Can be "RMSE", "SMRE", or "MAE". Default = "RMSE"}
  \item{exclude}{Integer or character vector. Specifies which chromosomes to exclude for model fitting}
  \item{penalty}{Numeric. Penalizes fits at lower cellularities. Suggested values between 0 and 1. Default = 0 (no penalty)}
  \item{penploidy}{Numeric. Penalizes fits that diverge from 2N with the formula (1+abs(ploidy-2))^penploidy. Default = 0}
  \item{highlightminima}{Logical. Minima are highlighted in the matrixplot by a black dot. Default = TRUE}
}

\details{
Unlike other functionality of ACE, \code{squaremodel} does not use the "standard", but it fits all tested ploidies to "standard = 1". It is therefore necessary that segment values are normalized to 1 (which they are by default coming from QDNAseq). The penalty parameter is the same as in \code{singlemodel}. Additionally, it is possible to penalize fits at ploidies diverging from 2N using the penploidy parameter. For other details on the fitting algorithm, see \code{\link{singlemodel}}. Range of ploidies is set by parameters ptop and pbottom, and resolution is determined by prows. To create good contrast in the matrixplot, the color scale derives from the inverse of the error, and the opacity of the dots marking the minima is calculated as min(error)/error.
}
\value{
Returns a list, containing
  \item{method}{Applied error method}
  \item{penalty}{Applied penalty factor for low cellularities}
  \item{penploidy}{Applied penalty factor for diverging ploidies}
  \item{errormatrix}{Numeric matrix with errors of all combinations of ploidy and cellularity}
  \item{minimatrix}{Logical matrix indicating whether the combination of ploidy and cellularity represents a minimum}
  \item{errordf}{Data frame with columns ploidy, cellularity, error, and minimum}
  \item{minimadf}{Same as errordf, but only containing minima and sorted by error}
  \item{matrixplot}{ggplot2-graph of the relative errors calculated at each combination of ploidy and cellularity}
}

\author{
Jos B. Poell
}
\note{
squaremodel() only needs a data frame with columns named \code{chr} and \code{segments}. Every row should contain an individual genomic feature, i.e. a bin or a probe. If you have data with each row representing a segment, and the size of the segment given in a column (e.g. NumBins or NumProbes), you can create the data frame as follows (giving the correct variable names of course):
\code{template <- data.frame(chr = rep(Chromosome, NumProbes), segments = rep(SegmentMean, NumProbes))}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{objectsampletotemplate}}, \code{\link{squaremodel}}, \code{\link{singleplot}}
}
\examples{
## toy data assuming each chromosome comprises 100 bins
s <- jitter(c(1, 1, 0.8, 1.2, rep(1, 5), 1.4, rep(1, 13)), amount = 0)
n <- c(100, 100, 40, 60, rep(100, 5), 100, rep(100, 13))
df <- data.frame(chr = rep(1:22, each = 100), segments = rep(s, n))
squaremodel(df)$matrixplot
sm <- squaremodel(df, method = 'MAE', penalty = 0.5, penploidy = 0.5)
sm$matrixplot
mdf <- sm$minimadf
head(mdf[order(mdf$error,-mdf$cellularity),])

## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
sqm <- squaremodel(copyNumbersSegmented, QDNAseqobjectsample = 2, penalty = 0.5, penploidy = 0.5, 
  ptop = 4.3, pbottom = 1.8, prows = 250)
sqm$matrixplot
mdf <- sqm$minimadf
head(mdf[order(mdf$error,-mdf$cellularity),])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
