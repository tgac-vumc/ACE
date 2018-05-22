\name{singlemodel}
\alias{singlemodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate potential fits for a single sample
}
\description{
\code{singlemodel} performs the basic fitting algorithm of ACE on a single sample. Input can be either a template or a QDNAseq-object with the index of the sample specified. Returns a list with input parameters (ploidy, standard, and penalty) and model characteristics (calculated minima, the relative error corresponding with the minima, and the errors calculated at every cellularity). It also returns the plot associated with the error list. The minima represent cellularities, as can be seen in the plot.
}
\usage{
singlemodel(template, QDNAseqobjectsample = FALSE, ploidy = 2, standard, 
            method = 'RMSE', exclude = c(), penalty = 0, highlightminima = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{Object. Either a template data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to analyze from the QDNAseqobject. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{ploidy}{Integer. Calculate fits assuming the median of segments has this absolute copy number. Default = 2}
  \item{standard}{Numeric. Force the given ploidy to represent this raw value. When omitted, the standard will be calculated from the data}
  \item{method}{String character specifying which error method to use. For more documentation, consult the vignette. Can be "RMSE", "SMRE", or "MAE". Default = "RMSE"}
  \item{exclude}{Integer or character vector. Specifies which chromosomes to exclude for model fitting}
  \item{penalty}{Numeric value. Penalizes fits at lower cellularities. Suggested values between 0 and 1. Default = 0 (no penalty)}
  \item{highlightminima}{Logical. Minima are highlighted in the errorplot by a red color. Default = TRUE}
}
\details{
All ACE fitting algorithms work by calculating "expected values" of integer copies given a certain cellularity. It calculates these expected values for 1-12 copies at cellularities 0.05-1 (in increments of 0.01). First of all, this means that fits at cellularities below 0.05 are not calculated. These low-cellularity fits will not give very meaningful results, and only obscure more plausible fits. Second, it means that 0 copies and >12 copies are not "fitted". This prevents fits predicting many and/or large segments with 0 or >12 copies, which is biologically unlikely. More explanation is given in the vignette.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jos B. Poell
}
\note{
singlemodel() only needs a data frame with columns named \code{chr} and \code{segments}. Every row should contain an individual genomic feature, i.e. a bin or a probe. If you have data with each row representing a segment, and the size of the segment given in a column (e.g. NumBins or NumProbes), you can create the data frame as follows (giving the correct variable names of course):
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
singlemodel(df)
singlemodel(df, ploidy = 3)
singlemodel(df, method = 'MAE', penalty = 0.5)
singlemodel(df, exclude = 1:3)

## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
singlemodel(copyNumbersSegmented, QDNAseqobjectsample = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
