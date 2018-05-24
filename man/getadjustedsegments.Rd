\name{getadjustedsegments}
\alias{getadjustedsegments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a data frame with segment information corresponding to a model
}
\description{
\code{getadjustedsegments} applies model parameters to segment data and returns a data frame with information of the individual segments, scaled according to the model.
}
\usage{
getadjustedsegments(template, QDNAseqobjectsample = FALSE, 
                    cellularity = 1, ploidy = 2, standard, 
                    log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{Object. Either a data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to analyze from the QDNAseq-object. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{cellularity}{Numeric. Used for rescaling bin and segment values. Default = 1}
  \item{ploidy}{Integer. Assume the median of segments has this absolute copy number. Default = 2}
  \item{standard}{Numeric. Force the given ploidy to represent this raw value. When omitted, the standard will be calculated from the data. When using parameters obtained from \code{squaremodel}, specify standard = 1}
  \item{log}{Logical. When TRUE, log2-values are calculated straigth from raw data, unadjusted! Convenience function to resemble DNAcopy output as used for ABSOLUTE and others. Default = FALSE}
}
\details{
The output contains two columns for segment mean. The first is the adjusted segment value, the second (Segment_Mean2) is the mean of the adjusted copy number values. I do not know how the QDNAseq or DNAcopy calculates the segment mean, but there is always a very small difference between the two. The P_log10 is the 10-base log of the two-sided probability that the real segment mean is the integer closest to the segment mean. While this gives an indication of the chance that a segment is subclonal, it should be interpreted with care. Because segments usually comprise many bins, these values can easily be very low. A small bias in the normalization can cause "significant", but not necessarily relevant results.
}
\value{
Returns a data frame with segment information
}

\author{
Jos B. Poell
}

\seealso{
  \code{\link{analyzegenomiclocations}}, \code{\link{postanalysisloop}}
}
\examples{
## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
singlemodel(copyNumbersSegmented, QDNAseqobjectsample = 2)
getadjustedsegments(copyNumbersSegmented, QDNAseqobjectsample = 2, 
  cellularity = 0.39)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
