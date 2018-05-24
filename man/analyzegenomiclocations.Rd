\name{analyzegenomiclocations}
\alias{analyzegenomiclocations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve adjusted copy number information for specific genomic locations
}
\description{
\code{analyzegenomiclocations} searches an adjusted segment data frame for the specified genomic locations and reports the associated segment mean as "Copynumbers". If frequency of a variant is given (in percentage), mutant copies are also calculated. Make sure the specified cellularity is the same number as was used to create the adjusted segment data frame.
}
\usage{
analyzegenomiclocations(segmentdf, Chromosome, Position, 
                        Frequency, cellularity)
}

\arguments{
  \item{segmentdf}{Data frame. Output of \code{getadjustedsegments}}
  \item{Chromosome}{Vector. Specifies the chromosome(s)}
  \item{Position}{Numeric vector. Specifies the base position(s) of interest}
  \item{Frequency}{Numeric vector. Optional. Used (in conjuction with cellularity) to calculate mutant copies. Enter a percentage.}
  \item{cellularity}{Numeric. Only required when calculating mutant copies. Use same number as used with \code{getadjustedsegments}}
}
\details{
The formula that calculates mutant copies only works if the variant is not present in normal tissue. In other words, make sure variants are not of germline origin. If you are interested in a heterozygous germline variant, you can calculate the "mutant copies" by subtracting (1/cellularity - 1).
}
\value{
Returns a data frame with Chromosome, Position, Copynumbers. When \code{Frequency} is given as an argument, also returns Mutant_copies.
}

\author{
Jos B. Poell
}
\note{
Chromosome, Position, and Frequency can be single values or vectors. In the latter case, they need to be of equal length. Make sure the position you enter matches with the genome build used for alignment / binning of sequence reads.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getadjustedsegments}}, \code{\link{linkmutationdata}}
}
\examples{
## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
segmentdf <- getadjustedsegments(copyNumbersSegmented, 
  QDNAseqobjectsample = 2, cellularity = 0.39)
analyzegenomiclocations(segmentdf, Chromosome = 5, Position = 26365569)
analyzegenomiclocations(segmentdf, Chromosome = 5, Position = 26365569, 
  Frequency = 25, cellularity = 0.39)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
