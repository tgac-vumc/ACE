\name{copyNumbersSegmented}
\alias{copyNumbersSegmented}
\docType{data}
\title{
Segmented data of two tumor samples
}
\description{
A QDNAseq-object with segmented copy number data of two tumor samples:
'sample1' and 'sample2'
}
\usage{data("copyNumbersSegmented")}
\format{
  A QDNAseqCopyNumbers S4 object with 3113 features (bins), each 
  represented by copy number values and segment values of two samples
  \describe{
    \item{\code{copyNumbersSegmented@assayData$copynumber}}{copy number data saved in a matrix, rownames are bins (chr:start-end), colnames are sampleNames}
    \item{\code{copyNumbersSegmented@assayData$segmented}}{segment values saved in a matrix, rownames are bins (chr:start-end), colnames are sampleNames}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data("copyNumbersSegmented")
## ploidyplotloop(copyNumbersSegmented, currentdir = ".")
singlemodel(copyNumbersSegmented, which(sampleNames(copyNumbersSegmented)=='sample1'))
singleplot(copyNumbersSegmented, 1, cellularity = 0.79)
model2 <- singlemodel(copyNumbersSegmented, which(sampleNames(copyNumbersSegmented)=='sample2'))
singleplot(copyNumbersSegmented, 2, cellularity = with(model2, minima[which(rerror==min(rerror))]))
}
\keyword{datasets}
