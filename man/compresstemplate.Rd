\name{compresstemplate}
\alias{compresstemplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduce the size of a template data frame
}
\description{
\code{compresstemplate} combines bins to reduce the size of a template data frame. 
}
\usage{
compresstemplate(template, factor = 20, funtype = "median")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{Data frame as created by \code{\link{objectsampletotemplate}}}
  \item{factor}{Integer. Number of bins to combine into a single bin. Default = 20}
  \item{funtype}{Character string. Specifies the function used to combine bins. Either "mean" or "median". Default = "median"}
}
\details{
\code{compresstemplate} can be used when the number of features in a template is too large for downstream analyses. Bins are combined, resulting in larger bin sizes and a reduced number of features.  
}
\value{
Returns a template data frame with increased bin size and reduced number of features compared to the input data frame.
}

\author{
Jos B. Poell
}
\note{
Segmentation is slightly affected by this function due to the fact that the breakpoints will usually fall somewhere within a newly created bin. In case of funtype "mean" the segment value of this bin will be different from both adjacent bins and therefore create a 1-bin-segment. This "oversegmentation" is greatly reduced when using funtype "median" (especially when argument factor is and odd number).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{objectsampletotemplate}}, \code{\link{segmentstotemplate}}
}

\examples{
## segmented data from a QDNAseq-object
data("copyNumbersSegmented")
template <- objectsampletotemplate(copyNumbersSegmented, index = 1)
ct <- compresstemplate(template, factor = 10)
length(template$bin)
length(ct$bin)
}
