\name{forcesegmentsontemplate}
\alias{forcesegmentsontemplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Custom resegmentation with user-defined segment information 
}
\description{
\code{forcesegmentsontemplate} forces a user-defined segmentation pattern on a template. The input containing the segment information requires the chromosome name, start position, and end position of the segments. This can be provided manually as a data frame or as a data frame such as obtained by \code{\link{getadjustedsegments}}. This means you can also use this function to force the segment pattern of one sample onto another sample. It is possible to retain the break points of the original input template.
}
\usage{
forcesegmentsontemplate(segmentinput, template, QDNAseqobjectsample = FALSE, 
                        combinesegments = FALSE, funtype = 'mean')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segmentinput}{Data frame containing segment information. Requires columns with chromosome name, start position, and end position. The function looks for column names containing "chr", "start", and "end" to find the required information}
  \item{template}{Object. Either a data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object. Note that segmented data is not required for this function, unless combining segments}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to use as input from the QDNAseq-object. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{combinesegments}{Logical. When TRUE, segment break points from the input template are retained. Default = FALSE}
  \item{funtype}{Character string. Specifies which function to use to calculate segment values. An alternative would be 'median'. Default = 'mean'}
}

\value{
Returns a template data frame with the columns bin, chromosome, copynumbers, and segments
}

\author{
Jos B. Poell
}
\details{
This function only changes (or adds) the segments column of the input template. Copynumbers values are always retained, even if they fall out of the range of the segmentinput.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getadjustedsegments}}, \code{\link{twosamplecompare}}
}
\examples{
data("copyNumbersSegmented")
template <- objectsampletotemplate(copyNumbersSegmented, index = 1)
segments <- getadjustedsegments(copyNumbersSegmented, 2)[,1:3]
newtemplate <- forcesegmentsontemplate(segments, template)
first50M <- forcesegmentsontemplate(data.frame(chr = 1:22, 
  start = rep(1, 22), end = rep(50000000, 22)), template)
singleplot(first50M)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
