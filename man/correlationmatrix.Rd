\name{correlationmatrix}
\alias{correlationmatrix}
\alias{correlationmatrixadjusted}
\title{
Create a correlation matrix of all samples in a QDNAseq-object
}
\description{
\code{correlationmatrix} runs the R stats \code{cor} function on the segmented data in a QDNAseqobject and returns the correlation matrix. The adjusted version equalizes the segments for each individual comparison, either through \code{\link{twosamplecompare}} or \code{\link{templatefromequalsegments}}
}
\usage{
correlationmatrix(object, trncname = FALSE)

correlationmatrixadjusted(object, trncname = FALSE, equalsegments = FALSE, funtype = 'mean')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{QDNAseq-object}
  \item{trncname}{Logical. The name of the sample is retrieved from the object and used as title. If set to TRUE, \code{trncname} truncates the sample name from the first instance of "_" in the name. You can also specify the regular expression here, e.g. trncname = "-.*" truncates the name from the first dash. Default = FALSE}
  \item{equalsegments}{Logical or integer. If FALSE, the function will simply take all combined breakpoints and "resegment" both samples accordingly. When an integer is given, the function will create artificial segments that are all comprised of as many bins as entered with this argument, or 20 when set to TRUE. Default = FALSE}
  \item{funtype}{Character string. Specifies the function used to calculate new segment values. Alternative is "median". Default = "mean"}
}
\details{
Calculating and plotting correlation of segments between samples can be helpful to examine similarity / dissimilarity among samples. 
}
\value{
Returns a matrix with sample names defining both rows and columns and cells containing the pearson correlation of segment values of all bins of the intersecting samples. 
}

\author{
Jos B. Poell
}
\note{
It is possible to visualize a correlation matrix with the R stats \code{heatmap} function. Although it is undoubtedly possible to adjust the function to give an interpretable plot, there are packages that facilitate this greatly. An example is the \code{\link[corrplot]{corrplot}} function from the eponymous package. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{twosamplecompare}}, \code{\link{templatefromequalsegments}}
}
\examples{
## using segmented data from a QDNAseq-object
data("copyNumbersSegmented")
correlationmatrix(copyNumbersSegmented)
correlationmatrixadjusted(copyNumbersSegmented)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
