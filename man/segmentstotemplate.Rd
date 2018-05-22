\name{segmentstotemplate}
\alias{segmentstotemplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a template data frame from input that only provides segment information
}
\description{
\code{segmentstotemplate} "explodes" segment data and creates a data frame with information for all bins. This enables the use of data limited to segments in ACE-related functions. 
}
\usage{
segmentstotemplate(segmentdf, chrci = 1, startci = 2, endci = 3, binsci = 4, 
                   meanci = 5, seci, sdci, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segmentdf}{Data frame or file path. Data frame or file path to tab-delimited file with segment data}
  \item{chrci}{Integer. Specifies the index of the column containing chromosome information. Default = 1}
  \item{starti}{Integer. Specifies the index of the column containing start positions. Default = 2}
  \item{endci}{Integer. Specifies the index of the column containing end positions. Default = 3}
  \item{binsci}{Integer. Specifies the index of the column containing number of bins or number of probes. Default = 4}
  \item{meanci}{Integer. Specifies the index of the column containing the segment value. Default = 5}
  \item{seci}{Integer. Optional. Specifies the index of the column containing standard errors of the segment mean. See details below}
  \item{sdci}{Integer. Optional. Specifies the index of the column containing standard deviations of the segment mean. See details below}
  \item{log}{Logical or integer. If the data is given in log-values, the data is converted to linear scale using the given log base, or inverting the natural logarithm when set to TRUE. Default = FALSE}
}
\details{
Even though there is no copy number information for each bin, this needs to be provided for several ACE-related functions. By default, \code{segmentstotemplate} will just "copy" the segment data. If the segment data comes with either standard deviations or standard errors, it is possible to "simulate" copy number data. To do so, \code{segmentstotemplate} will use the R stats \code{rnorm} function.  
}
\value{
Returns a template data frame. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jos B. Poell
}
\note{
\itemize{
\item Keep in mind that copy number data is used in some functions, notably \code{\link{ACEcall}}. If information on the variance of copy number values within segments is missing, calls are only based on adjusted segment values, because the p-values are all 0. 
\item Also note that resulting templates should not be used for any of the functions that perform resegmentation (e.g. \code{\link{twosamplecompare}}), or the results should be interpreted with due caution. 
\item The functions \code{\link{singlemodel}} and \code{\link{squaremodel}} only use  segment values and can be used to analyze the resulting template data frames. \code{\link{singleplot}} can also be used without problem.
\item This function works fine with high resolution data, but the resulting template can be a bit much for other functions. For instance, the segmented data files provided by TCGA are the result of SNP arrays that have a resolution of roughly 1.5 kbp. The function \code{\link{compresstemplate}} can help out to make input more manageable.
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{compresstemplate}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
