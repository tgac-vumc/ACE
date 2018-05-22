\name{objectsampletotemplate}
\alias{objectsampletotemplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts data of a sample in a QDNAseq-object to a template for ACE functions
}
\description{
\code{objectsampletotemplate} parses QDNAseq-objects into the dataframe structure used 
by the \code{\link{singlemodel}}, \code{\link{squaremodel}}, and \code{\link{singleplot}} functions. 
These functions call \code{objectsampletotemplate} itself when necessary, but it might be convenient 
to make a template if you expect some repeated use of the functions or if you 
want to make manual adjustments to the template.
}
\usage{
objectsampletotemplate(copyNumbersSegmented, index = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copyNumbersSegmented}{QDNAseqCopyNumbers object with segmented data}
  \item{index}{Integer. Specifies the index of the sample}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
If segment values are missing, the function will still return a data frame, but only containing the copynumbers column, not the segments column
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{singlemodel}}, \code{\link{squaremodel}}, \code{\link{singleplot}}
}
\examples{
data("copyNumbersSegmented")
template <- objectsampletotemplate(copyNumbersSegmented, index = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
