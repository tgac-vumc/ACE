\name{templatefromequalsegments}
\alias{templatefromequalsegments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Create a template data frame with artificial segments of equal length
}
\description{
Divide chromosomes in artificial segments with the specified number of bins. You can provide a template without segment values (for instance, an unsegmented QDNAseq-object) as input for this function. You can also use this function to only resegment a subset of chromosomes.
}
\usage{
templatefromequalsegments(template, QDNAseqobjectsample = FALSE, 
  equalsegments = 20, funtype = 'mean', chrsubset, 
  onlyautosomes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{template}{Object. Either a data frame as created by \code{\link{objectsampletotemplate}}, or a QDNAseq-object}
  \item{QDNAseqobjectsample}{Integer. Specifies which sample to analyze from the QDNAseqobject. Required when using a QDNAseq-object as template. Default = FALSE}
  \item{equalsegments}{Integer. \code{templatefromequalsegments} "resegments" simply by creating segments containing as many bins as specified in this argument. Default = 20}
  \item{funtype}{Character string. Specifies the function used to calculate new segment values. Alternative is "median". Default = "mean"}
  \item{chrsubset}{Integer vector. Specify the chromosomes you want to resegment}
  \item{onlyautosomes}{Logical or integer. You can fill in an integer to specify how many autosomes your species has. When TRUE, \code{twosamplecompare} defaults to 22 (human) autosomes. When FALSE, \code{templatefromequalsegments} will also resegment whichever other chromosomes are specified in the template, e.g. "X", "Y", "MT"}
}
\details{
\code{templatefromequalsegments} resegments the input template chromosome by chromosome. It makes segments with the number of bins specified in the argument equalsegments. Any leftover bins are divided equally among all newly created segments. Only bins with values (!is.na) are considered. If there are fewer bins on a chromosome than two times the specified value, then all bins of this chromosome are put in a single segment.
}
\value{
Returns a template data frame
}

\author{
Jos B. Poell
}
\note{
The chrsubset argument works slightly differently compared to the plotting functions. It does not necessarily segment only contiguous chromosomes. For instance, specifying chrsubset = c(8,12) will only resegment chromosomes 8 and 12. For the other chromosomes, segment values will not be changed.  
}

\seealso{
\code{\link{twosamplecompare}}, \code{\link{correlationmatrixadjusted}}
}

\examples{
data("copyNumbersSegmented")
tfes <- templatefromequalsegments(copyNumbersSegmented, 
  QDNAseqobjectsample = 2)
}

