\name{linkvariants}
\alias{linkvariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Append columns with total genomic copies and mutant copies to a file with variant/mutation data
}
\description{
\code{linkvariants} combines copy number data, estimated tumor cell percentage, and variant allele frequency (e.g. mutation data) to calculate how many variant (mutant) copies the tumor genome harbors. It requires a data frame or tab-delimited file with variant data and a data frame or tab-delimited file with adjusted segment data as obtained with \code{getadjustedsegments}. Also make sure to provide the correct (estimated) cellularity. Output is a file with "_ACE" added to the original file name before the extension. It can either be a copy of the input with Copynumbers and Variant_copies appended as extra columns at the end, or a file with the columns Chromosome, Position, Frequency, Copynumbers and Mutant_copies. To perform batch analysis, use \code{\link{postanalysisloop}}. \code{linkvariants} can provide upper and lower bounds of a confidence interval if read depth is given.
}
\usage{
linkvariants(variantdf, segmentdf, cellularity = 1, hetSNPs = FALSE,
             chrindex=1,posindex=2,freqindex,altreadsindex,
             totalreadsindex,refreadsindex,confidencelevel=FALSE,
             append=TRUE, outputdir, sgc = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variantdf}{Data frame or character string. File path to tab-delimited text (either .tsv, .csv, .txt or .xls) containing variant data or the corresponding data frame. File must contain a header and columns for chromosome, position, and frequency of the mutation. Optionally a column with read depth information. If frequency is missing, altreads + totalreads or altreads + refreads is required.}
  \item{segmentdf}{Data frame or character string containing file path of tab-delimited text with segment data. Expects data in the format provided by getadjustedsegments with argument log=FALSE.}
  \item{cellularity}{Numeric. Used to infer variant copies from frequency and total copies. Default = 1}
  \item{hetSNPs}{Logical. If TRUE, half of the germline copies are assumed to be variant. Default = FALSE}
  \item{chrindex}{Integer. Column index in input file specifying the chromosome associated with the genomic location. Default = 1}
  \item{posindex}{Integer. Column index in input file specifying the position on the chromosome associated with the genomic location. Default = 2}
  \item{freqindex}{Integer. Column index in input file specifying the frequency (as a percentage) of the variant}
  \item{altreadsindex}{Integer. Column index in input file specifying the number of variant-supporting reads}
  \item{totalreadsindex}{Integer. Column index in input file specifying the read depth at the genomic location of the variant}
  \item{refreadsindex}{Integer. Column index in input file specifying the number of reference-supporting reads}
  \item{confidencelevel}{Numeric or logical. If read depth information is available, calculate the upper and lower bounds of this confidence level for the frequency and the number of variant copies of each variant. Will be skipped if FALSE. Default = FALSE}
  \item{append}{Logical. When TRUE, appends the output columns to the original mutation input file, but it still saves the result in a new file. When FALSE, the output file will only contain the columns "Chromosome", "Position", "Frequency", "Copynumbers", and "Mutant_copies". Default = TRUE}
  \item{outputdir}{Character string. Convenience function to save output into a custom directory}
  \item{sgc}{Integer or character vector. Specify which chromosomes occur with only a single copy in the germline}
}
\details{
The default formula that calculates mutant copies works if the variant is not present in normal tissue. If you are interested in heterozygous germline variant, you can set the argument hetSNPs to TRUE. The confidence intervals are calculated using the binom.test function. 
}
\value{
Prints output to a tab-delimited file, or returns a data frame with columns added for copies and mutant copies.
}

\author{
Jos B. Poell
}
\note{
Make sure the variant data matches with the genome build used for alignment / binning of sequence reads for copy number analysis. If the resulting Variant_copies are very low, the variant allele frequencies were probably provided as fraction, not percentage. Just multiply by 100. If your data contains sex chromosomes, then make sure to specify sgc = c("X", "Y") when analyzing data from a male individual.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getadjustedsegments}}, \code{\link{analyzegenomiclocations}}, \code{\link{postanalysisloop}}
}
\examples{
## using manually simulated mutation data
## see vignette for more practical uses
data("copyNumbersSegmented")
segmentdf <- getadjustedsegments(copyNumbersSegmented, 
  QDNAseqobjectsample = 2, cellularity = 0.38)
Gene <- c("CASP8", "CDKN2A", "TP53")
Chromosome <- c(2, 9, 17)
Position <- c(202149589, 21971186, 7574003)
Frequency <- c(47.46, 36.28, 43.48)
AltReads <- c(345, 198, 284)
variantdf <- data.frame(Gene, Chromosome, Position, Frequency, AltReads)
linkvariants(variantdf, segmentdf, cellularity = 0.38, 
             chrindex = 2, posindex = 3, freqindex = 4)
linkvariants(variantdf, segmentdf, cellularity = 0.38, 
             chrindex = 2, posindex = 3, freqindex = 4,
             altreadsindex = 5, confidencelevel = 0.9)
}

