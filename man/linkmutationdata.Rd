\name{linkmutationdata}
\alias{linkmutationdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Append columns with total genomic copies and mutant copies to a file with mutation data
}
\description{
\code{linkmutationdata} combines copy number data, estimated tumor cell percentage, and variant allele frequency (e.g. mutation data) to calculate how many variant (mutant) copies the tumor genome harbors. It requires a tab-delimited file with mutation data and a data frame with adjusted segment data as obtained with \code{getadjustedsegments}. Also make sure to provide the correct (estimated) cellularity. Output is a file with "_ACE" added to the original file name before the extension. It can either be a copy of the input with Copynumbers and Mutant_copies appended as extra columns at the end, or a file with the columns Chromosome, Position, Frequency, Copynumbers and Mutant_copies. To perform batch analysis, use \code{\link{postanalysisloop}}.
}
\usage{
linkmutationdata(mutationdf, segmentdf, cellularity = 1, 
                 chrindex = 1, posindex = 2, freqindex = 3, 
                 append = TRUE, outputdir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutationdf}{Data frame or character string. File path to tab-delimited text (either .tsv, .csv, .txt or .xls) containing mutation data or the corresponding data frame. File must contain a header and columns for chromosome, position, and frequency of the mutation. Additional columns are ignored.}
  \item{segmentdf}{Data frame or character string containing file path of tab-delimited text with segment data. Expects data in the format provided by getadjustedsegments with argument log=FALSE.}
  \item{cellularity}{Numeric. Used to infer mutant copies from frequency and total copies. Default = 1}
  \item{chrindex}{Integer. Column index in input file specifying the chromosome associated with the genomic location. Default = 1}
  \item{posindex}{Integer. Column index in input file specifying the position on the chromosome associated with the genomic location. Default = 2}
  \item{freqindex}{Integer. Column index in input file specifying the frequency (as a percentage) of the variant. Default = 3}
  \item{append}{Logical. When TRUE, appends the output columns to the original mutation input file, but it still saves the result in a new file. When FALSE, the output file will only contain the columns "Chromosome", "Position", "Frequency", "Copynumbers", and "Mutant_copies". Default = TRUE}
  \item{outputdir}{Character string. Convenience function to save output into a custom directory.}
}
\details{
The formula that calculates mutant copies only works if the variant is not present in normal tissue. In other words, make sure variants are not of germline origin. If you are interested in a heterozygous germline variant, you can calculate the "mutant copies" by subtracting (1/cellularity - 1).
}
\value{
Prints output to a tab-delimited file, or returns a data frame with columns added for copies and mutant copies.
}

\author{
Jos B. Poell
}
\note{
Make sure the mutation data matches with the genome build used for alignment / binning of sequence reads for copy number analysis. If the resulting Mutant_copies are very low, the variant allele frequencies were probably provided as fraction, not percentage. Just multiply by 100.
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
mutationdf <- data.frame(Gene, Chromosome, Position, Frequency)
linkmutationdata(mutationdf, segmentdf, cellularity = 0.38, 
                 chrindex = 2, posindex = 3, freqindex = 4)
}

