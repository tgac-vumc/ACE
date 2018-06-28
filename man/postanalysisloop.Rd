\name{postanalysisloop}
\alias{postanalysisloop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch analysis of samples in a QDNAseq-object for which models have been chosen
}
\description{
When models have been chosen for all (or just multiple) samples in a QDNAseq-object, this function can be used to perform a batch analysis on those samples. This encompasses printing segment data, printing copy number plots, and linking mutation data. 
}
\usage{
postanalysisloop(copyNumbersSegmented, modelsfile, variantdata, 
                 prefix="", postfix="", trncname=FALSE, inputdir=FALSE,
                 hetSNPs=FALSE, chrindex=1, posindex=2, freqindex=3,
                 altreadsindex, totalreadsindex, refreadsindex,
                 confidencelevel=FALSE, append=TRUE,
                 dontmatchnames=FALSE, printsegmentfiles=TRUE,
                 printnewplots=TRUE, imagetype='pdf',
                 outputdir="./", log=FALSE, segext='tsv')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copyNumbersSegmented}{QDNAseq-object with segmented data or file path of an rds-file containing a QDNAseq-object}
  \item{modelsfile}{Character string or data frame. When a character, it specifies the file path of a tab-delimited text containing model variables of samples. Expects columns with a header. It contains at least two columns: the first specifying the sample names and the second specifying the cellularity. The third column is the ploidy of the samples. When omitted, it is assumed to be 2. The fourth column is the standard of the samples. When omitted, it is calculated from the data in the object. The fitpicker.tsv file created by \code{\link{runACE}} can be used as modelsfile after cellularity of the likely fit is specified in the second column.}
  \item{variantdata}{Character string. Specifies directory containing variant data of samples. Optional. When argument inputdir is used, the function will first see if this argument is specified, if not it will check if the directory inputdir/variantdata exists, if not it will look for the variant files in inputdir itself. When inputdir is not used and this argument is omitted, the function will not link variant data. Mutation files need to have the same file extensions, which can be either .csv, .tsv, .txt, or .xls.}
  \item{prefix}{Character string. Used when a uniform character string precedes the sample name in the file name. E.g. "mutations_sample1.csv" has prefix "mutations_". Default = ""}
  \item{postfix}{Character string. As \code{prefix}, but then after the sample name. E.g. "sample1_somatics.csv" has postfix "_somatics". Default = ""}
  \item{trncname}{Logical. When TRUE, truncates the sample names of the QDNAseq-object starting from the first "_", or specify a character string with your regular expression of choice (\code{trncname} uses the \code{gsub} function). NOTE: use only when this will provide matches with the mutation files and the sample names in the modelsfile. Default = FALSE}
  \item{inputdir}{Character string. Specifies the directory which contains the files to be analyzed. Convenience function. Reduces the amount of arguments required when all data is available in the same directory: the QDNAseq-object, a file named "models.tsv" with the model parameters, and the mutation data (either in the inputdir itself or in a subdirectory "mutationdata"). Specifying the first arguments (copyNumbersSegmented, modelsfile, mutationdata) will take priority. When missing it will look in inputdir. When multiple rds-files are present in the inputdir, it will try the first one. Note: the path specified has no consequences for the location of the output. Default = FALSE}
  \item{hetSNPs}{Logical. If TRUE, half of the germline copies are assumed to be variant. Default = FALSE}
  \item{chrindex}{Integer. Column index in input file specifying the chromosome associated with the genomic location. Default = 1}
  \item{posindex}{Integer. Column index in input file specifying the position on the chromosome associated with the genomic location. Default = 2}
  \item{freqindex}{Integer. Column index in input file specifying the frequency (as a percentage) of the variant. Default = 3}
  \item{altreadsindex}{Integer. Column index in input file specifying the number of variant-supporting reads}
  \item{totalreadsindex}{Integer. Column index in input file specifying the read depth at the genomic location of the variant}
  \item{refreadsindex}{Integer. Column index in input file specifying the number of reference-supporting reads}
  \item{confidencelevel}{Numeric or logical. If read depth information is available, calculate the upper and lower bounds of this confidence level for the frequency and the number of variant copies of each variant. Will be skipped if FALSE. Default = FALSE}
  \item{append}{Logical. When TRUE, appends the output columns to the original mutation input file, but it still saves the result in a new file. When FALSE, the output file will only contain the columns "Chromosome", "Position", "Frequency", "Copynumbers", and "Mutant_copies" (and including the upper and lower bounds of the frequency and variant copies confidence interval, when applicable). Default = TRUE}
  \item{dontmatchnames}{Logical. When TRUE, the model variables are called by the index of the sample in the QDNAseq-object. This will only work if the order of samples in the object exactly matches the order of samples in the modelsfile. Use with caution! This is somewhat of an emergency option if for some reason the name matching is not working. I recommend trying to get the name matching to work. Default = FALSE}
  \item{printsegmentfiles}{Logical. When TRUE, prints a tab-delimited text file for each sample into a "segmentfiles" folder. Default = TRUE}
  \item{printnewplots}{Logical. When TRUE, prints plots into a "newplots" folder in the specified image type. Default = TRUE}
  \item{imagetype}{Character string specifying the image type graphics device. Default = "pdf"}
  \item{outputdir}{Character string. Save output into this custom directory. Default = "./"}
  \item{log}{Logical or integer. Use log conversion for creating segments output. Default = FALSE}
  \item{segext}{Character string specifying the extension for the segments output. Default = "tsv"}
}
\details{
If your input is tailored for this function, you could run it without any arguments! Most arguments help with matching sample names in the QDNAseq-object, the modelsfile, and the names and columns of the files containing variant data. You can "trim" the name of the file with variant data using the \code{prefix} (everything before the name) and \code{postfix} (everything after the name, but before the file extension) arguments to match your sample names. \code{trncname} migth help trimming the name in the QDNAseq-object, but be sure it still matches the sample names in the modelsfile (and mutation data file names when applicable). 
}
\value{
Prints the specified output to an indicated directory. Returns a list of copy number plots.
}

\author{
Jos B. Poell
}
\note{
The use of \code{inputdir} and \code{outputdir} should be fairly robust. However, using irregular file paths might cause problems. If you suspect problems with file paths, try setting the working directory to the intended \code{inputdir}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getadjustedsegments}}, \code{\link{linkvariants}}, \code{\link{runACE}}
}
\examples{
## see the vignette for examples
\dontrun{
  data("copyNumbersSegmented")
  postanalysisloop(copyNumbersSegmented, "models.tsv", "variantdata", 
  outputdir = "loop_output")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
