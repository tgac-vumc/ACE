\name{runACE}
\alias{runACE}
\alias{ploidyplotloop}
\title{
  Absolute Copy number Estimation
}
\description{
ACE scales copy number data to fit with integer copy numbers, providing an 
estimate for tumor cell percentage in the process. \code{runACE} uses segmented data from 
the QDNAseq package. A folder with either bam-files (aligned sequencing 
data) or rds-files of QDNAseq-objects can be used as input. 
Model fitting and production of all output files (except "parameters.tsv") is 
executed by \code{ploidyplotloop}, which handles one QDNAseq-object at a time. 
}
\usage{
runACE(inputdir = "./", outputdir, filetype = 'rds', binsizes, 
  ploidies = 2, imagetype = 'pdf', method = 'RMSE', penalty = 0, 
  cap = 12, bottom = 0, trncname = FALSE, printsummaries = TRUE, 
  autopick = FALSE)
  
ploidyplotloop(copyNumbersSegmented, currentdir, ploidies = 2, 
  imagetype = 'pdf', method = 'RMSE', penalty = 0, cap = 12, 
  bottom = 0, trncname = FALSE, printsummaries = TRUE, 
  autopick = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputdir}{Character string specifying the directory containing the files you want analyzed. Note: will analyze ALL rds-files or bam-files in the given directory. Default = "./"}
  \item{copyNumbersSegmented}{QDNAseq-object with segmented data}
  \item{outputdir, currentdir}{Character string specifying the directory to which ACE should write the output. When missing, ACE will try to write to inputdir. For ploidyplotloop, currentdir is required.}
  \item{filetype}{Character string specifying the file type of your input, either "bam" or "rds". Default = "rds"}
  \item{binsizes}{Numeric vector, specifying which binsizes (in kbp) to analyze. Possible values are 1, 5, 10, 15, 30, 50, 100, 500, and 1000. Default = c(30,100,500,1000)}
  \item{ploidies}{Numeric vector, specifying which ploidies (N) to analyze. Use positive natural numbers. Default = 2}
  \item{imagetype}{Character string specifying the image type graphics device, default = "pdf"}
  \item{method}{Character string specifying what error method to use. See also section "Error methods". Default = "RMSE"}
  \item{penalty}{Numeric value. Penalizes fits at lower cellularities. Suggested values between 0 and 1. Default = 0 (no penalty)}
  \item{cap, bottom}{Integer. Influences your output copy number graphs. The upper and lower limits of the y-axis are set at these values. Bins and segments that exceed/subceed the cap/bottom are represented by a special mark. Default = 12 and 0 respectively}
  \item{trncname}{Logical. Convenience functionality. If all your samples have a certain extension to their name, you can use this to truncate this extension and be left with the actual sample name. When TRUE, the regular expression is "_.*". That means it will chop off everything from the sample name starting with the first underscore. Instead of a logical, you can specify a character string to match your regular expression of choice. You can test whether it will work with the gsub function, since this is what ACE uses to truncate names. Default = FALSE}
  \item{printsummaries}{Logical. If you do not want the big summary files, you can set this argument to FALSE. If you still want the summary files containing only error plots, you can set this to 2. Default = TRUE}
  \item{autopick}{Logical. When set to TRUE, ACE will fill in the cellularity of the best fit in the column likely_fit of the fitpicker file(s). Default = FALSE}
}
\details{
Since this is the core functionality of ACE, extensive documentation is available in the vignette.
}
\value{
\code{runACE} and \code{ploidyplotloop} do not return any values, they print all their output to the indicated location. The output comprises
\itemize{
    \item the file "parameters.tsv" which simply reports the used parameters
    \item rds-files (only in case you had bam-files as input)
    \item for each ploidy a "fitpicker.tsv" file which can be used for selecting the most likely fits
    \item a summary file of likely fits and error plot of each sample (if printsummaries is set to TRUE)
    \item a summary file of all error plots (if printsummaries is set to TRUE or 2)
    \item a directory with copy number plots of the likely fits of all samples
    \item a directory for each sample, containing the error plot, a summary file with all fits of the sample, and individual copy number plots of all fits in a subdirectory
    }
}


\author{
Jos B. Poell
}
\note{
You can use the example data for testing: see Examples
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{singlemodel}}, \code{\link{squaremodel}}, \code{\link{singleplot}}
}
\examples{
\dontrun{
  runACE("./bam/", outputdir = "./results", penalty = 0.5, 
    binsizes = c(100, 1000), imagetype = 'png')
  data("copyNumbersSegmented")
  ploidyplotloop(copyNumbersSegmented, ".", ploidies = c(2,3))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
