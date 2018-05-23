\name{loopsquaremodel}
\alias{loopsquaremodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create \code{squaremodel} summaries for all samples in a QDNAseq-object
}
\description{
\code{loopsquaremodel} performs a \code{\link{squaremodel}} fitting on all samples in a QDNAseq-object. It prints the summaries of the squaremodels to the specified output directory. It also returns a list with squaremodels which can be saved to a variable.
}
\usage{
loopsquaremodel(object, ptop = 5, pbottom = 1, prows = 100, 
                method = 'RMSE', penalty = 0, penploidy = 0,
                outputdir, imagetype = 'pdf', trncname = FALSE, 
                printplots = TRUE, printobjectsummary = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{QDNAseq-object}
  \item{ptop}{Numeric. Sets the highest ploidy at which to start testing fits. Default = 5}
  \item{pbottom}{Numeric. Sets the lowest ploidy to be tested. Default = 1}
  \item{prows}{Integer. Sets the resolution of the ploidy-axis. Determines how many decrements are used to get from ptop to pbottom (see below). Therefore, the actual number of rows is actually prows + 1. Default = 100}
  \item{method}{Character string specifying which error method to use. For more documentation, consult the vignette. Can be "RMSE", "SMRE", or "MAE". Default = "RMSE"}
  \item{penalty}{Numeric. Penalizes fits at lower cellularities. Suggested values between 0 and 1. Default = 0 (no penalty)}
  \item{penploidy}{Numeric. Penalizes fits that diverge from 2N with the formula (1+abs(ploidy-2))^penploidy. Default = 0}
  \item{outputdir}{Character string. Print the plots to this directory}
  \item{imagetype}{Character string. Plots are printed to file using this graphics device. "pdf" will result in pdf-files containing 8 pages with individual plots, while the other devices print 2x4 mosaics per sample. Default = "pdf"}
  \item{trncname}{Logical. If set to TRUE, \code{trncname} truncates the sample name from the first instance of "_" in the name. You can also specify the regular expression here, e.g. trncname = "-.*" truncates the name from the first dash. Default = FALSE}
  \item{printplots}{Logical. Print the plots to file. Note that the object summary is not affected by this argument. Default = TRUE}
  \item{printobjectsummary}{Logical. Print the object summary to file. Default = TRUE}
}
\details{
This function is basically the squaremodel equivalent of \code{\link{ploidyplotloop}}. One key difference is the output. \code{loopsquaremodel} returns a list with the squaremodels of each sample. The squaremodels are supplemented with the sample names. The output printed to file is optional (though by default enabled). 
}
\value{
Returns a list with squaremodels for each sample in a QDNAseq-object. If \code{printplots} is TRUE, the plots will be printed to file. If \code{printobjectsummary} is TRUE, an object summary will be printed, containing the matrixplot and the copy number plot of the best fit of each sample.  
}

\author{
Jos B. Poell
}
\note{
In case of large numbers of samples, you may have to set \code{printobjectsummary} to FALSE to prevent file size or memory issues.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{squaremodelsummary}}, \code{\link{squaremodel}}, \code{\link{ploidyplotloop}}
}
\examples{
\dontrun{
  data("copyNumbersSegmented")
  loopsquaremodel(copyNumbersSegmented, penalty = 0.5, penploidy = 0.5)
  }
lsm <- loopsquaremodel(copyNumbersSegmented, printplots = FALSE, printobjectsummary = FALSE, 
                       penalty = 0.5, penploidy = 0.5)
ls(lsm[[1]])
lsm[[1]]$samplename
lsm[[1]]$matrixplot
}

