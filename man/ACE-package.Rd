\name{ACE-package}
\alias{ACE-package}
\alias{ACE}
\docType{package}
\title{ACE package}

\description{
ACE is developed to analyze (low-coverage) whole genome sequencing data, infer absolute copies of chromosomal segments, estimate tumor cell percentage, and create visually appealing and highly interpretable copy number profiles. For preprocessing it utilizes the \code{\link[QDNAseq:QDNAseq-package]{QDNAseq package}}, but in principle any segmented copy number data can be used. The core function, \code{\link{runACE}}, performs the whole pipeline from mapped reads to suggested copy number profiles. Functions such as \code{\link{singlemodel}} and \code{\link{singleplot}} enable further inspection and customization of individual samples. For further instructions, please consult the vignette for a walkthrough, or the individual function documentation for detailed function specifics.
}

\section{Details}{
ACE basically starts from segmented copy number data. To arrive at this point from mapped sequencing reads, it runs a set of default QDNAseq functions to obtain \code{\link[QDNAseq]{QDNAseqCopyNumbers}} objects with segmented data. For copy number plots, ACE uses \code{\link[ggplot2:ggplot2-package]{ggplot2}}. The following functions comprise the main functionality of ACE:
\describe{
\item{\code{\link{runACE}}}{Perform model fitting on all samples in a directory}
\item{\code{\link{singlemodel}}}{Perform model fitting on a single sample}
\item{\code{\link{singleplot}}}{Create an absolute copy number plot of a single sample}
\item{\code{\link{squaremodel}}}{Perform model fitting on a single sample using both cellularity and ploidy as variable}
\item{\code{\link{getadjustedsegments}}}{Use obtained model parameters to calculate adjusted segment values}
\item{\code{\link{linkmutationdata}}}{Link mutation data to copy number data to estimate the number of mutant copies}
\item{\code{\link{analyzegenomiclocations}}}{Look up adjusted segment values of specific genomic locations}
\item{\code{\link{postanalysisloop}}}{Use obtained model parameters for all samples in a QDNAseq-object to print final output, e.g. copy number plots, adjusted segment data, mutation data linked to adjusted copy number data}
}
The following functions provide some more "niche" functionality
\describe{
\item{\code{\link{ACEcall}}}{Call losses and gains for all segments and visualize in copy number profile}
\item{\code{\link{twosamplecompare}}}{Compare segments and copy number profile of two samples and plot overlay}
\item{\code{\link{squaremodelsummary}}}{Return or print a summary of a squaremodel}
\item{\code{\link{loopsquaremodel}}}{Perform squaremodel fitting on all samples in a QDNAseq-object}
}
}

\section{License}{
This package is licensed under GPL.
}

\author{Jos B. Poell}

