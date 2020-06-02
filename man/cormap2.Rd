\name{cormap2}
\alias{cormap2}
\title{Draw correlation maps from a large dataset.}
\usage{
  cormap2(eset, cormat = NULL, lab = NULL, biomart = TRUE,
          cluster_correlations = TRUE, main = "", cex = 0.7, minfrac = 0.1,
          cor.thr = NULL, cor.mar = 0.8, cut.thr = NULL, cut.size = 1,
          labelheight = 0.2, labelwidth = 0.2, add.sig = FALSE,
          genes2highl = NULL, doPlot = TRUE, updateProgress = NULL)
}
\arguments{
  \item{eset}{ExpressionSet object or numeric matrix.}

  \item{cormat}{\code{matrix}. Numerical matrix of correlation values.}

  \item{lab}{\code{data.frame}. Optional data frame with one column for ID conversion. The column should have the target IDs and row names are the input IDs.}

  \item{biomart}{\code{logical}. Should Biomart be used for ID conversion?}

  \item{cluster_correlations}{\code{logical}. Should the correlation matrix be clustered before drawing the heatmap? If \code{TRUE} the correlation matrix is converted to a distance matrix and clustered.}

  \item{main}{\code{character}. Plot title.}

  \item{cex}{\code{numeric}. Factor to apply to default text size.}

  \item{minfrac}{\code{numeric}. Percent of NAs allowed per row. Only rows and columns with less than the given fraction of missing values are retained; defaults to 10 per cent.}

  \item{cor.thr}{\code{numeric}. Correlation threshold. This value is sign-sensitive: a negative threshold will retain rows and columns of the correlation matrix with correlation values between -1 and the threshold, a positive value will retain rows and columns with values between the threshold and 1. Zero (0) is treated as positive. Default is \code{NULL}, disabling filtering.}

  \item{cor.mar}{\code{numeric}. Margin of values needed meeting \code{cor.thr}. Only rows (and columns) with at least (100 * \code{cor.mar}) per cent of columns (or rows, respectively) meeting \code{cor.thr} and vice versa are retained. Defaults to 80 per cent.}

  \item{cut.thr}{\code{numeric}. Threshold for tree cutting, i.e., the height at which branches are to be cut. Passed to \code{WGCNA::cutreeStatic}.}

  \item{cut.size}{\code{numeric}. The minimum number of object on a branch to be considered a cluster. Passed to \code{WGCNA::cutreeStatic}.}

  \item{labelheight}{\code{numeric}. Height of lower plot margin. Can be increased to fit long gene names.}

  \item{labelwidth}{\code{numeric}. Width of right plot margin. Can be increased to fit long gene names.}

  \item{add.sig}{\code{logical}. Should significance stars be added to the plot? See details section.}

  \item{genes2highl}{\code{character}. Vector of gene symbols that should be highlighted in the heatmap.}

  \item{doPlot}{\code{logical}. Should the heatmap be produced? Set to \code{FALSE}, if only the correlation matrix should be returned.}

  \item{updateProgress}{\code{function}. Optional function to update a progress bar in a Shiny app.}
}
\value{
  The correlation matrix.
}
\description{
  Draw correlation maps from a large dataset.
}
\details{
  P-Values are calcluated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
  where r is the correlation coefficient, n is the number of samples with no missing values for each gene (row-wise \code{ncol(eset)}
                                                                                                           minus the number of columns that have an NA). (Source: http://stats.stackexchange.com/questions/120199/calculate-p-value-for-the-correlation-coefficient)
  P-Values are the calculated using \code{\link{pt}} and corrected account for the two-tailed nature of the test, i.e., the possibility of
  positive as well as negative correlation.
}
\seealso{
  \code{\link{pt}}

  \code{\link{tcrossprod}}
}

