\name{eset_cor}
\alias{eset_cor}
\title{Helper function to calculate the correlation matrix.}
\usage{
  eset_cor(eset, cormat = NULL, with.pvalues = FALSE)
}
\arguments{
  \item{eset}{\code{ExpressionSet}. Object or numerical matrix.}

  \item{cormat}{\code{matrix}. Numerical matrix of correlation values.}

  \item{with.pvalues}{\code{logical}. Should P-Values be calculated from the correlation values? If \code{TRUE} P-Values will be depicted in the heatmap by significance stars.}
}
\description{
  Helper function to calculate the correlation matrix.
}

