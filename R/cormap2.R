#' @title Draw Correlation Heatmaps
#' @docType package
#' @name coreheat
#' @description Create correlation heatmaps from a numeric matrix. Ensembl Gene ID row names can be converted to Gene Symbols
#'     using, e.g., BioMart. Optionally, data can be clustered and filtered by correlation, tree cutting and/or number
#'     of missing values. Genes of interest can be highlighted in the plot and correlation significance be indicated by
#'     asterisks encoding corresponding P-Values. Plot dimensions and label measures are adjusted automatically by default.
#'     The plot features rely on the heatmap.n2() function in the 'heatmapFlex' package.
#' @author Vidal Fey <vidal.fey@gmail.com>, Henri Sara <henri.sara@gmail.com>
#' Maintainer: Vidal Fey <vidal.fey@gmail.com>
#' @details \tabular{ll}{
#' Package: \tab coreheat\cr
#' Type: \tab Package\cr
#' Initial version: \tab 0.1.0\cr
#' Created: \tab 2016-08-11\cr
#' License: \tab GPL-3\cr
#' The main function to be called by end users is \command{cormap2} which is wrapper performing all necessary steps to create a heatmap.
#' }
#'
#' @keywords package
#' @import Biobase
#' @import heatmapFlex
#' @import convertid
#' @import stats
#' @import graphics
#' @import grDevices
#' @importFrom methods is
NULL
#' Draw correlation maps from large datasets.
#' @description \command{cormap2()} generates pair-wise correlations from an input ExpressionSet object, a \code{data.frame} or a
#'     numerical \code{matrix}. With the default options it also produces a heatmap.
#' @param x (\code{ExpressionSet}, \code{data.frame} or \code{numeric}). A numeric data frame, matrix or an ExpressionSet object.
#' @param cormat (\code{numeric}). A correlation matrix. If this not \code{NULL} then \option{x} is ignored. Defaults to \code{NULL}.
#' @param lab (\code{character}). Optional row/column labels for the heatmap. Defaults to NULL meaning the row names of the input data
#'     are used. Note that the order of the labels must match the order of the row names of the input data!
#' @param convert (\code{logical}). Should an attempt be made to convert IDs provided as row names of the input or in \option{lab}?
#'     Defaults to \code{TRUE}. Conversion will be done using BioMart or an annotation package, depending on \option{biomart}.
#' @param biomart (\code{logical}). Should BioMart (or an annotation package) be used to convert IDs? If \code{TRUE}
#'     the \code{todisp2} function in package \code{convertid} attempts to access the BioMart API to convert ENSG IDs to Gene Symbols
#'     Defaults to \code{FALSE} which will use the traditional \code{AnnotationDbi} Bimap interface.
#' @param cluster_correlations (\code{logical}). Should the correlation matrix be clustered before plotting? Defaults to \code{TRUE}.
#' @param main (\code{character}). The main title of the plot. Defaults to \code{""}.
#' @param postfix (\code{character} of \code{logical}). A plot sub-title. Will be printed below the main title. Defaults to \code{NULL}.
#' @param cex (\code{numeric}). Font size. Defaults to \code{0.5} if \code{autoadj} is \code{FALSE}. See 'Details'.
#' @param na.frac (\code{numeric}). Fraction of missing values allowed per row of the input matrix. Defaults to \code{0.1} which
#'     means LESS than 10 per cent of the values in one row are allowed to be NAs.
#' @param cor.cluster (\code{numeric}). The correlation cluster along the diagonal 'line' in the heatmap that should be
#'     zoomed into. A sliding window of size \code{cor.window} will be moved along the diagonal of the correlation
#'     matrix to find the cluster with the most corelation values meeting \code{core.thr}. Defaults to \code{1}.
#' @param cor.window (\code{numeric}). The size of the sliding window (see \code{cor.cluster}). Defaults to \code{NULL}.
#'     Note that this works only for positive correlations.
#' @param cor.thr (\code{numeric}). Correlation threshold to filter the correlation matrix for plotting. Defaults to \code{NULL} meaning
#'     no filtering. Note that this value will be applied to margin \option{cor.mar} of the values per row.
#' @param cor.mar (\code{numeric}). Margin of the values per row of the correlation matrix the \option{cor.thr} filter needs to
#'     meet. Defaults to \code{0.5} meaning at least 50 per cent of the values in a row need to meet the threshold in order to keep the row.
#' @param cut.thr (\code{numeric}). Threshold at which dendrogram branches are to be cut. Passed on to argument \code{cutHeight} in
#'     \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{NULL} meaning no cutting.
#' @param cut.size (\code{numeric}). Minimum number of objects on a dendrogram branch considered a cluster. Passed on to argument
#'     \code{minSize} in \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{5}.
#' @param autoadj (\code{logical}). Should plot measures be adjusted automatically? Defaults to \code{TRUE}.
#' @param labelheight (\code{numeric} or \code{lcm(numeric)}). Relative or absolute height (using \code{\link[graphics]{lcm}},
#'     see \code{\link[graphics]{layout}}) of the labels. Defaults to \code{0.2} if \code{autoadj} is \code{FALSE}. See 'Details'.
#' @param labelwidth (\code{numeric} or \code{lcm(numeric)}). Relative or absolute width (using \code{\link[graphics]{lcm}},
#'     see \code{\link[graphics]{layout}}) of the labels. Defaults to \code{0.2} if \code{autoadj} is \code{FALSE}. See 'Details'.
#' @param add.sig (\code{logical}). Should significance asterisks be drawn? If \code{TRUE} P-Values for correlation significance
#'     are calculated and encoded as asterisks. See 'Details'.
#' @param genes2highl (\code{character}). Vector of gene symbols (or whatever labels are used) to be highlighted.
#'     If not \code{NULL} will draw a semi-transparent rectangle around the labels and rows or columns in the heatmap
#'     labels.
#' @param order.list (\code{logical}). Should the order of the correlation matrix, i.e. the 'list' of labels be reversed?
#'     Meaningful if the order of input variables should be preserved because \code{\link[graphics]{image}} turns the input
#'     matrix. Defaults to \code{TRUE}.
#' @param doPlot (\code{logical}). Draw the plot? Defaults to \code{TRUE}.
#' @param updateProgress (\code{function}). Function for updating a progress bar in a Shiny web application. This was added here
#'     for the \strong{BioCPR} application.
#' @param verbose (\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.
#' @details
#'     P-Values are calculated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
#'     where r is the correlation coefficient, n is the number of samples with no missing values for each gene (row-wise
#'     \code{ncol(eset)} minus the number of columns that have an NA). P-Values are the calculated using \code{\link[stats]{pt}} and
#'     corrected account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation.
#'     The approach to calculate correlation significance was adopted from \emph{Miles, J., & Banyard, P. (2007)} on
#'     "Calculating the exact significance of a Pearson correlation in MS Excel".
#'
#'     The asterisks encode significance as follows:
#'     \tabular{ll}{
#'     \tab P < 0.05:  *\cr
#'     \tab P < 0.01:  **\cr
#'     \tab P < 0.001: ***\cr
#'     }
#'
#'     The label measures (\code{labelheight}, \code{labelwidth} and \code{cex}) are adjusted automatically by default
#'     with argument \code{autoadj=TRUE} and have default values which are hard coded into the helper function
#'     \code{heatmap.cor}. The values calculated by the helper function \code{plotAdjust} can be overridden by setting
#'     any of those arguments to a valid \code{numeric} or \code{lcm(numeric)} value.
#' @references
#'     Miles, J., & Banyard, P. (2007). \emph{Understanding and using statistics in psychology: A practical introduction.}
#'     Sage Publications Ltd. \url{https://psycnet.apa.org/record/2007-06525-000}.
#' @seealso \code{\link[stats]{pt}}
#' @seealso \code{\link[base]{tcrossprod}}
#' @return Invisibly returns the correlation matrix, though the function is mainly called for its side-effect of producing
#'     a heatmap (if \code{doPlot = TRUE} which is the default).
#' @examples
#' # 1. Generate a random 10x10 matrix with two distinct sets and plot it with
#' # default settings without ID conversion since the IDs are made up:
#' set.seed(1234)
#' mat <- matrix(c(rnorm(100, mean = 1), rnorm(100, mean = -1)), nrow = 20)
#' rownames(mat) <- paste0("gene-", 1:20)
#' colnames(mat) <- paste0(c("A", "B"), rep(1:5, 2))
#' cormap2(mat, convert=FALSE, main="Random matrix")
#'
#' # 2. Use a real-world dataset from TCGA (see README file in inst/extdata directory).
#' # Package 'convertid' is used to convert Ensembl Gene IDs to HGNC Symbols
#' ## Read data and prepare input data frame
#' fl <- system.file("extdata", "PrCaTCGASample.txt", package = "coreheat", mustWork = TRUE)
#' dat0 <- read.delim(fl, stringsAsFactors=FALSE)
#' dat1 <- data.frame(dat0[, grep("TCGA", names(dat0))], row.names=dat0$ensembl_gene_id)
#' cormap2(dat1, main="TCGA data frame + ID conversion")
#'
#' # 3. Use separately supplied IDs with a matrix created from the data frame of the
#' # previous example and highlight genes of interest
#' dat2 <- as.matrix(dat0[, grep("TCGA", names(dat0))])
#' sym <- dat0$hgnc_symbol
#' cormap2(dat1, convert=FALSE, lab=sym, genes2highl=c("GNAS","NCOR1","AR", "ATM"),
#' main="TCGA matrix + custom labels")
#'
#' # 4. Use an ExpressionSet object and add significance asterisks
#' ## For simplicity reasons we create the ExpressionSet from a matrix created
#' ## from the data frame in the second example
#' expr <- Biobase::ExpressionSet(as.matrix(dat1))
#' cormap2(expr, add.sig=TRUE, main="TCGA ExpressionSet object + ID conversion")
#'
#' # More examples can be found in the vignette.
#' @export
cormap2 <- function(x, cormat = NULL, lab = NULL, convert = TRUE, biomart = FALSE, cluster_correlations = TRUE, main = "",
                    postfix = NULL, cex = NULL, na.frac = 0.1, cor.cluster = 1, cor.window = NULL, cor.thr = NULL,
                    cor.mar = 0.5, cut.thr = NULL, cut.size = 5, autoadj = TRUE, labelheight= NULL, labelwidth = NULL,
                    add.sig = FALSE, genes2highl = NULL, order.list = TRUE, doPlot = TRUE, updateProgress = NULL,
                    verbose = FALSE)
{
  if (missing(x) && is.null(cormat)) stop("Need one of 'x' or 'cormat'!")
  list.output <- TRUE
  if (!is.null(cormat)) {
    if (is.matrix(cormat)) {
      if (verbose) message("Input is a correlation matrix.")
      add.sig <- FALSE
    }
    x <- cormat
  }
  if (is.data.frame(x) || is(x, "ExpressionSet") || is.numeric(x)) {
    if (verbose) message("@ Calculating correlation matrix...")
    if (is.function(updateProgress)) {
      updateProgress(detail = "Calculating correlation matrix")
    }
    if (!is.null(lab)) {
      if (verbose) message("    Working with user-provided labels...")
      if (length(lab) != nrow(x))
        stop("Number of labels does not match number of rows in input matrix!")
      rownames(x) <- lab
    }
    cormat <- eset_cor(x, with.pvalues = TRUE, order.list = order.list, verbose = verbose)
  } else if (is.list(x)) {
    if (is.null(cormat)) {
      stop("'x' needs to be an ExpressionSet, a data frame or a numeric matrix!")
    } else {
      if (all(c("cormat", "pvalues") %in% names(cormat)) && !is.null(cormat$cormat) && !is.null(cormat$pvalues)) {
        if (verbose) message("Input is a complete list with correlation matrix and P-Values.")
        if (order.list) {
          cormat$pvalues <- cormat$pvalues[nrow(cormat$pvalues):1, ]
          cormat$cormat <- cormat$cormat[nrow(cormat$cormat):1, ]
        }
      } else {
        stop("If 'cormap' is provided and a list it needs at least non-empty slots 'cormat' and 'pvalues'!")
      }
    }
  }
  if (is.function(updateProgress)) {
    updateProgress(detail = "Clustering and filtering correlation matrix")
  }
  if (cluster_correlations) {
    if (verbose) message("@ Clustering and filtering correlation matrix...")
    cormat <- clust_cormap(cormat, na.frac=na.frac, cor.cluster=cor.cluster, cor.window=cor.window,
                           cor.thr=cor.thr, cor.mar=cor.mar, cut.thr=cut.thr, cut.size=cut.size,
                           list.output = list.output, verbose = verbose)
    if (!is.null(cut.thr) || !is.null(cor.thr)) {
      if (is.logical(postfix) && !postfix) {
        postfix <- NULL
      } else {
        if (is.list(cormat)) {
          cmp <- cormat$cormat
        } else {
          cmp <- cormat
        }
        postfix <- paste(postfix, paste0("(", nrow(cmp), " rows after filtering)"))
      }
    }
  }
  if (list.output) {
    cormat.p <- cormat$pvalues[nrow(cormat$cormat):1, ]
    corm <- cormat <- cormat$cormat[nrow(cormat$cormat):1, ]
  } else {
    cormat.p <- NULL
    corm <- cormat <- cormat[nrow(cormat):1, ]
  }
  if (verbose) message("@ Formatting output matrix...")
  if (is.function(updateProgress)) {
    updateProgress(detail = "Formatting matrix")
  }
  if (verbose) message("  Checking input IDs:")
  if (convert) {
    if (verbose) message("    Attempting to convert input matrix row names...")
    xl <- xlab <- convertid::todisp2(rownames(cormat), biomart = biomart, verbose = verbose)
    yl <- ylab <- rev(xl)
  } else {
    if (verbose) message("    Using input matrix row names...")
    xl <- xlab <- rownames(cormat)
    yl <- ylab <- colnames(cormat)
  }

  if (any(duplicated(xl))) {
    if (verbose) message("  Fixing duplicated row names...")
    if (is.function(updateProgress)) {
      updateProgress(detail = "Fixing duplicated row names")
    }
    xl[xl %in% xl[duplicated(xl)]] <- paste(xl[xl %in% xl[duplicated(xl)]], rownames(cormat)[xl %in% xl[duplicated(xl)]], sep=".")
  }
  if (any(duplicated(yl))) {
    if (verbose) message("  Fixing duplicated column names...")
    if (is.function(updateProgress)) {
      updateProgress(detail = "Fixing duplicated row names")
    }
    yl[yl %in% yl[duplicated(yl)]] <- paste(yl[yl %in% yl[duplicated(yl)]], rownames(cormat)[yl %in% yl[duplicated(yl)]], sep=".")
  }
  rownames(corm) <- xl
  colnames(corm) <- yl
  if (doPlot) {
    if (verbose) message("@ Generating heatmap...")
    if (is.function(updateProgress)) {
      updateProgress(detail = "Generating heatmap")
    }
    heatmap.cor(cormat, order.list=order.list, main=main, main_postfix=postfix, x.labels = xlab, y.labels = ylab,
                autoadj = autoadj, labelheight=labelheight, labelwidth=labelwidth, cex=cex, add.sig=add.sig,
                pv=cormat.p, genes2highl=genes2highl, verbose = verbose)
  }
  if (list.output) {
    return(invisible(list(cormat=cormat, pvalues=cormat.p)))
  } else {
    return(invisible(corm))
  }
}

#' Helper function to calculate the correlation matrix.
#' @param x (\code{ExpressionSet}, \code{data.frame} or \code{numeric}). A numeric data frame, matrix or an ExpressionSet object.
#' @param with.pvalues (\code{logical}). Should P-Values be calculated for the correlations? If \code{TRUE} P-Values will be
#'     depicted in the heatmap by significance asterisks. See 'Details'. Defaults to \code{FALSE}.
#' @param order.list (\code{logical}). Is the input matrix column order reversed? Only applicable if input is correlation matrix.
#'     Defaults to \code{TRUE}.
#' @param verbose (\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.
#' @details
#'     P-Values are calculated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
#'     where r is the correlation coefficient, n is the number of samples with no missing values for each gene (row-wise
#'     \code{ncol(eset)} minus the number of columns that have an NA). P-Values are then calculated using \code{\link[stats]{pt}} and
#'     corrected account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation.
#'     The approach to calculate correlation significance was adopted from \emph{Miles, J., & Banyard, P. (2007)} on
#'     "Calculating the exact significance of a Pearson correlation in MS Excel".
#' @references
#'     Miles, J., & Banyard, P. (2007). \emph{Understanding and using statistics in psychology: A practical introduction.}
#'     Sage Publications Ltd. \url{https://psycnet.apa.org/record/2007-06525-000}.
#' @return A correlation \code{matrix} or a \code{list} with three slots: the correlation matrix, the number of samples with no
#'     missing value for each gene and the P-Values matrix.
eset_cor <- function(x, with.pvalues=TRUE, order.list = TRUE, verbose = FALSE) {
  if (is(x, "ExpressionSet")) {
    expr <- exprs(x)
  } else {
    expr <- as.matrix(x)
  }
  if (!is.numeric(expr)) {
    stop("'x' needs to be a numeric matrix, data frame or an ExpressionSet object!")
  }

  if (order.list) {
    oe <- expr[nrow(expr):1, ]
  } else {
    oe <- expr
  }
  if (all(diag(oe)==1) && identical(colSums(oe), rowSums(oe))) {
    if (verbose) message("  Input is a correlation matrix. Doing nothing...")
    if (order.list) {
      if (verbose) message("Returning ordered row names")
      return(expr[nrow(expr):1, ])
    } else {
      if (verbose) message("Returning input unsorted")
      return(expr)
    }
  }
  # remove rows and columns with all values missing
  expr <- expr[apply(expr, 1, function(x) !all(is.na(x))), ]
  expr <- expr[, apply(expr, 2, function(x) !all(is.na(x)))]
  if (verbose) message("  Computing correlations...")
  cormat <- cor(t(expr), method="pearson", use="pairwise.complete")

  if (with.pvalues) {
    if (verbose) message("  Calculating P-Values...")
    if (verbose) message("   Getting NAs...")
    nas <- !is.na(expr) # get a logical matrix with all non-NAs as 'TRUE' to effectively have a matrix of '0's and '1's for the next step
    if (verbose) message("   Getting counts...")
    counts <- nas %*% t(nas) # get number of samples with no missing value for each gene, i.e., ncol(expr) minus the cols that have an NA for that gene in a matrix with the same dimensions as cormat by calculating the matrix cross-product (tcrossprod)
    df <- counts - 2 # use the above calculated matrix to obtain the degrees of freedom: N-2 since we always have a 2-dimensional correlation matrix
    t <- sqrt(df) * cormat / sqrt(1 - cormat^2) # calculate the t-test value (from "Understanding and Using Statistics in Psychology --- A Practical Introduction" on "Calculating the exact significance of a Pearson correlation in MS Excel",
    # found at http://stats.stackexchange.com/questions/120199/calculate-p-value-for-the-correlation-coefficient)
    p <- pt(t, df) # convert the test value to it's P-Value
    if (verbose) message("   Getting P-Values...")
    pvalues <- 2 * pmin(p, 1-p) # account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation
    list(cormat=cormat, pvalues=pvalues)
  } else {
    cormat
  }
}

#' Cluster a correlation matrix and return the sorted matrix for plotting.
#' @description Helper function to cluster the correlation matrix and return the sorted matrix for plotting.
#' @param cormat (\code{numeric} or \code{list}). The correlation matrix or a list containing matrices of correlation values
#'     and P-Values as generated by \code{eset_cor}. If a list then the correlation matrix is expected to be in slot
#'     'cor' and the matrix of P-Values in slot 'pvalues'.
#' @param na.frac (\code{numeric}). Fraction of missing values allowed per row of the input matrix. Defaults to \code{0.1} which
#'     means LESS than 10per centof the values in one row are allowed to be NAs.
#' @param distfn (\code{function}). Function to calculate the dissimilarity matrix for clustering. Defaults to \code{function(cm) (1-cm)}.
#' @param method (\code{character}). The agglomeration method used for clustering. See help for \code{\link[stats]{hclust}}.
#'     Defaults to "complete".
#' @param cor.cluster (\code{numeric}). The correlation cluster along the diagonal 'line' in the heatmap that should be
#'     zoomed into. A sliding window of size \code{cor.window} will be moved along the diagonal of the correlation
#'     matrix to find the cluster with the most corelation values meeting \code{core.thr}. Defaults to \code{1}.
#' @param cor.window (\code{numeric}). The size of the sliding window (see \code{cor.cluster}). Defaults to \code{NULL}.
#'     Note that this works only for positive correlations.
#' @param cor.thr (\code{numeric}). Correlation threshold to filter the correlation matrix for plotting. Defaults to \code{NULL} meaning
#'     no filtering. See also \option{cor.mar}. This value is sign-sensitive: a negative threshold will retain rows and columns of
#'     the correlation matrix with correlation values between -1 and the threshold, a positive value will retain rows and columns
#'     with values between the threshold and 1. Zero (0) is treated as positive.
#' @param cor.mar (\code{numeric}). Margin of the values per row of the correlation matrix the \option{cor.thr} filter needs to
#'     meet. Defaults to \code{0.5} meaning at least 50 per cent of the values in a row need to meet the threshold in order to keep the row.
#' @param cut.thr (\code{numeric}). Threshold at which dendrogram branches are to be cut. Passed on to argument \code{cutHeight} in
#'     \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{NULL} meaning no cutting.
#' @param cut.size (\code{numeric}). Minimum number of objects on a dendrogram branch considered a cluster. Passed on to argument
#'     \code{minSize} in \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{5}.
#' @param list.output (\code{logical}). Should the output be a list of different object created in the call? Depends also on input
#'     type. If \code{FALSE} only the (filtered) correlation matrix is returned.
#' @param verbose (\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.
#' @return A correlation \code{matrix} or a \code{list} or the matrix and other values needed for plotting.
clust_cormap <- function(cormat, na.frac=0.1, distfn=function(cm) (1-cm), method="complete",
                         cor.cluster = 1, cor.window = NULL, cor.thr=0.8, cor.mar=0.05, cut.thr=0.9, cut.size=5,
                         list.output=FALSE, verbose = FALSE) {
  ori.pv <- FALSE
  ordrd <- FALSE
  if (is.list(cormat)) {
    if (all(c("cormat", "pvalues", "hclust") %in% names(cormat)) && !is.null(cormat$cormat) && !is.null(cormat$pvalues) && !is.null(cormat$hclust)) {
      if (verbose) message("Input is a complete list with correlation matrix, dendrogram and P-Values.")
      ori.pv <- TRUE
      pval <- cormat$pvalues
      cormat <- cormat$cormat
      ordrd <- TRUE
    } else if (all(c("cormat", "hclust") %in% names(cormat)) && !is.null(cormat$cormat) && !is.null(cormat$hclust)) {
      if (verbose) message("Input is a list with correlation matrix and dendrogram.")
      cormat <- cormat$cormat
      ordrd <- TRUE
    } else if (all(c("cormat", "pvalues") %in% names(cormat)) && !is.null(cormat$cormat) && !is.null(cormat$pvalues)) {
      if (verbose) message("Input is a list with correlation matrix and pvalues. Assuming non-clustered matrix.")
      ori.pv <- TRUE
      pval <- cormat$pvalues
      cormat <- cormat$cormat
    } else {
      stop("If 'cormap' is a list it needs at least non-empty slots 'cormat' and/or 'hclust' and/or 'pvalues'!")
    }
  }

  if (!ordrd) {
    if (verbose) message("Filtering NAs...")
    # filtering: keep only rows and columns with less than the given fraction of missing values; defaults to 10 per cent
    filt <- apply(cormat, 1, function(x) sum(is.na(x)) < na.frac*length(x))
    cormat2 <- cormat[filt, filt]
    if (ori.pv) {
      pval <- pval[filt, filt]
    }
  } else {
    if (verbose) message("Input is correlation matrix. No filtering is done...")
    cormat2 <- cormat
  }
  if (verbose) message("Getting distance matrix and hierarchical clustering...")
  dissimilarity <- distfn(cormat2)
  distance <- as.dist(dissimilarity)
  cl <- hclust(distance, method=method)

  # cut tree at given dissimilarity threshold
  if (!is.null(cut.thr)) {
    if (verbose) message("Cutting the dendrogram at ", cut.thr, "...")
    clust <- WGCNA::cutreeStatic(cl, cutHeight = cut.thr, minSize = cut.size)
    keep <- (clust==1)
    if (verbose) message("  Isolated tree section with ", length(which(keep)), " genes")
    cormat2 <- cormat2[keep, keep]
    if (ori.pv) {
      pval <- pval[keep, keep]
    }
    if (verbose) message("  Recalculating distance matrix and hierarchical clustering...")
    dissimilarity <- distfn(cormat2)
    distance <- as.dist(dissimilarity)
    cl <- hclust(distance, method=method)
  }

  if (verbose) message("Ordering according to hierarchical clustering...")
  cl_ord <- cl$order
  cormat3 <- cormat2[cl_ord, cl_ord]
  if (ori.pv) {
    pval <- pval[cl_ord, cl_ord]
  }

  if (!is.null(cor.thr)) {
    if (verbose) message("Filtering by correlation...")
    if (!is.null(cor.window)) {
      if (cor.window < 2) {
        if (verbose) message("  Correlation cluster window too small. Setting to 2.")
        cor.window <- 2
      }
      if (verbose) message("  Filtering by correlation threshold in a window of ", cor.window)
      if (sign(cor.thr) == -1) {
        stop("The sliding window approach does not work for negative correlations.")
      }
      steps <- nrow(cormat3) - cor.window + 1
      ki <- 1
      cfi <- cfi1 <- 1
      cf <- list()
      for (i in 1:steps) {
        cw <- cormat3[i:(i+cor.window-1), i:(i+cor.window-1)]
        nm <- sum(cw > cor.thr)
        if (nm == cor.window^2) {
          if (ki == cor.cluster) {
            if (verbose) message(paste0("Cycle ", i, "  -> Found ", ki, ". cluster of correlating values (cor. threshold ", cor.thr, ")\n   -> Expanding window..."))
          }
          k1 <- list()
          for (k in 1:(nrow(cormat3)-i-cor.window+1)) {
            cw1 <- cormat3[i:(i+cor.window-1+k), i:(i+cor.window-1+k)]
            nm1 <- sum(cw1 > cor.thr)
            if (nm1 < (cor.window+k)^2) {
              k1[[k]] <- cor.window+k-1
              break
            } else {
              k1[[k]] <- NA_integer_
            }
          }
          k1 <- na.omit(unlist(k1))
          if (length(k1) == 0 && ki != cor.cluster) {
            if (verbose) message(paste0("  ! No ", cor.cluster, ". region found. Returning latest cluster!"))
            cf[[i]] <- i:(i + k + 1)
            break
          } else if (length(k1) == 0 && ki == cor.cluster) {
            cf[[i]] <- i:(i + k + 1)
            break
          } else {
            cf[[i]] <- i:(i + k1[1] - 1)
          }
          if (ki == cor.cluster) {
            if (verbose) message("   -> done.")
            break
          } else if (k1 == cor.window) {
            ki <- ki + 1
          }
        } else {
          cf[[i]] <- NA_integer_
        }
      }
      cf <- cf[[length(cf)]]
    } else {
      if (verbose) message("  Filtering by correlation threshold with a margin of ", cor.mar*100, "%")
      cf <- apply(cormat3, 1, function(x) {
        cond <- sum(x >= cor.thr)
        if (sign(cor.thr) == -1) {
          cond <- sum(x <= cor.thr)
        }
        cond >= cor.mar*length(x)
      })
      cf <- which(cf)
    }
    if (length(cf) < 2) {
      if (verbose) message("NOTE: Filtering left no values to plot. Returning unfiltered matrix.")
    } else {
      if (verbose) message("Filtered ", nrow(cormat3)-length(cf), paste0(" rows not meeting correlation threshold (", cor.thr, " in ", cor.mar*100, "% of the columns)"))
      cormat3 <- cormat3[cf, cf]
      if (ori.pv) {
        pval <- pval[cf, cf]
      }
    }
  }

  if (list.output) {
    if (ori.pv) {
      list(cormat=cormat3, hclust=cl, pvalues=pval)
    } else {
      list(cormat=cormat3, hclust=cl)
    }
  } else {
    cormat3
  }
}

#' Helper function to draw the heatmap. Wrapper around 'heatmapFlex::heatmap.n2'.
#' @noRd
heatmap.cor <- function(cormat, order.list = TRUE, main="", main_postfix="", x.labels = NULL, y.labels = NULL,
                        autoadj = TRUE, labelheight=NULL, labelwidth=NULL, cex=NULL, add.sig=FALSE, pv=NULL,
                        genes2highl=NULL, verbose = FALSE) {
  main.map <- paste(main, main_postfix, sep="")
  col <- colorRampPalette(c("blue", "white", "red"))(50)
  if(order.list){
    if (verbose) message("  Ordering Y-axis...")
    cormat <- cormat[nrow(cormat):1, ]
    if (add.sig) pv <- pv[nrow(pv):1, ]
    x.labels <- rev(x.labels)
  }
  if (autoadj) {
    adj.l <- plotAdjust(cormat, cormap = TRUE)
    if (!is.null(labelheight)) adj.l$labelheight <- labelheight
    if (!is.null(labelwidth)) adj.l$labelwidth <- labelwidth
    if (!is.null(cex)) adj.l$r.cex <- adj.l$c.cex <- cex
  } else {
    if (is.null(labelheight)) labelheight <- 0.2
    if (is.null(labelwidth)) labelwidth <- 0.2
    if (is.null(cex)) cex <- 0.5
    adj.l <- list(labelwidth = labelwidth, labelheight = labelheight, r.cex = cex, c.cex = cex)
  }

  heatmapFlex::heatmap.n2(cormat, order_list = order.list, reorder=c(FALSE, FALSE), labRow = x.labels,
                          labCol = y.labels, labelheight=adj.l$labelheight, labelwidth=adj.l$labelwidth,
                          col=col, main=main.map, r.cex=adj.l$r.cex, c.cex=adj.l$r.cex, add.sig=add.sig,
                          pv=pv, genes2highl=genes2highl)
}

#' Helper function to extract leaves from a dendrogram. This is adapted from stats:::cut.dendrogram.
#' @noRd
dendro2leaves <- function(x) {
  leaves <- list()
  getNodes <- function(subtree) {
    if (!is.leaf(subtree)) {
      if (!(K <- length(subtree)))
        stop("non-leaf subtree of length 0")
      for (k in 1L:K) {
        getNodes(subtree[[k]])
      }
    } else {
      leaves <<- append(leaves, subtree)
    }
  }
  getNodes(x)
  as.numeric(leaves)
}

#' Helper function for signal metric for a correlation p-value matrix: -log10(median(pval))
#' @noRd
signal_metric <- function(x) {
  -log10(median(pmax(c(x, recursive=TRUE), 1e-38), na.rm=TRUE))
}

#' Automatically split clusters based on noise level and hierarchy
#' @description \command{cormat_filt} splits (cuts) the dendrogram at a given threshold dividing it into larger or
#'     smaller "sub-clusters". Correlation P-Values (see \code{\link[coreheat]{eset_cor}}) are converted to represent
#'     significance as a sub-cluster-wise signal metric used for filtering. Optionally, up to 3 plots are produced,
#'     the third one being a filtered heatmap based on significance and three height cutting.
#'
#' @param x (\code{ExpressionSet}, \code{data.frame} or \code{numeric}). A numeric data frame, matrix or an ExpressionSet object.
#' @param na.frac (\code{numeric}). Fraction of missing values allowed per row of the input matrix. Defaults to \code{0.1} which
#'     means LESS than 10 per cent of the values in one row are allowed to be NAs.
#' @param method (\code{character}). The agglomeration method used for clustering. See help for \code{\link[stats]{hclust}}.
#'     Defaults to "ward.D".
#' @param do.abs (\code{logical}). Should the distances for clustering be calculated based on the absolute correlation values?
#'     In other words, should the sign of the correlation be ignored in favor of its strength?
#' @param main (\code{character}). The main title of the plot. Defaults to \code{""}.
#' @param postfix (\code{character} of \code{logical}). A plot sub-title. Will be printed below the main title. Defaults to \code{NULL}.
#' @param p.thr (\code{numeric}). P-Value threshold for filtering sub-clusterd with significant correlations. Defaults to \code{0.01}.
#' @param cex (\code{numeric}). Font size for the heatmap of the unfiltered correlation matrix. Defaults to \code{0.2}.
#' @param cex.clust (\code{numeric}). Font size for the dendrogram plot of the unfiltered correlation matrix clusters.
#'     Defaults to \code{cex}.
#' @param cex.filt (\code{numeric}). Font size for the heatmap of the filtered correlation matrix. Defaults to \code{cex}.
#' @param cut.thr (\code{numeric}). Threshold at which dendrogram branches are to be cut. Passed on to argument \code{h} in
#'     \code{\link[stats]{cut.dendrogram}}. Defaults to \code{NULL} meaning no cutting.
#' @param cor.thr (\code{numeric}). Correlation threshold to filter the correlation matrix for plotting. Defaults to \code{NULL} meaning
#'     no filtering. Note that this value will be applied to margin \option{cor.mar} of the values per row.
#' @param cor.cluster (\code{numeric}). The correlation cluster along the diagonal 'line' in the heatmap that should be
#'     zoomed into. A sliding window of size \code{cor.window} will be moved along the diagonal of the correlation
#'     matrix to find the cluster with the most corelation values meeting \code{core.thr}. Defaults to \code{1}.
#' @param cor.window (\code{numeric}). The size of the sliding window (see \code{cor.cluster}). Defaults to \code{NULL}.
#'     Note that this works only for positive correlations.
#' @param do.plots (\code{character}). The plots to be produced. A character vector containing one or more of "dend"
#'     to produce the dendrogram plot, "full.heat" to produce the heatmap of the unfiltered correlation matrix, and
#'     "filt.heat" to produce the heatmap of the filtered correlation matrix. Defaults to all three plots.
#' @param genes2highl (\code{character}). Vector of gene symbols (or whatever labels are used) to be highlighted.
#'     If not \code{NULL} will draw a semi-transparent rectangle around the labels and rows or columns in the heatmap
#'     labels.
#' @param order.list (\code{logical}). Should the order of the correlation matrix, i.e. the 'list' of labels be reversed?
#'     Meaningful if the order of input variables should be preserved because \code{\link[graphics]{image}} turns the input
#'     matrix. Defaults to \code{TRUE}.
#' @param convert (\code{logical}). Should an attempt be made to convert IDs provided as row names of the input or in \option{lab}?
#'     Defaults to \code{TRUE}. Conversion will be done using BioMart or an annotation package, depending on \option{biomart}.
#' @param biomart (\code{logical}). Should BioMart (or an annotation package) be used to convert IDs? If \code{TRUE}
#'     the \code{todisp2} function in package \code{convertid} attempts to access the BioMart API to convert ENSG IDs to Gene Symbols
#'     Defaults to \code{FALSE} which will use the traditional \code{AnnotationDbi} Bimap interface.
#' @param add.sig (\code{logical}). Should significance asterisks be drawn? If \code{TRUE} P-Values for correlation significance
#'     are calculated and encoded as asterisks. See 'Details'.
#' @param verbose (\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.
#' @details P-Values are calculated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
#'     where r is the correlation coefficient, n is the number of samples with no missing values for each gene (row-wise
#'     \code{ncol(eset)} minus the number of columns that have an NA). P-Values are then calculated using \code{\link[stats]{pt}} and
#'     corrected account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation.
#'     The approach to calculate correlation significance was adopted from \emph{Miles, J., & Banyard, P. (2007)} on
#'     "Calculating the exact significance of a Pearson correlation in MS Excel".
#'
#'     To obtain a suitable metric for isolating significant sub-clusters, P-Values are represented as \eqn{-log10(median(pval))}
#'     where \code{pval} is the \emph{median of the parallel maximum of all P-Values belonging to the sub-cluster and
#'     \code{1e-38}} to avoid values of zero (0).
#' @return A \code{list}. If the dendrogram is being cut, i.e., \code{cut.thr} is not \code{NULL}, a list of
#' \tabular{ll}{
#' \tab clusters: the list of cluster labels from \code{lower} component of the \code{cut.dendrogram} output which
#'     is list with the branches obtained from cutting the tree\cr
#' \tab filt: the index of the cluster labels passing the signal metrics threshold\cr
#' \tab filt_cluster: the list of the filtered cluster labels\cr
#' \tab h: the cut threshold\cr
#' \tab p.thr: the P-Value threshold for filtering sub-clusters\cr
#' \tab metric: the signal metrics for all sub-clusters\cr
#' \tab cormat: the clustered (ordered) correlation matrix\cr
#' \tab hclust: a list of hierarchical clustering metrics (output of \code{\link[stats]{hclust}})\cr
#' \tab pvalues: the correlation P-Value matrix\cr
#'}
#'
#' If no tree cutting is applied, a list of
#' \tabular{ll}{
#' \tab cormat: the clustered (ordered) correlation matrix\cr
#' \tab hclust: a list of hierarchical clustering metrics (output of \code{\link[stats]{hclust}})\cr
#' \tab pvalues: the correlation P-Value matrix\cr
#' }
#'
#' @export
cormap_filt <- function(x, na.frac=0.1, method="ward.D", do.abs=TRUE, main="correlation map",
                        postfix = NULL, p.thr=0.01, cex=0.2, cex.clust=cex, cex.filt=cex, cut.thr=NULL,
                        cor.thr = NULL, cor.cluster = 1, cor.window = NULL, do.plots=c("dend", "full.heat", "filt.heat"),
                        genes2highl=NULL, order.list=TRUE, convert = TRUE, biomart = FALSE, add.sig = FALSE,
                        verbose = FALSE) {
  if (dev.interactive()) {
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar["mfrow"]))
    par(mfrow=c(1,1))
  }

  if (verbose) message("@ Computing and clustering correlations...")
  distfn <- if (do.abs) function(cm) { 1-abs(cm) } else function(cm) { 1-cm }
  cormat <- eset_cor(x, with.pvalues = TRUE, order.list = order.list, verbose = verbose)
  cormat.cl <- clust_cormap(cormat, cor.window = NULL, cor.thr=NULL, cut.thr = NULL, na.frac=na.frac, distfn=distfn,
                            method=method, verbose = verbose, list.output=TRUE)
  if (is.logical(postfix) && !postfix) {
    postfix <- NULL
  } else if (is.null(postfix)) {
    postfix <- if (do.abs) " Dissimilarity = 1 - abs(Correlation)" else "Dissimilarity = 1 - Correlation"
  }

  if ("dend" %in% do.plots) {
    if (verbose) message(">> Plotting dendrogram...")
    ccl <- cormat.cl
    if (convert) {
      ccl$hclust$labels <- convertid::todisp2(ccl$hclust$labels, biomart = biomart, verbose = verbose)
    }
    plot(ccl$hclust, cex=cex.clust)
    if (!is.null(cut.thr)) abline(h=cut.thr)
  }
  if ("full.heat" %in% do.plots) {
    if (verbose) message(">> Plotting full correlation map...")
    if (order.list) {
      cormat.cl$pvalues <- cormat.cl$pvalues[nrow(cormat.cl$pvalues):1, ]
      cormat.cl$cormat <- cormat.cl$cormat[nrow(cormat.cl$cormat):1, ]
    }
    xl <- yl <- NA
    if (convert) {
      xl <- convertid::todisp2(rownames(cormat.cl$cormat), biomart = biomart, verbose = verbose)
      yl <- convertid::todisp2(colnames(cormat.cl$cormat), biomart = biomart, verbose = verbose)
    }
    heatmap.cor(cormat.cl$cormat, order.list = order.list, main=main, main_postfix=postfix, x.labels = xl,
                      y.labels = yl, cex=cex, genes2highl=genes2highl, add.sig = add.sig, pv = cormat.cl$pvalues,
                      verbose = verbose)
  }

  if (!is.null(cut.thr)) {
    # compare the two above to find interesting clusters
    cut_cluster <- cut(as.dendrogram(cormat.cl$hclust), h=cut.thr)
    cluster_labels <- lapply(cut_cluster$lower, function(x) cormat.cl$hclust$labels[dendro2leaves(x)])

    cluster_signal_metrics <- sapply(cluster_labels, function(ids) {
      signal_metric(cormat.cl$pvalues[ids, ids])
    })
    # threshold is converted to -log10 to match signal metric above
    p.post <- p.thr <- -log10(p.thr)
    filt <- cluster_signal_metrics > p.thr
    ids <- lunion(cluster_labels[filt])
    cm.filt <- list(cormat=cormat.cl$cormat[ids, ids], pvalues=cormat.cl$pvalues[ids, ids])
    # filtered cormap
    cormat.cl <- clust_cormap(cm.filt, na.frac=1, distfn=distfn, method=method, cor.thr = cor.thr, cor.cluster = cor.cluster,
                              cor.window = cor.window, cut.thr = NULL, list.output=TRUE, verbose = verbose)
    if (order.list) {
      cormat.cl$pvalues <- cormat.cl$pvalues[nrow(cormat.cl$pvalues):1, ]
      cormat.cl$cormat <- cormat.cl$cormat[nrow(cormat.cl$cormat):1, ]
    }
    if ("filt.heat" %in% do.plots) {
      if (verbose) message(">> Plotting filtered correlation map... (Tree height ", cut.thr, ", P-Value ", p.post, ")")
      if (is.logical(postfix) && !postfix) {
        postfix <- NULL
      } else {
        postfix <- paste0("(Cut ", cut.thr, ", PVal ", round(p.post, 2), ", ", length(ids), " rows after filtering)")
      }
      xl <- yl <- NA
      if (convert) {
        xl <- convertid::todisp2(rownames(cormat.cl$cormat), biomart = biomart, verbose = verbose)
        yl <- convertid::todisp2(colnames(cormat.cl$cormat), biomart = biomart, verbose = verbose)
      }
      heatmap.cor(cormat.cl$cormat, order.list = order.list, main=main, main_postfix=postfix, x.labels = xl,
                        y.labels = yl, cex=cex.filt, add.sig = add.sig, pv = cormat.cl$pvalues,
                        genes2highl=genes2highl, verbose = verbose)
    }
    invisible(append(list(clusters=cluster_labels, filt=filt, filt_clusters=cluster_labels[filt], h=cut.thr, p.thr=p.thr,
                          metric=cluster_signal_metrics), cormat.cl))
  } else {
    invisible(cormat.cl)
  }
}
