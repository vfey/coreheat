#' @title Draw Correlation Heatmaps
#' @docType package
#' @name coreheat
#' @description Create correlation heatmaps from a numeric matrix. Optionally, data can be annotated using, e.g., BioMart,
#'     clustered and filtered by correlation and/or tree cutting.
#' @author Vidal Fey <vidal.fey@gmail.com>, Henri Sara <henri.sara@gmail.com>
#' Maintainer: Vidal Fey <vidal.fey@gmail.com>
#' @details \tabular{ll}{
#' Package: \tab coreheat\cr
#' Type: \tab Package\cr
#' Initial version: \tab 0.1-0\cr
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
#'     are used.
#' @param convert (\code{logical}). Should an attempt be made to convert IDs provided as row names of the input or in \option{lab}?
#'     Defaults to \code{TRUE}. Conversion will be done using BioMart or an annotation package, depending on \option{biomart}.
#' @param biomart (\code{logical}). Should BioMart (or an annotation package) be used to convert IDs? If \code{TRUE} (the default)
#'     the \code{todisp2} function in package \code{convertid} attempts to access the BioMart API to convert ENSG IDs to Gene Symbols
#'     but if that fails the user can set this option to \code{FALSE} to use the traditional \code{AnnotationDbi} Bimap interface.
#' @param cluster_correlations (\code{logical}). Should the correlation matrix be clustered before plotting? Defaults to \code{TRUE}.
#' @param main (\code{character}). The main title of the plot. Defaults to \code{""}.
#' @param cex (\code{numeric}). Font size. Defaults to \code{0.7}.
#' @param na.frac (\code{numeric}). Fraction of missing values allowed per row of the input matrix. Defaults to \code{0.1} which
#'     means LESS than 10 per cent of the values in one row are allowed to be NAs.
#' @param cor.thr (\code{numeric}). Correlation threshold to filter the correlation matrix for plotting. Defaults to \code{NULL} meaning
#'     no filtering. Note that this value will be applied to margin \option{cor.mar} of the values per row.
#' @param cor.mar (\code{numeric}). Margin of the values per row of the correlation matrix the \option{cor.thr} filter needs to
#'     meet. Defaults to \code{0.5} meaning at least 50 per cent of the values in a row need to meet the threshold in order to keep the row.
#' @param cut.thr (\code{numeric}). Threshold at which dendrogram branches are to be cut. Passed on to argument \code{cutHeight} in
#'     \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{NULL} meaning no cutting.
#' @param cut.size (\code{numeric}). Minimum number of objects on a dendrogram branch considered a cluster. Passed on to argument
#'     \code{minSize} in \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{1}.
#' @param labelheight (\code{numeric} or \code{lcm(numeric)}). Relative or absolute height (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the labels.
#' @param labelwidth (\code{numeric} or \code{lcm(numeric)}). Relative or absolute width (using \code{\link[graphics]{lcm}}, see \code{\link[graphics]{layout}}) of the labels.
#' @param add.sig (\code{logical}). Should significance asterisks be drawn? If \code{TRUE} P-Values for correlation significance
#'     are calculated and encoded as asterisks. See 'Details'.
#' @param genes2highl (\code{character}). Vector of gene symbols (or whatever labels are used) to be highlighted.
#'     If not \code{NULL} will draw a semi-transparent rectangle around the labels and rows or columns in the heatmap
#'     labels.
#' @param order.list (\code{logical}). Should the order of the correlation matrix, i.e. the 'list' of labels be reversed?
#'     Meaningful if the order of input variables should be preserved because \code{\link[graphics]{image}} turns the input matrix.
#' @param doPlot (\code{logical}). Draw the plot? Defaults to \code{TRUE}.
#' @param updateProgress (\code{function}). Function for updating a progress bar in a Shiny web application. This was added here
#'     for the \strong{BioCPR} application.
#' @details
#'     P-Values are calcluated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
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
#' @references
#'     Miles, J., & Banyard, P. (2007). \emph{Understanding and using statistics in psychology: A practical introduction.}
#'     Sage Publications Ltd. \url{https://psycnet.apa.org/record/2007-06525-000}.
#' @seealso \code{\link[stats]{pt}}
#' @seealso \code{\link[base]{tcrossprod}}
#' @return Invisibly returns the correlation matrix, though the function is mainly called for its side-effect of producing
#'     a heatmap (if \code{doPlot = TRUE} which is the default).
#' @examples
#' # Generate a random 10x10 matrix with two distinct sets and plot it with
#' # default settings but without biomart since the IDs are made up:
#' set.seed(1234)
#' mat <- matrix(c(rnorm(100, mean = 1), rnorm(100, mean = -1)), nrow = 20)
#' rownames(mat) <- paste0("gene-", 1:20)
#' colnames(mat) <- paste0(c("A", "B"), rep(1:5, 2))
#' cormap2(mat, biomart=FALSE)
#'
#' # Use a real-world dataset from TCGA (see README file in inst/extdata directory).
#' # BioMart is used to convert Ensembl Gene IDs to HGNC Symbols
#' ## Read data and prepare input data frame
#' fl <- system.file("extdata", "PrCaTCGASample.txt", package = "coreheat", mustWork = TRUE)
#' dat <- read.delim(fl)
#' dat <- data.frame(dat[, grep("TCGA", names(dat))], row.names=dat$ensembl_gene_id)
#' ## Plot correlation map
#' cormap2(dat, biomart=FALSE)
#'
#' @export
cormap2 <- function(x, cormat = NULL, lab = NULL, convert = TRUE, biomart = TRUE, cluster_correlations = TRUE, main = "",
                    cex = 0.7, na.frac = 0.1, cor.thr = NULL, cor.mar = 0.5, cut.thr = NULL, cut.size = 1, labelheight= 0.2,
                    labelwidth = 0.2, add.sig = FALSE, genes2highl = NULL, order.list = TRUE, doPlot = TRUE,
                    updateProgress = NULL) {
	cat("@ Calculating correlation matrix...\n")
	if (is.function(updateProgress)) {
		updateProgress(detail = "Calculating correlation matrix")
	}
  if (!is.null(cormat)) {
    x <- cormat
  }
	if (add.sig) {
		cormat <- eset_cor(x, with.pvalues = TRUE)
		list.output <- TRUE
	} else {
		cormat <- eset_cor(x)
		list.output <- FALSE
	}
	if (is.function(updateProgress)) {
		updateProgress(detail = "Clustering and filtering correlation matrix")
	}
	if (cluster_correlations) {
		cat("@ Clustering and filtering correlation matrix...\n")
		cormat <- clust_cormap(cormat, na.frac=na.frac, cor.thr=cor.thr, cor.mar=cor.mar, cut.thr=cut.thr, cut.size=cut.size, list.output = list.output)
	}
	if (list.output) {
		cormat.p <- cormat$pval[nrow(cormat$cormat):1, ]
		corm <- cormat <- cormat$cormat[nrow(cormat$cormat):1, ]
	} else {
		cormat.p <- NULL
		corm <- cormat <- cormat[nrow(cormat):1, ]
	}
	if (is.null(lab)) {
		lab.df <- data.frame(ID=rownames(cormat), stringsAsFactors=F)
		rownames(lab.df) <- rownames(cormat)
	} else if (is.character(lab)) {
	  lab.df <- data.frame(ID=lab, stringsAsFactors=F)
	  rownames(lab.df) <- lab
	}
	cat("@ Formatting output matrix...\n")
	if (is.function(updateProgress)) {
		updateProgress(detail = "Formatting matrix")
	}
	cat("  Row names:\n")
	if (convert) {
	  if (!is.null(lab)) {
	    xl <- xlab <- convertid::todisp2(lab, biomart = biomart)
	  } else {
	    xl <- xlab <- convertid::todisp2(rownames(cormat), biomart = biomart)
	  }
	}
	cat("  Column names:\n")
	yl <- ylab <- convertid::todisp2(colnames(cormat), lab.df, biomart)
	if (any(duplicated(xl))) {
		cat("  Fixing duplicated row names...\n")
		if (is.function(updateProgress)) {
			updateProgress(detail = "Fixing duplicated row names")
		}
		xl[xl %in% xl[duplicated(xl)]] <- paste(xl[xl %in% xl[duplicated(xl)]], rownames(cormat)[xl %in% xl[duplicated(xl)]], sep=".")
	}
	if (any(duplicated(yl))) {
		cat("  Fixing duplicated column names...\n")
		if (is.function(updateProgress)) {
			updateProgress(detail = "Fixing duplicated row names")
		}
		yl[yl %in% yl[duplicated(yl)]] <- paste(yl[yl %in% yl[duplicated(yl)]], rownames(cormat)[yl %in% yl[duplicated(yl)]], sep=".")
	}
	rownames(corm) <- xl
	colnames(corm) <- yl
	if (doPlot) {
		cat("@ Generating heatmap...\n")
		if (is.function(updateProgress)) {
			updateProgress(detail = "Generating heatmap")
		}
		heatmap.cor(cormat, order.list=order.list, main=main, x.labels = xlab, y.labels = ylab, labelheight=labelheight, labelwidth=labelwidth, cex=cex, add.sig=add.sig, pv=cormat.p, genes2highl=genes2highl)
	}
	return(invisible(corm))
}

#' Helper function to calculate the correlation matrix.
#' @param x (\code{ExpressionSet}, \code{data.frame} or \code{numeric}). A numeric data frame, matrix or an ExpressionSet object.
#' @param with.pvalues (\code{logical}). Should P-Values be calculated for the correlations? If \code{TRUE} P-Values will be
#'     depicted in the heatmap by significance asterisks. See 'Details'. Defaults to \code{FALSE}.
#' @details
#'     P-Values are calcluated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
#'     where r is the correlation coefficient, n is the number of samples with no missing values for each gene (row-wise
#'     \code{ncol(eset)} minus the number of columns that have an NA). P-Values are the calculated using \code{\link[stats]{pt}} and
#'     corrected account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation.
#'     The approach to calculate correlation significance was adopted from \emph{Miles, J., & Banyard, P. (2007)} on
#'     "Calculating the exact significance of a Pearson correlation in MS Excel".
#' @references
#'     Miles, J., & Banyard, P. (2007). \emph{Understanding and using statistics in psychology: A practical introduction.}
#'     Sage Publications Ltd. \url{https://psycnet.apa.org/record/2007-06525-000}.
#' @return A correlation \code{matrix} or a \code{list} with three slots: the correlation matrix, the number of samples with no
#'     missing value for each gene and the P-Values matrix.
eset_cor <- function(x, with.pvalues=FALSE) {
	if (is(x, "ExpressionSet")) {
		expr <- exprs(x)
	} else {
		expr <- as.matrix(x)
	}
	if (!is.numeric(expr)) {
		stop("'x' needs to be a numeric matrix, data frame or an ExpressionSet object!")
	}

	if (all(diag(expr)==1) && identical(colSums(expr), rowSums(expr))) {
	  cat("  Input is a correlation matrix. Doing nothing...\n")
	  cormat <- expr
	} else {
	  # remove rows and columns with all values missing
	  expr <- expr[apply(expr, 1, function(x) !all(is.na(x))), ]
	  expr <- expr[, apply(expr, 2, function(x) !all(is.na(x)))]
	  cat("  Computing correlations...\n")
	  cormat <- cor(t(expr), method="pearson", use="pairwise.complete")
	}

	if (with.pvalues) {
		cat("  Calculating P-Values...\n")
		cat("   Getting NAs...\n")
		nas <- !is.na(expr) # get a logical matrix with all non-NAs as 'TRUE' to effectively have a matrix of '0's and '1's for the next step
		cat("   Getting counts...\n")
		counts <- nas %*% t(nas) # get number of samples with no missing value for each gene, i.e., ncol(expr) minus the cols that have an NA for that gene in a matrix with the same dimensions as cormat by calculating the matrix cross-product (tcrossprod)
		df <- counts - 2 # use the above calculated matrix to obtain the degrees of freedom: N-2 since we always have a 2-dimensional correlation matrix
		t <- sqrt(df) * cormat / sqrt(1 - cormat^2) # calculate the t-test value (from "Understanding and Using Statistics in Psychology --- A Practical Introduction" on "Calculating the exact significance of a Pearson correlation in MS Excel",
		# found at http://stats.stackexchange.com/questions/120199/calculate-p-value-for-the-correlation-coefficient)
		p <- pt(t, df) # convert the test value to it's P-Value
		cat("   Getting P-Values...\n")
		pvalues <- 2 * pmin(p, 1-p) # account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation
		list(cormat=cormat, counts=counts, pvalues=pvalues)
	} else {
		cormat
	}
}

#' Cluster a correlation matrix and return the sorted matrix for plotting.
#' @description Helper function to cluster the correlation matrix and return the sorted matrix for plotting.
#' @param cormat (\code{numeric} or \code{list}). The correlation matrix or a list containing matrices of correlation values,
#'     counts and P-Values as generated by \code{eset_cor}. If a list then the correlation matrix is expected to be in slot
#'     'cor', the counts matrix in slot 'counts' and the matrix of P-Values in slot 'pvalues'.
#' @param na.frac (\code{numeric}). Fraction of missing values allowed per row of the input matrix. Defaults to \code{0.1} which
#'     means LESS than 10per centof the values in one row are allowed to be NAs.
#' @param distfn (\code{function}). Function to calculate the dissimilarity matrix for clustering. Defaults to \code{function(cm) (1-cm)}.
#' @param method (\code{character}). The agglomeration method used for clustering. See help for \code{\link[stats]{hclust}}.
#'     Defaults to "complete".
#' @param cor.thr (\code{numeric}). Correlation threshold to filter the correlation matrix for plotting. Defaults to \code{NULL} meaning
#'     no filtering. See also \option{cor.mar}. This value is sign-sensitive: a negative threshold will retain rows and columns of
#'     the correlation matrix with correlation values between -1 and the threshold, a positive value will retain rows and columns
#'     with values between the threshold and 1. Zero (0) is treated as positive.
#' @param cor.mar (\code{numeric}). Margin of the values per row of the correlation matrix the \option{cor.thr} filter needs to
#'     meet. Defaults to \code{0.5} meaning at least 50 per cent of the values in a row need to meet the threshold in order to keep the row.
#' @param cut.thr (\code{numeric}). Threshold at which dendrogram branches are to be cut. Passed on to argument \code{cutHeight} in
#'     \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{NULL} meaning no cutting.
#' @param cut.size (\code{numeric}). Minimum number of objects on a dendrogram branch considered a cluster. Passed on to argument
#'     \code{minSize} in \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{1}.
#' @param list.output (\code{logical}). Should the output be a list of different object created in the call? Depends also on input
#'     type. If \code{FALSE} only the (filtered) correlation matrix is returned.
#' @return A correlation \code{matrix} or a \code{list} or the matrix and other values needed for plotting.
clust_cormap <- function(cormat, na.frac=0.1, distfn=function(cm) (1-cm), method="complete", cor.thr=0.8, cor.mar=0.05, cut.thr=0.9, cut.size=1, list.output=FALSE) {
	cormat.orig <- cormat
	if (is.list(cormat)) {
		cormat <- cormat.orig$cormat
		counts <- cormat.orig$counts
		pval <- cormat.orig$pvalues
	}

	# filtering: keep only rows and columns with less than the given fraction of missing values; defaults to 10 per cent
	filt <- apply(cormat, 1, function(x) sum(is.na(x)) < na.frac*length(x))

	cormat2 <- cormat[filt, filt]

	dissimilarity <- distfn(cormat2)
	distance <- as.dist(dissimilarity)
	cl <- hclust(distance, method=method)

	# cut tree at given dissimilarity threshold
	if (!is.null(cut.thr)) {
		clust <- WGCNA::cutreeStatic(cl, cutHeight = cut.thr, minSize = cut.size)
		keepSamples <- (clust==1)
		cat("Isolated tree section with", length(which(keepSamples)), "genes\n")
		cormat2 <- cormat2[keepSamples, keepSamples]
		dissimilarity <- distfn(cormat2)
		distance <- as.dist(dissimilarity)
		cl <- hclust(distance, method=method)
	}

	cl_ord <- cl$order
	cormat3 <- cormat2[cl_ord, cl_ord]

	if (!is.null(cor.thr)) {
		cf <- apply(cormat3, 1, function(x) {
					cond <- sum(x >= cor.thr)
					if (sign(cor.thr) == -1) {
						cond <- sum(x <= cor.thr)
					}
					cond >= cor.mar*length(x)
		})
		if (length(which(cf)) < 2) {
		  cat("NOTE: Filtering left no values to plot. Returning unfiltered matrix.\n")
		} else {
		  cat("Filtered", nrow(cormat3)-length(which(cf)), paste0("rows not meeting correlation threshold (", cor.thr, " in ", cor.mar*100, "% of the columns)"), "\n")
		  cormat3 <- cormat3[cf, cf]
		}
	}

	if (list.output) {
		if (is.list(cormat.orig)) {
			list(cormat=cormat3, hclust=cl, n=counts[filt, filt][cl_ord, cl_ord], pval=pval[filt, filt][cl_ord, cl_ord])
		} else {
			list(cormat=cormat3, hclust=cl)
		}
	} else {
		cormat3
	}
}

#' Helper function to draw the heatmap. Wrapper around 'heatmapFlex::heatmap.n2'.
#' @noRd
heatmap.cor <- function(cormat, order.list = TRUE, main="", main_postfix="Dissimilarity = 1 - Correlation", x.labels = NULL, y.labels = NULL, labelheight=0.2, labelwidth=0.2, cex = 0.5, add.sig=FALSE, pv=NULL, genes2highl=NULL) {
	main.map <- paste(main, main_postfix, sep="\n")
	col <- colorRampPalette(c("blue", "white", "red"))(50)
	if(order.list){
		cormat <- cormat[nrow(cormat):1, ]
		if (add.sig) pv <- pv[nrow(pv):1, ]
		x.labels <- rev(x.labels)
	}
	heatmapFlex::heatmap.n2(cormat, order_list = order.list, reorder=c(FALSE, FALSE), labRow = x.labels, labCol = y.labels, labelheight=labelheight, labelwidth=labelwidth, col=col, main=main.map, r.cex=cex, c.cex=cex, add.sig=add.sig, pv=pv, genes2highl=genes2highl)
}


#' Helper function to test for the difference in strengths of two correlations
#' returns the "z-scores" (approximately standard normal distributed) or their p-values
#' works on matrices, vectors, single values
#' @noRd
compare_cor <- function(cor1, counts1, cor2, counts2, pvalues=TRUE) {
	# according to http://www.fon.hum.uva.nl/Service/Statistics/Two_Correlations.html
	# from Papoulis, Athanasios "Probability and Statistics"
	#      Prentence-Hall International Editions
	#      ISBN: 0 13 711730 2, 1990
	# Zf = 1/2 * ln( (1+R) / (1-R) )
	# z = (Zf1 - Zf2) / SQRT( 1/(N1-3) + 1/(N2-3) )
	Zf1 <- log( (1+cor1) / (1-cor1) ) / 2
	Zf2 <- log( (1+cor2) / (1-cor2) ) / 2
	z <- (Zf1 - Zf2) / sqrt( 1/(counts1-3) + 1/(counts2-3) )
	if (pvalues) {
		p <- pnorm(z)
		pvalues <- 2 * pmin(p, 1-p)
		pvalues
	} else {
		z
	}
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
	# try negative mean of the log10 of the significances of the correlations or similar - log10(0) can be problematic
	-log10(median(pmax(c(x, recursive=TRUE), 1e-38), na.rm=TRUE))
}

#' Function to automatically split clusters based on noise level and hierarchy, produce a list of vectors of genes, ...
#' @noRd
cormap3 <- function(x, na.frac=0.1, method="ward", do.abs=TRUE, main="correlation map", threshold=2, cex=0.2, h=NULL, cex.clust=cex, cex.filt=cex, do.plots=1:3, genes2highl=NULL) {
	if (dev.interactive()) par(mfrow=c(1,1))

	distfn <- if (do.abs) function(cm) { 1-abs(cm) } else function(cm) { 1-cm }
	cormat <- eset_cor(x, with.pvalues = TRUE)
	cormat.cl <- clust_cormap(cormat, na.frac=na.frac, distfn=distfn, method=method, list.output=TRUE)
	postfix <- if (do.abs) "Dissimilarity = 1 - abs(Correlation)" else "Dissimilarity = 1 - Correlation"

	if (1 %in% do.plots) {
		plot(cormat.cl$hclust, cex=cex.clust, labels=todisp2(rownames(cormat.cl$cormat)))
		if (!is.null(h)) abline(h=h)
	}
	if (2 %in% do.plots) {
		heatmap.cor(cormat.cl$cormat, main=main, main_postfix=postfix, x.labels = todisp2(rownames(cormat.cl$cormat)), y.labels = todisp2(colnames(cormat.cl$cormat)), cex=cex, genes2highl=genes2highl)
	}

	if (!is.null(h)) {
		# compare the two above to find interesting clusters
		cut_cluster <- cut(as.dendrogram(cormat.cl$hclust), h=h)
		cluster_labels <- lapply(cut_cluster$lower, function(x) cormat.cl$hclust$labels[dendro2leaves(x)])
		cluster_signal_metrics <- sapply(cluster_labels, function(ids) signal_metric(cormat.cl$pval[ids, ids]) )
		# threshold of two corresponds to a median p-value of 0.01
		filt <- cluster_signal_metrics > threshold
		ids <- lunion(cluster_labels[filt])
		# filtered cormap
		cormat.cl <- clust_cormap(cormat$cor[ids, ids], na.frac=na.frac, distfn=distfn, method=method, list.output=TRUE)
		if (3 %in% do.plots) {
			heatmap.cor(cormat.cl$cormat, main=main, main_postfix=postfix, x.labels = todisp2(rownames(cormat.cl$cormat)), y.labels = todisp2(colnames(cormat.cl$cormat)), cex=cex.filt, genes2highl=genes2highl)
		}
		invisible(append(list(clusters=cluster_labels, filt=filt, filt_clusters=cluster_labels[filt], h=h, threshold=threshold, metric=cluster_signal_metrics), cormat.cl))
	} else {
		invisible(cormat.cl)
	}
}
