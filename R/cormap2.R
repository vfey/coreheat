cormap2 <- function(eset, cormat=NULL, lab=NULL, biomart=TRUE, cluster_correlations=TRUE, main="", cex=0.7, minfrac=0.1, cor.thr=NULL, cor.mar=0.8, cut.thr=NULL, cut.size=1, labelheight=0.2,
		labelwidth=0.2, add.sig=FALSE, genes2highl=NULL, doPlot=TRUE, updateProgress=NULL) {
	cat("@ Calculating correlation matrix...\n")
	if (is.function(updateProgress)) {
		updateProgress(detail = "Calculating correlation matrix")
	}
	if (add.sig) {
		cormat <- eset_cor(eset, cormat=cormat, with.pvalues = TRUE)
		list.output <- TRUE
	} else {
		cormat <- eset_cor(eset, cormat=cormat)
		list.output <- FALSE
	}
	if (is.function(updateProgress)) {
		updateProgress(detail = "Clustering and filtering correlation matrix")
	}
	if (cluster_correlations) {
		cat("@ Clustering and filtering correlation matrix...\n")
		cormat <- clust_cormap(cormat, minfrac=minfrac, cor.thr=cor.thr, cor.mar=cor.mar, cut.thr=cut.thr, cut.size=cut.size, list.output = list.output)
	}
	if (list.output) {
		cormat.p <- cormat$pval[nrow(cormat$cormat):1, ]
		corm <- cormat <- cormat$cormat[nrow(cormat$cormat):1, ]
	} else {
		cormat.p <- NULL
		corm <- cormat <- cormat[nrow(cormat):1, ]
	}
	if (is.null(lab)) {
		lab <- data.frame(ID=rownames(cormat), stringsAsFactors=F)
		rownames(lab) <- rownames(cormat)
	}
	cat("@ Formatting output matrix...\n")
	if (is.function(updateProgress)) {
		updateProgress(detail = "Formatting matrix")
	}
	cat("  Row names:\n")
	xl <- xlab <- todisp2(rownames(cormat), lab, biomart)
	cat("  Column names:\n")
	yl <- ylab <- todisp2(colnames(cormat), lab, biomart)
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
		heatmap.cor(cormat, main=main, x.labels = xlab, y.labels = ylab, labelheight=labelheight, labelwidth=labelwidth, cex=cex, add.sig=add.sig, pv=cormat.p, genes2highl=genes2highl)
	}
	return(corm)
}

# Helper function to calculate the correlation matrix.

eset_cor <- function(eset, cormat=NULL, with.pvalues=FALSE) {
	if (is(eset, "ExpressionSet")) {
		expr <- exprs(eset)
	} else{
		expr <- as.matrix(eset)
	}
	if (!is.numeric(expr)) {
		stop("'eset' needs to be a numerical matrix or an ExpressionSet object!")
	}

	# remove rows and columns with all values missing
	expr <- expr[apply(expr, 1, function(x) !all(is.na(x))), ]
	expr <- expr[, apply(expr, 2, function(x) !all(is.na(x)))]

	if (!is.null(cormat)) {
		if (!all(diag(cormat)==1)) {
			stop("Input is not a correlation matrix!")
		}
		cat("  Input is correlation matrix. Doing nothing...\n")
	} else {
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

# Helper function to cluster the correlation matrix and return the sorted matrix for plotting.

clust_cormap <- function(cormat, minfrac=0.1, distfn=function(cm) (1-cm), method="complete", cor.thr=0.8, cor.mar=0.05, cut.thr=0.9, cut.size=1, list.output=FALSE) {
	cormat.orig <- cormat
	if (is.list(cormat)) {
		cormat <- cormat.orig$cormat
		counts <- cormat.orig$counts
		pval <- cormat.orig$pvalues
	}

	# filtering: keep only rows and columns with less than the given fraction of missing values; defaults to 10%
	filt <- apply(cormat, 1, function(x) sum(is.na(x)) < minfrac*length(x))

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
		cat("Filtered", nrow(cormat3)-length(which(cf)), paste0("rows not meeting correlation threshold (", cor.thr, " in ", cor.mar*100, "% of the columns)"), "\n")
		cormat3 <- cormat3[cf, cf]
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

# Helper function to draw the heatmap. Wrapper around 'heatmap.n2' from the 'heatmapGen2' package.

heatmap.cor <- function(cormat, order_list = TRUE, main="", main_postfix="Dissimilarity = 1 - Correlation", x.labels = NULL, y.labels = NULL, labelheight=0.2, labelwidth=0.2, cex = 0.5, add.sig=FALSE, pv=NULL, genes2highl=NULL) {
	main.map <- paste(main, main_postfix, sep="\n")
	col <- colorRampPalette(c("blue", "white", "red"))(50)
	if(order_list){
		cormat <- cormat[nrow(cormat):1, ]
		if (add.sig) pv <- pv[nrow(pv):1, ]
		x.labels <- rev(x.labels)
	}
	heatmap.n2(cormat, order_list = order_list, reorder=c(FALSE, FALSE), labRow = x.labels, labCol = y.labels, labelheight=labelheight, labelwidth=labelwidth, col=col, main=main.map, r.cex=cex, c.cex=cex, add.sig=add.sig, pv=pv, genes2highl=genes2highl)
}


#' Helper function test for the difference in strengths of two correlations
# returns the "z-scores" (approximately standard normal distributed) or their p-values
# works on matrices, vectors, single values
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

#' Wrapper function to draw up to three different heatmaps from three different data sets and in addition, draw two heatmaps
#' showing the difference between the correlation maps of the first tow data sets.
#'
cor_maps <- function(full.eset.first,full.eset.second, full.eset.third=NULL, heading=NULL,clust=TRUE, groupnames=c("cancer", "normal", "cell lines"), cex=0.5, add.sig=FALSE, genes2highl=NULL){
	#First map cormat
	cormat <- eset_cor(full.eset.first, with.pvalues=TRUE)
	#Second map cormat
	cormat2 <- eset_cor(full.eset.second, with.pvalues=TRUE)
	# Third map cormat
	if (!is.null(full.eset.third)) cormat3 <- eset_cor(full.eset.third, with.pvalues=TRUE)

	shared_genes <- intersect(rownames(cormat$cor), rownames(cormat2$cor))
	row_ok1 <- rownames(cormat$cor) %in% shared_genes
	row_ok2 <- rownames(cormat2$cor) %in% shared_genes
	if (!is.null(full.eset.third)) row_ok3 <- rownames(cormat3$cor) %in% shared_genes
	cormat.cor <- cormat$cor[row_ok1,row_ok1]
	cormat.N <- cormat$counts[row_ok1,row_ok1]
	cormat2.cor <- cormat2$cor[row_ok2,row_ok2]
	cormat2.N <- cormat2$counts[row_ok2,row_ok2]
	if (!is.null(full.eset.third)) cormat3.cor <- cormat3$cor[row_ok3,row_ok3]
	if (!is.null(full.eset.third)) cormat3.N <- cormat3$counts[row_ok3,row_ok3]
	# cormat p-values
	cormat.p1 <- cormat$pvalues[row_ok1,row_ok1]
	cormat.p2 <- cormat2$pvalues[row_ok2,row_ok2]
	if (!is.null(full.eset.third)) cormat.p3 <- cormat3$pvalues[row_ok3,row_ok3]

	#Using hclust to order the map
	if(clust){
		Na_values <- apply(cormat.cor,1,function(x)any(is.na(x)))
		cormat.cor <- cormat.cor[!Na_values,!Na_values]
		tmp_hclust_cormat <- clust_cormap(cormat.cor)
		ordered_genelist <- colnames(tmp_hclust_cormat)
	} else {
		#using chromosome position to order the map
		ensg_order <- rownames(cormat.cor)
		chr_pos <- ga$ga$gene_chrom_start[match(ensg_order, ga$ga$ensg_id)]
		ordered_genelist <- colnames(cormat.cor)[order(chr_pos)]
	}
	cormat.cor <- cormat.cor[ordered_genelist, ordered_genelist]
	cormat.N <- cormat.N[ordered_genelist, ordered_genelist]
	cormat2.cor <- cormat2.cor[ordered_genelist, ordered_genelist]
	cormat2.N <- cormat2.N[ordered_genelist, ordered_genelist]
	if (!is.null(full.eset.third)) cormat3.cor <- cormat3.cor[ordered_genelist, ordered_genelist]
	if (!is.null(full.eset.third)) cormat3.N <- cormat3.N[ordered_genelist, ordered_genelist]
	cormat.p1 <- cormat.p1[ordered_genelist, ordered_genelist]
	cormat.p2 <- cormat.p2[ordered_genelist, ordered_genelist]
	if (!is.null(full.eset.third)) cormat.p3 <- cormat.p3[ordered_genelist, ordered_genelist]

	#First (cancer) heatmap
	heatmap.cor(cormat.cor[nrow(cormat.cor):1, ],order_list = FALSE,main=heading,main_postfix=paste("correlation map ", groupnames[1], ", Dissimilarity = 1 - Correlation", sep=""),x.labels = todisp(rev(ordered_genelist)), y.labels = todisp(ordered_genelist), cex=cex, add.sig=add.sig, pv=cormat.p1[nrow(cormat.cor):1, ], genes2highl=genes2highl)
	#Second (normal) heatmap
	heatmap.cor(cormat2.cor,main=heading, main_postfix=paste("genes correlating map ",groupnames[2], ", Dissimilarity = 1 - Correlation", sep=""),x.labels = todisp(rownames(cormat2.cor)), y.labels = todisp(colnames(cormat2.cor)), cex=cex, add.sig=add.sig, pv=cormat.p2, genes2highl=genes2highl)
	#Third (cell lines) heatmap
	if (!is.null(full.eset.third)) heatmap.cor(cormat3.cor,main=heading, main_postfix=paste("genes correlating map ",groupnames[3], ", Dissimilarity = 1 - Correlation", sep=""),x.labels = todisp(rownames(cormat2.cor)), y.labels = todisp(colnames(cormat2.cor)), cex=cex, add.sig=add.sig, pv=cormat.p3, genes2highl=genes2highl)

	#Correlation matrix comparison using significance level

	sig.cormat <- compare_cor(cormat.cor,cormat.N,cormat2.cor,cormat2.N,pvalues=TRUE)
	#Removing Nas and Inf values in p-value matrix
	diag(sig.cormat) <- 1
	# non-ideal way to prevent infinite values if p-value == 0 . The /10 is not really correct, but sets 0 p-values to a small non-zero value.
	sig.cormat <- -log10(pmax(sig.cormat, min(sig.cormat[sig.cormat > 0],na.rm = TRUE) / 10))

	sig.cormat[which(sig.cormat == Inf)] <- NA
	sig.cormat[which(sig.cormat == -Inf)] <- NA
	heatmap.cor(sig.cormat,main="Plot of significance (-log10(p-value)) of difference between two correlation maps",x.labels = todisp(rownames(sig.cormat)), y.labels = todisp(colnames(sig.cormat)), cex=cex, genes2highl=genes2highl)

	#Correlation matrix comparison using z-values
	sig.cormat <- compare_cor(cormat.cor,cormat.N,cormat2.cor,cormat2.N,pvalues=FALSE)
	sig.cormat[which(sig.cormat == Inf)] <- NA
	sig.cormat[which(sig.cormat == -Inf)] <- NA

	heatmap.cor(sig.cormat,main="Plot of z-values showing differences between two correlation maps",main_postfix="Dissimilarity = 1 - Correlation", x.labels = todisp(rownames(sig.cormat)), y.labels = todisp(colnames(sig.cormat)), cex=cex, genes2highl=genes2highl)

}

#' Helper function to organise correlation tables.
organised_cor_table <- function(eset){
	x <- eset_cor(eset,with.pvalues=TRUE)
	flat_table <- data.frame()
	for(cormat11 in x){
		flatcor <- do.call("rbind", lapply(cormat11, function(x) x))
		chrom1 <- ga$ga$chrom_name[match(as.character(rownames(cormat11)), ga$ga$ensg_id)]
		chrom2 <- ga$ga$chrom_name[match(as.character(colnames(cormat11)), ga$ga$ensg_id)]
		if(length(flat_table)==0){
			flatcor_genes <- data.frame(gene1=rep(todisp(rownames(cormat11))), chrom1=rep(chrom1,ncol(cormat11)), gene2=rep(todisp(colnames(cormat11)), each=nrow(cormat11)),chrom1=rep(chrom1,each=ncol(cormat11)))
			flat_table <- data.frame(flatcor_genes,flatcor)
		}else{
			flat_table <- data.frame(flat_table,flatcor)
		}
	}
	table_names <- c("gene1","chrom1","gene2","chrom2","cor","counts","pval")
	colnames(flat_table) <- table_names
	flat_table
}

#' Helper function to extract leaves from a dendrogram. This is adapted from stats:::cut.dendrogram.
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
signal_metric <- function(x) {
	# try negative mean of the log10 of the significances of the correlations or similar - log10(0) can be problematic
	-log10(median(pmax(c(x, recursive=TRUE), 1e-38), na.rm=TRUE))
}

#' Function to automatically split clusters based on noise level and hierarchy, produce a list of vectors of genes, ...
cormap3 <- function(eset, minfrac=0.1, method="ward", do.abs=TRUE, main="correlation map", threshold=2, cex=0.2, h=NULL, cex.clust=cex, cex.filt=cex, do.plots=1:3, genes2highl=NULL) {
	if (dev.interactive()) par(mfrow=c(1,1))

	distfn <- if (do.abs) function(cm) { 1-abs(cm) } else function(cm) { 1-cm }
	cormat <- eset_cor(eset, with.pval=TRUE)
	cormat.cl <- clust_cormap(cormat, minfrac=minfrac, distfn=distfn, method=method, list.output=TRUE)
	postfix <- if (do.abs) "Dissimilarity = 1 - abs(Correlation)" else "Dissimilarity = 1 - Correlation"

	if (1 %in% do.plots) {
		plot(cormat.cl$hclust, cex=cex.clust, labels=todisp(rownames(cormat.cl$cormat)))
		if (!is.null(h)) abline(h=h)
	}
	if (2 %in% do.plots) {
		heatmap.cor(cormat.cl$cormat, main=main, main_postfix=postfix, x.labels = todisp(rownames(cormat.cl$cormat)), y.labels = todisp(colnames(cormat.cl$cormat)), cex=cex, genes2highl=genes2highl)
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
		cormat.cl <- clust_cormap(cormat$cor[ids, ids], minfrac=minfrac, distfn=distfn, method=method, list.output=TRUE)
		if (3 %in% do.plots) {
			heatmap.cor(cormat.cl$cormat, main=main, main_postfix=postfix, x.labels = todisp(rownames(cormat.cl$cormat)), y.labels = todisp(colnames(cormat.cl$cormat)), cex=cex.filt, genes2highl=genes2highl)
		}
		invisible(append(list(clusters=cluster_labels, filt=filt, filt_clusters=cluster_labels[filt], h=h, threshold=threshold, metric=cluster_signal_metrics), cormat.cl))
	} else {
		invisible(cormat.cl)
	}
}
