#' @noRd
plotAdjust <-
function (dm, cormap = FALSE)
{

  adj.l <- vector(mode="list")
  adj.l$r.cex <- round(1.05/log(sqrt(nrow(dm))), 1)
  if (nrow(dm) <= 10) {
    adj.l$r.cex <- 0.8
  }
  adj.l$c.cex <- round(1.4/log(sqrt(ncol(dm))), 1)
  if (ncol(dm) <= 10) {
    adj.l$c.cex <- 0.8
  }
  adj.l$pdf.width <- ceiling(sqrt(log(ncol(dm))) * ceiling(sqrt(ncol(dm)))) * 1.8
  adj.l$pdf.height <- ceiling(sqrt(nrow(dm)))
  if(adj.l$pdf.height < 5) {
    adj.l$pdf.height <- 5
  }
  mcc <- max(nchar(colnames(dm)))
  adj.l$labelheight <- round(mcc/100/1.15, 2) / ceiling(adj.l$pdf.height/9)
  mc <- max(nchar(rownames(dm)))
  adj.l$labelwidth <- round(log10(sqrt(mc)) / (sqrt(mc) * mc ^ (-mc / 100)), 2) * 1.2 * adj.l$r.cex
  if (cormap) {
    if (nrow(dm) > 190) {
      adj.l$r.cex <- adj.l$c.cex <- round(adj.l$r.cex / log(nrow(dm)), 1)
    }
    if (adj.l$r.cex < 0.1)
      adj.l$r.cex <- adj.l$c.cex <- 0.1
  }
  return(adj.l)

}
