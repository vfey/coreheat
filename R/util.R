lintersect <- function(lst) {
        res <- lst[[1]]
        if (length(lst) > 1) {
                for (i in 2:length(lst)) {
                        res <- intersect(res, lst[[i]])
                }
        }
        res
}

lunion <- function(lst) {
        res <- lst[[1]]
        if (length(lst) > 1) {
                for (i in 2:length(lst)) {
                        res <- union(res, lst[[i]])
                }
        }
        res
}


id <- function(x) {
	x
}

