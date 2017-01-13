#' Use tree metrics to compare multiple inferred ancestries (ARGs) over a genomic region
#'
#' Runs genome.trees.dist to compare multiple ancestry estimates against a known original
#' See ?genome.trees.dist for more information.
#' @param treeseq.base The base multiPhylo object, or a path to a .nex file
#' @param treeseq.multi A list of multiPhylo objects (or list of paths to .nex files).
#' If the list has names which can be all be converted to numbers, then these are treated as 
#' weights, and the returned value is a single average set of distance measures over all
#' the different tree sequences.
#' @param acceptable.length.diff.pct How much difference in sequence length is allows between the 2 trees? (Default: 0.1 percent)
#' @param variant.positions A list of positions of each variant (not implemented)
#' @param randomly.resolve.polytomies Some distance metrics only operate on binary trees. Set this to TRUE to force trees to be binary
#' by randomly resolving polytomies where necessary. If a number, it is passed to set.seed as a RNG seed.
#' @export
#' @examples
#' genome.trees.dist.multi()

genome.trees.dist.multi <- function(treeseq.base, treeseq.multi, acceptable.length.diff.pct = 0.1, variant.positions=NULL, randomly.resolve.polytomies=FALSE) { 
    if (class(treeseq.multi) == "multiPhylo") {
        stop("treeseq.multi should contain a *list* of multiPhylo objects, not simply a single multiPhylo object.")
    }
    metrics <- do.call(rbind,lapply(treeseq.multi, 
                                    genome.trees.dist, 
                                    treeseq.base,
                                    acceptable.length.diff.pct = acceptable.length.diff.pct,
                                    variant.positions = variant.positions,
                                    randomly.resolve.polytomies = randomly.resolve.polytomies))
    weights <- suppressWarnings(as.numeric(rownames(metrics)))
    if (any(is.na(weights))) {
        return(metrics)
    } else {
        return(colSums(metrics * weights)/sum(weights))
    }
}
