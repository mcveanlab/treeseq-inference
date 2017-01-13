#' Use tree metrics to compare multiple inferred ancestries (ARGs) over a genomic region
#'
#' Runs genome.trees.dist to compare multiple ancestry estimates against a known original
#' Tips labels in the two objects should correspond (and are assumed to be 1..N) 
#' Trees within the multiPhylo object should be in order along the genome, and each tree
#' must have a name which can be converted to a number. The name (number) should give
#' the (non-inclusive) rightmost genome position covered by that tree.
#' @param treeseq.base The base multiPhylo object, or a path to a .nex file
#' @param treeseq.multi A list of multiPhylo objects (or list of paths to .nex files).
#' If the list has names which can be all be converted to numbers, then these are treated as 
#' weights, and the returned value is a single average set of distance measures over all
#' the different tree sequences.
#' @param acceptable.length.diff.pct How much difference in sequence length is allows between the 2 trees? (Default: 0.1 percent)
#' @param variant.positions A list of positions of each variant (not implemented)
#' @export
#' @examples
#' genome.trees.dist()

genome.trees.dist.multi <- function(treeseq.base, treeseq.multi, acceptable.length.diff.pct = 0.1, variant.positions=NULL) { 
    if (class(treeseq.multi) == "multiPhylo") {
        stop("treeseq.multi should contain a *list* of multiPhylo objects, not simply a single multiPhylo object.")
    }
    metrics <- lapply(treeseq.multi, genome.trees.dist, treeseq.base, acceptable.length.diff.pct, variant.positions)
    weights <- suppressWarnings(as.numeric(names(metrics))
    if (any(is.na(weights))) {
        return(metrics)
    } else {
        
    }
}
