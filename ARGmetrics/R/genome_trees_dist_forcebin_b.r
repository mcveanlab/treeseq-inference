#' Use binary tree metrics to compare two inferred ancestries (ARGs) over a genomic region
#'
#' Runs genome.trees.dist(..., randomly.resolve.polytomies=TRUE) multiple times on the same
#' two tree sequences, breaking polytomies in the second tree sequence differentially (at random)
#' each time.
#' See ?genome.trees.dist for more information.
#' @param treeseq.a The base multiPhylo object, or a path to a .nex file
#' @param treeseq.b The multiPhylo object (or path to .nex files) containing a nonbinary tree 
#' to compare to treeseq.a.
#' @param replicates The number of times to run genome.trees.dist()
#' @param seed The random seed to give to set.seed before starting random polytomy breaking
#' @param acceptable.length.diff.pct How much difference in sequence length is allowed between the 2 trees? (Default: 0.1 percent)
#' @param variant.positions A list of positions of each variant (not implemented)
#' @return an average of the metrics over each of the replicates 
#' @export
#' @examples
#' genome.trees.dist.forcebin()

genome.trees.dist.forcebin.b <- function(treeseq.a, treeseq.b, replicates=1, seed=NA, acceptable.length.diff.pct = 0.1, variant.positions=NULL, threads=1) { 
   require(ape)
   require(parallel) #always available, as this is in R base
   if (class(treeseq.a) != "multiPhylo") {
        a <- read.nexus(treeseq.a, force.multi=TRUE)
    } else {
        a <- treeseq.a
    }
    if (class(treeseq.b) != "multiPhylo") {
         b <- read.nexus(treeseq.b, force.multi=TRUE)
    } else {
         b <- treeseq.b
    }
    
    if (!is.na(seed)) {
        set.seed(seed)
        seeds <- sample.int(2^31-1, replicates) #pick an integer seed - see ?sample.int
    } else {
        seeds <- TRUE
    }
    return(colMeans(do.call(rbind,mcmapply(
        function(r, seed) genome.trees.dist(treeseq.a, treeseq.b, randomly.resolve.b=seed,
                             acceptable.length.diff.pct = acceptable.length.diff.pct,
                             variant.positions = variant.positions), 
        1:replicates,
        seeds, 
        SIMPLIFY = FALSE, mc.cores=threads))))
}
