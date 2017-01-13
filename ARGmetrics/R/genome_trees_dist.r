#' Compare multiple trees across a genomic region
#'
#' Requires multiPhylo objects representing sequential trees along a genomic region.
#'  The tree names should be numbers, which represent the rightmost position along a
#'  genome (non-inclusive) covered by each tree.
#' @param a The first multiPhylo object, or a path to a .nex file
#' @param b The second multiPhylo object, or a path to a .nex file
#' @param output.full.table Output tree metrics for each overlapping region, rather than simply a weighted summary.
#' @param acceptable.length.diff.pct How much difference in sequence length is allows between the 2 trees? (Default: 0.1 percent)
#' @param variant.positions A list of positions of each variant (not implemented)
#' @export
#' @examples
#' genome.trees.dist()

genome.trees.dist <- function(a, b, output.full.table = FALSE, acceptable.length.diff.pct = 0.1, variant.positions=NULL) { 
    #a and b should be multiPhylo objects containing multiple trees
    #if variant.positions is given, it should be a vector of genome positions for each variants, and the metric will be compared per variant site
    #rather than per genomic region
    require(phangorn) #to use the various treedist metrics
    if (class(a) != "multiPhylo") {
         a <- new.read.nexus(a)
    }
    if (class(b) != "multiPhylo") {
         b <- new.read.nexus(b)
    }
    brk.a <- as.numeric(names(a))
    if (is.unsorted(brk.a))
        stop("Tree names should correspond to numerical breakpoints or variant totals, sorted from low to high, but tree names in the first trees object are not numbers in ascending order.")
    brk.b <- as.numeric(names(b))
    if (is.unsorted(brk.b))
        stop("Tree names should correspond to numerical breakpoints or variant totals, sorted from low to high, but tree names in the second trees object are not numbers in ascending order.")
    if ((max(brk.a) * (100+ acceptable.length.diff.pct)/100 < max(brk.b)) || 
        (max(brk.b) * (100+ acceptable.length.diff.pct)/100 < max(brk.a)))
        warning("The sequence lengths of the two trees files differ markedly: ", max(brk.a), " vs. ", max(brk.b), immediate. = TRUE)
    breaks.table <- stack(list('a'=brk.a,'b'=brk.b))
    breaks.table <- by(breaks.table, breaks.table$values, function(x) x) #put identical breakpoints together
    tree.index.counter=c(a=1, b=1) 
    lft <- 0
    results=data.frame(lft=numeric(), rgt=numeric(), RF.rooted=numeric(), RF.unrooted=numeric(),
        wRF.rooted=numeric(), wRF.unrooted=numeric(), SPR.unrooted=numeric(), path.unrooted=numeric())
    for (o in order(as.numeric(names(breaks.table)))) {
        brk = breaks.table[[o]]
        if (any(tree.index.counter[brk$ind] > c(length(a),length(b))[brk$ind])) {
            warning("Reached the end of the trees with ", max(brk.a, brk.b)-lft, " left to go")
            break
        }
        rgt <- brk$values[0:1]
        RF.rooted <- RF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=TRUE)
        RF.unrooted <- RF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=FALSE)
        wRF.rooted <- wRF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=TRUE)
        wRF.unrooted <- wRF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=FALSE)
        SPR.unrooted <- SPR.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]])
        path.unrooted <- path.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]])
        results[nrow(results)+1,] <- c(lft,rgt,RF.rooted, RF.unrooted, wRF.rooted, wRF.unrooted, SPR.unrooted, path.unrooted)
        lft <- rgt
        tree.index.counter[brk$ind] <- tree.index.counter[brk$ind] + 1 #NB, brk$ind is a factor with levels (m1,m2), so we hope that m1==1 and m2==2
    }
    if (output.full.table) {
        return(results)
    } else {
        return(c(RF.rooted=weighted.mean(results$RF.rooted, (results$rgt-results$lft)),
                 RF.unrooted=weighted.mean(results$RF.unrooted, (results$rgt-results$lft)),
                 wRF.rooted=weighted.mean(results$wRF.rooted, (results$rgt-results$lft)),
                 wRF.unrooted=weighted.mean(results$wRF.unrooted, (results$rgt-results$lft)),
                 SPR.unrooted=weighted.mean(results$SPR.unrooted, (results$rgt-results$lft)),
                 path.unrooted=weighted.mean(results$path.unrooted, (results$rgt-results$lft))))
    }
}
