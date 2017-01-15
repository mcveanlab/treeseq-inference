#' Use tree metrics to compare two inferred ancestries (ARGs) over a genomic region
#'
#' Requires two multiPhylo objects which contain sequential trees (tree sequences) along a genomic region.
#' Tips labels in the two tree sequences should correspond (and are assumed to be 1..N) 
#' Trees within the multiPhylo object should be in order along the genome, and each tree
#' must have a name which can be converted to a number. The name (number) should give
#' the (non-inclusive) rightmost genome position covered by that tree.
#' @param treeseq.a The first multiPhylo object, or a path to a .nex file
#' @param treeseq.b The second multiPhylo object, or a path to a .nex file
#' @param output.full.table Output tree metrics for each overlapping region, rather than simply a weighted summary.
#' @param acceptable.length.diff.pct How much difference in sequence length is allows between the 2 trees? (Default: 0.1 percent)
#' @param variant.positions A vector of genome positions of each variant. If given the metric will be 
#'  calculated for each variant site and averaged over all sites, rather than averaged over every point on the genome (not yet implemented)
#' @param randomly.resolve.a Some distance metrics only operate on binary trees. Set this to TRUE to 
#' force trees in treeseq.a to be binary by randomly resolving polytomies where necessary, using the
#' multi2di routine in the 'ape' library. If this is a number, it is additionally used as an RNG seed.
#' @param randomly.resolve.b Some distance metrics only operate on binary trees. Set this to TRUE to 
#' force trees in treeseq.b to be binary by randomly resolving polytomies where necessary, using the
#' multi2di routine in the 'ape' library. If this is a number, it is additionally used as an RNG seed.
#' @export
#' @examples
#' genome.trees.dist()

genome.trees.dist <- function(treeseq.a, treeseq.b, output.full.table = FALSE, acceptable.length.diff.pct = 0.1, variant.positions=NULL, randomly.resolve.a=FALSE, randomly.resolve.b=FALSE) { 
    require(ape)
    require(phangorn) #to use the various treedist metrics

    if (randomly.resolve.a) {
        if (is.numeric(randomly.resolve.polytomies))
            set.seed(randomly.resolve.polytomies)
        process.a = multi2di
    } else {
        process.a = identity
    }

    if (randomly.resolve.b) {
        if (is.numeric(randomly.resolve.polytomies))
            set.seed(randomly.resolve.polytomies)
        process.b = multi2di
    } else {
        process.b = identity
    }

    #check if trees or filenames
    if (class(treeseq.a) != "multiPhylo") {
        a <- process.a(new.read.nexus(treeseq.a, force.multi=TRUE))
    } else {
        a <- process.a(treeseq.a)
    }
    if (class(treeseq.b) != "multiPhylo") {
         b <- process.b(new.read.nexus(treeseq.b, force.multi=TRUE))
    } else {
         b <- process.b(treeseq.b)
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
    results=data.frame(lft=numeric(), rgt=numeric(), RFrooted=numeric(), RFunrooted=numeric(),
        wRFrooted=numeric(), wRFunrooted=numeric(), SPRunrooted=numeric(), pathunrooted=numeric())
    for (o in order(as.numeric(names(breaks.table)))) {
        brk = breaks.table[[o]]
        if (any(tree.index.counter[brk$ind] > c(length(a),length(b))[brk$ind])) {
            warning("Reached the end of the trees with ", max(brk.a, brk.b)-lft, " left to go")
            break
        }
        rgt <- brk$values[0:1]
        RFrooted <- RFunrooted <- wRFrooted <- wRFunrooted <- SPRunrooted  <- pathunrooted <- NA
        tryCatch({RFrooted <- RF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=TRUE)}, 
                 error=function(e){})
        tryCatch({RFunrooted <- RF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=FALSE)},
                 error=function(e){})
        tryCatch({wRFrooted <- wRF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=TRUE)},
                 error=function(e){})
        tryCatch({wRFunrooted <- wRF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=FALSE)},
                 error=function(e){})
        tryCatch({SPRunrooted <- SPR.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]])},
                 error=function(e){})
        tryCatch({pathunrooted <- path.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]])},
                 error=function(e){})
        results[nrow(results)+1,] <- c(lft,rgt,RFrooted, RFunrooted, wRFrooted, wRFunrooted, SPRunrooted, pathunrooted)
        lft <- rgt
        tree.index.counter[brk$ind] <- tree.index.counter[brk$ind] + 1 #NB, brk$ind is a factor with levels (m1,m2), so we hope that m1==1 and m2==2
    }
    if (output.full.table) {
        return(results)
    } else {
        return(c(RFrooted=weighted.mean(results$RFrooted, (results$rgt-results$lft), na.rm = TRUE),
                 RFunrooted=weighted.mean(results$RFunrooted, (results$rgt-results$lft), na.rm = TRUE),
                 wRFrooted=weighted.mean(results$wRFrooted, (results$rgt-results$lft), na.rm = TRUE),
                 wRFunrooted=weighted.mean(results$wRFunrooted, (results$rgt-results$lft), na.rm = TRUE),
                 SPRunrooted=weighted.mean(results$SPRunrooted, (results$rgt-results$lft), na.rm = TRUE),
                 pathunrooted=weighted.mean(results$pathunrooted, (results$rgt-results$lft), na.rm = TRUE)))
    }
}
