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
#' @param acceptable.length.diff.pct How much difference in sequence length is allowed between the 2 trees? (Default: 0.1 percent)
#' @param variant.positions A vector of genome positions of each variant. If given the metric will be 
#'  calculated for each variant site and averaged over all sites, rather than averaged over every point on the genome
#' @param randomly.resolve.a Some distance metrics only operate on binary trees. Set this to TRUE to 
#' force trees in treeseq.a to be binary by randomly resolving polytomies where necessary, using the
#' multi2di routine in the 'ape' library. If this is a number, it is additionally used as an RNG seed.
#' @param randomly.resolve.b Some distance metrics only operate on binary trees. Set this to TRUE to 
#' force trees in treeseq.b to be binary by randomly resolving polytomies where necessary, using the
#' multi2di routine in the 'ape' library. If this is a number, it is additionally used as an RNG seed.
#' @param force.rooted Force input trees to be treated as rooted, even if there is a polytomy at the root.
#' @export
#' @examples
#' genome.trees.dist()

genome.trees.dist <- function(treeseq.a=NA, treeseq.b=NA, output.full.table = FALSE, acceptable.length.diff.pct = 0.1, variant.positions=NULL, randomly.resolve.a=FALSE, randomly.resolve.b=FALSE, force.rooted=TRUE) { 
    results=data.frame(unchanged.tree=numeric(), lft=numeric(), rgt=numeric(), RFrooted=numeric(), RFunrooted=numeric(),
        wRFrooted=numeric(), wRFunrooted=numeric(), SPRunrooted=numeric(), pathunrooted=numeric(), KCrooted=numeric())
    
    if (is.na(treeseq.a) || is.na(treeseq.b)) {
        #one or other is null, so simply return the column names for reference
        if (output.full.table)
            return(results)
        else
            return(results[-1:-2])
    } else {
        if (identical(randomly.resolve.a,FALSE)) {
            process.a = identity
        } else {
            if (is.numeric(randomly.resolve.a))
                set.seed(randomly.resolve.a)
            process.a = multi2di
        }

    
        if (identical(randomly.resolve.b, FALSE)) {
            process.b = identity
        } else {
            if (is.numeric(randomly.resolve.b))
                set.seed(randomly.resolve.b)
            process.b = multi2di
        }
    
        #check if trees or filenames
        if (class(treeseq.a) != "multiPhylo") {
            multiphy <- read.nexus(treeseq.a, force.multi=TRUE)
            if (force.rooted)
                for (i in seq_along(multiphy)) 
                    multiphy[[i]]$root.edge <- 0
            a <- process.a(multiphy)
        } else {
            a <- process.a(treeseq.a)
        }
        if (class(treeseq.b) != "multiPhylo") {
            multiphy <- read.nexus(treeseq.b, force.multi=TRUE)
            if (force.rooted)
                for (i in seq_along(multiphy)) 
                    multiphy[[i]]$root.edge <- 0
            b <- process.b(multiphy)
        } else {
            b <- process.b(treeseq.b)
        }
        catchTreeDistErrors <- function(f, metric="", rooted=FALSE) {
            #allow processing to continue even if there are errors
            tryCatch(withCallingHandlers(f, 
                message=function(cond){
                    cond$message <- paste("For", metric, "metric (rooted = ", rooted,"): ", cond$message, "\n")
                    message(cond)
                }),
                error=function(cond){
                    cond$message <- paste("For", metric, "metric (rooted = ", rooted,"): ", cond$message, "\n")
                    message(cond)
                })
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
        breaks.table <- stack(list('1'=brk.a,'2'=brk.b)) #set a:index=1 and b:index=2
        breaks.table$ind <- as.numeric(as.character(breaks.table$ind)) #convert to numeric
        breaks.table <- by(breaks.table, breaks.table$values, function(x) x) #put identical breakpoints together
        tree.index.ctr=c(1, 1)
        lft <- 0
        for (o in order(as.numeric(names(breaks.table)))) {
            brk = breaks.table[[o]]
            if (any(tree.index.ctr[brk$ind] > c(length(a),length(b))[brk$ind])) {
                warning("Reached the end of the trees with ", max(brk.a, brk.b)-lft, " left to go")
                break
            }
            rgt <- brk$values[0:1];
            RFrooted <- RFunrooted <- wRFrooted <- wRFunrooted <- SPRunrooted  <- pathunrooted <- KCrooted  <- NA
            catchTreeDistErrors({RFrooted <- phangorn::RF.dist(a[[tree.index.ctr[1]]], b[[tree.index.ctr[2]]], rooted=TRUE)},
                'RF', rooted=TRUE)
            catchTreeDistErrors({RFunrooted <- phangorn::RF.dist(a[[tree.index.ctr[1]]], b[[tree.index.ctr[2]]], rooted=FALSE)},
                'RF')
            catchTreeDistErrors({wRFrooted <- phangorn::wRF.dist(a[[tree.index.ctr[1]]], b[[tree.index.ctr[2]]], rooted=TRUE)},
                'weighted RF', rooted=TRUE)
            catchTreeDistErrors({wRFunrooted <- phangorn::wRF.dist(a[[tree.index.ctr[1]]], b[[tree.index.ctr[2]]], rooted=FALSE)},
                'weighted RF')
            catchTreeDistErrors({SPRunrooted <- phangorn::SPR.dist(a[[tree.index.ctr[1]]], b[[tree.index.ctr[2]]])},
                'subtree prune & regraft')
            catchTreeDistErrors({pathunrooted <- phangorn::path.dist(a[[tree.index.ctr[1]]], b[[tree.index.ctr[2]]])},
                'path distance')
            catchTreeDistErrors({KCrooted <- kc.dist(a[[tree.index.ctr[1]]], b[[tree.index.ctr[2]]])},
                'Kendall-Colijn', rooted=TRUE)
            results[nrow(results)+1,] <- c(ifelse(length(brk$ind)>1,NA,setdiff(1:2,brk$ind)),lft,rgt,RFrooted, RFunrooted, wRFrooted, wRFunrooted, SPRunrooted, pathunrooted, KCrooted)
            lft <- rgt
            tree.index.ctr[brk$ind] <- tree.index.ctr[brk$ind] + 1 #NB, brk$ind is a factor with levels (m1,m2), so we hope that m1==1 and m2==2
        }
    }
    if (output.full.table) {
        return(results)
    } else {
        if (is.null(variant.positions))
            weights <- results$rgt-results$lft
        else
            weights <- tabulate(findInterval(variant.positions, results$rgt, rightmost.closed=TRUE)+1, nbins=nrow(results))
        return(c(RFrooted=weighted.mean(results$RFrooted,        weights, na.rm = TRUE),
                 RFunrooted=weighted.mean(results$RFunrooted,    weights, na.rm = TRUE),
                 wRFrooted=weighted.mean(results$wRFrooted,      weights, na.rm = TRUE),
                 wRFunrooted=weighted.mean(results$wRFunrooted,  weights, na.rm = TRUE),
                 SPRunrooted=weighted.mean(results$SPRunrooted,  weights, na.rm = TRUE),
                 pathunrooted=weighted.mean(results$pathunrooted,weights, na.rm = TRUE),
                 KCrooted=weighted.mean(results$KCrooted,        weights, na.rm = TRUE)))
    }
}
