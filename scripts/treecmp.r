#script to take two multiphylo objects, whose names can be converted to floating
# point integers, denoting the end position on a genome of each tree in the file
# the
setwd("/Users/yan/Documents/Research/Wellcome/treeseq-inference/test_files")
library(ape) #for read.nexus
base = "500-1000000"

orig_trees <- read.nexus(paste(base,".nex", sep=""))
fastarg_trees <- read.nexus(paste(base,".fa.nex", sep=""))
#argweaver_trees <- read.nexus(paste(base,".aw.nwk", sep=""))
#msprime_trees <- read.nexus(paste(base,".mp.nwk", sep=""))


genome.trees.dist <- function(a, b, output_full_table = FALSE, acceptable_length_diff_pct = 0.1) { #a and b should be multiPhylo objects containing multiple trees
    library(phangorn) #to use the various treedist metrics
    brk.a <- as.numeric(names(a))
    if (is.unsorted(brk.a))
        stop("Tree names should correspond to numerical breakpoints or variant totals, sorted from low to high, but tree names in the first trees object are not numbers in ascending order.")
    brk.b <- as.numeric(names(b))
    if (is.unsorted(brk.b))
        stop("Tree names should correspond to numerical breakpoints or variant totals, sorted from low to high, but tree names in the second trees object are not numbers in ascending order.")
    if ((max(brk.a) * (100+ acceptable_length_diff_pct)/100 < max(brk.b)) || 
        (max(brk.b) * (100+ acceptable_length_diff_pct)/100 < max(brk.a)))
        warning("The sequence lengths of the two trees files differ markedly: ", max(brk.a), " vs. ", max(brk.b), immediate. = TRUE)
    breaks.table <- stack(list('a'=brk.a,'b'=brk.b))
    breaks.table <- by(breaks.table, breaks.table$values, function(x) x) #put identical breakpoints together
    tree.index.counter=c(a=1, b=1) 
    lft <- 0
    results=data.frame(lft=numeric(), rgt=numeric(), RF=numeric(), wRF=numeric())
    for (o in order(as.numeric(names(breaks.table)))) {
        brk = breaks.table[[o]]
        rgt <- brk$values[0:1]
        RF <- RF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=TRUE)
        wRF <- wRF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=TRUE)
        results[nrow(results)+1,] <- c(lft,rgt,RF,wRF)
        lft <- rgt
        tree.index.counter[brk$ind] <- tree.index.counter[brk$ind] + 1 #NB, brk$ind is a factor with levels (m1,m2), so we hope that m1==1 and m2==2
        if (any(tree.index.counter[brk$ind] > c(length(a),length(b))[brk$ind])) {
            warning("Reached the end of the trees with ", max(brk.a, brk.b)-lft, " left to go")
            break
        }
    }
    if (output_full_table) {
        return(results)
    } else {
        return(c(RF=weighted.mean(results$RF, (results$rgt-results$lft)), wRF=weighted.mean(results$wRF, (results$rgt-results$lft))))
    }
}

test <- function() {
    #here we should take some trees of know topology distance and place them in the measure.genome.trees format, to check we are calculating correctly.
}

fa.res <- genome.trees.dist(orig_trees, fastarg_trees)
weighted.mean(fa.res$RF, (fa.res$rgt-fa.res$lft))


aw.res <- genome.trees.dist(orig_trees, argweaver_trees)