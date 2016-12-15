#script to take two multiphylo objects, whose names can be converted to floating
# point integers, denoting the end position on a genome of each tree in the file
# the
library(phangorn) #to use the various treedist metrics
setwd("/Users/yan/Documents/Research/Wellcome/treeseq-inference/test_files")

base = "500-1000000"

orig_trees <- read.nexus(paste(base,".nex", sep=""))
fastarg_trees <- read.nexus(paste(base,".fa.nex", sep=""))
#argweaver_trees <- read.nexus(paste(base,".aw.nwk", sep=""))
#msprime_trees <- read.nexus(paste(base,".mp.nwk", sep=""))


tree.measure <- function(a, b, acceptable_length_diff_pct = 0.1) { #a and b should be multiPhylo objects containing multiple trees
    brk.a <- as.numeric(names(a))
    if (is.unsorted(brk.a))
        stop("Tree names should correspond to numerical breakpoints, sorted from low to high, but trees names in the first trees object are not sorted.")
    brk.b <- as.numeric(names(b))
    if (is.unsorted(brk.b))
        stop("Tree names should correspond to numerical breakpoints, sorted from low to high, but trees names in the second trees object are not sorted.")
    if ((max(brk.a) * (100+ acceptable_length_diff_pct)/100 > max(brk.b)) || 
        (max(brk.b) * (100+ acceptable_length_diff_pct)/100 > max(brk.a)))
        warning("The sequence lengths of the two trees files differ markedly: ", max(brk.a), " vs. ", max(brk.b), immediate. = TRUE)
    breaks.table <- stack(list('a'=brk.a,'b'=brk.b))
    tree.index=c(a=1, b=1) 
    lft <- 0
    results=data.frame(lft=numeric(), rgt=numeric(), RF=numeric(), wRF=numeric())
    for (o in order(breaks.table$values)) {
        brk = breaks.table[o,]
        rgt <- brk$values
        RF <- RF.dist(a[[tree.index['a']]], b[[tree.index['b']]], rooted=TRUE)
        wRF <- wRF.dist(a[[tree.index['a']]], b[[tree.index['b']]], rooted=TRUE)
        results[nrow(results)+1,] <- c(lft,rgt,RF,wRF)
        lft <- rgt
        tree.index[brk$ind] <- tree.index[brk$ind] + 1 #NB, brk$ind is a factor with levels (m1,m2), so we hope that m1==1 and m2==2
        if (tree.index[brk$ind] > c(length(a),length(b))[brk$ind]) {
            warning("Reached the end of the trees with ", max(brk.a, brk.b)-lft, " left to go")
            break
        }
    }
    return(results)
}

fa.res <- tree.measure(orig_trees, fastarg_trees)
aw.res <- tree.measure(orig_trees, argweaver_trees)