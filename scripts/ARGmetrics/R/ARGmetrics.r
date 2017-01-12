#' Replacement for read.nexus in the ape library
#'
#' This should be identical to read.nexus, but with a parameter force.multi
#' to force the return of a multiPhylo object, even if there is only one
#' tree in the nexus file
#' @param file A path to the nexus file.
#' @param tree.names Override or provide names for the trees in the file.
#' @param force.multi Return a multiPhylo object, even if there is only one tree in the file.
#' @export
#' @examples
#' See ?read.nexus

new.read.nexus <- function (file, tree.names = NULL, force.multi=FALSE) 
{
    require(ape) #for read.nexus
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) {
        w <- LEFT == RIGHT
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        }
        w <- !w
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        }
    }
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- if (length(i2) == 1 && i2 > i1) 
        TRUE
    else FALSE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end]
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        n <- dim(TRANS)[1]
    }
    start <- if (translation) 
        semico[semico > i2][1] + 1
    else i1 + 1
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    rm(X)
    tree <- tree[tree != ""]
    semico <- grep(";", tree)
    Ntree <- length(semico)
    if (Ntree == 1 && length(tree) > 1) 
        STRING <- paste(tree, collapse = "")
    else {
        if (any(diff(semico) != 1)) {
            STRING <- character(Ntree)
            s <- c(1, semico[-Ntree] + 1)
            j <- mapply(":", s, semico)
            if (is.list(j)) {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], 
                  collapse = "")
            }
            else {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[, 
                  i]], collapse = "")
            }
        }
        else STRING <- tree
    }
    rm(tree)
    STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
    Ntree <- length(STRING)
    nms.trees <- sub(" *= *.*", "", STRING)
    nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", 
        nms.trees, ignore.case = TRUE)
    STRING <- sub("^.*= *", "", STRING)
    STRING <- gsub(" ", "", STRING)
    colon <- grep(":", STRING)
    if (!length(colon)) {
        trees <- lapply(STRING, clado.build)
    }
    else if (length(colon) == Ntree) {
        trees <- if (translation) 
            lapply(STRING, ape:::.treeBuildWithTokens)
        else lapply(STRING, tree.build)
    }
    else {
        trees <- vector("list", Ntree)
        trees[colon] <- lapply(STRING[colon], tree.build)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        trees[nocolon] <- lapply(STRING[nocolon], clado.build)
        if (translation) {
            for (i in 1:Ntree) {
                tr <- trees[[i]]
                for (j in 1:n) {
                  ind <- which(tr$tip.label[j] == TRANS[, 1])
                  tr$tip.label[j] <- TRANS[ind, 2]
                }
                if (!is.null(tr$node.label)) {
                  for (j in 1:length(tr$node.label)) {
                    ind <- which(tr$node.label[j] == TRANS[, 
                      1])
                    tr$node.label[j] <- TRANS[ind, 2]
                  }
                }
                trees[[i]] <- tr
            }
            translation <- FALSE
        }
    }
    for (i in 1:Ntree) {
        tr <- trees[[i]]
        if (!translation) 
            n <- length(tr$tip.label)
        ROOT <- n + 1
        if (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] > 
            1) {
            stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading NEXUS file aborted at tree no.", 
                i, sep = ""))
        }
    }
    if ((Ntree == 1) && (force.multi == FALSE)) {
        trees <- trees[[1]]
        if (translation) {
            trees$tip.label <- if (length(colon)) 
                TRANS[, 2]
            else TRANS[, 2][as.numeric(trees$tip.label)]
        }
    }
    else {
        if (!is.null(tree.names)) 
            names(trees) <- tree.names
        if (translation) {
            if (length(colon) == Ntree) 
                attr(trees, "TipLabel") <- TRANS[, 2]
            else {
                for (i in 1:Ntree) trees[[i]]$tip.label <- TRANS[, 
                  2][as.numeric(trees[[i]]$tip.label)]
                trees <- .compressTipLabel(trees)
            }
        }
        class(trees) <- "multiPhylo"
        if (!all(nms.trees == "")) 
            names(trees) <- nms.trees
    }
    trees
}

#' Compare multiple trees across a genomic region
#'
#' Requires multiPhylo objects representing sequential trees along a genomic region.
#' The tree names should be numbers, which represent the rightmost position along a
#' genome (non-inclusive) covered by each tree.
#' @param a The first multiPhylo object
#' @param b The second multiPhylo object
#' @param output.full.table Output tree metrics for each overlapping region, rather than simply a weighted summary.
#' @param acceptable.length.diff.pct How much difference in sequence length is allows between the 2 trees? (Default: 0.1%)
#' @param variant.positions A list of positions of each variant (not implemented)
#' @export
#' @examples
#' genome.trees.dist()
genome.trees.dist <- function(a, b, output.full.table = FALSE, acceptable.length.diff.pct = 0.1, variant.positions=NULL) { 
    #a and b should be multiPhylo objects containing multiple trees
    #if variant.positions is given, it should be a vector of genome positions for each variants, and the metric will be compared per variant site
    #rather than per genomic region
    require(phangorn) #to use the various treedist metrics
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
    if (output_full_table) {
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
