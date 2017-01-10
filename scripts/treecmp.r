#script to take two multiphylo objects, whose names can be converted to floating
# point integers, denoting the end position on a genome of each tree in the file



setwd("/Users/yan/Documents/Research/Wellcome/treeseq-inference/test_files")
new.read.nexus <- function (file, tree.names = NULL, force.multi=FALSE) 
{
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
library(ape) #for read.nexus
unlockBinding("read.nexus", as.environment("package:ape"))
assign("read.nexus", new.read.nexus, pos=as.environment("package:ape"))
lockBinding("read.nexus", as.environment("package:ape"))



genome.trees.dist <- function(a, b, output_full_table = FALSE, acceptable_length_diff_pct = 0.1, rooted=FALSE, variant.positions=NULL) { 
    #a and b should be multiPhylo objects containing multiple trees
    #if variant.positions is given, it should be a vector of genome positions for each variants, and the metric will be compared per variant site
    #rather than per genomic region
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
    results=data.frame(lft=numeric(), rgt=numeric(), RF=numeric(), wRF=numeric(), SPR=numeric(), path=numeric())
    for (o in order(as.numeric(names(breaks.table)))) {
        brk = breaks.table[[o]]
        if (any(tree.index.counter[brk$ind] > c(length(a),length(b))[brk$ind])) {
            warning("Reached the end of the trees with ", max(brk.a, brk.b)-lft, " left to go")
            break
        }
        rgt <- brk$values[0:1]
        RF <- RF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=rooted)
        wRF <- wRF.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]], rooted=rooted)
        if (rooted==FALSE) {
            SPR <- SPR.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]])
            path <- path.dist(a[[tree.index.counter['a']]], b[[tree.index.counter['b']]])
        } else {
            SPR <- path <- NA
        }
        results[nrow(results)+1,] <- c(lft,rgt,RF,wRF, SPR, path)
        lft <- rgt
        tree.index.counter[brk$ind] <- tree.index.counter[brk$ind] + 1 #NB, brk$ind is a factor with levels (m1,m2), so we hope that m1==1 and m2==2
    }
    if (output_full_table) {
        return(results)
    } else {
        return(c(RF=weighted.mean(results$RF, (results$rgt-results$lft)),
                 wRF=weighted.mean(results$wRF, (results$rgt-results$lft)),
                 SPR=weighted.mean(results$SPR, (results$rgt-results$lft)),
                 path=weighted.mean(results$path, (results$rgt-results$lft))))
    }
}

test <- function() {
    #here we should take some trees of known topology distance and place them in the measure.genome.trees format, to check we are calculating correctly.
}

sim.params.extract <- function(filename) {
    #extracts a string like "-n100_Ne10000_l100000_rc0.00000002_mu0.00000002-" from the filename
    #extract the last part of the filename delimited by '-'
    p <- sub("^[^-]*-([^-]+)-.*", "\\1", filename)
    p <- strsplit(p, "_")[[1]]
    params <- as.numeric(gsub("[^0-9.]","",p))
    names(params) <- gsub("[0-9.]","",p)
    return(params)
}

sim.seeds.extract <- function(filename) {
    p <- sub("^[^-]*-[^-]+-([0-9_]+).*", "\\1", filename)
    seeds <- as.numeric(strsplit(p, "_")[[1]])
    names(seeds) <- paste0('simseed', 1:length(seeds))
    return(seeds)
}

fa.name.extract <- function(filename) {
    p <- sub("^[^+]+[+][^+]+[+](.*)$", "\\1", filename)
    vals <- as.numeric(p)
    #fastarg has a single seed number at the end
    return(c("fastarg_seed"=vals))
}

aw.samples.extract <- function(filename) {
    p <- sub("^[^+]+[+][^+]+[+](.*)$", "\\1", filename, perl=TRUE)
    #argweaver has seed.iter
    vals <- as.numeric(strsplit(p, ".")[[1]])
    #fastarg has a single number as the seed
    names(vals) <- c("argweaver_seed", "iteration")
    return(vals)
}

fa.aw.compare <- function(original.simulation.files="tmp/nexus_files/msprime*.nex", save.tables=TRUE) {
    fa.gmetrics <- data.frame()
    fa.vmetrics <- data.frame()
    dr <- dirname(original.simulation.files)
    dr <- ifelse(dr=='', dr, paste0(dr,"/"))
    found.files=FALSE
    for (path1 in Sys.glob(original.simulation.files)) {
        found.files=TRUE
        #each of these should have an equivalent set of fastarg-XXXXX files
        fn <- sub("\\.\\w+$", "", basename(path1), perl=TRUE)
        orig.nex <- read.nexus(path1, force.multi=TRUE)
        sim.params <- sim.params.extract(fn)
        sim.seeds <- sim.seeds.extract(fn)
        vpos.file = sub("\\.\\w+$", ".vpos", path1, perl=TRUE)
        if (file.exists(vpos.file)) {
            variant.positions <- scan(vpos.file)
        } else {
            variant.positions <- NULL
        }
        #look at fastarg files first
        for (path2 in Sys.glob(paste0(dr, 'fastarg+',fn,'*.nex'))) {
            fastarg.file <- sub("\\.\\w+$", "", basename(path2), perl=TRUE)
            print(fastarg.file)
            #inference.files should contain all the results for inferences based on these msprime sims
            fa.nex <- read.nexus(path2, force.multi=TRUE)
           
            #there may also be a .stats file
            metrics <- genome.trees.dist(orig.nex, fa.nex)
            fa.gmetrics <- rbind(fa.gmetrics, as.list(c(sim.params, sim.seeds, fa.name.extract(basename(fastarg.file)), metrics)))
            if (length(variant.positions)) {
                metrics <- genome.trees.dist(orig.nex, fa.nex, variant.positions=variant.positions)
                fa.vmetrics <- rbind(fa.gmetrics, c(sim.params, sim.seeds, fa.name.extract(basename(fastarg.file)), metrics))
           }
        }
    }
    if (found.files==FALSE) {
        warning(paste0("No files found in ", getwd(), "/", original.simulation.files))
    } else {
        fa.gmetrics=fa.gmetrics[order(fa.gmetrics$simseed1, fa.gmetrics$mu)
        fa.vmetrics=fa.vmetrics[order(fa.vmetrics$simseed1, fa.vmetrics$mu)
        if (save.tables) {
            write.table(fa.gmetrics, file=paste0(dr,'fastarg_metrics_by_genome.txt'))
            write.table(fa.vmetrics, file=paste0(dr,'fastarg_metrics_by_variant.txt'))
        }
    }
    return(list(fa.gmetrics=fa.gmetrics, fa.vmetrics=fa.vmetrics))
}

metrics <- fa.aw.compare()

equipment

#plot(as.numeric(sub("muts\\d+", "", rownames(gpos_metrics))), gpos_metrics$RF)
#points(as.numeric(sub("muts\\d+", "", rownames(mpos_metrics))), mpos_metrics$RF, colors="red")

#aw.res <- genome.trees.dist(orig_trees, argweaver_trees)