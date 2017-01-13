#script to take two multiphylo objects, whose names can be converted to floating
# point integers, denoting the end position on a genome of each tree in the file



setwd("/Users/yan/Documents/Research/Wellcome/treeseq-inference/test_files")
library("ARGmetrics")
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
        orig.nex <- new.read.nexus(path1, force.multi=TRUE)
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
            fa.nex <- new.read.nexus(path2, force.multi=TRUE)
           
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
        #fa.gmetrics=fa.gmetrics[order(fa.gmetrics$simseed1, fa.gmetrics$mu)
        #fa.vmetrics=fa.vmetrics[order(fa.vmetrics$simseed1, fa.vmetrics$mu)
        if (save.tables) {
            write.table(fa.gmetrics, file=paste0(dr,'fastarg_metrics_by_genome.txt'))
            write.table(fa.vmetrics, file=paste0(dr,'fastarg_metrics_by_variant.txt'))
        }
    }
    return(list(fa.gmetrics=fa.gmetrics, fa.vmetrics=fa.vmetrics))
}

metrics <- fa.aw.compare()


#plot(as.numeric(sub("muts\\d+", "", rownames(gpos_metrics))), gpos_metrics$RF)
#points(as.numeric(sub("muts\\d+", "", rownames(mpos_metrics))), mpos_metrics$RF, colors="red")

#aw.res <- genome.trees.dist(orig_trees, argweaver_trees)