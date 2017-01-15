import logging
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
importr("ape")
importr("phangorn")
ARGmetrics = importr("ARGmetrics")
#NB the above requires your version of R to have the bespoke ARGmetrics package installed
#Open R and set the working dir to treeseq-inference, then do
# > install("ARGmetrics")

def get_ARG_metrics(true_nexus_fn, **inferred_nexus_fns):
    """
    Pass in a set of inferred nexus files (each of which could be an array)
    e.g. ARG_metrics(true.nex, msinfer='ms.nex', argWeaver=['aw1.nex', 'aw2.nex'])

    or (more complex):
     ARG_metrics(true.nex, msinfer={'nexus':'fa.nex', 'make_bin_seed'=123, 'reps':10}

    or (even worse)
     ARG_metrics(true.nex, aw ={'nexus':['aw1.nex', 'aw2.nex'], 'weights':[1.2,1.3]})

    
    Returns a dictionary of one-line Data Frames, one for each tool passed in
    """
    #make sure that the
    # load the true_nexus into the R session
    orig_tree = ARGmetrics.new_read_nexus(true_nexus_fn)
    weights = 1
    metrics ={}
    for tool, metric_params in inferred_nexus_fns.items():
        try:
            nexus_files = metric_params['nexus']
            weights = metric_params('weights') or 1
        except:
            nexus_files = metric_params
        if isinstance(nexus_files, str):
            try:
                break_binary_reps = int(metric_params['reps'])
                seed = metric_params.get('make_bin_seed')
                logging.debug("calculating tree metrics to compare binarised '{}' vs '{}'.".format(
                    true_nexus_fn, nexus_files))
                if seed:
                    m = ARGmetrics.genome.trees.dist.forcebin.b(orig_tree, nexus_files, 
                         replicates = break_binary_reps, seed = int(seed))
                else:
                    m = ARGmetrics.genome.trees.dist.forcebin.b(orig_tree, nexus_files, 
                         replicates = break_binary_reps)
            except:
                logging.debug("calculating tree metrics to compare '{}' vs '{}'.".format(
                    true_nexus_fn, nexus_files))
                m = ARGmetrics.genome_trees_dist(orig_tree, nexus_files)
        else:
            logging.debug("calculating tree metrics to compare '{}' vs '{}'.".format(
                true_nexus_fn, nexus_files))
            m = ARGmetrics.genome_trees_dist_multi(orig_tree, nexus_files, weights=weights)
        metrics[tool]=dict(zip(m.names, tuple(m)))
    return metrics