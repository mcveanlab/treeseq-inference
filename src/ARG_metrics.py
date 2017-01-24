import logging
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

ape=importr("ape")
assert robjects.r('packageVersion("ape") >= "4.0.0.2"')[0], \
    "You must install an 'ape' version in R > 4.0.0.2, e.g. by running" + \
    'install.packages("ape", repos = "http://ape-package.ird.fr/", type="source")'
importr("phangorn")
ARGmetrics = importr("ARGmetrics")
#NB the above requires your version of R to have the bespoke ARGmetrics package installed
#Open R and set the working dir to treeseq-inference, then do
# > install("ARGmetrics")

def get_ARG_metrics(true_nexus_fn=None, threads=1, **inferred_nexus_fns):
    """
    The complicated thing here is the inferred_nexus_fns parameters.
    Each should contain a dictionary with a 'nexus' item pointing to 
    one or more nexus files. Optionally, this dictionary can also have 
    keys named 'make_bin_seed' and 'reps', which specify the number
    of times to replicate the metric and the starting seed for replicates
    
    e.g. ARG_metrics(true.nex, tsinfer={'nexus':'ms.nex'}, argWeaver={'nexus':['aw1.nex', 'aw2.nex']})

    or (more complex):
    ARG_metrics(true.nex, 
                tsinfer={'nexus':'fa.nex', 'make_bin_seed'=123, 'reps':10}, 
                Aweaver={'nexus':['aw1.nex', 'aw2.nex']})
    
    Returns a Data Frame with rows labelled by the parameters (tool names) passed in
    ('tsinfer', 'ARGweaver', etc), and one column for each metric.
    
    If called with no params, simply returns a dummy data frame with the right 
    column names (currently this has a row of NaNs as the data)
    
    threads gives the number of subprocesses to fork() when doing extra random replicates
    during polytomy resolution. This needs testing to see if it provides a speedup
    
    """
    if true_nexus_fn is None:
        return pd.DataFrame(columns = ARGmetrics.genome_trees_dist().names)
    logging.debug("running get_ARG_metrics({},{},{})".format(true_nexus_fn,
        threads, inferred_nexus_fns))
    # load the true_nexus into the R session (but don't convert it to a python obj)
    orig_tree = ape.read_nexus(true_nexus_fn, force_multi=True)
    weights = 1
    metrics = None
    for tool, metric_params in inferred_nexus_fns.items():
        logging.info(" * calculating ARG metrics for {} on R process {}.".format(tool,
            robjects.r("Sys.getpid()")[0]))
        nexus_files = metric_params['nexus']
        if isinstance(nexus_files, str):
            try:
                break_binary_reps = int(metric_params['reps'])
                seed = metric_params.get('make_bin_seed')
                logging.debug("calculating tree metrics to compare binarised '{}' vs '{}'.".format(
                    true_nexus_fn, nexus_files))
                if seed:
                
                    m = ARGmetrics.genome_trees_dist_forcebin_b(orig_tree, nexus_files, 
                         replicates = break_binary_reps, seed = int(seed), threads=threads)
                else:
                    m = ARGmetrics.genome_trees_dist_forcebin_b(orig_tree, nexus_files, 
                         replicates = break_binary_reps, threads=threads)
            except:
                logging.debug("calculating tree metrics to compare '{}' vs '{}'.".format(
                    true_nexus_fn, nexus_files))
                m = ARGmetrics.genome_trees_dist(orig_tree, nexus_files)
        else:
            if len(nexus_files) == 0:
                logging.debug("No inference files give to compare")
                m = ARGmetrics.genome_trees_dist()
            else:
                logging.debug("calculating tree metrics to compare '{}' vs '{}'.".format(
                    true_nexus_fn, nexus_files))
                m = ARGmetrics.genome_trees_dist_multi(orig_tree, nexus_files, weights=1)
        if metrics is None:
            metrics = pd.DataFrame(columns = m.names)
        metrics.loc[tool, tuple(m.names)] = tuple(m)
    return metrics