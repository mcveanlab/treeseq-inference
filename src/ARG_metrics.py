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
    e.g. ARG_metrics(true.nex, fastARG='fa.nex', argWeaver=['aw1.nex', 'aw2.nex'])

    Where there are multiple nex files for a given method, we should also allow 
    different weights to be passed in, to provide a weighted average of metrics
    over the files. Yet to be coded***
    
    Returns a dictionary of Data Frames, one for each nexus file passed in
    """
    #make sure that the
    # load the true_nexus into the R session
    orig_tree = ARGmetrics.new_read_nexus(true_nexus_fn)
    metrics={}
    for tool, nexus_files in inferred_nexus_fns.items():
        if isinstance(nexus_files, str):
            m = ARGmetrics.genome_trees_dist(orig_tree, nexus_files)
        else:
            m = ARGmetrics.genome_trees_dist_multi(orig_tree, nexus_files, weights=1)
        metrics[tool]=pd.DataFrame(data=[tuple(m)], columns=m.names)
