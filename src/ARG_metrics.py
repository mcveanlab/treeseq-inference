import logging
import pandas as pd
import warnings
import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects.packages import importr

# Suppress noisy warnings from R.
warnings.simplefilter("ignore", rinterface.RRuntimeWarning)

ape=importr("ape")
assert robjects.r('packageVersion("ape") >= "4.0.0.2"')[0], \
    "You must install an 'ape' version in R > 4.0.0.2, e.g. by running" + \
    'install.packages("ape", repos = "http://ape-package.ird.fr/", type="source")'
importr("phangorn")
ARGmetrics = importr("ARGmetrics")
#NB the above requires your version of R to have the bespoke ARGmetrics package installed
#Open R and set the working dir to treeseq-inference, then do
# > install("ARGmetrics")


def get_metric_names():
    """
    Returns the list of the names of computed metrics.
    """
    return list(pd.DataFrame(columns=ARGmetrics.genome_trees_dist().names))


def get_ARG_metrics(true_nexus_fn=None, threads=1, variant_positions = None, **inferred_nexus_fns):
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

    threads gives the number of subprocesses to fork() when doing extra random replicates
    during polytomy resolution. This needs testing to see if it provides a speedup

    """
    logging.debug(
        "get_ARG_metrics() is comparing {} against inferred trees for the "
        "following tools (with file numbers) {}".format(
        true_nexus_fn, threads, {k: len(v['nexus']) for k,v in inferred_nexus_fns.items()}))
    # load the true_nexus into the R session (but don't convert it to a python obj)
    orig_tree = ape.read_nexus(true_nexus_fn, force_multi=True)
    var_pos = variant_positions or rinterface.NULL
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
                    m = ARGmetrics.genome_trees_dist_forcebin_b(orig_tree, nexus_files, seed = int(seed),
                         variant_positions=var_pos, replicates = break_binary_reps, threads=threads)
                else:
                    m = ARGmetrics.genome_trees_dist_forcebin_b(orig_tree, nexus_files,
                         variant_positions=var_pos, replicates = break_binary_reps, threads=threads)
            except:
                logging.debug("calculating tree metrics to compare '{}' vs '{}'.".format(
                    true_nexus_fn, nexus_files))
                m = ARGmetrics.genome_trees_dist(orig_tree, nexus_files, variant_positions=var_pos)
        else:
            if len(nexus_files) == 0:
                logging.debug("No inference files give to compare")
                m = ARGmetrics.genome_trees_dist()
            else:
                logging.debug("calculating tree metrics to compare '{}' vs {} other files.".format(
                    true_nexus_fn, len(nexus_files)))
                m = ARGmetrics.genome_trees_dist_multi(orig_tree, nexus_files,
                    variant_positions=var_pos, weights=1)
        if metrics is None:
            metrics = pd.DataFrame(columns = m.names)
        metrics.loc[tool, tuple(m.names)] = tuple(m)
    return metrics

if __name__ == "__main__":
    """Test whether we get sensible metrics for a tree sequence"""
    import tempfile

    ts = [{3:"(1,(2,3));", 4:"((1,2),3);"},
          {4:"((1,2),3);"}]
    #this should give an RF distance of 2 for the 1st tree and 0 for the second
    with tempfile.NamedTemporaryFile("w+") as nex1, \
        tempfile.NamedTemporaryFile("w+") as nex2:
        for fh, trees in zip([nex1, nex2], ts):
            print("#NEXUS\nBEGIN TREES;",file=fh)
            print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i,i) for i in [1,2,3]])), file = fh)
            for endpoint in sorted(trees.keys()):
                print("TREE " + str(endpoint) + " = " + trees[endpoint], file=fh)
            print("END;", file = fh)
            fh.flush()
        #Genome average should be (2*3 + 0*1) /4 = 6/4 = 1.5
        print(get_ARG_metrics(nex1.name, test={'nexus':nex2.name}))
        #Positional average for should be (2*2 + 0*2) /5 = 4/5 = 0.8
        print(get_ARG_metrics(nex1.name, variant_positions=[0, 2.9999,3,3.5, 4],test={'nexus':nex2.name}))
