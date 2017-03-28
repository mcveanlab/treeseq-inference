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
# > library(devtools)
# > install("ARGmetrics")


def get_metric_names():
    """
    Returns the list of the names of computed metrics.
    """
    # We could do it with :
    # return list(pd.DataFrame(columns=ARGmetrics.genome_trees_dist().names))
    # but it's extremely slow. Just return the list of strings instead.
    return [
        'RFrooted', 'RFunrooted', 'wRFrooted', 'wRFunrooted', 'SPRunrooted',
        'pathunrooted', 'KCrooted']


def get_metrics(true_nexus_fn, inferred_nexus_fn):
    """
    Returns a dictionary of metrics for the specified pair of nexus files.
    """
    logging.debug("get_ARG_metrics() is comparing {} against {}".format(
        true_nexus_fn, inferred_nexus_fn))
    # load the true_nexus into the R session (but don't convert it to a python obj)
    orig_tree = ape.read_nexus(true_nexus_fn, force_multi=True)
    var_pos = rinterface.NULL
    weights = 1
    m = ARGmetrics.genome_trees_dist_multi(
            orig_tree, [inferred_nexus_fn], variant_positions=var_pos, weights=1)
    return dict(m.items())


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
