import os
import logging
import warnings
import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects.packages import importr

# Suppress noisy warnings from R.
warnings.simplefilter("ignore", rinterface.RRuntimeWarning)

ape=importr("ape")
importr("phangorn")
ARGmetrics = importr("ARGmetrics")
assert robjects.r('packageVersion("ARGmetrics") >= "0.0.1.0"')[0], (
    "Your version of ARGmetrics installed in R is too old (requires >= 0.0.1.0). "
    'Install the latest version from within R by syncing with git and typing e.g.\n'
    '> library(devtools)\n'
    '> setwd("{}")\n'
    '> install("ARGmetrics")\n'
    .format(os.path.join(os.path.basename(__file__),'..')))
#NB the above requires your version of R to have the bespoke ARGmetrics package installed
#Open R and set the working dir to treeseq-inference, then do
# > library(devtools)
# > install("ARGmetrics")


def get_metric_names():
    """
    Returns the list of the names of the computed metrics.
    """
    # We could do it with :
    # return list(pandas.DataFrame(columns=ARGmetrics.genome_trees_dist().names))
    # but it's extremely slow. Just return the list of strings instead.
    return list(ARGmetrics.genome_trees_dist().names)


def get_metrics(true_nexus_fn, inferred_nexus_fns, variant_positions = None, randomly_resolve_inferred=False):
    """
    Returns a dictionary of metrics for the specified pair of nexus files.
    :param str true_nexus_fn: The path to the file to write the SVG. If None, do not
        write to file.
    :param str (or array of str) inferred_nexus_fns: Other filenames.
    :param list variant_positions: The height of the image in pixels.
    :param int randomly_resolve_inferred: If Falsey do not randomly resolve polytomies.
    :return: A dictionary giving the average tree stats, whose keys are the method
        names as returned by get_metric_names().
    :rtype: dict

    """
    logging.debug("get_ARG_metrics() is comparing {} against {}{}".format(
        true_nexus_fn, inferred_nexus_fns,
        ' randomly breaking {} polytomies before comparison'.format(
        randomly_resolve_inferred if isinstance(randomly_resolve_inferred, (int,)) else ""
        ) if randomly_resolve_inferred else ''))
    if variant_positions is None:
        variant_positions  = rinterface.NULL
    if isinstance(inferred_nexus_fns, str):
        orig_tree = ape.read_nexus(true_nexus_fn, force_multi=True)
        inferred_tree = ape.read_nexus(inferred_nexus_fns, force_multi=True)
        m = ARGmetrics.genome_trees_dist(
            orig_tree, inferred_tree, variant_positions=variant_positions)
    else:
        # load the true_nexus into the R session (but don't convert it to a python obj)
        orig_tree = ape.read_nexus(true_nexus_fn, force_multi=True)
        m = ARGmetrics.genome_trees_dist_multi(
                orig_tree, inferred_nexus_fns, variant_positions=variant_positions, 
                weights=1, randomly_resolve_multi = randomly_resolve_inferred)
    return dict(m.items())

def get_full_metrics(true_nexus_fn, inferred_nexus_fn, variant_positions = rinterface.NULL):
    """
    Returns the full metric array (several metrics for each tree) for the specified pair of nexus files
    """
    logging.debug("get_ARG_metrics() is comparing {} against {}".format(
        true_nexus_fn, inferred_nexus_fn))
    if variant_positions is None:
        variant_positions  = rinterface.NULL
    # load the true_nexus into the R session (but don't convert it to a python obj)
    orig_trees = ape.read_nexus(true_nexus_fn, force_multi=True)
    infer_trees = ape.read_nexus(inferred_nexus_fn, force_multi=True)
    m = ARGmetrics.genome_trees_dist(
            orig_trees, infer_trees, output_full_table=True, variant_positions=variant_positions)
    return m


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
