#!/usr/bin/env python3
description = """
Use native methods implemented in msprime by Peter Ralph et al to display measures of selection across a genomic region, using the tree-based methods available in msprime.
"""
import os
import sys
import logging
from fractions import Fraction

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot


#set up for local imports
curr_dir = os.path.dirname(os.path.abspath(__file__))
# import the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(curr_dir,'..','msprime'))
import msprime


def H(n):
    """Returns the n-th harmonic number (or an approximation for large n)

       http://en.wikipedia.org/wiki/Harmonic_number
    """
    if n<100:
        return sum(Fraction(1, d) for d in xrange(1, n + 1))
    else:
        return np.euler_gamma + np.log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)


def Theta_W(sfs, n_samples):
    """
    Calculate the (mutationless) Watterson estimator
    """
    return sum(sfs) / H(n_samples)

def Theta_pi(sfs, n_samples):
    """
    Calculate the (mutationless) nucleotide diversity
    2 sum(i(n−i)*ξi) / n / (n-1)
    """
    assert len(sfs)==n_samples, (sfs, len(sfs), n_samples)
    return 2 * sum([i*(n_samples-i)*sfs[i] for i in range(len(sfs))]) / n_samples / (n_samples-1)


def plot_tajimas_D(tree_sequences):
    """
    pass this function an interator over a list of tree sequences, and it will average over all of them
    """
    def TajimasD(sfs, n_samples):
        print(sfs)
        return Theta_pi(sfs, n_samples)-Theta_W(sfs, n_samples)
    tD = None # will contain the DatFrame summarizing Tajima's D values across the genome 
    for ts in tree_sequences:
        if D is None:
            D = pd.DataFrame(0, index=np.arange(len(data)), columns=("n","sum","sum_sq"))
        elif ts.get_sequence_length() != D.shape[0]:
            logging.warning("Tree sequences have different lengths")
        branch_tsc = msprime.stats.BranchLengthStatCalculator(ts)
        all_samples = ts.samples()
        branch_tsc.tree_stat_vector(
            sample_sets=[all_samples],
            weight_fun =lambda x: TajimasD(x,len(all_samples)), 
            windows=list(range(int(ts.get_sequence_length())+1))) #calculate for each BP
        
        
def main()
    import argparse
    
    from matplotlib.backends.backend_pdf import PdfPages
    
    parser = argparse.ArgumentParser(
        description=description
    )

    parser.add_argument(
        'files', nargs='+',
        help='Any number of hdf5 files containing tree sequences to be averaged'
    )

    parser.add_argument(
        '--focal-file-number', '-ff', type=int, default=1,
        help=''
    )

    parser.add_argument(
        '--output-file', '-o', type=argparse.FileType('w'), default=sys.stdout,
        help=''
    )

    parser.add_argument('--verbosity', '-v', action='count', default=0)

    args = parser.parse_args()

    log_level = logging.WARNING
    if args.verbosity == 1:
        log_level = logging.INFO
    if args.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(
        format='%(asctime)s %(message)s', level=log_level, stream=sys.stderr)


    with PdfPages(args.output_file) as pdf:
        plot()
        pdf.savefig()
        #plot
        #pdf.savefig()
        #plot
        #pdf.savefig()
        

if __name__ == "__main__":
    main()