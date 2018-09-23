
import sys
import simplebgen
import bgen_reader
import numpy as np


def check():
    filename = sys.argv[1]

    bgen = bgen_reader.read_bgen(filename, verbose=False)
    num_variants = len(bgen["variants"])
    num_samples = len(bgen["samples"])

    br = simplebgen.BgenReader(filename)
    assert br.num_variants == num_variants
    assert br.num_samples == num_samples

    for j in range(num_variants):
        p1  = br.get_probabilities(j)
        p2 = np.array(bgen["genotype"][j].compute())
        np.testing.assert_equal(p1, p2)

    print("all good")


def memory():
    
    filename = sys.argv[1]
    for j in range(1000):
        print("Iter ", j)
        br = simplebgen.BgenReader(filename)
        for j in range(br.num_variants):
            p1  = br.get_probabilities(j)


memory()

# check()
