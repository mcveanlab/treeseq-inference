"""
Example client code for the tsinf module.
"""
import tsinf
import msprime
import numpy as np


def main():
    rho = 2
    ts = msprime.simulate(10, mutation_rate=10, recombination_rate=rho)
    S = np.zeros((ts.sample_size, ts.num_mutations), dtype="u1")
    for variant in ts.variants():
        S[:,variant.index] = variant.genotypes
    panel = tsinf.ReferencePanel(S)
    P = panel.infer_paths(rho, num_workers=1)
    ts_new = panel.convert_records(P)
    ts_simplified = ts_new.simplify()

    # Check that we get the same variants.
    assert ts.num_mutations == ts_simplified.num_mutations
    for v1, v2 in zip(ts.variants(), ts_simplified.variants()):
        assert np.all(v1.genotypes == v2.genotypes)


if __name__ == "__main__":
    main()
