"""
Implemenation of the Li and Stephens algorithm for inferring a tree
sequence.

Python 3 only.
"""

import sys
import time
import collections
import threading
import tempfile
import resource
import subprocess
import multiprocessing

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import pandas as pd
import seaborn as sns

import msprime
import _msprime
import _tsinf

if sys.version_info[0] < 3:
    raise Exception("Python 3 you idiot!")

def li_and_stephens_simulation(H, rho, theta):
    """
    Simulate new haplotypes from the Li and Stephens model with
    the specified haplotypes, recombination and mutation rate.
    """
    np.random.seed(1)
    print("rho = ", rho)
    print("theta = ", theta)
    print(H)
    n, m = H.shape
    P_mut =  0.5 * theta / (n + theta)
    P_recomb = (1 - np.exp(-rho / n)) / n
    while True:
        hp = np.zeros(m, dtype="u1")
        x = np.random.randint(0, n)
        for j in range(m):
            hp[j] = H[x][j]
            if np.random.random() < P_mut:
                hp[j] = (hp[j] + 1) % 2
                print("mut")
            if np.random.random() < P_recomb:
                # TODO: is this with or without replacement?
                x = np.random.randint(0, n)
                print("switch")
        print(hp)

def simple_infer_ancestors(H):
    """
    Given a set of input haplotypes infer their ancestors using frequency
    information. This is a simple unvectorised version suitable for implementation
    in C.
    """
    n, m = H.shape
    sites = []
    for l in range(m):
        samples = []
        for j in range(n):
            if H[j, l] == 1:
                samples.append(j)
        if len(samples) > 1:
            sites.append((len(samples), l, samples))
    sites.sort(key=lambda x: (-x[0], x[1]))
    A = np.zeros((len(sites) + 1, m), dtype="uint8")
    for j in range(len(sites)):
        freq, site, samples = sites[j]
        # For each older site, find the consensus among these samples.
        for k in range(j):
            _, other_site, _ = sites[k]
            f = 0
            for s in samples:
                f += H[s, other_site]
            if f >= len(samples) / 2:
                A[j + 1, other_site] = 1
        A[j + 1, site] = 1
    return A

def infer_ancestors(H):
    n, m = H.shape
    freq = np.sum(H, axis=0)
    non_singletons = np.sum(freq > 1)
    sites = sorted(
        [(freq[site], site) for site in range(m) if freq[site] > 1],
        key=lambda x: (-x[0], x[1]))
    mask = np.zeros(m, dtype=int)
    A = np.zeros((non_singletons + 1, m), dtype="uint8")
    for j, (_, site) in enumerate(sites):
        assert freq[site] > 1
        mask[site] = 1
        # Find the set of haplotypes that have a 1 at this site
        # and mask out any younger mutations.
        Hp = H[H[:,site] == 1] & mask
        # Find the majority allele at each locus among these haplotypes
        A[j + 1:] = np.sum(Hp, axis=0) >= Hp.shape[0] / 2
    return A


def thread(H, haplotype_index, n, rho):
    """
    Thread the haplotype with the specified index through the
    reference panel of the n oldest haplotypes.
    """
    h = H[haplotype_index]
    N, m = H.shape
    r = 1 - np.exp(-rho / n)
    recomb_proba = np.log(r / n)
    no_recomb_proba = np.log(1 - r + r / n)

    V = np.zeros(N)
    T = np.zeros((N, m), dtype=int)
    # The truth table for the emission probabilities. Index
    # (0, 1) corresponds to zero in the reference haplotype
    # and 1 in the sequence we are threading. The probablility
    # of mutation must be less than the probability of many
    # recombinations, or we will have multiple mutations.
    mutation_proba = recomb_proba * m;
    E = np.array([
        [0, mutation_proba],
        [-sys.float_info.max, 0]])

    max_V_index = N - n
    for j in range(N - n, N):
        V[j] = np.log(1 / n) + E[H[j, 0], h[0]]
        if V[j] > V[max_V_index]:
            max_V_index = j
    for l in range(1, m):
        Vn = np.zeros(N)
        max_Vn_index = N - n
        for j in range(N - n, N):
            x = V[j] + no_recomb_proba
            y = V[max_V_index] + recomb_proba
            if x >= y:
                max_p = x
                max_k = j
            else:
                max_k = max_V_index
                max_p = y
            Vn[j] = max_p + E[H[j, l], h[l]]
            T[j, l] = max_k
            if Vn[j] > Vn[max_Vn_index]:
                max_Vn_index = j
        V = Vn
        max_V_index = max_Vn_index
    P = np.zeros(m, dtype=int)
    P[m - 1] = max_V_index
    for l in range(m - 1, 0, -1):
        P[l - 1] = T[P[l], l]
    return P

def infer_paths_threads(reference_panel, rho, num_workers=None):

    N = reference_panel.num_haplotypes
    n = reference_panel.num_samples
    m = reference_panel.num_sites
    P = np.zeros((N, m), dtype=np.uint32)
    P[-1,:] = -1

    work = []
    for j in range(N - 2, n - 1, -1):
        work.append((j, N - j - 1))
    for j in range(n):
        work.append((j, N - n))
    work_index = 0
    lock = threading.Lock()

    def worker():
        nonlocal work_index
        threader = _tsinf.Threader(reference_panel)
        while True:
            with lock:
                if work_index >= len(work):
                    break
                haplotype_index, panel_size = work[work_index]
                work_index += 1
            threader.run(haplotype_index, panel_size, rho, P[haplotype_index])

    threads = []
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    for _ in range(num_workers):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    return P

def infer_paths(reference_panel, rho):
    N = reference_panel.num_haplotypes
    n = reference_panel.num_samples
    m = reference_panel.num_sites
    threader = _tsinf.Threader(reference_panel)
    P = np.zeros((N, m), dtype=np.uint32)
    p = np.zeros(m, dtype=np.uint32)
    P[-1,:] = -1
    for j in range(N - 2, n - 1, -1):
        # p1 = thread(H, j, num_haplotypes - j - 1, rho)
        threader.run(j, N - j - 1, rho, p)
        # if np.any(p1 != p):
        #     print("ERROR", j)
        #     print(p1)
        #     print(p)
        #     mismatches = np.nonzero(p1 != p)
        #     print(mismatches)
        #     for k in mismatches[0]:
        #         print("\t", k, p[k], p1[k])
        P[j] = p
    for j in range(n):
        threader.run(j, N - n, rho, p)
        # p = thread(H, j, num_haplotypes - n, rho)
        P[j] = p
    return P

def check_paths(H, P, n):
    """
    Checks that the specified paths through the haplotypes make sense.
    """
    N, m = P.shape
    assert H.shape == P.shape
    i = np.arange(m)
    for j in range(N - 1):
        p = P[j]
        h = H[j]
        # hp = np.array([H[p[l], l] for l in range(n)])
        hp = H[p, i]
        # print("path:")
        # print(p)
        # print(h)
        # print(hp)
        for l in np.where(hp != h)[0]:
            assert hp[l] == 0

def convert_records(H, P, n):
    N, m = P.shape
    assert H.shape == P.shape
    C = [[[] for _ in range(m)] for _ in range(N)]

    for j in range(N - 1):
        for l in range(m):
            C[P[j][l]][l].append(j)

    # print("children = ")
    # for j in range(N):
    #     print(j, "\t", C[j])

    mutations = []
    records = []
    # Create the mutations by finding the oldest 1 in each locus.
    # print("mutations = ")
    for l in range(m):
        u = np.where(H[:,l] == 1)[0][0]
        # u is a sample with this mutations. Follow its path upwards until
        # we find the oldest node with the mutation.
        while H[u, l] == 1:
            v = u
            u = P[u][l]
        mutations.append((l, v))
    for u in range(n, N):
        row = C[u]
        last_c = row[0]
        left = 0
        for l in range(1, m):
            if row[l] != last_c:
                if len(last_c) > 0:
                    records.append(msprime.CoalescenceRecord(
                        left=left, right=l, node=u, children=tuple(last_c),
                        time=u, population=0))
                left = l
                last_c = row[l]
        if len(last_c) > 0:
            records.append(msprime.CoalescenceRecord(
                left=left, right=m, node=u, children=tuple(last_c),
                time=u, population=0))
    records.sort(key=lambda r: r.time)
    # print("records = ")
    # for r in records:
    #     print(r)
    # for m in mutations:
    #     print(m)
    ll_ts = _msprime.TreeSequence()
    ll_ts.load_records(records)
    ll_ts.set_mutations(mutations)
    ts = msprime.TreeSequence(ll_ts)
    return ts


def plot_paths(n, P):

    P = (-1 * (P - P.shape[1] + 1)).T
    x = 200
    label_col = np.zeros((P.shape[0], 1))
    label_col[:,0] = np.arange(P.shape[0])
    # label_col[,P.shape[0] - n:] = 10**6
    P = np.hstack((label_col[::-1], P))
    print(P)
    fig, ax = pyplot.subplots() #figsize=(x, x))
    fig.tight_layout(pad=0)
    ax.imshow(P, interpolation="none")
    ax.axhline(P.shape[0] - n - 0.5, color="black")
    ax.axvline(0.5, color="black")
    pyplot.savefig("paths.png")

def verify_variants(V1, ts2):
    V2 = np.empty((ts2.num_mutations, ts2.sample_size), dtype=np.uint8)
    for v in ts2.variants():
        V2[v.index] = v.genotypes

    # print("Comparing:")
    # print(V1)
    # print(V2)
    # print("Equal", np.all(V1 == V2))
    # print(V1 == V2)
    # for mut in ts2.mutations():
    #     print(mut)
    # if not np.all(V1 == V2):
    #     print("err")
    #     for j in range(V1.shape[0]):
    #         if not np.all(V1[j] == V2[j]):
    #             mismatches = np.nonzero(V1[j] != V2[j])
    #             print("Mismatch at ", j, mismatches)
    #             for k in mismatches[0]:
    #                 print("\t", k, V1[j, k], V2[j, k])
    #                 print("\t\t", V1[j])
    #                 print("\t\t", V2[j])

    assert np.all(V1 == V2)

def make_errors(v, p):
    """
    Flip each bit with probability p.
    """
    m = v.shape[0]
    mask = np.random.random(m) < p
    return np.logical_xor(v.astype(bool), mask).astype(int)

def get_samples(ts, error_p):
    """
    Returns samples with a bits flipped with a specified probability.

    Rejects any variants that result in a fixed column.
    """
    S = np.zeros((ts.sample_size, ts.num_mutations), dtype="u1")
    for variant in ts.variants():
        done = False
        # Reject any columns that have no 1s or no zeros
        while not done:
            S[:,variant.index] = make_errors(variant.genotypes, error_p)
            s = np.sum(S[:, variant.index])
            done = 0 < s < ts.sample_size
    return S

def run_checks():
    np.random.seed(10)
    np.set_printoptions(linewidth=1000)
    theta = 1.5
    rho = 2.5
    # n = 1000
    # n = 500
    # length = 500
    n = 50
    length = 100
    error_p = 0.5
    for seed in range(1, 100000):
    # for seed in [1696]:
        ts = msprime.simulate(
            n, Ne=0.5, mutation_rate=theta, length=length,
            recombination_rate=rho, random_seed=seed)
        m = ts.get_num_mutations()
        assert m > 0
        S = get_samples(ts, error_p)
        V = S.T
        # np.savetxt("tmp__NOBACKUP__/samples.txt", S, fmt="%d", delimiter="")
        before = time.clock()
        panel = _tsinf.ReferencePanel(S)
        duration = time.clock() - before
        print("Inferred panel with {} haplotypes and {} sites in {}s".format(
            panel.num_haplotypes, panel.num_sites, duration))
        before = time.clock()
        P = infer_paths_threads(panel, rho)
        # P1 = infer_paths(panel, rho)
        # assert np.all(P1 == P)
        duration = time.clock() - before
        print("Inferred paths in {}s".format(duration))
        H = panel.get_haplotypes()
        # check_paths(H, P, n)
        before = time.clock()
        ts_new = convert_records(H, P, n)
        ts_subset = ts_new.simplify()
        assert ts_new.num_mutations == m
        assert ts_subset.num_mutations == m
        duration = time.clock() - before
        print("Converted records in ", duration, "seconds")
        before = time.clock()
        verify_variants(V, ts_new)
        verify_variants(V, ts_subset)
        duration = time.clock() - before
        print("Verified records in ", duration, "seconds")



def run_replicates(rho, replicates, num_replicates, error):
    num_sites = np.zeros(num_replicates)
    num_samples = np.zeros(num_replicates)
    num_source_records = np.zeros(num_replicates)
    num_raw_records = np.zeros(num_replicates)
    num_simplified_records = np.zeros(num_replicates)
    for k, ts in enumerate(replicates):
        S = get_samples(ts, error)
        panel = _tsinf.ReferencePanel(S)
        P = infer_paths_threads(panel, rho)
        ts_new = convert_records(panel.get_haplotypes(), P, ts.sample_size)
        ts_simplified = ts_new.simplify()
        # Record results
        num_samples[k] = ts.sample_size
        num_sites[k] = ts.num_mutations
        num_source_records[k] = ts.num_records
        num_raw_records[k] = ts_new.num_records
        num_simplified_records[k] = ts_simplified.num_records
    return (
        np.mean(num_samples), np.mean(num_sites), num_replicates, np.mean(num_source_records),
        np.mean(num_raw_records), np.mean(num_simplified_records))


def run_sample_size_plot(output_file, error=0):
    theta = 1.5
    rho = 2.5
    length = 50
    sample_sizes = np.linspace(10, 500, num=10).astype(int)
    num_replicates = 10
    print("Running samples with error rate", error)
    df = pd.DataFrame(columns=(
        "samples", "sites", "replicates", "source_records", "raw_records",
        "simplified_records"))
    for j, n in enumerate(sample_sizes):
        replicates = msprime.simulate(
            n, Ne=0.5, mutation_rate=theta, length=length,
            recombination_rate=rho, num_replicates=num_replicates)
        row = run_replicates(rho, replicates, num_replicates, error)
        df.loc[j] = row
        print("Done for ", n, output_file)
        print(df)
        df.to_csv(output_file)

def run_length_plot(output_file, error=0):
    theta = 1.5
    rho = 2.5
    lengths = np.linspace(10, 500, num=10).astype(int)
    num_replicates = 10
    sample_size = 100
    df = pd.DataFrame(columns=(
        "samples", "sites", "replicates", "source_records", "raw_records",
        "simplified_records"))
    for j, length in enumerate(lengths):
        replicates = msprime.simulate(
            sample_size, Ne=0.5, mutation_rate=theta, length=length,
            recombination_rate=rho, num_replicates=num_replicates)
        row = run_replicates(rho, replicates, num_replicates, error)
        df.loc[j] = row
        print("Done for ", length)
        print(df)
        df.to_csv(output_file)

def run_sample_size():
    data_file = "results/small_example_sample_size.csv"
    run_sample_size_plot(data_file, 0)

def run_length():
    data_file = "results/small_example_length.csv"
    run_length_plot(data_file, 0)


def make_illustration(n, H, P, filename, focal_individual, focal_site):
    N, M = H.shape
    out = tempfile.NamedTemporaryFile("w", prefix="ls_fig_")
    matrix_gap = 2
    haplotype_colours = ["0.95 * white", "0.95 * white"]
    palette = sns.color_palette("Dark2", N - n)
    copy_colours = {-1: "0.95 * white"}
    for j, (r, g, b) in enumerate(palette):
        copy_colours[n + j] = "rgb({}, {}, {})".format(r, g, b)
    print('size(11cm);', file=out)
    print('path cellbox = scale(0.95) * unitsquare;', file=out)
    print('path focal_individual_marker = scale(0.25) * unitcircle;', file=out)
    print('pen cellpen = fontsize(5);', file=out)
    print('defaultpen(fontsize(8));', file=out)
    print('frame f;', file=out)
    print('label(f, "Samples", W);', file=out)
    y = 1
    print('label("Haplotypes", ({}, {}), N);'.format(M / 2, y), file=out)
    print('label("Parent matrix", ({}, {}), N);'.format(M + matrix_gap + M / 2, y), file=out)
    y = -n / 2
    x = -3
    print('add(rotate(90) * f, ({}, {}));'.format(x, y), file=out)
    print('frame f;', file=out)
    print('label(f, "Ancestors", W);', file=out)
    y = -n - (N - n) / 2
    print('add(rotate(90) * f, ({}, {}));'.format(x, y), file=out)
    ancestors = list(range(n)) + list(range(max(focal_individual, n), N))
    for j in ancestors:
        x = -1.5
        y = -j
        print('label("{}", ({}, {}), cellpen);'.format(j, x, y), file=out)
        for k in range(M):
            x = k
            y = -j
            colour = haplotype_colours[H[j, k]]
            if P[focal_individual][k] == j:
                colour = copy_colours[j]
            print('fill(shift(({}, {})) * cellbox, {});'.format(
                x - 0.5, y - 0.5, colour), file=out)
            print('label("{}", ({}, {}), cellpen);'.format(H[j, k], x, y), file=out)
    for j in range(N - 1, focal_individual - 1, -1):
        x = 2 * M + matrix_gap + 0.5
        y = -j
        print('label("{}", ({}, {}), cellpen);'.format(j, x, y), file=out)
        for k in range(M):
            x = k + M + matrix_gap
            y = -j
            print('fill(shift(({}, {})) * cellbox, {});'.format(
                x - 0.5, y - 0.5, copy_colours[P[j, k]]), file=out)
            print('label("{}", ({}, {}), cellpen);'.format(P[j, k], x, y), file=out)
    # Highlight the focal individual and site
    print('filldraw(shift(({}, {})) * focal_individual_marker, black);'.format(
        M + matrix_gap / 2 - 0.5, -focal_individual), file=out)
    if focal_site == -1:
        pen = "invisible"
    else:
        pen = "dotted"
    x, y = focal_site - 0.5, 1
    y_bottom = -N
    print('draw({}--{}--{}--{}--cycle, {});'.format(
        (x, y), (x + 1, y), (x + 1, y_bottom), (x, y_bottom), pen), file=out)
    # Draw the frames around the matrices.
    y_middle = -n + 0.5
    y_bottom = -N + 0.5
    y_top = 0.5
    for x_min in [-0.5, M + matrix_gap - 0.5]:
        x_max = x_min + M
        print('draw(({}, {})--({}, {}));'.format(x_min, y_top, x_max, y_top), file=out)
        print('draw(({}, {})--({}, {}));'.format(x_min, y_middle, x_max, y_middle), file=out)
        print('draw(({}, {})--({}, {}));'.format(x_min, y_bottom, x_max, y_bottom), file=out)
        print('draw(({}, {})--({}, {}));'.format(x_min, y_bottom, x_min, y_top), file=out)
        print('draw(({}, {})--({}, {}));'.format(x_max, y_bottom, x_max, y_top), file=out)

    out.flush()
    subprocess.check_call(["asy", "-f", "pdf", out.name, "-o", filename])
    out.close()

def make_tree(n, H, P, filename, t_source, t_simplified):
    N, M = H.shape

    def get_tree_coords(t, x_min, x_max):
        coords = {}
        interval = {}
        interval[t.root] = x_min, x_max
        # Start at the root and traverse downwards
        for u in t.nodes():
            x1, x2 = interval[u]
            x = x1 + (x2 - x1) / 2
            if u < n:
                y = -u
                coords[u] = (x, y)
            else:
                # Awkward special case for u=n
                time = t.get_time(u) if u > n else n
                delta = M / (M - n)
                y = -time
                coords[u] = (x, y)
            # Now get the interval for the children
            if len(t.children(u)) > 0:
                num_leaves = []
                for v in t.children(u):
                    num_leaves.append(
                        len([w for w in t.nodes(v) if len(t.children(w)) == 0]))
                total = sum(num_leaves)
                delta = (x2 - x1) / total
                x = x1
                for v, leaves in zip(t.children(u), num_leaves):
                    d = delta * leaves
                    interval[v] = x, x + d
                    x += d
        return coords

    def draw_tree(t, coords):
        for j in t.nodes(t.root):
            x, y = coords[j]
            if t.parent(j) != -1:
                dest_x, dest_y = coords[t.parent(j)]
                mid_y = dest_y + 0.75
                print(
                    'draw({}--{}--{}--{});'.format(
                        (x, y), (x, mid_y), (dest_x, mid_y), (dest_x, dest_y)), file=out)
        for u in t.nodes(t.root):
            # We do some weird tricks here to map the labels and colours back to
            # the original tree for the simplified tree.
            if u < n:
                j = u
            elif u == n:
                # work around awkward special case for n
                j = n
            else:
                j = int(t.get_time(u))
            x, y = coords[u]
            colour = copy_colours[j]
            print('label("{}", ({}, {}), cellpen);'.format(j, x, y), file=out)
            print('fill(shift(({}, {})) * cellbox, {});'.format(
                x - 0.5, y - 0.5, colour), file=out)

    group = list(range(*map(int, t_source.interval)))
    out = tempfile.NamedTemporaryFile("w", prefix="ls_fig_")
    matrix_gap = 2.5
    haplotype_colours = ["0.95 * white", "0.95 * white"]
    palette = sns.color_palette("Dark2", N - n)
    copy_colours = {}
    for j in range(-1, n):
        copy_colours[j] = "0.95 * white"
    for j, (r, g, b) in enumerate(palette):
        copy_colours[n + j] = "rgb({}, {}, {})".format(r, g, b)
    print('size(11cm);', file=out)
    print('path cellbox = scale(0.95) * unitsquare;', file=out)
    print('path focal_marker = scale(0.25) * unitcircle;', file=out)
    print('pen cellpen = fontsize(5);', file=out)
    print('pen arrowpen = fontsize(12);', file=out)
    print('pen coverpen = 0.9 * white + opacity(0.8);', file=out)
    print('defaultpen(fontsize(8));', file=out)
    print('frame f;', file=out)
    # print('label(f, "Samples", W);', file=out)
    y = 1
    print('label("Parent matrix", ({}, {}), N);'.format(M / 2, y), file=out)
    print('label("Tree", ({}, {}), N);'.format(M + matrix_gap + M / 4, y), file=out)
    print('label("Simplified", ({}, {}), N);'.format(M + 2 * matrix_gap + 3 * M / 4, y), file=out)
    y = -n / 2
    x = -3
    print('add(rotate(90) * f, ({}, {}));'.format(x, y), file=out)
    print('frame f;', file=out)
    y = -n - (N - n) / 2
    print('add(rotate(90) * f, ({}, {}));'.format(x, y), file=out)

    for j in range(N - 1, - 1, -1):
        x = M + 0.1
        y = -j
        print('label("{}", ({}, {}), cellpen);'.format(j, x, y), file=out)
        for k in range(M):
            x = k
            y = -j
            print('fill(shift(({}, {})) * cellbox, {});'.format(
                x - 0.5, y - 0.5, copy_colours[P[j, k]]), file=out)
            print('label("{}", ({}, {}), cellpen);'.format(P[j, k], x, y), file=out)
            if k not in group:
                print('fill(shift(({}, {})) * cellbox, coverpen);'.format(
                    x - 0.5, y - 0.5), file=out)
    # Draw the frames around the matrices.
    y_middle = -n + 0.5
    y_bottom = -N + 0.5
    y_top = 0.5
    x_min = -0.5
    x_max = x_min + M
    print('draw(({}, {})--({}, {}));'.format(x_min, y_top, x_max, y_top), file=out)
    print('draw(({}, {})--({}, {}));'.format(x_min, y_middle, x_max, y_middle), file=out)
    print('draw(({}, {})--({}, {}));'.format(x_min, y_bottom, x_max, y_bottom), file=out)
    print('draw(({}, {})--({}, {}));'.format(x_min, y_bottom, x_min, y_top), file=out)
    print('draw(({}, {})--({}, {}));'.format(x_max, y_bottom, x_max, y_top), file=out)

    left = M + matrix_gap
    width = M / 2
    coords = get_tree_coords(t_source, left, left + width)
    draw_tree(t_source, coords)
    print('label("$\\rightarrow$", {}, arrowpen);'.format((left - 1, -N / 2)), file=out);
    left = left + width + matrix_gap
    coords = get_tree_coords(t_simplified, left, left + width)
    draw_tree(t_simplified, coords)
    print('label("$\\rightarrow$", {}, arrowpen);'.format((left - 1, -N / 2)), file=out);

    site = group[0] + (group[-1] - group[0]) / 2
    print('filldraw(shift(({}, {})) * focal_marker, black);'.format(
        site, 1), file=out)
    out.flush()
    subprocess.check_call(["asy", "-f", "pdf", out.name, "-o", filename])
    out.close()

def run_figures_simulation():
    seed = 1
    n = 5
    theta = 2.5
    rho = 2.5
    length = 2
    ts = msprime.simulate(
        n, Ne=0.5, mutation_rate=theta, length=length,
        recombination_rate=rho, random_seed=seed)
    S = get_samples(ts, 0)
    panel = _tsinf.ReferencePanel(S)
    P = infer_paths(panel, rho)
    P = np.array(P, dtype=int)
    P[-1] = -1
    H = panel.get_haplotypes()
    return ts, P, H

def run_ls_figures():
    ts, P, H = run_figures_simulation()
    n = ts.sample_size
    for j in range(H.shape[0]):
        filename = "figures/ls-{}.pdf".format(j)
        print("Writing", filename)
        if j == H.shape[0] - 1 or j < n:
            focal_site = -1
        else:
            for l in range(ts.num_mutations):
                if H[j, l] == 1 and H[j + 1, l] == 0:
                    focal_site = l
        make_illustration(n, H, P, filename, j, focal_site)

def run_ls_tree_figures():
    ts, P, H = run_figures_simulation()
    n = ts.sample_size
    ts_new = convert_records(H, P, n)
    ts_simplified = ts_new.simplify()
    ts_simplified_trees = ts_simplified.trees()
    t_simplified = next(ts_simplified_trees)
    for t_source in ts_new.trees():
        filename = "figures/ls-tree-{}.pdf".format(t_source.index)
        print("Writing", filename)
        if t_simplified.interval[1] <= t_source.interval[0]:
            t_simplified = next(ts_simplified_trees)
        make_tree(n, H, P, filename, t_source, t_simplified)

def main():
    cmd = sys.argv[1]
    if cmd == "run-sample-size":
        run_sample_size()
    elif cmd == "run-length":
        run_length()
    elif cmd == "ls-figures":
        run_ls_figures()
    elif cmd == "ls-tree-figures":
        run_ls_tree_figures()
    else:
        print("BAD COMMAND")
        sys.exit(1)

if __name__ == "__main__":
    main()
