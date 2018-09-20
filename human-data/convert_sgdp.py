"""
Import the SGDP VCFs.
"""
import argparse
import subprocess
import os
import sys

import numpy as np
import tsinfer
import attr
import cyvcf2
import tqdm


@attr.s()
class Site(object):
    position = attr.ib(None)
    alleles = attr.ib(None)
    genotypes = attr.ib(None)
    metadata = attr.ib({})



def filter_duplicates(vcf):
    """
    Returns the list of variants from this specified VCF with duplicate sites filtered
    out. If any site appears more than once, throw all variants away.
    """
    # TODO this had not been tested properly.
    row = next(vcf, None)
    bad_pos = -1
    for next_row in vcf:
        if bad_pos == -1 and next_row.POS != row.POS:
            yield row
        else:
            if bad_pos == -1:
                bad_pos = row.POS
            elif bad_pos != next_row.POS:
                bad_pos = -1
        row = next_row
    if row is not None and bad_pos != -1:
        yield row


def vcf_num_rows(vcf_path):
    """
    Return the number of rows in the VCF. Requires an index to be present.
    """
    output = subprocess.check_output(["bcftools", "index", "--nrecords", vcf_path])
    return int(output)


def variants(vcf_path, ancestral_states, show_progress=False, max_sites=None):
    """
    Yield a tuple of position, alleles, genotypes, metadata. Ancestral_states is a
    dictionary mapping RSID to the ancestral allele.
    """
    tot_sites = vcf_num_rows(vcf_path)
    progress = tqdm.tqdm(
        total=tot_sites, desc="Read genotypes", disable=not show_progress)

    vcf = cyvcf2.VCF(vcf_path)

    sites_used = 0
    non_biallelic = 0
    no_ancestral_state = 0
    missing_data = 0
    indels = 0
    unphased = 0
    invariant = 0

    num_diploids = len(vcf.samples)
    num_samples = 2 * num_diploids
    j = 0
    for row in filter_duplicates(vcf):
        progress.update()
        key = str(row.POS) + ":" + row.REF
        # only use biallelic sites with data for all samples and where ancestral state is known
        if len(row.ALT) != 1:
            non_biallelic += 1
        elif len(row.ALT[0]) != 1 or len(row.REF) != 1:
            indels += 1
        elif row.num_called != num_diploids:
            missing_data += 1
        elif key not in ancestral_states:
            no_ancestral_state += 1
        else:
            ancestral_state = ancestral_states[key]
            a = np.zeros(num_samples, dtype=np.uint8)
            all_alleles = set([ancestral_state])
            # Fill in a with genotypes.
            bases = np.array(row.gt_bases)
            for j in range(num_diploids):
                alleles = bases[j].split("|")
                if len(alleles) != 2:
                    unphased += 1
                    break
                for allele in alleles:
                    if allele == ".":
                        missing_data += 1
                        break
                    all_alleles.add(allele)
                a[2 * j] = alleles[0] != ancestral_state
                a[2 * j + 1] = alleles[1] != ancestral_state
            else:
                # The loop above exited without breaking, so we have valid data.
                if len(all_alleles) == 1:
                    invariant += 1
                elif len(all_alleles) > 2:
                    non_biallelic += 1
                else:
                    progress.set_postfix(used=sites_used)
                    all_alleles.remove(ancestral_state)
                    alleles = [ancestral_state, all_alleles.pop()]
                    metadata = {"ID": row.ID, "REF": row.REF}
                    sites_used += 1
                    yield Site(position=row.POS, alleles=alleles, genotypes=a, metadata=metadata)
                    if max_sites == sites_used:
                        break
    progress.close()
    vcf.close()
    print(
        "Used {} out of {} sites. {} non biallelic, {} without ancestral state, "
        "{} indels, {} with missing data, {} invariant, and {} unphased".format(
            sites_used, tot_sites, non_biallelic, no_ancestral_state, indels,
            missing_data, invariant, unphased))


def add_samples(ped_file, individual_names, sample_data):
    """
    Reads the specified PED file to get information about the samples.
    Individuals are to be added in the order in the specified list.
    """
    columns = next(ped_file).lstrip("#").split("\t")
    sane_names = [col.lower().strip() for col in columns]
    j = sane_names.index("sample_id(aliases)")
    sane_names[j] = "aliases"
    for j, name in enumerate(sane_names):
        if name.startswith("\"sgdp-lite category"):
            # There's a very long key that doesn't impart any information here. Remove it.
            sane_names[j] = "DELETE"
    rows = {}
    populations = {}
    population_id_map = {}
    locations = {}
    for line in ped_file:
        metadata = dict(zip(sane_names, line.strip().split("\t")))
        del metadata["DELETE"]
        name = metadata["sgdp_id"]
        population_name = metadata.pop("population_id")
        if population_name not in population_id_map:
            sample_data.add_population(metadata={"name": population_name})
            population_id_map[population_name] = len(population_id_map)
        populations[name] = population_id_map[population_name]
        rows[name] = metadata
        location = [float(metadata.pop("latitude")), float(metadata.pop("longitude"))]
        locations[name] = location
        if metadata["town"] == "?":
            metadata["town"] = None

    # Add in the metadata rows in the order of the VCF.
    for name in individual_names:
        if name not in rows:
            sample_data.add_individual(ploidy=2)
        else:
            metadata = rows[name]
            sample_data.add_individual(
                metadata=metadata, location=locations[name], ploidy=2,
                population=populations[name])


def read_ancestral_states(ancestor_file, show_progress=False):
    states = {}
    sites = vcf_num_rows(ancestor_file)
    vcf = cyvcf2.VCF(ancestor_file)
    iterator = tqdm.tqdm(
        vcf, total=sites, desc="Read ancestral state", disable=not show_progress)
    for site in iterator:
        key = str(site.POS) + ":" + site.REF
        try:
            aa = site.INFO["AA"]
        except KeyError:
            aa = None
        if aa is not None:
            # Format: for indels = AA|REF|ALT|IndelType; for SNPs = AA
            splits = aa.split("|")
            # We're only interest in SNPs
            if len(splits[0]) == 1:
                base = splits[0].upper()
                if base in "ACTG":
                    states[key] = base
    vcf.close()
    print("Read {} ancestral alleles out of {} sites from {}".format(
        len(states), sites, ancestor_file))
    return states


def convert(
        vcf_file, ancestral_states_file, pedigree_file, output_file, max_variants=None,
        show_progress=False):

    git_hash = subprocess.check_output(["git", "rev-parse", "HEAD"])
    git_provenance = {
        "repo": "git@github.com:mcveanlab/treeseq-inference.git",
        "hash": git_hash.decode().strip(),
        "dir": "human-data",
        "notes:": (
            "Use the Makefile to download and process the upstream data files")}

    with tsinfer.SampleData(path=output_file, num_flush_threads=2) as sample_data:

        vcf = cyvcf2.VCF(vcf_file)
        individual_names = list(vcf.samples)
        vcf.close()

        # The file contains some non UTF-8 codepoints for a contributors name.
        with open(pedigree_file, "r", encoding="ISO-8859-1") as ped_file:
            add_samples(ped_file, individual_names, sample_data)

        ancestral_states = read_ancestral_states(ancestral_states_file, show_progress)

        iterator = variants(vcf_file, ancestral_states, show_progress, max_variants)
        for site in iterator:
            sample_data.add_site(
                position=site.position, genotypes=site.genotypes,
                alleles=site.alleles, metadata=site.metadata)
        sample_data.record_provenance(
            command=sys.argv[0], args=sys.argv[1:], git=git_provenance)



def main():
    parser = argparse.ArgumentParser(
        description="Script to convert VCF files into tsinfer input.")
    parser.add_argument(
        "vcf_file", help="The input VCF file pattern.")
    parser.add_argument(
        "ancestral_states_file",
        help="A vcf file containing ancestral allele states. ")
    parser.add_argument(
        "pedigree_file",
        help="The pedigree file containing population and sample data")
    parser.add_argument(
        "output_file",
        help="The tsinfer output file")
    parser.add_argument(
        "-n", "--max-variants", default=None, type=int,
        help="Keep only the first n variants")
    parser.add_argument(
        "-p", "--progress", action="store_true",
        help="Show progress bars and output extra information when done")

    args = parser.parse_args()

    if not os.path.exists(args.vcf_file):
        raise ValueError("{} does not exist".format(args.vcf_file))
    if not os.path.exists(args.pedigree_file):
        raise ValueError("{} does not exist".format(args.pedigree_file))

    convert(
        args.vcf_file, args.ancestral_states_file, args.pedigree_file, args.output_file,
        args.max_variants, show_progress=args.progress)


if __name__ == "__main__":
    main()
