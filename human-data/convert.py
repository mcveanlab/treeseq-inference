"""
Convert input data from various sources to samples format.
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


class Converter(object):
    """
    Superclass of converters.
    """
    def __init__(self, data_file, ancestral_states_file, samples):
        self.data_file = data_file
        self.ancestral_states_file = ancestral_states_file
        self.samples = samples
        self.num_samples = -1

    def process_metadata(self, metadata_file):
        pass


def get_ancestral_state(row):
    try:
        aa = row.INFO["AA"]
    except KeyError:
        aa = None
    if aa is not None:
        # Format: for indels = AA|REF|ALT|IndelType; for SNPs = AA
        splits = aa.split("|")
        # We're only interest in SNPs
        if len(splits[0]) == 1:
            base = splits[0].upper()
            if base in "ACTG":
                aa = base
    return aa

class VcfConverter(Converter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.num_no_ancestral_state = 0
        self.num_unphased = 0
        self.num_missing_data = 0
        self.num_invariant = 0
        self.num_non_biallelic = 0

    def report(self):
        print("no_ancestral_state  :", self.num_no_ancestral_state)
        print("unphased            :", self.num_unphased)
        print("missing_data        :", self.num_missing_data)
        print("invariant           :", self.num_invariant)
        print("non_biallelic       :", self.num_non_biallelic)

    def convert_genotypes(self, row, ancestral_state):
        ret = None
        num_diploids = self.num_samples // 2
        a = np.zeros(self.num_samples, dtype=np.uint8)
        all_alleles = set([ancestral_state])
        # Fill in a with genotypes.
        bases = np.array(row.gt_bases)
        for j in range(num_diploids):
            alleles = bases[j].split("|")
            if len(alleles) != 2:
                self.num_unphased += 1
                break
            for allele in alleles:
                if allele == ".":
                    self.num_missing_data += 1
                    break
                all_alleles.add(allele)
            a[2 * j] = alleles[0] != ancestral_state
            a[2 * j + 1] = alleles[1] != ancestral_state
        else:
            # The loop above exited without breaking, so we have valid data.
            if len(all_alleles) == 1:
                self.num_invariant += 1
            elif len(all_alleles) > 2:
                self.num_non_biallelic += 1
            else:
                all_alleles.remove(ancestral_state)
                alleles = [ancestral_state, all_alleles.pop()]
                metadata = {"ID": row.ID, "REF": row.REF}
                ret = Site(position=row.POS, alleles=alleles, genotypes=a, metadata=metadata)
        return ret

    def process_sites(self, show_progress=False, max_sites=None):
        num_ancestral_sites = int(subprocess.check_output(
            ["bcftools", "index", "--nrecords", self.ancestral_states_file]))
        # We tie the iterator to the ancestral states file so it must be updated
        # each time we advance that iterator.
        progress = tqdm.tqdm(
            total=num_ancestral_sites, disable=not show_progress)

        num_sites = 0
        vcf_a = cyvcf2.VCF(self.ancestral_states_file)
        vcf_d = cyvcf2.VCF(self.data_file)
        row_a = next(vcf_a, None)
        row_d = next(vcf_d, None)
        while row_a is not None and row_d is not None:

            # Filtering by POS and REF should guarantee that any sites with equal positions
            # as we are excluding indels from the ancestral states explicitly and we're also
            # dropping any sites that are not single letter SNPs.
            if row_a.POS == row_d.POS and row_a.REF == row_d.REF:
                ancestral_state = get_ancestral_state(row_a)
                if ancestral_state is not None:
                    site = self.convert_genotypes(row_d, ancestral_state)
                    if site is not None:
                        self.samples.add_site(
                            position=site.position, genotypes=site.genotypes,
                            alleles=site.alleles, metadata=site.metadata)
                        progress.set_postfix(used=str(num_sites))
                        num_sites += 1
                        if num_sites == max_sites:
                            break
                else:
                    self.num_no_ancestral_state += 1
            if row_a.POS == row_d.POS:
                # Advance both iterators.
                row_a = next(vcf_a, None)
                row_d = next(vcf_d, None)
                progress.update()
            elif row_a.POS < row_d.POS:
                row_a = next(vcf_a, None)
                progress.update()
            elif row_d.POS < row_a.POS:
                row_d = next(vcf_d, None)

        progress.close()
        self.report()


class ThousandGenomesConverter(VcfConverter):
    """
    Converts data for the 1000 Genomes.
    """

    def process_metadata(self, metadata_file):
        """
        Adds the 1000 genomes populations metadata.
        """
        # Based on
        # http://www.internationalgenome.org/faq/which-populations-are-part-your-study/
        populations = [
            ["CHB", "Han Chinese in Beijing, China", "EAS"],
            ["JPT", "Japanese in Tokyo, Japan", "EAS"],
            ["CHS", "Southern Han Chinese", "EAS"],
            ["CDX", "Chinese Dai in Xishuangbanna, China", "EAS"],
            ["KHV", "Kinh in Ho Chi Minh City, Vietnam", "EAS"],
            ["CEU", "Utah Residents (CEPH) with Northern and Western European Ancestry",
                "EUR"],
            ["TSI", "Toscani in Italia", "EUR"],
            ["FIN", "Finnish in Finland", "EUR"],
            ["GBR", "British in England and Scotland", "EUR"],
            ["IBS", "Iberian Population in Spain", "EUR"],
            ["YRI", "Yoruba in Ibadan, Nigeria", "AFR"],
            ["LWK", "Luhya in Webuye, Kenya", "AFR"],
            ["GWD", "Gambian in Western Divisions in the Gambia", "AFR"],
            ["MSL", "Mende in Sierra Leone", "AFR"],
            ["ESN", "Esan in Nigeria", "AFR"],
            ["ASW", "Americans of African Ancestry in SW USA", "AFR"],
            ["ACB", "African Caribbeans in Barbados", "AFR"],
            ["MXL", "Mexican Ancestry from Los Angeles USA", "AMR"],
            ["PUR", "Puerto Ricans from Puerto Rico", "AMR"],
            ["CLM", "Colombians from Medellin, Colombia", "AMR"],
            ["PEL", "Peruvians from Lima, Peru", "AMR"],
            ["GIH", "Gujarati Indian from Houston, Texas", "SAS"],
            ["PJL", "Punjabi from Lahore, Pakistan", "SAS"],
            ["BEB", "Bengali from Bangladesh", "SAS"],
            ["STU", "Sri Lankan Tamil from the UK", "SAS"],
            ["ITU", "Indian Telugu from the UK", "SAS"],
        ]
        population_id_map = {}
        for pop in populations:
            pop_id = self.samples.add_population(
                dict(zip(["name", "description", "super_population"], pop)))
            population_id_map[pop[0]] = pop_id

        with open(metadata_file, "r") as ped_file:
            # Parse the individual metadata out of the ped file.
            columns = next(ped_file).split("\t")
            sane_names = [col.replace(" ", "_").lower().strip() for col in columns]
            metadata = {}
            populations = {}
            for line in ped_file:
                row = dict(zip(sane_names, line.strip().split("\t")))
                name = row["individual_id"]
                population_name = row.pop("population")
                populations[name] = population_id_map[population_name]
                # The value '0' seems to be used to encode missing, so insert None
                # instead to be more useful.
                nulled = {}
                for key, value in row.items():
                    if value == "0":
                        value = None
                    nulled[key] = value
                metadata[name] = nulled

        vcf = cyvcf2.VCF(self.data_file)
        individual_names = list(vcf.samples)
        vcf.close()
        self.num_samples = len(individual_names) * 2
        # Add in the metadata rows in the order of the VCF.
        for name in individual_names:
            self.samples.add_individual(
                metadata=metadata[name], population=populations[name], ploidy=2)




def main():
    parser = argparse.ArgumentParser(
        description="Script to convert VCF files into tsinfer input.")
    parser.add_argument(
        "source", choices=["1kg"],
        help="The source of the input data.")
    parser.add_argument(
        "data_file", help="The input data file pattern.")
    parser.add_argument(
        "ancestral_states_file",
        help="A vcf file containing ancestral allele states. ")
    parser.add_argument(
        "metadata_file",
        help="The metadata file containing population and sample data")
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

    with tsinfer.SampleData(path=args.output_file, num_flush_threads=2) as samples:
        if args.source == "1kg":
            converter = ThousandGenomesConverter(
                args.data_file, args.ancestral_states_file, samples)

        converter.process_metadata(args.metadata_file)
        converter.process_sites(args.progress, args.max_variants)
    print(samples)

    # convert(
    #     args.vcf_file, args.ancestral_states_file, args.pedigree_file, args.output_file,
    #     args.max_variants, show_progress=args.progress)


if __name__ == "__main__":
    main()
