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
try:
    import bgen_reader
    # Local module used to work around slow genotype access in bgen_reader
    import simplebgen
except ImportError:
    # bgen-reader isn't available for Python 3.4.
    print("WARNING: Cannot import bgen reader")


@attr.s()
class Site(object):
    position = attr.ib(None)
    alleles = attr.ib(None)
    genotypes = attr.ib(None)
    metadata = attr.ib({})


def filter_duplicates(vcf):
    """
    Returns the variants from this VCF with duplicate sites filtered
    out. If any site position appears more than once, throw all variants away.
    """
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


class VcfConverter(Converter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.num_no_ancestral_state = 0
        self.num_unphased = 0
        self.num_missing_data = 0
        self.num_invariant = 0
        self.num_non_biallelic = 0
        # Rows that don't overlap
        self.num_ancestral_state_no_data = 0
        self.num_data_no_ancestral_state = 0

    def report(self):
        print("no_ancestral_state  :", self.num_no_ancestral_state)
        print("unphased            :", self.num_unphased)
        print("missing_data        :", self.num_missing_data)
        print("invariant           :", self.num_invariant)
        print("non_biallelic       :", self.num_non_biallelic)
        print("ancestral_state_no_data :", self.num_ancestral_state_no_data)
        print("data_no_ancestral_state :", self.num_data_no_ancestral_state)

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
                ret = Site(
                    position=row.POS, alleles=alleles, genotypes=a, metadata=metadata)
        return ret

    def process_sites(self, show_progress=False, max_sites=None):
        num_data_sites = int(subprocess.check_output(
            ["bcftools", "index", "--nrecords", self.data_file]))
        # We tie the iterator to the ancestral states file so it must be updated
        # each time we advance that iterator.
        progress = tqdm.tqdm(total=num_data_sites, disable=not show_progress)

        num_sites = 0
        vcf_a = filter_duplicates(cyvcf2.VCF(self.ancestral_states_file))
        vcf_d = filter_duplicates(cyvcf2.VCF(self.data_file))
        row_a = next(vcf_a, None)
        row_d = next(vcf_d, None)
        while row_a is not None and row_d is not None:

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
                self.num_ancestral_state_no_data += 1
            elif row_d.POS < row_a.POS:
                progress.update()
                self.num_data_no_ancestral_state += 1
                row_d = next(vcf_d, None)

        progress.close()
        self.report()


class ThousandGenomesConverter(VcfConverter):
    """
    Converts data for the 1000 Genomes.
    """

    def process_metadata(self, metadata_file, show_progress=False):
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


class SgdpConverter(VcfConverter):
    """
    Converts data for the Simons Genome Diversity project data.
    """

    def process_metadata(self, metadata_file, show_progress=False):
        """
        Adds the SGDP populations metadata.
        """
        # All populations in SGDP mapped to their regions.
        region_map = {
            "Abkhasian": "WestEurasia",
            "Adygei": "WestEurasia",
            "Albanian": "WestEurasia",
            "Aleut": "CentralAsiaSiberia",
            "Altaian": "CentralAsiaSiberia",
            "Ami": "EastAsia",
            "Armenian": "WestEurasia",
            "Atayal": "EastAsia",
            "Australian": "Oceania",
            "Balochi": "SouthAsia",
            "BantuHerero": "Africa",
            "BantuKenya": "Africa",
            "BantuTswana": "Africa",
            "Basque": "WestEurasia",
            "BedouinB": "WestEurasia",
            "Bengali": "SouthAsia",
            "Bergamo": "WestEurasia",
            "Biaka": "Africa",
            "Bougainville": "Oceania",
            "Brahmin": "SouthAsia",
            "Brahui": "SouthAsia",
            "Bulgarian": "WestEurasia",
            "Burmese": "EastAsia",
            "Burusho": "SouthAsia",
            "Cambodian": "EastAsia",
            "Chane": "America",
            "Chechen": "WestEurasia",
            "Chipewyan": "America",
            "Chukchi": "CentralAsiaSiberia",
            "Cree": "America",
            "Crete": "WestEurasia",
            "Czech": "WestEurasia",
            "Dai": "EastAsia",
            "Daur": "EastAsia",
            "Dinka": "Africa",
            "Druze": "WestEurasia",
            "Dusun": "Oceania",
            "English": "WestEurasia",
            "Esan": "Africa",
            "Eskimo_Chaplin": "CentralAsiaSiberia",
            "Eskimo_Naukan": "CentralAsiaSiberia",
            "Eskimo_Sireniki": "CentralAsiaSiberia",
            "Estonian": "WestEurasia",
            "Even": "CentralAsiaSiberia",
            "Finnish": "WestEurasia",
            "French": "WestEurasia",
            "Gambian": "Africa",
            "Georgian": "WestEurasia",
            "Greek": "WestEurasia",
            "Han": "EastAsia",
            "Hawaiian": "Oceania",
            "Hazara": "SouthAsia",
            "Hezhen": "EastAsia",
            "Hungarian": "WestEurasia",
            "Icelandic": "WestEurasia",
            "Igbo": "Africa",
            "Igorot": "Oceania",
            "Iranian": "WestEurasia",
            "Iraqi_Jew": "WestEurasia",
            "Irula": "SouthAsia",
            "Itelman": "CentralAsiaSiberia",
            "Japanese": "EastAsia",
            "Jordanian": "WestEurasia",
            "Ju_hoan_North": "Africa",
            "Kalash": "SouthAsia",
            "Kapu": "SouthAsia",
            "Karitiana": "America",
            "Kashmiri_Pandit": "SouthAsia",
            "Kharia": "SouthAsia",
            "Khomani_San": "Africa",
            "Khonda_Dora": "SouthAsia",
            "Kinh": "EastAsia",
            "Kongo": "Africa",
            "Korean": "EastAsia",
            "Kurumba": "SouthAsia",
            "Kusunda": "SouthAsia",
            "Kyrgyz": "CentralAsiaSiberia",
            "Lahu": "EastAsia",
            "Lemande": "Africa",
            "Lezgin": "WestEurasia",
            "Luhya": "Africa",
            "Luo": "Africa",
            "Madiga": "SouthAsia",
            "Makrani": "SouthAsia",
            "Mala": "SouthAsia",
            "Mandenka": "Africa",
            "Mansi": "CentralAsiaSiberia",
            "Maori": "Oceania",
            "Masai": "Africa",
            "Mayan": "America",
            "Mbuti": "Africa",
            "Mende": "Africa",
            "Miao": "EastAsia",
            "Mixe": "America",
            "Mixtec": "America",
            "Mongola": "CentralAsiaSiberia",
            "Mozabite": "Africa",
            "Nahua": "America",
            "Naxi": "EastAsia",
            "North_Ossetian": "WestEurasia",
            "Norwegian": "WestEurasia",
            "Onge": "SouthAsia",
            "Orcadian": "WestEurasia",
            "Oroqen": "EastAsia",
            "Palestinian": "WestEurasia",
            "Papuan": "Oceania",
            "Pathan": "SouthAsia",
            "Piapoco": "America",
            "Pima": "America",
            "Polish": "WestEurasia",
            "Punjabi": "SouthAsia",
            "Quechua": "America",
            "Relli": "SouthAsia",
            "Russian": "WestEurasia",
            "Saami": "WestEurasia",
            "Saharawi": "Africa",
            "Samaritan": "WestEurasia",
            "Sardinian": "WestEurasia",
            "She": "EastAsia",
            "Sherpa": "SouthAsia",
            "Sindhi": "SouthAsia",
            "Somali": "Africa",
            "Spanish": "WestEurasia",
            "Surui": "America",
            "Tajik": "WestEurasia",
            "Thai": "EastAsia",
            "Tibetan": "SouthAsia",
            "Tlingit": "CentralAsiaSiberia",
            "Tubalar": "CentralAsiaSiberia",
            "Tu": "EastAsia",
            "Tujia": "EastAsia",
            "Turkish": "WestEurasia",
            "Tuscan": "WestEurasia",
            "Ulchi": "CentralAsiaSiberia",
            "Uygur": "EastAsia",
            "Xibo": "EastAsia",
            "Yadava": "SouthAsia",
            "Yakut": "CentralAsiaSiberia",
            "Yemenite_Jew": "WestEurasia",
            "Yi": "EastAsia",
            "Yoruba": "Africa",
            "Zapotec": "America",
        }
        population_id_map = {}
        for name in sorted(region_map.keys()):
            pop_id = self.samples.add_population(
                {"name": name, "region": region_map[name]})
            population_id_map[name] = pop_id

        # The file contains some non UTF-8 codepoints for a contributors name.
        with open(metadata_file, "r", encoding="ISO-8859-1") as md_file:
            columns = next(md_file).lstrip("#").split("\t")
            sane_names = [col.lower().strip() for col in columns]
            j = sane_names.index("sample_id(aliases)")
            sane_names[j] = "aliases"
            for j, name in enumerate(sane_names):
                if name.startswith("\"sgdp-lite category"):
                    # There's a very long key that doesn't impart any information here.
                    # Remove it.
                    sane_names[j] = "DELETE"
            rows = {}
            populations = {}
            locations = {}
            for line in md_file:
                metadata = dict(zip(sane_names, line.strip().split("\t")))
                del metadata["DELETE"]
                name = metadata["sgdp_id"]
                population_name = metadata.pop("population_id")
                populations[name] = population_id_map[population_name]
                rows[name] = metadata
                location = [
                    float(metadata.pop("latitude")), float(metadata.pop("longitude"))]
                locations[name] = location
                if metadata["town"] == "?":
                    metadata["town"] = None

        vcf = cyvcf2.VCF(self.data_file)
        individual_names = list(vcf.samples)
        vcf.close()
        self.num_samples = 2 * len(individual_names)

        # Add in the metadata rows in the order of the VCF.
        for name in individual_names:
            if name == "S_Naxi-2":
                # This individual is missing from the metadata. Skip metadata and location,
                # but assume the population is the same as S_Naxi-1
                self.samples.add_individual(population=populations["S_Naxi-1"], ploidy=2)
            else:
                metadata = rows[name]
                self.samples.add_individual(
                    metadata=metadata, location=locations[name], ploidy=2,
                    population=populations[name])


class UkbbConverter(Converter):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.num_no_ancestral_state = 0
        self.num_missing_ancestral_state = 0
        self.num_rsid_mismatches = 0
        self.num_non_biallelic = 0

    def report(self):
        print("no_ancestral_state        :", self.num_no_ancestral_state)
        print("missing_ancestral_states  :", self.num_missing_ancestral_state)
        print("non_biallelic             :", self.num_non_biallelic)
        print("rsid_mismatches           :", self.num_rsid_mismatches)

    def process_metadata(self, metadata_file, show_progress=False):
        bgen = bgen_reader.read_bgen(self.data_file, verbose=False)
        sample_df = bgen['samples']
        num_individuals = len(sample_df)
        self.num_samples = 2 * num_individuals
        for _ in tqdm.tqdm(range(num_individuals), disable=not show_progress):
            self.samples.add_individual(ploidy=2)

    def process_sites(self, show_progress=False, max_sites=None):

        bgen = bgen_reader.read_bgen(self.data_file, verbose=False)
        num_alleles = np.array(bgen["variants"]["nalleles"])
        position = np.array(bgen["variants"]["pos"])
        rsid = np.array(bgen["variants"]["rsid"])
        allele_id = np.array(bgen["variants"]["allele_ids"])
        del bgen
        vcf_a = cyvcf2.VCF(self.ancestral_states_file)
        row_a = next(vcf_a, None)

        bg = simplebgen.BgenReader(self.data_file)
        N = 2 * bg.num_samples
        for j in tqdm.tqdm(range(bg.num_variants)):
            while row_a is not None and row_a.POS < position[j]:
                row_a = next(vcf_a, None)
            if row_a is None:
                break

            if row_a.POS == position[j]:
                ancestral_state = get_ancestral_state(row_a)
                alleles = allele_id[j].split(",")
                if row_a.ID != rsid[j]:
                    self.num_rsid_mismatches += 1
                elif num_alleles[j] != 2:
                    self.num_non_biallelic += 1
                elif ancestral_state not in alleles:
                    self.num_missing_ancestral_state += 1
                else:
                    P = bg.get_probabilities(j).astype(np.int8).reshape((N, 2))
                    # The probabilities for each site is a (num_diploids, 4) array,
                    # in the form (n0_a0, n0_a1, n1_a0, n1_a1). These are always zero
                    # or one for the different alleles. We first flatten this array so
                    # that it's (N, 2) and then generate the genotypes based on that.
                    genotypes = np.zeros(N, dtype=np.int8)
                    if ancestral_state == alleles[0]:
                        genotypes[P[:, 1] == 1] = 1
                    else:
                        genotypes[P[:, 0] == 1] = 1
                    metadata = {"ID": row_a.ID, "REF": row_a.REF}
                    self.samples.add_site(
                        position=float(position[j]), genotypes=genotypes,
                        alleles=alleles, metadata=metadata)
            else:
                self.num_no_ancestral_state += 1
            if j == max_sites:
                break
        self.report()


def main():
    parser = argparse.ArgumentParser(
        description="Script to convert VCF files into tsinfer input.")
    parser.add_argument(
        "source", choices=["1kg", "sgdp", "ukbb"],
        help="The source of the input data.")
    parser.add_argument(
        "data_file", help="The input data file pattern.")
    parser.add_argument(
        "ancestral_states_file",
        help="A vcf file containing ancestral allele states. ")
    parser.add_argument(
        "output_file",
        help="The tsinfer output file")
    parser.add_argument(
        "-m", "--metadata_file", default=None,
        help="The metadata file containing population and sample data")
    parser.add_argument(
        "-n", "--max-variants", default=None, type=int,
        help="Keep only the first n variants")
    parser.add_argument(
        "-p", "--progress", action="store_true",
        help="Show progress bars and output extra information when done")
    parser.add_argument(
        "--ancestral-states-url", default=None,
        help="The source of ancestral state information for provenance.")
    parser.add_argument(
        "--reference-name", default=None,
        help="The name of the reference for provenance.")

    args = parser.parse_args()

    git_hash = subprocess.check_output(["git", "rev-parse", "HEAD"])
    git_provenance = {
        "repo": "git@github.com:mcveanlab/treeseq-inference.git",
        "hash": git_hash.decode().strip(),
        "dir": "human-data",
        "notes:": (
            "Use the Makefile to download and process the upstream data files")}
    data_provenance = {
        "ancestral_states_url": args.ancestral_states_url,
        "reference_name": args.reference_name
    }

    try:
        with tsinfer.SampleData(path=args.output_file, num_flush_threads=2) as samples:
            if args.source == "1kg":
                converter = ThousandGenomesConverter(
                    args.data_file, args.ancestral_states_file, samples)
            if args.source == "sgdp":
                converter = SgdpConverter(
                    args.data_file, args.ancestral_states_file, samples)
            if args.source == "ukbb":
                converter = UkbbConverter(
                    args.data_file, args.ancestral_states_file, samples)
            converter.process_metadata(args.metadata_file, args.progress)
            converter.process_sites(args.progress, args.max_variants)
            samples.record_provenance(
                command=sys.argv[0], args=sys.argv[1:], git=git_provenance,
                data=data_provenance)
    except Exception as e:
        os.unlink(args.output_file)
        raise e
    print(samples)


if __name__ == "__main__":
    main()
