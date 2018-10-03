# Generates all the actual figures. Run like
# python3 src/plot.py PLOT_NAME

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import pandas as pd
import numpy as np
import humanize


class Figure(object):
    """
    Superclass of figures for the paper. Each figure is a concrete subclass.
    """
    def __init__(self):
        datafile_name = "data/{}.csv".format(self.name)
        self.data = pd.read_csv(datafile_name)

    def save(self):
        pyplot.savefig("figures/{}.pdf".format(self.name))



class StoringEveryone(Figure):
    """
    Figure showing how tree sequences can store the entire human population
    worth of variation data.
    """
    name = "storing_everyone"

    def plot(self):
        df = self.data
        df = df[df.sample_size > 100]

        fig = pyplot.figure()
        ax1 = fig.add_subplot(111)
        xytext = (18, 0)
        GB = 1024**3
        largest_n = np.array(df.sample_size)[-1]

        index = df.vcf > 0
        line, = ax1.loglog(df.sample_size[index], df.vcf[index], "^", label="VCF")
        ax1.loglog(df.sample_size, df.vcf_fit, "--", color=line.get_color())
        largest_value = np.array(df.vcf_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(df.sample_size[index], df.bcf[index], "s", label="BCF")
        ax1.loglog(df.sample_size, df.bcf_fit, "--", color=line.get_color())
        largest_value = np.array(df.bcf_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.uncompressed, "o", label=".trees")
        ax1.loglog(df.sample_size, df.tsk_fit, "--", color=line.get_color())
        largest_value = np.array(df.tsk_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.compressed, "*", label=".trees.gz")
        ax1.loglog(df.sample_size, df.tskz_fit, "--", color=line.get_color())
        largest_value = np.array(df.tskz_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        ax1.set_xlabel("Number of chromosomes")
        ax1.set_ylabel("File size (GiB)")
        pyplot.legend()
        # pyplot.tight_layout()
        self.save()


def main():
    figures = [StoringEveryone]

    name_map = {fig.name: fig for fig in figures}

    parser = argparse.ArgumentParser(description="Make the plots for specific figures.")
    parser.add_argument(
        "name", type=str, help="figure name", choices=list(name_map.keys()))

    args = parser.parse_args()
    fig = name_map[args.name]()
    fig.plot()


if __name__ == "__main__":
    main()
