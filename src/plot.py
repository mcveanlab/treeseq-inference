# Generates all the actual figures. Run like
# python3 src/plot.py PLOT_NAME

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

    def save(self, figure_name=None):
        if figure_name is None:
            figure_name = self.name
        print("Saving figure ", figure_name)
        plt.savefig("figures/{}.pdf".format(figure_name))
        plt.clf()


class StoringEveryone(Figure):
    """
    Figure showing how tree sequences can store the entire human population
    worth of variation data.
    """
    name = "storing_everyone"

    def plot(self):
        df = self.data
        df = df[df.sample_size > 100]

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        xytext = (18, 0)
        GB = 1024**3
        largest_n = np.array(df.sample_size)[-1]

        index = df.vcf > 0
        line, = ax1.loglog(df.sample_size[index], df.vcf[index], "^", label=".vcf")
        ax1.loglog(df.sample_size, df.vcf_fit, "--", color=line.get_color(), label="")
        largest_value = np.array(df.vcf_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(df.sample_size[index], df.vcfz[index], "s", label=".vcf.gz")
        ax1.loglog(df.sample_size, df.vcfz_fit, "--", color=line.get_color(), label="")
        largest_value = np.array(df.vcfz_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.uncompressed, "o", label=".trees")
        ax1.loglog(df.sample_size, df.tsk_fit, "--", color=line.get_color(), label="")
        largest_value = np.array(df.tsk_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.compressed, "*", label=".trees.gz")
        ax1.loglog(df.sample_size, df.tskz_fit, "--", color=line.get_color(), label="")
        largest_value = np.array(df.tskz_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        ax1.set_xlabel("Number of chromosomes")
        ax1.set_ylabel("File size (GiB)")
        plt.legend()
        # plt.tight_layout()
        self.save()


class SampleEdges(Figure):
    name = "sample_edges"

    def plot_region(self, df, dataset, region):
        fig = plt.figure(figsize=(14,6))
        ax = fig.add_subplot(111)
        ax.plot(df.sample_edges.values, "o")
        breakpoints = np.where(df.population.values[1:] != df.population.values[:-1])[0]
        breakpoints = list(breakpoints) + [len(df)]
        last = 0
        y = df.sample_edges.min()
        for bp in breakpoints:
            label = df.population[bp - 1]
            x = last + (bp - last) / 2
            ax.annotate(label, xy=(x, y), horizontalalignment='centre', verticalalignment='top', rotation=90)
            last = bp
        ax.set_xticks(np.array(breakpoints) + 0.5)
        ax.set_xticklabels([])
        ax.set_xlim(-0.5, len(df) - 0.5)
        ax.set_title("{}:{}".format(dataset.upper(), region))
        ax.grid(axis="x")
        self.save("{}_{}_{}".format(self.name, dataset, region))

    def plot(self):
        full_df = self.data

        fig, axes = plt.subplots(2, figsize=(14,6))
        for ax, dataset in zip(axes, ["1kg", "sgdp"]):
            df = full_df[full_df.dataset == dataset]
            df = df.sort_values(by=["region", "population", "sample", "strand"])
            df = df.reset_index()

            ax.plot(df.sample_edges.values)
            breakpoints = np.where(df.region.values[1:] != df.region.values[:-1])[0]
            for bp in breakpoints:
                ax.axvline(x=bp, ls="--", color="black")

            last = 0
            for j, bp in enumerate(list(breakpoints) + [len(df)]):
                x = last + (bp - last) / 2
                ax.annotate(df.region[bp - 1], xy=(x, 200), horizontalalignment='center')
                last = bp

            breakpoints = np.where(df.population.values[1:] != df.population.values[:-1])[0]
            breakpoints = list(breakpoints) + [len(df)]
            ax.set_xticks(breakpoints)
            ax.set_xticklabels([])
            ax.grid(axis="x")
            ax.set_xlim(0, len(df))

            if dataset == "1kg":
                last = 0
                for bp in breakpoints:
                    x = last + (bp - last) / 2
                    last = bp
                    ax.annotate(
                        df.population[int(x)], xy=(x, 0), horizontalalignment='right',
                        verticalalignment='top', rotation=270)

        axes[0].set_ylim(0, 1500)
        axes[1].set_ylim(0, 3500)
        self.save()

        # Also plot each region
        for dataset, region in set(zip(full_df.dataset, full_df.region)):
            df = full_df[(full_df.dataset == dataset) & (full_df.region == region)]
            df = df.sort_values(by=["population", "sample", "strand"])
            df = df.reset_index()
            self.plot_region(df, dataset, region)


def main():
    figures = [StoringEveryone, SampleEdges]

    name_map = {fig.name: fig for fig in figures}

    parser = argparse.ArgumentParser(description="Make the plots for specific figures.")
    parser.add_argument(
        "name", type=str, help="figure name", choices=list(name_map.keys()))

    args = parser.parse_args()
    fig = name_map[args.name]()
    fig.plot()


if __name__ == "__main__":
    main()
