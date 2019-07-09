#!/usr/bin/env python3
"""
Generates all the actual figures. Run like
 python3 src/plot.py PLOT_NAME
"""
import argparse
import collections

import pandas as pd
import numpy as np
import humanize
import matplotlib
matplotlib.use('Agg')  # NOQA
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition, inset_axes
import seaborn as sns
import scipy

tgp_region_pop = {
    'AMR': ['CLM', 'MXL', 'PUR', 'PEL'],
    'AFR': ['LWK', 'ASW', 'GWD', 'MSL', 'YRI', 'ACB', 'ESN'],
    'EAS': ['CHS', 'KHV', 'JPT', 'CHB', 'CDX'],
    'SAS': ['BEB', 'STU', 'GIH', 'PJL', 'ITU'],
    'EUR': ['FIN', 'GBR', 'IBS', 'CEU', 'TSI']
}

# Standard order.
tgp_populations = [
    'CHB', 'JPT', 'CHS', 'CDX', 'KHV',
    'CEU', 'TSI', 'FIN', 'GBR', 'IBS',
    'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB',
    'MXL', 'PUR', 'CLM', 'PEL',
    'GIH', 'PJL', 'BEB', 'STU', 'ITU']


def get_tgp_region_colours():
    return {
        "EAS": sns.color_palette("Greens", 2)[1],
        "EUR": sns.color_palette("Blues", 1)[0],
        "AFR": sns.color_palette("Wistia", 3)[0],
        "AMR": sns.color_palette("Reds", 2)[1],
        "SAS": sns.color_palette("Purples", 2)[1],
    }


def get_sgdp_region_colours():
    cols = get_tgp_region_colours()
    return {
        'Africa': cols["AFR"],
        'America': cols["AMR"],
        'EastAsia': cols["EAS"],
        'SouthAsia': cols["SAS"],
        'Oceania': "brown",
        'WestEurasia': cols["EUR"],
        'CentralAsiaSiberia': "pink"
     }


def get_tgp_colours():
    # TODO add option to give shades for the different pops.
    region_colours = get_tgp_region_colours()
    pop_colour_map = {}
    for region, pops in tgp_region_pop.items():
        for pop in pops:
            pop_colour_map[pop] = region_colours[region]
    return pop_colour_map


class Figure(object):
    """
    Superclass of figures for the paper. Each figure is a concrete subclass.
    """
    name = None

    def __init__(self):
        datafile_name = "data/{}.csv".format(self.name)
        self.data = pd.read_csv(datafile_name)

    def save(self, figure_name=None, bbox_inches="tight"):
        if figure_name is None:
            figure_name = self.name
        print("Saving figure '{}'".format(figure_name))
        plt.savefig("figures/{}.pdf".format(figure_name), bbox_inches='tight', dpi=400)
        plt.savefig("figures/{}.png".format(figure_name), bbox_inches='tight', dpi=400)
        plt.close()

    def error_label(self, error, label_for_no_error = "No genotyping error"):
        """
        Make a nice label for an error parameter
        """
        try:
            error = float(error)
            return "Error rate = {}".format(error) if error else label_for_no_error
        except (ValueError, TypeError):
            try: # make a simplified label
                if "Empirical" in error:
                    error = "With genotyping"
            except:
                pass
            return "{} error".format(error) if error else label_for_no_error


class StoringEveryone(Figure):
    """
    Figure showing how tree sequences can store the entire human population
    worth of variation data.
    """
    name = "storing_everyone"

    def plot(self):
        df = self.data
        df = df[df.sample_size > 10]

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        xytext = (18, 0)
        GB = 1024**3
        largest_n = np.array(df.sample_size)[-1]

        index = df.vcf > 0
        line, = ax1.loglog(df.sample_size, df.vcf_fit, "-", color="tab:pink", label="")
        ax1.loglog(
            df.sample_size[index], df.vcf[index], "d", label="VCF",
            color=line.get_color())
        largest_value = np.array(df.vcf_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.vcfz_fit, ":", label="", color=line.get_color())
        ax1.loglog(
            df.sample_size[index], df.vcfz[index], "d", label="Compressed VCF",
            color=line.get_color(), markerfacecolor='w')
        largest_value = np.array(df.vcfz_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        index = df.pbwt > 0
        line, = ax1.loglog(
            df.sample_size[index], df.pbwt[index], "-", color="tab:orange", label="")
        ax1.loglog(
            df.sample_size[index], df.pbwt[index], "s", label="pbwt",
            color=line.get_color())
        line, = ax1.loglog(
            df.sample_size[index], df.pbwtz[index], ":", label="", color=line.get_color())
        ax1.loglog(
            df.sample_size[index], df.pbwtz[index], "s", label="Compressed pbwt",
            color=line.get_color(), markerfacecolor='w')

        line, = ax1.loglog(
            df.sample_size, df.uncompressed, "o", label="Trees", color="b")
        ax1.loglog(df.sample_size, df.tsk_fit, "-", color=line.get_color(), label="")
        largest_value = np.array(df.tsk_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.tskz_fit, ":", label="", color=line.get_color())
        ax1.loglog(df.sample_size, df.compressed, "o", label="Compressed trees",
            color=line.get_color(), markerfacecolor='w')
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

        plt.clf()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        largest_n = 10**7
        index = df.sample_size <= largest_n
        # Rescale to MiB
        uncompressed = np.array(df.uncompressed[index] * 1024)
        compressed = np.array(df.compressed[index] * 1024)
        pbwt = np.array(df.pbwt[index] * 1024)
        pbwtz = np.array(df.pbwtz[index] * 1024)
        ax1.loglog(df.sample_size[index], uncompressed, "-o", label="trees")
        largest_value = uncompressed[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * 1024**2, binary=True, format="%.1f"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        ax1.loglog(df.sample_size[index], compressed, "-o", label="Compressed trees")
        largest_value = compressed[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * 1024**2, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value - 10), xycoords="data")

        ax1.loglog(df.sample_size[index], pbwt, "-s", label="pbwt")
        largest_value = pbwt[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * 1024**2, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        ax1.loglog(df.sample_size[index], pbwtz, "-s", label="Compressed pbwt")
        largest_value = pbwtz[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * 1024**2, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value + 10), xycoords="data")

        ax1.set_ylabel("File size (MiB)")
        ax1.set_xlabel("Number of chromosomes")
        plt.legend()
        self.save("storing_everyone_pbwt")


class SampleEdges(Figure):
    name = "sample_edges"

    def plot_region(self, df, ax, rotate_labels):
        ax.plot(df.sample_edges.values, "o")
        breakpoints = np.where(df.population.values[1:] != df.population.values[:-1])[0]
        breakpoints = np.array([-1] + list(breakpoints) + [len(df)-1])+0.5
        x_labels = []
        x_pos = []
        last = -0.5
        for bp in breakpoints[1:]:
            x_labels.append(df.population[int(bp - 1.5)])
            x_pos.append(last + (bp - last) / 2)
            last = bp
        # use major ticks for labels, so they are not cut off
        ax.tick_params(axis="x", which="major", length=0)
        ax.set_xticks(x_pos)
        if rotate_labels:
            ax.set_xticklabels(x_labels, rotation=30)
        else:
            ax.set_xticklabels(x_labels)
        ax.tick_params(axis="x", which="minor", length=0)
        ax.set_xticks(breakpoints, minor=True)
        ax.set_xticklabels([], minor=True)
        ax.set_xlim(-0.5, len(df) - 0.5)
        ax.set_title(df.region.unique()[0])
        ax.grid(which="minor", axis="x")


    def plot_summary(self):
        full_df = self.data
        fig, axes = plt.subplots(2, 1, figsize=(14, 6))
        plt.subplots_adjust(hspace=0.5)

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
                y = -400
                if dataset == "1kg":
                    y = -200
                ax.annotate(
                    df.region[bp - 1], xy=(x, y), horizontalalignment='center',
                    annotation_clip=False)
                last = bp

            breakpoints = np.where(
                df.population.values[1:] != df.population.values[:-1])[0]
            breakpoints = list(breakpoints) + [len(df)]
            ax.set_xticks(breakpoints)
            ax.set_xticklabels([])
            ax.set_ylabel("Sample Edges")
            ax.grid(axis="x")
            ax.set_xlim(0, len(df))
            ax.xaxis.set_ticks_position('none')

            title = "SGDP"
            if dataset == "1kg":
                title = "TGP"
                last = 0
                for bp in breakpoints:
                    x = last + (bp - last) / 2
                    last = bp
                    ax.annotate(
                        df.population[int(x)], xy=(x, 100), horizontalalignment='center',
                        annotation_clip=False)

                outliers = ["NA20289", "HG02789"]
                for outlier in outliers:
                    tmp_df = df[df["sample"] == outlier]
                    x = tmp_df.index.values[0]
                    ax.annotate(
                        outlier, xy=(x, 1550), horizontalalignment='center',
                        annotation_clip=False, style='italic')
            ax.set_title(title + " individuals")
        axes[0].set_ylim(0, 1500)
        axes[1].set_ylim(0, 3500)
        self.save("sample_edges_summary")

    def plot(self):
        self.plot_summary()
        full_df = self.data
        for ds in ["1kg", "sgdp"]:

            df_ds = full_df[full_df.dataset == ds]
            print(ds, "Overall", df_ds.sample_edges.mean(), sep="\t")
            fig, axes = plt.subplots(5, 1, figsize=(14, 16))
            plt.subplots_adjust(hspace=0.5)
            if ds == "1kg":
                plt.title("TGP sample edges per population")
            else:
                plt.title("SGDP sample edges per population")
            for ax, region  in zip(axes, df_ds.region.unique()):
                df = df_ds[df_ds.region == region]
                df = df.sort_values(by=["population", "sample", "strand"])
                df = df.reset_index()
                print(ds, region, df.sample_edges.mean(), sep="\t")
                self.plot_region(df, ax, rotate_labels=ds == "sgdp")

            self.save("{}_{}".format(self.name, ds))


class FrequencyDistanceAccuracy(Figure):
    """
    Plot accuracy of frequency ordering pairs of mutations vs distance between mutations
    The csv file is created by running
        python3 ./src/freq_dist_simulations.py
    or, if you have, say 40 processors available, you can run it in parallel like
        python3 -p 40 ./src/freq_dist_simulations.py

    """
    name = "frequency_distance_accuracy_singletons"

    def plot(self):
        df = self.data
        plt.plot((df.SeparationDistanceStart + df.SeparationDistanceEnd)/2/1e3,
            df.Agree/df.Total,label=self.error_label(None),
            color="k", linestyle="-")

        plt.plot((df.SeparationDistanceStart + df.SeparationDistanceEnd)/2/1e3,
            df.ErrorAgree/df.Total,label=self.error_label("EmpiricalError"),
            color="k", linestyle="-.")

        plt.xlabel("Distance between variants (kb)")
        plt.ylabel("Proportion of mutation pairs correctly ordered")
        plt.legend()

        self.save()


class AncestorAccuracy(Figure):
    """
    Compare lengths of real vs reconstructed ancestors, using 2 csv files generated by
    TSINFER_DIR=../tsinfer #set to your tsinfer directory
    python3 ${TSINFER_DIR}/evaluation.py aq -l 5 -d data -C -s 321 -e 0
    python3 ${TSINFER_DIR}/evaluation.py aq -l 5 -d data -C -s 321 -e data/EmpiricalErrorPlatinum1000G.csv
    cd data
    cat anc-qual_n=100_Ne=5000_L=5.0_mu=1e-08_rho=1e-08_err=data_EmpiricalErrorPlatinum1000G_error_data.csv > ancestor_accuracy.csv
    tail +2 anc-qual_n=100_Ne=5000_L=5.0_mu=1e-08_rho=1e-08_err=0.0_error_data.csv >> ancestor_accuracy.csv

    """ # noqa
    name = "ancestor_accuracy"

    def __init__(self):
        super().__init__()
        # rescale length to kb
        self.data["Real length"] /= 1e3
        self.data["Estim length"] /= 1e3
        # put high inaccuracy first
        self.data = self.data.sort_values("Inaccuracy")

    def plot(self):
        n_bins=50
        max_length = max(np.max(self.data["Real length"]), np.max(self.data["Estim length"]))* 1.1
        min_length = min(np.min(self.data["Real length"]), np.min(self.data["Estim length"])) * 0.9
        fig = plt.figure(figsize=(20, 8))
        gs = matplotlib.gridspec.GridSpec(1, 4, width_ratios=[5,5,0.5,2.5])
        ax0 = fig.add_subplot(gs[0])
        axes = [ax0, fig.add_subplot(gs[1], sharex=ax0, sharey=ax0, yticklabels=[])]
        c_ax = fig.add_subplot(gs[2])
        h_ax = fig.add_subplot(gs[3], sharey=c_ax)
        for ax, error in zip(axes, sorted(self.data.seq_error.unique())):
            df = self.data.query("seq_error == @error")
            ls = "-" if ax == axes[0] else "-."
            im = ax.scatter(df["Real length"], df["Estim length"], c=1-df["Inaccuracy"],
                s=20, cmap=matplotlib.cm.viridis)
            ax.plot([0, max_length], [0, max_length], '-',
                color='grey', zorder=-1, linestyle=ls)
            ax.set_title(self.error_label(error))
            ax.set_xscale('log')
            ax.set_yscale('log')
            n_greater_eq = sum(df["Estim length"]/df["Real length"] >= 1)
            n_less = sum(df["Estim length"]/df["Real length"] < 1)
            ax.text(min_length*1.1, min_length*2,
                "{} haplotypes $\geq$ true length".format(n_greater_eq),
                rotation=45, va='bottom', ha='left', color="#2ca02c")
            ax.text(min_length*2, min_length*1.1,
                "{} haplotypes $<$ true length".format(n_less),
                rotation=45, va='bottom', ha='left', color="#d62728")
            ax.set_aspect(1)
            ax.set_xlim(min_length, max_length)
            ax.set_ylim(min_length, max_length)
            ax.set_xlabel("True ancestral haplotype length (kb)")
            if ax == axes[0]:
                ax.set_ylabel("Inferred ancestral haplotype length (kb)")
            n, bins, patches = h_ax.hist(1-df["Inaccuracy"],
                bins=n_bins, orientation="horizontal", alpha=0.5,
                edgecolor='black', linewidth=1, linestyle=ls);
            norm = matplotlib.colors.Normalize(bins.min(), bins.max())
            # set a color for every bar (patch) according
            # to bin value from normalized min-max interval
            for bin, patch in zip(bins, patches):
                color = matplotlib.cm.viridis(norm(bin))
                patch.set_facecolor(color)
        c_ax.set_axes_locator(InsetPosition(axes[1], [1.05,0,0.05,1]))
        cbar = fig.colorbar(im, cax=c_ax)
        cbar.set_label("Accuracy", rotation=270, va="center")
        h_ax.set_axes_locator(InsetPosition(c_ax, [3.5,0,7,1]))
        h_ax.set_title("Accuracy distribution")
        h_ax.axis('off')
        self.save()


class ToolsFigure(Figure):
    """
    Superclass of all figures where different tools (e.g. ARGweaver, fastarg) are compared
    """

    # Colours taken from Matplotlib default color wheel.
    # https://matplotlib.org/users/dflt_style_changes.html
    tools_format = collections.OrderedDict([
        ("ARGweaver", {"mark":"*", "col":"#d62728"}),
        ("RentPlus",  {"mark":"d", "col":"#2ca02c"}),
        ("fastARG",   {"mark":"^", "col":"#ff7f0e"}),
        ("tsinfer",   {"mark":"o", "col":"#1f77b4"}),
    ])

    error_bars = True


class CputimeAllToolsBySampleSizeFigure(ToolsFigure):
    """
    Compare cpu times for tsinfer vs other tools. We can only really get the CPU times
    for all four methods in the same scale for tiny examples.
    We can show that ARGWeaver and RentPlus are much slower than tsinfer
    and FastARG here and compare tsinfer and FastARG more thoroughly
    in a dedicated figure.
    """
    name = "cputime_all_tools_by_sample_size"

    def plot(self):
        df = self.data
        # Scale time to hours
        time_scale = 3600
        df.cputime_mean /= time_scale
        df.cputime_se /= time_scale
        sample_sizes = df.sample_size.unique()
        fig, (ax_hi, ax_lo) = plt.subplots(2, 1, sharex=True)
        lengths = df.length.unique()
        # check these have fixed lengths
        assert len(lengths) == 1
        max_non_AW = 0
        for tool in df.tool.unique():
            line_data = df.query("tool == @tool")
            if tool != 'ARGweaver':
                max_non_AW = max(max_non_AW, max(line_data.cputime_mean+line_data.cputime_se))
            for ax in (ax_lo, ax_hi):
                ax.errorbar(
                    line_data.sample_size,
                    line_data.cputime_mean,
                    yerr=line_data.cputime_se,
                    color=self.tools_format[tool]["col"],
                    marker=self.tools_format[tool]['mark'],
                    elinewidth=1,
                    label=tool)
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax_hi.transAxes, color='k', clip_on=False)
        ax_hi.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax_hi.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
        ax_lo.set_xlabel("Sample Size")
        ax_hi.set_ylabel("CPU time (hours)")
        #ax_lo.set_xlim(sample_sizes.min(), sample_sizes.max())

        # zoom-in / limit the view to different portions of the data
        ax_hi.set_ylim(bottom = max_non_AW*40)  # outliers only
        ax_lo.set_ylim(bottom = 0-max_non_AW/20, top=max_non_AW+max_non_AW/20)  # most of the data
        #ax_hi.set_ylim(0.01, 3)  # outliers only
        #ax_lo.set_ylim(0, 0.002)  # most of the data

        # hide the spines between ax and ax2
        ax_hi.spines['bottom'].set_visible(False)
        ax_lo.spines['top'].set_visible(False)
        ax_hi.xaxis.tick_top()
        ax_hi.tick_params(labeltop=False)  # don't put tick labels at the top
        ax_lo.xaxis.tick_bottom()

        kwargs.update(transform=ax_lo.transAxes)  # switch to the bottom axes
        ax_lo.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax_lo.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        ax_hi.legend(loc="lower right")
        self.save()


class MemTimeFastargTsinferFigure(ToolsFigure):
    name = "mem_time_fastarg_tsinfer"
    def __init__(self):
        super().__init__()
        # Rescale the length to Mb
        length_scale = 10**6
        self.data.length /= length_scale
        # Scale time to hours
        time_scale = 3600
        cpu_cols = [c for c in self.data.columns if c.startswith("cputime")]
        self.data[cpu_cols] /= time_scale
        # Scale memory to GiB
        mem_cols = [c for c in self.data.columns if c.startswith("memory")]
        self.data[mem_cols] /= 1024 * 1024 * 1024
        length_sample_size_combos = self.data[["length", "sample_size"]].drop_duplicates()
        self.fixed_length = length_sample_size_combos['length'].value_counts().idxmax()
        self.fixed_sample_size = length_sample_size_combos['sample_size'].value_counts().idxmax()

    def plot(self):
        fig, axes = plt.subplots(2, 2, sharey="row", sharex="col", figsize=(8, 5.5), constrained_layout=True)
        xticks = [[0,2,4,6,8,10], [0,5000,10000,15000,20000]] #hardcode a hack here
        ylim_inset = [0.32, 0.35] #hardcode a hack here
        for i, (plotted_column, y_label) in enumerate(
                zip(["cputime", "memory"], ["CPU time (hours)", "Memory (GiB)"])):
            df = self.data.query("sample_size == @self.fixed_sample_size")
            df_inset = df.query("tool == 'tsinfer'")

            for tool in df.tool.unique():
                line_data = df.query("tool == @tool")
                axes[i][0].errorbar(
                    line_data.length,
                    line_data[plotted_column+"_mean"],
                    yerr=line_data[plotted_column+"_se"],
                    label=tool,
                    color=self.tools_format[tool]["col"],
                    marker=self.tools_format[tool]['mark'],
                    elinewidth=1)
            axes[i][0].set_xlim(0)
            axes[i][0].get_yaxis().set_label_coords(-0.08,0.5)
            axes[i][0].set_ylabel(y_label)

            tool = "tsinfer"
            axins1 = inset_axes(axes[i][0], width="50%", height="40%", loc=2, borderpad=1)
            axins1.errorbar(
                df_inset.length, df_inset[plotted_column+"_mean"],
                yerr=df_inset[plotted_column + "_se"], label=tool,
                color=self.tools_format[tool]["col"], marker=self.tools_format[tool]['mark'],
                linewidth=1, elinewidth=1, markersize=3)


            df = self.data.query("length == @self.fixed_length")
            df_inset = df.query("tool == 'tsinfer'")

            for tool in df.tool.unique():
                line_data = df.query("tool == @tool")
                axes[i][1].errorbar(
                    line_data.sample_size,
                    line_data[plotted_column+"_mean"],
                    yerr=line_data[plotted_column+"_se"],
                    label=tool,
                    color=self.tools_format[tool]["col"],
                    marker=self.tools_format[tool]['mark'],
                    elinewidth=1)

            tool = "tsinfer"
            axins2 = inset_axes(axes[i][1], width="50%", height="40%", loc=2, borderpad=1)
            axins2.errorbar(
                df_inset.sample_size, df_inset[plotted_column+"_mean"],
                yerr=df_inset[plotted_column+"_se"], label=tool,
                color=self.tools_format[tool]["col"], marker=self.tools_format[tool]['mark'],
                linewidth=1, elinewidth=1, markersize=3)

            axins1.tick_params(axis='both', which='major', labelsize=7)
            axins1.set_ylim(-0.01, ylim_inset[i])
            axins1.set_xticks(xticks[0])
            axins1.yaxis.tick_right()
            axins2.tick_params(axis='both', which='major', labelsize=7)
            axins2.set_ylim(-0.01, ylim_inset[i])
            axins2.set_xticks(xticks[1])
            axins2.yaxis.tick_right()
            #axins2.set_yticklabels(["{:.0f}".format(s) for s in yticks[i]])
        axes[0][0].legend(loc="upper right", numpoints=1, fontsize="small")
        axes[1][0].set_xlabel("Length (Mb) for fixed sample size of {}".format(
            self.fixed_sample_size))
        axes[1][1].set_xlabel("Sample size for fixed length of {:g} Mb".format(
            self.fixed_length))
        axes[1][0].set_xticks(xticks[0])
        axes[1][1].set_xticks(xticks[1])
        self.save()


class PerformanceLengthSamplesFigure(ToolsFigure):
    """
    The performance metrics figures. Each of these figures
    has two panels; one for scaling by sequence length and the other
    for scaling by sample size. Different lines are given for each
    of the different combinations of tsinfer parameters
    """
    y_name = "plotted_column"
    y_axis_label = None

    def __init__(self):
        super().__init__()
        # Rescale the length to Mb
        length_scale = 10**6
        self.data.length /= length_scale
        length_sample_size_combos = self.data[["length", "sample_size"]].drop_duplicates()
        self.fixed_length = length_sample_size_combos['length'].value_counts().idxmax()
        self.fixed_sample_size = length_sample_size_combos['sample_size'].value_counts().idxmax()

    def plot(self):
        df = self.data
        recombination_linestyles = [':', '-', '--']
        recombination_rates = np.sort(df.recombination_rate.unique())
        mutation_rates = df.mutation_rate.unique()
        tools = df.tool.unique()
        assert len(recombination_linestyles) >= len(recombination_rates)
        assert len(mutation_rates) == len(tools) == 1
        mu = mutation_rates[0]
        tool = tools[0]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
        ax1.set_title("Fixed number of chromosomes ({})".format(self.fixed_sample_size))
        ax1.set_xlabel("Sequence length (MB)")
        ax1.set_ylabel(self.y_axis_label)
        for linestyle, rho in zip(recombination_linestyles, recombination_rates):
            line_data = df.query("(sample_size == @self.fixed_sample_size) and (recombination_rate == @rho)")
            ax1.errorbar(
                line_data.length, line_data[self.plotted_column + "_mean"],
                yerr=None, # line_data[self.plotted_column + "_se"],
                linestyle=linestyle,
                color=self.tools_format[tool]["col"],
                #marker=self.tools_format[tool]['mark'],
                #elinewidth=1
                )


        ax2.set_title("Fixed sequence length ({:g} Mb)".format(self.fixed_length))
        ax2.set_xlabel("Sample size")
        for linestyle, rho in zip(recombination_linestyles, recombination_rates):
            line_data = df.query("(length == @self.fixed_length) and (recombination_rate == @rho)")
            ax2.errorbar(
                line_data.sample_size, line_data[self.plotted_column + "_mean"],
                yerr=None, # line_data[self.plotted_column + "_se"],
                linestyle=linestyle,
                color=self.tools_format[tool]["col"],
                #marker=self.tools_format[tool]['mark'],
                #elinewidth=1
                )
        params = [
            plt.Line2D((0,0),(0,0), color=self.tools_format[tool]["col"],
                linestyle=linestyle, linewidth=2)
            for linestyle, rho in zip(recombination_linestyles, recombination_rates)]
        ax1.legend(
            params,
            [r"$\rho$ = {}".format("$\mu$" if rho==mu else
                r"{:g}$\mu$ (more recombination)".format(rho/mu) if rho>mu else
                r"$\mu$/{:g} (less recombination)".format(mu/rho))
                for rho_index, rho in enumerate(recombination_rates)],
            loc="upper right", fontsize=10, title="Relative rate of recombination")

        self.save()


class TSCompressionFigure(PerformanceLengthSamplesFigure):
    name = "tsinfer_ts_filesize_ln"
    plotted_column = "ts_relative_filesize"
    y_axis_label = "File size relative to simulated tree sequence"


class VCFCompressionFigure(PerformanceLengthSamplesFigure):
    name = "tsinfer_vcf_compression_ln"
    plotted_column = "vcf_compression_factor"
    y_axis_label = "Compression factor relative to vcf.gz"


class NodesWithMutations(PerformanceLengthSamplesFigure):
    name = "mutation_ancestors_ln"
    plotted_column = "prop_nodes_with_mutation"
    y_axis_label = "Proportion of internal nodes associated with at least one mutation"


class TreeMetricsFigure(ToolsFigure):

    metric_titles = {
        "wRF": "weighted Robinson-Foulds metric",
        "RF": "Robinson-Foulds metric",
        "SPR": "estimated SPR difference",
        "path": "Path difference",
        "KC": "Kendall-Colijn metric",
    }

    polytomy_and_averaging_format = collections.OrderedDict([
        ("broken", {
            "per site":    {"linestyle":"--"},
            "per variant": {"linestyle":":"}}),
        ("retained", {
            "per site":    {"linestyle":"-"},
            "per variant": {"linestyle":"-."}})
    ])

    sample_size_format = [
        {'fillstyle':'full'}, #smaller ss
        {'fillstyle':'none'} #assume only max 2 sample sizes per plot
        ]

    length_format = {'tsinfer':[{'col':'k'}, {'col':'#1f77b4'}, {'col':'#17becf'}]}

    def single_metric_plot(self, df, x_variable, ax, av_method,
        rho = None, markers = True, x_jitter = None):
        """
        A single plot on an ax. This requires plotting separate lines, e.g. for each tool
        If rho is give, plot an x=rho vertical line, assuming x is the mutation_rate.
        x_jitter can be None, 'log' or 'linear'
        """
        v_cols = ['length', 'sample_size', 'tool', 'polytomies']
        v_order = df[v_cols].drop_duplicates() # find unique combinations
        # sort for display
        v_order = v_order.sort_values(v_cols, ascending=[False, True, True, False])
        ss_order = {v:k for k,v in enumerate(v_order.sample_size.unique())}
        l_order = {v:k for k,v in enumerate(v_order.length.unique())}
        for i, r in enumerate(v_order.itertuples()):
            query = []
            query.append("length == @r.length")
            query.append("sample_size == @r.sample_size")
            query.append("tool == @r.tool")
            query.append("polytomies == @r.polytomies")
            line_data = df.query("(" + ") and (".join(query) + ")")
            if not line_data.empty:
                if len(v_order.length.unique()) > 1:
                    # all tsinfer tools: use colours for length for polytomy format
                    colour = self.length_format[r.tool][l_order[r.length]]["col"]
                else:
                    # no variable lengths: use standard tool colours
                    colour = self.tools_format[r.tool]["col"]
                x = line_data[x_variable]
                if x_jitter:
                    if x_jitter == 'log':
                        x *= 1 + (2*i/len(v_order)-1) * (max(x)/min(x))/5000
                    else:
                        x += (2 * i - 1) * (max(x)-min(x))/400
                ax.errorbar(
                    x, line_data.treedist_mean,
                    yerr=line_data.treedist_se if self.error_bars else None,
                    linestyle=self.polytomy_and_averaging_format[r.polytomies][av_method]["linestyle"],
                    fillstyle=self.sample_size_format[ss_order[r.sample_size]]['fillstyle'],
                    color=colour,
                    marker=self.tools_format[r.tool]['mark'] if markers else None,
                    elinewidth=1)
        if rho is not None:
            ax.axvline(x=rho, color = 'gray', zorder=-1, linestyle=":", linewidth=1)
            ax.text(rho, ax.get_ylim()[1]/40,  r'$\mu=\rho$',
                va="bottom",  ha="right", color='gray', rotation=90)
        return v_order

class MetricsAllToolsFigure(TreeMetricsFigure):
    """
    Simple figure that shows all the metrics at the same time.
    Assumes at most 2 sample sizes
    """
    name = "metrics_all_tools"

    def plot(self):
        averaging_method = self.data.averaging.unique()
        eff_sizes = self.data.Ne.unique()
        rhos = self.data.recombination_rate.unique()
        lengths = self.data.length.unique()
        assert len(averaging_method) == len(eff_sizes) == len(rhos) == 1
        rho = rhos[0]
        method = averaging_method[0]

        sample_sizes = self.data.sample_size.unique()
        # x-direction is different error rates
        seq_error_params = self.data.error_param.unique()
        # y-direction is the permutations of metric + whether it is rooted
        metric_and_rooting = self.data.groupby(["metric", "rooting"]).groups
        metric_and_rooting = collections.OrderedDict( # make a consistent order
            sorted(metric_and_rooting.items(), key=lambda x: x[0]))
        # sort this so that metrics come out in a set order (TO DO)
        fig, axes = plt.subplots(len(metric_and_rooting), len(seq_error_params),
            squeeze=False, sharey='row',
            figsize=(6*len(seq_error_params), 15))
        for j, ((metric, root), rows) in enumerate(metric_and_rooting.items()):
            for k, error in enumerate(seq_error_params):
                # we are in the j,k th subplot
                ax = axes[j][k]
                ax.set_xscale('log')
                display_order = self.single_metric_plot(
                    self.data.loc[rows].query("error_param == @error"), "mutation_rate",
                    ax, method, rho, markers = (len(sample_sizes)!=1))
                # Use integers for labels
                ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
                # Make the labels on each Y axis line up properly
                ax.get_yaxis().set_label_coords(-0.08,0.5)
                if j == 0:
                    ax.set_title(self.error_label(error))
                if j == len(metric_and_rooting) - 1:
                    ax.set_xlabel("Mutation rate")
                if k == 0:
                    ax.set_ylim(getattr(self,'ylim', 0))
                    rooting_suffix = " (unrooted)" if root=="unrooted" else ""
                    ylab = getattr(self, 'y_axis_label', self.metric_titles[metric] + rooting_suffix)
                    ax.set_ylabel(ylab)

        artists = [
            plt.Line2D((0,1),(0,0), linewidth=2,
                color=self.tools_format[d.tool]["col"],
                linestyle=self.polytomy_and_averaging_format[d.polytomies][method]["linestyle"],
                marker = None if len(sample_sizes)==1 else self.tools_format[d.tool]['mark'])
            for d in display_order[['tool', 'polytomies']].drop_duplicates().itertuples()]
        tool_labels = [d.tool + ("" if d.polytomies == "retained" else (" (polytomies " + d.polytomies + ")"))
            for d in display_order[['tool', 'polytomies']].drop_duplicates().itertuples()]
        axes[0][0].legend(
            artists, tool_labels, numpoints=1, labelspacing=0.1)

        fig.tight_layout()
        self.save()


class MetricsAllToolsAccuracyFigure(MetricsAllToolsFigure):
    """
    Show the metrics tending to 0 as mutation rate increases
    """
    name = "metrics_all_tools_accuracy"


class MetricAllToolsFigure(TreeMetricsFigure):
    """
    Plot each metric in a different pdf file.
    For the publication: make symbols small and solid
    """
    plot_height = 4.5
    name = "metric_all_tools"
    #y_axis_label="Average distance from true trees"
    hide_polytomy_breaking = True
    output_metrics = [("KC","rooted"), ("RF", "rooted")] #can add extras in here if necessary

    def plot(self):
        if getattr(self,"hide_polytomy_breaking", None):
            df = self.data.query("polytomies != 'broken'")
        else:
            df = self.data
        for metric, rooting in self.output_metrics:
            query = ["metric == @metric", "rooting == @rooting"]
            averaging_method = df.averaging.unique()
            eff_sizes = df.Ne.unique()
            rhos = df.recombination_rate.unique()
            lengths = df.length.unique()
            assert len(averaging_method) == len(eff_sizes) == len(rhos) == 1
            rho = rhos[0]
            method = averaging_method[0]

            # x-direction is different sequencing error rates
            try:
                seq_error_params = df.error_param.unique()
            except AttributeError:
                seq_error_params = [0]

            # y-direction is different ancestral allele error rates (if present)
            try:
                aa_error_params = self.data.ancestral_state_error_param.unique()
            except AttributeError:
                aa_error_params = [0]

            fig, axes = plt.subplots(len(aa_error_params), len(seq_error_params),
                squeeze=False, sharey=True,
                figsize=getattr(self,'figsize',(6*len(seq_error_params), self.plot_height)))
                # Used for final version of Figure 3.
                # figsize=(7, 3.5))
            for j, aa_error in enumerate(aa_error_params):
                for k, seq_error in enumerate(seq_error_params):
                    ax = axes[j][k]
                    subquery = query
                    if len(seq_error_params) > 1:
                        subquery.append("error_param == @seq_error")
                    if len(aa_error_params) > 1:
                        subquery.append("ancestral_state_error_param == @aa_error")
                    display_order = self.single_metric_plot(
                        df.query("(" + ") and (".join(subquery) + ")"),
                        "mutation_rate", ax, method, rho)
                    ax.set_title(self.error_label(seq_error))
                    ax.set_xlabel("Mutation rate")
                    ax.set_xscale('log')
                    if k == 0:
                        ax.set_ylim(getattr(self,'ylim', 0))
                        rooting_suffix = " (unrooted)" if rooting=="unrooted" else ""
                        ylab = getattr(self, 'y_axis_label', None)
                        if ylab is None:
                            ylab = metric + rooting_suffix + " metric"
                        if len(aa_error_params)>1 or aa_error != 0:
                            ylab += ", {:g}% ancestral state error".format(aa_error*100)
                        ax.set_ylabel(ylab)

            # Create legends from custom artists
            artists = [
                plt.Line2D((0,1),(0,0),
                    color=self.tools_format[d.tool]["col"],
                    linestyle=self.polytomy_and_averaging_format[d.polytomies][method]["linestyle"],
                    marker = self.tools_format[d.tool]['mark'])
                for d in display_order[['tool', 'polytomies']].drop_duplicates().itertuples()]
            tool_labels = [d.tool + ("" if d.polytomies == "retained" else (" (polytomies " + d.polytomies + ")"))
                for d in display_order[['tool', 'polytomies']].drop_duplicates().itertuples()]
            first_legend = axes[0][0].legend(
                artists, tool_labels, numpoints=1, labelspacing=0.1, loc="upper right")
            fig.tight_layout()
            if len(self.output_metrics)==1:
                self.save()
            else:
                self.save("_".join([self.name, metric, rooting]))


class MetricAllToolsAccuracyBadAncestorsSummary(MetricAllToolsFigure):
    name = "metric_all_tools_accuracy_bad_ancestors"
    output_metrics = [("KC","rooted"), ("RF", "rooted")]
    y_axis_label = None
    hide_polytomy_breaking = False
    plot_height = 10.5


class MetricAllToolsAccuracyDemographyFigure(MetricAllToolsFigure):
    """
    Simple figure that shows an ARG metrics for a genome under a more complex demographic
    model (the Gutenkunst Out Of Africa model), as mutation rate increases to high values
    """
    plot_height = 3.5 # To more-or-less match to other supplementary figures
    name = "metric_all_tools_accuracy_demography"
    hide_polytomy_breaking = True


class MetricAllToolsAccuracySweepFigure(TreeMetricsFigure):
    """
    Figure for simulations with selection.
    Each page should be a single figure for a particular metric, with error on the
    """
    name = "metric_all_tools_accuracy_sweep"
    error_bars = True
    hide_polytomy_breaking = True
    output_metrics = [("KC","rooted")] #can add extras in here if necessary

    def plot(self):
        if getattr(self,"hide_polytomy_breaking", None):
            df = self.data.query("polytomies != 'broken'")
        else:
            df = self.data
        for metric, rooting in self.output_metrics:
            df = df.query("(metric == @metric) and (rooting == @rooting)")
            output_freqs = df[['output_frequency', 'output_after_generations']].drop_duplicates()
            averaging_method = df.averaging.unique()
            eff_sizes = df.Ne.unique()
            rhos = df.recombination_rate.unique()
            lengths = df.length.unique()
            assert len(averaging_method) == len(eff_sizes) == len(rhos) == 1
            rho = rhos[0]
            method = averaging_method[0]
            # x-direction is different error rates
            seq_error_params = df.error_param.unique()
            fig, axes = plt.subplots(len(output_freqs), len(seq_error_params),
                figsize=getattr(self,'figsize',(6*len(seq_error_params), 2.5*len(output_freqs))),
                squeeze=False, sharey=True)
            for j, output_data in enumerate(output_freqs.itertuples()):
                for k, error in enumerate(seq_error_params):
                    ax = axes[j][k]
                    freq = output_data.output_frequency
                    gens = output_data.output_after_generations
                    query = ["error_param == @error"]
                    query.append("output_frequency == @freq")
                    query.append("output_after_generations == @gens")
                    display_order = self.single_metric_plot(
                        df.query("(" + ") and (".join(query) + ")"),
                        "mutation_rate", ax, method, rho)
                    ax.set_xscale('log')
                    if j == 0:
                        ax.set_title(self.error_label(error))
                    if j == len(output_freqs) - 1:
                        ax.set_xlabel("Neutral mutation rate")
                    if k == 0:
                        ax.set_ylabel(getattr(self, 'y_axis_label', metric + " metric") +
                            " @ {}{}".format(
                                "fixation " if np.isclose(freq, 1.0) else "freq {}".format(freq),
                                "+{} gens".format(int(gens)) if gens else ""))
                    if np.isclose(freq, 1.0) and not gens:
                        # This is *at* fixation - set the plot background colour
                        ax.set_facecolor('0.9')
            # Create legends from custom artists
            artists = [
                plt.Line2D((0,1),(0,0),
                    color=self.tools_format[d.tool]["col"],
                    linestyle=self.polytomy_and_averaging_format[d.polytomies][method]["linestyle"],
                    marker = self.tools_format[d.tool]['mark'])
                for d in display_order[['tool', 'polytomies']].drop_duplicates().itertuples()]
            tool_labels = [d.tool + ("" if d.polytomies == "retained" else (" (polytomies " + d.polytomies + ")"))
                for d in display_order[['tool', 'polytomies']].drop_duplicates().itertuples()]
            first_legend = axes[0][0].legend(
                artists, tool_labels, numpoints=1, labelspacing=0.1, loc="upper right")
            fig.tight_layout()
            if len(self.output_metrics)==1:
                self.save()
            else:
                self.save("_".join([self.name, metric, rooting]))

class MetricSubsamplingFigure(TreeMetricsFigure):
    """
    Figure that shows whether increasing sample size helps with the accuracy of
    reconstructing the ARG for a fixed subsample. We only use tsinfer for this.
    """
    name = "metric_subsampling"
    hide_polytomy_breaking = True
    output_metrics = [("KC","rooted")] #can add extras in here if necessary

    def plot(self):
        self.polytomy_and_averaging_format['retained']['per variant']['linestyle'] = "-"
        for metric, rooting in self.output_metrics:
            query = ["metric == @metric", "rooting == @rooting"]
            if getattr(self,"hide_polytomy_breaking", None):
                query.append("polytomies != 'broken'")
            df = self.data.query("(" + ") and (".join(query) + ")")
            subsample_size = df.subsample_size.unique()
            # all should have the same similarion sample size (but inference occurs on
            # different subsample sizes, and tree comparisons on a fixed small tip #.
            averaging_method = self.data.averaging.unique()
            sample_sizes = df.sample_size.unique()
            all_tree_tips = df.restrict_sample_size_comparison.unique()
            mutation_rates = df.mutation_rate.unique()
            assert len(all_tree_tips) == len(mutation_rates) == len(sample_sizes) == len(averaging_method) == 1
            tree_tips = all_tree_tips[0]
            method = averaging_method[0]
            lengths = df.length.unique()
            seq_error_params = df.error_param.unique()
            fig, axes = plt.subplots(1, len(seq_error_params),
                figsize=(12, 6), squeeze=False, sharey=True)
            for k, error in enumerate(seq_error_params):
                ax = axes[0][k]
                display_order = self.single_metric_plot(
                    df.query("error_param == @error"), "subsample_size",
                    ax, method, rho = None, markers = False, x_jitter = 'log')
                ax.set_title(self.error_label(error))
                if k == 0:
                    ylab = getattr(self, 'y_axis_label', self.metric_titles[metric])
                    ax.set_ylabel("Average " + ylab + " for trees reduced to first {} samples".format(tree_tips))
                ax.set_xlabel("Sample size used for inference")
                ax.set_xscale('log')
            if len(display_order)>1:
                l_order = {v:k for k,v in enumerate(display_order.length.unique())}
                artists = [
                    plt.Line2D((0,1),(0,0),
                        color=self.length_format[d.tool][l_order[d.length]]["col"],
                        linestyle=self.polytomy_and_averaging_format[d.polytomies][method]["linestyle"],
                        marker = False)
                    for d in display_order[['length', 'tool', 'polytomies']].drop_duplicates().itertuples()]
                labels = ["{} kb".format(d.length//1000)
                    for d in display_order[['length']].drop_duplicates().itertuples()]
                first_legend = axes[0][0].legend(
                    artists, labels, numpoints=1, labelspacing=0.1, loc="upper right")
            fig.tight_layout()
            if len(self.output_metrics)==1:
                self.save()
            else:
                self.save("_".join([self.name, metric, rooting]))




def rotate_linkage(linkage, index):
    x, y = linkage[index][0:2]
    linkage[index][0] = y
    linkage[index][1] = x


class UkbbStructureFigure(Figure):
    """
    Figure showing the structure for UKBB using heatmaps.
    """
    name = "ukbb_structure"
    def __init__(self):
        # We don't have a CSV called this, name.csv so skip loading.
        pass

    def plot_ukbb_region_clustermap(self):
        df = pd.read_csv("data/ukbb_ukbb_british_centre.csv").set_index("centre")

        # Zscore normalise
        for col in list(df):
            df[col] = scipy.stats.zscore(df[col])

        row_linkage = scipy.cluster.hierarchy.linkage(df, method="average", optimal_ordering=True)

        # Tweaks to the clade rotation order
        # Flip the north vs north-west cities
        rotate_linkage(row_linkage, -3)

        # Put Welsh towns next to Birmingham
        rotate_linkage(row_linkage, -8)

        # Do Leeds - Sheffield - Nottingham, not Nottingham - Sheffield - Leeds
        # (this simply helps the top-to-bottom labelling scheme
        rotate_linkage(row_linkage, -9)
        rotate_linkage(row_linkage, -11)

        # Bristol near Welsh towns
        rotate_linkage(row_linkage, -12)

        # Swansea then Cardiff (highlights association between Cardiff & Bristol)
        rotate_linkage(row_linkage, -15)

        # The south/south-east centres are a bit visually messy - try putting east at end
        # Oxford near Bristol
        rotate_linkage(row_linkage, -16)
        # Reading nearer Oxford
        rotate_linkage(row_linkage, -20)
        # Push Croydon (furthest east) to the end
        rotate_linkage(row_linkage, -21)

        order = scipy.cluster.hierarchy.leaves_list(row_linkage)
        x_pop = df.index.values[order]

        cg = sns.clustermap(df[x_pop], row_linkage=row_linkage, col_cluster=False, rasterized=True)
        cg.ax_heatmap.set_ylabel("")
        for tick in cg.ax_heatmap.get_xticklabels():
            tick.set_rotation(-45)
            tick.set_ha('left')
            tick.set_rotation_mode("anchor")
        self.save("ukbb_ukbb_clustermap_british")

    def plot_1kg_ukbb_clustermap(self):

        df = pd.read_csv("data/1kg_ukbb_ethnicity.csv").set_index("ethnicity")

        colour_map = get_tgp_colours()
        colours = [colour_map[pop] for pop in tgp_populations]
        df = df[tgp_populations]

        row_linkage = scipy.cluster.hierarchy.linkage(df, method="average")

        cg = sns.clustermap(
            df, row_linkage=row_linkage, col_cluster=False, col_colors=colours, rasterized=True)
        cg.ax_heatmap.set_ylabel("")
        for region, col in get_tgp_region_colours().items():
            cg.ax_col_dendrogram.bar(0, 0, color=col, label=region, linewidth=0)
        cg.ax_col_dendrogram.legend(bbox_to_anchor=(1.2, 0.8))

        self.save("ukbb_1kg_clustermap_ethnicity")


    def plot_1kg_ukbb_british_centre_clustermap(self):

        df = pd.read_csv("data/1kg_ukbb_british_centre.csv").set_index("centre")

        colour_map = get_tgp_colours()
        colours = [colour_map[pop] for pop in tgp_populations]
        df = df[tgp_populations]

        # Zscore normalise
        for col in list(df):
            df[col] = scipy.stats.zscore(df[col])

        row_linkage = scipy.cluster.hierarchy.linkage(df, method="average")

        cg = sns.clustermap(
            df, row_linkage=row_linkage, col_cluster=False, col_colors=colours, rasterized=True)
        cg.ax_heatmap.set_ylabel("")
        for region, col in get_tgp_region_colours().items():
            cg.ax_col_dendrogram.bar(0, 0, color=col, label=region, linewidth=0)
        cg.ax_col_dendrogram.legend(bbox_to_anchor=(1.2, 0.8))

        self.save("ukbb_1kg_clustermap_british_centre")


    def plot(self):
        self.plot_ukbb_region_clustermap()
        self.plot_1kg_ukbb_clustermap()
        self.plot_1kg_ukbb_british_centre_clustermap()


class GlobalStructureFigure(Figure):
    name = "global_structure"

    def __init__(self):
        # We don't have a CSV called this, name.csv so skip loading.
        pass


    def plot_clustermap(self, df, pop_colours, region_colours, figsize):

        dfg = df.groupby("population").mean()
        # Zscore normalise
        for col in list(dfg):
            dfg[col] = scipy.stats.zscore(dfg[col])

        row_linkage = scipy.cluster.hierarchy.linkage(dfg, method="average")
        order = scipy.cluster.hierarchy.leaves_list(row_linkage)
        x_pop = dfg.index.values[order]

        colours = pd.Series(pop_colours)
        cg = sns.clustermap(
            dfg[x_pop], row_linkage=row_linkage, col_cluster=False,
            row_colors=colours, figsize=figsize, rasterized=True)
        cg.ax_heatmap.set_ylabel("")

        for region, col in region_colours.items():
            cg.ax_col_dendrogram.bar(0, 0, color=col, label=region, linewidth=0)
        return cg

    def plot_1kg_clustermap(self):
        df = pd.read_csv("data/1kg_gnn.csv")
        cg = self.plot_clustermap(
            df, get_tgp_colours(), get_tgp_region_colours(), figsize=(10, 10))
        cg.ax_col_dendrogram.legend(ncol=5)
        self.save("1kg_gnn_clustermap")

    def plot_sgdp_clustermap(self):
        df = pd.read_csv("data/sgdp_gnn.csv")
        colours = {}
        region_colours = get_sgdp_region_colours()
        for pop, region in df.groupby(["population", "region"]).size().index:
            colours[pop] = region_colours[region]

        cg = self.plot_clustermap(df, colours, region_colours, figsize=(10, 10))
        cg.ax_col_dendrogram.legend(ncol=2)
        cg.ax_heatmap.set_yticks([])
        cg.ax_heatmap.set_xticks([])
        self.save("sgdp_gnn_clustermap")
        cg = self.plot_clustermap(df, colours, region_colours, figsize=(30, 30))
        cg.ax_col_dendrogram.legend(ncol=2)
        self.save("sgdp_gnn_clustermap_full_size")

    def plot_pel_population(self):
        # Same as plot_composite but we don't include the fancy annotations.

        colours = get_tgp_region_colours()
        region_order = ['EUR', 'EAS', 'SAS', 'AFR', 'AMR']
        full_df = pd.read_csv("data/1kg_gnn.csv")
        df = full_df[full_df.population == "PEL"].reset_index()
        A = np.zeros((len(tgp_region_pop), len(df)))

        regions = ['EUR', 'EAS', 'SAS', 'AFR', 'AMR']
        for j, region in enumerate(regions):
            A[j, :] = np.sum([df[pop].values for pop in tgp_region_pop[region]], axis=0)

        focal_ind = "HG01933"
        index = np.argsort(A[0])[::-1]
        inds = df.individual.values

        inds = list(inds[index])
        x1 = inds.index(focal_ind)
        x2 = inds.index(focal_ind, x1 + 1)

        A = A[:, index]

        fig, ax = plt.subplots(1, 1, figsize=(17, 3))
        x = np.arange(len(df))
        for j, region in enumerate(regions):
            ax.bar(
                x, A[j], bottom=np.sum(A[:j, :], axis=0), label=region,
                width=1, color=colours[region], align="edge")
        ax.set_xlim(0, len(df) - 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        self.save("1kg_gnn_pel")

        for j in range(2):
            df = pd.read_csv("data/HG01933_local_gnn_{}.csv".format((j + 1) % 2), index_col = 0)
            left = df["left"]
            width = np.diff(np.append(left, df.right.max()))
            total = np.zeros_like(width)
            fig, ax = plt.subplots(1, 1, figsize=(17, 1.5))
            for region in region_order:
                ax.bar(
                    left, df[region], bottom=total, width=width, align="edge",
                    label=region, color=colours[region])
                total += df[region]
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(0, df.right.max())
            ax.set_ylim(0, 1)
            ax.axis('off')
            self.save("1kg_gnn_HG01933_{}".format(j))

    def plot(self):
        self.plot_pel_population()
        self.plot_1kg_clustermap()
        self.plot_sgdp_clustermap()
        self.plot_composite()

    def plot_composite(self):

        colours = get_tgp_region_colours()
        gs = matplotlib.gridspec.GridSpec(2, 2, height_ratios=[4, 1], hspace=0.6)
        fig = plt.figure(figsize=(14, 4))

        full_df = pd.read_csv("data/1kg_gnn.csv")
        df = full_df[full_df.population == "PEL"].reset_index()
        A = np.zeros((len(tgp_region_pop), len(df)))

        regions = ['EUR', 'EAS', 'SAS', 'AFR', 'AMR']
        for j, region in enumerate(regions):
            A[j, :] = np.sum([df[pop].values for pop in tgp_region_pop[region]], axis=0)

        focal_ind = "HG01933"
        index = np.argsort(A[0])[::-1]
        inds = df.individual.values

        inds = list(inds[index])
        x1 = inds.index(focal_ind)
        x2 = inds.index(focal_ind, x1 + 1)

        A = A[:, index]

        ax = plt.subplot(gs[0,:])
        x = np.arange(len(df))
        for j, region in enumerate(regions):
            ax.bar(
                x, A[j], bottom=np.sum(A[:j, :], axis=0), label=region, width=1,
                color=colours[region], align="edge")
        ax.set_xlim(0, len(df) - 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_ylabel("GNN Fraction")
        ax.set_title("TGP PEL population haplotypes sorted by GNN")
        for x in [x1, x2]:
            p = matplotlib.patches.Rectangle(
                (x, 0), width=1, height=1, fill=False, linestyle="--", color="grey")
            ax.add_patch(p)

        ax.legend(bbox_to_anchor=(1.02, 0.76))
        # ax.annotate(
        #     "(C)".format(label), xy=(-0.1, 0.5), xycoords="axes fraction", fontsize=15)
        ax_pop = ax
        ax_left = plt.subplot(gs[1, 0])
        ax_right = plt.subplot(gs[1, 1])
        for j, ax in enumerate([ax_left, ax_right]):
            df = pd.read_csv("data/HG01933_local_gnn_{}.csv".format((j + 1) % 2))
            left = df.left
            width = df.right - left
            total = np.zeros_like(width)
            for region in regions:
                ax.bar(
                    left, df[region].values, bottom=total, width=width, align="edge",
                    label=region, color=colours[region])
                total += df[region].values
            ax.set_title("HG01933 haplotype ({})".format(j + 1))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(0, df.right.max())
            ax.set_ylim(0, 1)

            centre_index = np.where(df.left.values[1:] != df.right.values[:-1])[0][0]
            start = df.right.iloc[centre_index]
            end = df.left.iloc[centre_index + 1]

            L = df.right.max()
            size = L / 18
            paths = [
                [(0, 0), (-size, 0.5), (0, 1), (0, 0)],
                [(start, 0), (start + size, 0.5), (start, 1), (start, 0)],
                [(end, 0), (end - size, 0.5), (end, 1), (end, 0)],
                [(L, 0), (L + size, 0.5), (L, 1), (L, 0)],
            ]
            for path in paths:
                Path = matplotlib.path.Path
                pp = matplotlib.patches.PathPatch(
                    Path(path, [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CLOSEPOLY]),
                    transform=ax.transData, clip_on=False,
                    facecolor='grey', alpha=0.5)
                ax.add_patch(pp)

            ax.axis('off')
        # Note: we have to manually tweak the label on this axis left a bit
        # ax_left.annotate(
        #     "(D)".format(label), xy=(-0.225, 0.5), xycoords="axes fraction", fontsize=15)

        L = df.right.max()
        transFigure = fig.transFigure.inverted()
        endpoints = [
            (ax_pop.transData.transform([x1, 0]), ax_left.transData.transform([0, 1])),
            (ax_pop.transData.transform([x1 + 1, 0]), ax_left.transData.transform([L, 1])),
            (ax_pop.transData.transform([x2, 0]), ax_right.transData.transform([0, 1])),
            (ax_pop.transData.transform([x2 + 1, 0]), ax_right.transData.transform([L, 1])),
        ]
        for (a, b) in endpoints:
            coord1 = transFigure.transform(a)
            coord2 = transFigure.transform(b)
            line = matplotlib.lines.Line2D(
                (coord1[0], coord2[0]),(coord1[1], coord2[1]), transform=fig.transFigure,
                linestyle="--", color="grey")
            fig.lines.append(line)

        self.save()


######################################
#
# Helper functions
#
######################################


def get_subclasses(cls):
    for subclass in cls.__subclasses__():
        yield from get_subclasses(subclass)
        yield subclass

def latex_float(f):
    """
    Return an exponential number in nice LaTeX form.
    In titles etc for plots this needs to be encased in $ $ signs, and r'' strings used
    """
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

######################################
#
# Main
#
######################################

def main():
    figures = list(get_subclasses(Figure))

    name_map = {fig.name: fig for fig in figures if fig.name is not None}

    parser = argparse.ArgumentParser(description="Make the plots for specific figures.")
    parser.add_argument(
        "name", type=str, help="figure name",
        choices=sorted(list(name_map.keys()) + ['all']))

    args = parser.parse_args()
    if args.name == 'all':
        for name, fig in name_map.items():
            if fig in figures:
                fig().plot()
    else:
        fig = name_map[args.name]()
        fig.plot()


if __name__ == "__main__":
    main()
