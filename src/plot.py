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
import seaborn as sns


class Figure(object):
    """
    Superclass of figures for the paper. Each figure is a concrete subclass.
    """
    name = None

    def __init__(self):
        datafile_name = "data/{}.csv".format(self.name)
        self.data = pd.read_csv(datafile_name)

    def save(self, figure_name=None):
        if figure_name is None:
            figure_name = self.name
        print("Saving figure '{}'".format(figure_name))
        plt.savefig("figures/{}.pdf".format(figure_name))
        plt.savefig("figures/{}.png".format(figure_name))
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
        line, = ax1.loglog(df.sample_size[index], df.vcf[index], "^", label="vcf")
        ax1.loglog(df.sample_size, df.vcf_fit, "--", color=line.get_color(), label="")
        largest_value = np.array(df.vcf_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(df.sample_size[index], df.vcfz[index], "s", label="vcf.gz")
        ax1.loglog(df.sample_size, df.vcfz_fit, "--", color=line.get_color(), label="")
        largest_value = np.array(df.vcfz_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.uncompressed, "o", label="trees")
        ax1.loglog(df.sample_size, df.tsk_fit, "--", color=line.get_color(), label="")
        largest_value = np.array(df.tsk_fit)[-1]
        ax1.annotate(
            humanize.naturalsize(largest_value * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(largest_n, largest_value), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.compressed, "*", label="trees.gz")
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
        fig = plt.figure(figsize=(14, 6))
        ax = fig.add_subplot(111)
        ax.plot(df.sample_edges.values, "o")
        breakpoints = np.where(df.population.values[1:] != df.population.values[:-1])[0]
        breakpoints = list(breakpoints) + [len(df)]
        last = 0
        y = df.sample_edges.min()
        for bp in breakpoints:
            label = df.population[bp - 1]
            x = last + (bp - last) / 2
            ax.annotate(
                label, xy=(x, y), horizontalalignment='centre', verticalalignment='top',
                rotation=90)
            last = bp
        ax.set_xticks(np.array(breakpoints) + 0.5)
        ax.set_xticklabels([])
        ax.set_xlim(-0.5, len(df) - 0.5)
        ax.set_title("{}:{}".format(dataset.upper(), region))
        ax.grid(axis="x")
        self.save("{}_{}_{}".format(self.name, dataset, region))

    def plot(self):
        full_df = self.data

        fig, axes = plt.subplots(2, figsize=(14, 6))
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

            breakpoints = np.where(
                df.population.values[1:] != df.population.values[:-1])[0]
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
        plt.plot(df["Agree"]/df["Total"],label="No Error")

        plt.plot(df["ErrorAgree"]/df["Total"],label="Error")
        # plt.xlabel("Distance Separating Alleles (bp)")

        plt.xlabel("Kb separating Alleles")
        plt.ylabel("Proportion of Mutation Pairs Correctly Ordered")
        plt.legend()
        
        self.save()


class AncestorAccuracy(Figure):
    """
    Compare lengths of real vs reconstructed ancestors, using 2 csv files generated by
    TSINFER_DIR=../tsinfer #set to your tsinfer directory
    python3 ${TSINFER_DIR}/evaluation.py aq -l 2 -d data -C -s 123 -e 0
    python3 ${TSINFER_DIR}/evaluation.py aq -l 2 -d data -C -s 123 -e data/EmpiricalErrorPlatinum1000G.csv
    cd data
    cat anc-qual_n=100_L=2.0_mu=1e-08_rho=1e-08_err=data_EmpiricalErrorPlatinum1000G_error_data.csv > ancestor_accuracy.csv
    tail +2 anc-qual_n=100_L=2.0_mu=1e-08_rho=1e-08_err=0.0_error_data.csv >> ancestor_accuracy.csv

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
        max_length = max(np.max(self.data["Real length"]), np.max(self.data["Estim length"]))* 1.1
        min_length = min(np.min(self.data["Real length"]), np.min(self.data["Estim length"])) * 0.9
        fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12, 5.5))
        for ax, error in zip(axes, self.data.seq_error.unique()):
            df = self.data.query("seq_error == @error")
            im = ax.scatter(df["Real length"], df["Estim length"], c=df["Inaccuracy"], s=20)
            ax.plot([0, max_length], [0, max_length], '-', color='lightgrey', zorder=-1)
            #print(np.mean(df["Inaccuracy"]), error)
            ax.set_title(self.error_label(error))
            ax.set_xlabel("True ancestral haplotype length (kb)")
            if ax == axes[0]:
                ax.set_ylabel("Inferred ancestral haplotype length (kb)")
            ax.set_xscale('log')
            ax.set_yscale('log')
        cbar = fig.colorbar(im, ax=axes.ravel().tolist())
        cbar.set_label("Inaccuracy", labelpad=6, rotation=270, va="center")
        plt.xlim(min_length, max_length)
        plt.ylim(min_length, max_length)
        self.save()


class ToolsFigure(Figure):
    """
    Superclass of all figures where different tools (e.g. ARGweaver, fastarg) are compared
    """

    # Colours taken from Matplotlib default color wheel.
    # https://matplotlib.org/users/dflt_style_changes.html
    tools_format = collections.OrderedDict([
        ("ARGweaver", {"mark":"o", "col":"#d62728"}),
        ("RentPlus",  {"mark":"^", "col":"#2ca02c"}),
        ("fastARG",   {"mark":"s", "col":"#ff7f0e"}),
        ("tsinfer",   {"mark":"*", "col":"#1f77b4"}),
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
        length_sample_size_combos = self.data[["length", "sample_size"]].drop_duplicates()
        self.fixed_length = length_sample_size_combos['length'].value_counts().idxmax()
        self.fixed_sample_size = length_sample_size_combos['sample_size'].value_counts().idxmax()
        # Scale time to hours
        time_scale = 3600
        cpu_cols = [c for c in self.data.columns if c.startswith("cputime")]
        self.data[cpu_cols] /= time_scale
        # Scale memory to GiB
        mem_cols = [c for c in self.data.columns if c.startswith("memory")]
        self.data[mem_cols] /= 1024 * 1024 * 1024

    def plot(self):
        fig, axes = plt.subplots(2, 2, sharey="row", sharex="col", figsize=(8, 5.5))
        for i, (plotted_column, y_label) in enumerate(
                zip(["cputime", "memory"], ["CPU time (hours)", "Memory (GiB)"])):
            df = self.data.query("sample_size == @self.fixed_sample_size")
            for tool in df.tool.unique():
                line_data = df.query("tool == @tool")
                axes[i][0].errorbar(
                    line_data.length,
                    line_data[plotted_column+"_mean"],
                    yerr=line_data[plotted_column+"_se"],
                    color=self.tools_format[tool]["col"],
                    marker=self.tools_format[tool]['mark'],
                    elinewidth=1,
                    label=tool)
            axes[i][0].get_yaxis().set_label_coords(-0.08,0.5)
            axes[i][0].set_ylabel(y_label)

            df = self.data.query("length == @self.fixed_length")
            for tool in df.tool.unique():
                line_data = df.query("tool == @tool")
                axes[i][1].errorbar(
                    line_data.sample_size,
                    line_data[plotted_column+"_mean"],
                    yerr=line_data[plotted_column+"_se"],
                    color=self.tools_format[tool]["col"],
                    marker=self.tools_format[tool]['mark'],
                    elinewidth=1,
                    label=tool)

        axes[0][0].legend(
            loc="upper right", numpoints=1, fontsize="small")
        axes[1][0].set_xlabel("Length (Mb) for fixed sample size of {}".format(
            self.fixed_sample_size))
        axes[1][1].set_xlabel("Sample size for fixed length of {:g} Mb".format(
            self.fixed_length))
        fig.tight_layout()

        self.save()

class PerformanceLengthSamplesFigure(ToolsFigure):
    """
    Superclass for the performance metrics figures. Each of these figures
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
        recombination_rates = df.recombination_rate.unique()
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
            params, [r"$\rho$ = {}".format("$\mu$" if rho==mu else r"{:g}$\mu$".format(rho/mu) if rho>mu else r"$\mu$/{:g}".format(mu/rho))
                for rho_index, rho in enumerate(recombination_rates)],
            loc="upper right", fontsize=10, title="Relative rate of\nrecombination")

        self.save()


class TSCompressionFigure(PerformanceLengthSamplesFigure):
    name = "tsinfer_ts_filesize_ln"
    plotted_column = "ts_relative_filesize"
    y_axis_label = "File size relative to simulated tree sequence"


class VCFCompressionFigure(PerformanceLengthSamplesFigure):
    name = "tsinfer_vcf_compression_ln"
    plotted_column = "vcf_compression_factor"
    y_axis_label = "Compression factor relative to vcf.gz"


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
        error_params = self.data.error_param.unique()
        # y-direction is the permutations of metric + whether it is rooted
        metric_and_rooting = self.data.groupby(["metric", "rooting"]).groups
        # sort this so that metrics come out in a set order (TO DO)
        fig, axes = plt.subplots(len(metric_and_rooting),
            len(error_params), figsize=(6*len(error_params), 15), sharey='row')
        for j, ((metric, root), rows) in enumerate(metric_and_rooting.items()):
            for k, error in enumerate(error_params):
                # we are in the j,k th subplot
                ax = axes[j][k] if len(error_params)>1 else axes[j]
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
    name = "metric_all_tools"
    figsize = (9,4.5)
    y_axis_label="Average distance from true trees"
    hide_polytomy_breaking = True
    output_metrics = [("KC","rooted")] #can add extras in here if necessary

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

            # x-direction is different error rates
            error_params = df.error_param.unique()

            fig, axes = plt.subplots(1, len(error_params),
                figsize=getattr(self,'figsize',(12, 6)), sharey=True)
            for k, error in enumerate(error_params):
                ax = axes[k] if len(error_params)>1 else axes
                display_order = self.single_metric_plot(
                    df.query("(" + ") and (".join(query + ["error_param == @error"]) + ")"),
                    "mutation_rate", ax, method, rho)
                ax.set_title(self.error_label(error))
                ax.set_xlabel("Mutation rate")
                ax.set_xscale('log')
                if k == 0:
                    ax.set_ylim(getattr(self,'ylim', 0))
                    rooting_suffix = " (unrooted)" if rooting=="unrooted" else ""
                    ylab = getattr(self, 'y_axis_label', self.metric_titles[metric] + rooting_suffix)
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
            first_legend = axes[0].legend(
                artists, tool_labels, numpoints=1, labelspacing=0.1, loc="upper right")
            if len(self.output_metrics)==1:
                self.save()
            else:
                self.save("_".join([self.name, metric, rooting]))


class MetricAllToolsAccuracyDemographyFigure(MetricAllToolsFigure):
    """
    Simple figure that shows an ARG metrics for a genome under a more complex demographic
    model (the Gutenkunst Out Of Africa model), as mutation rate increases to high values
    """
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
            error_params = df.error_param.unique()
            fig, axes = plt.subplots(len(output_freqs), len(error_params),
                figsize=getattr(self,'figsize',(6*len(error_params), 2.5*len(output_freqs))),
                sharey=True)
            for j, output_data in enumerate(output_freqs.itertuples()):
                for k, error in enumerate(error_params):
                    ax = axes[j][k] if len(error_params)>1 else axes[j]
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
            tree_tips = df.restrict_sample_size_comparison.unique()
            mutation_rates = df.mutation_rate.unique()
            assert len(tree_tips) == len(mutation_rates) == len(sample_sizes) == len(averaging_method) == 1
            method = averaging_method[0]
            lengths = df.length.unique()
            error_params = df.error_param.unique()
            fig, axes = plt.subplots(1, len(error_params),
                figsize=(12, 6), sharey=True)
            for k, error in enumerate(error_params):
                ax = axes[k]
                display_order = self.single_metric_plot(
                    df.query("error_param == @error"), "subsample_size",
                    ax, method, rho = None, markers = False, x_jitter = 'log')
                ax.set_title(self.error_label(error))
                if k == 0:
                    ylab = getattr(self, 'y_axis_label', self.metric_titles[metric])
                    ax.set_ylabel(ylab)
                ax.set_xlabel("Original sample size")
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
                first_legend = axes[0].legend(
                    artists, labels, numpoints=1, labelspacing=0.1, loc="upper right")
            if len(self.output_metrics)==1:
                self.save()
            else:
                self.save("_".join([self.name, metric, rooting]))


class UkbbStructureFigure(Figure):
    """
    Figure showing the structure for UKBB using heatmaps.
    """
    name = "ukbb_structure"

    def plot(self):
        # Current document width is 4.7747 inches
        dfs = [
            pd.read_csv("data/ukbb_1kg_ethnicity.csv").set_index("Ethnicity"),
            pd.read_csv("data/ukbb_1kg_british_centre.csv").set_index("CentreName"),
            self.data.set_index("CentreName")
        ]
        # print(self.data)

        fig, axes = plt.subplots(1, 3, figsize=(18, 8))
        axes[0].set_title("(A)")
        axes[1].set_title("(B)")
        axes[2].set_title("(C)")
        cbar_ax = fig.add_axes([.92, .3, .03, .4])
        # fig.tight_layout(rect=[0, 0, .9, 1])

        vmax = max([df.values.max() for df in dfs])
        vmin = min([df.values.min() for df in dfs])

        # FIXME this should use hierarchical clustering or something to
        # show the correct order.
        df = dfs[0]
        row_order = [
            'Chinese',
            'British',
            'Irish',
            'White',
            'Any other white background',
            'Prefer not to answer',
            'Do not know',
            'African',
            'Caribbean',
            'Any other Black background',
            'Black or Black British',
            'White and Black African',
            'White and Black Caribbean',
            'Mixed',
            'Any other mixed background',
            'Other ethnic group',
            'White and Asian',
            'Pakistani',
            'Indian',
            'Bangladeshi',
            'Asian or Asian British',
            'Any other Asian background',
        ]
        row_labels = df.index.unique()
        index = np.searchsorted(row_labels, row_order)
        assert len(row_labels) == len(row_order)
        V = df.values
        sns.heatmap(
            V[index], xticklabels=list(df), yticklabels=row_order,
            ax=axes[0], vmax=vmax, vmin=vmin, cbar=True,
            cbar_ax=cbar_ax)

        df = dfs[1]
        V = df.values
        sns.heatmap(df, ax=axes[1], vmax=vmax, vmin=vmin, cbar=False)

        df = dfs[2]
        sns.heatmap(df[df.index], ax=axes[2], vmax=vmax, vmin=vmin, cbar=False)

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
