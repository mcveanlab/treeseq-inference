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
        plt.clf()


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

    def plot(self):
        raise NotImplementedError()


class ToolsFigure(Figure):
    """
    Superclass of all figures where different tools (e.g. ARGweaver, fastarg) are compared
    """

    tools_format = collections.OrderedDict([
        ("ARGweaver", {"mark":"o", "col":"green"}),
        ("RentPlus",  {"mark":"^", "col":"magenta"}),
        ("fastARG",   {"mark":"s", "col":"red"}),
        ("tsinfer",   {"mark":"*", "col":"blue"}),
    ])
    
    error_bars = True

    def error_label(self, error, label_for_no_error = "No error"):
        """
        Make a nice label for an error parameter
        """
        try:
            error = float(error)
            return "Error rate = {}".format(error) if error else label_for_no_error
        except (ValueError, TypeError):
            try: # make a simplified label
                if error.startswith("Empirical"):
                    error = "Empirical"
            except:
                pass            
            return "{} error".format(error) if error else label_for_no_error

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
        ax_hi.set_ylim(bottom = max_non_AW*20)  # outliers only
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

class FastargTsinferComparisonFigure(ToolsFigure):
    """
    Superclass for the program comparison figures (comparing tsinfer with fastarg)
    Each figure has two panels; one for scaling by sequence length and the other
    for scaling by sample size.
    """

    def __init__(self):
        super().__init__()
        # Rescale the length to Mb
        length_scale = 10**6
        self.data.length /= length_scale
        length_sample_size_combos = self.data[["length", "sample_size"]].drop_duplicates()
        self.fixed_length = length_sample_size_combos['length'].value_counts().idxmax()
        self.fixed_sample_size = length_sample_size_combos['sample_size'].value_counts().idxmax()
        

    def plot(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(8, 5.5))
        df = self.data.query("sample_size == @self.fixed_sample_size")
        for tool in df.tool.unique():
            line_data = df.query("tool == @tool")
            ax1.errorbar(
                line_data.length,
                line_data[self.plotted_column+"_mean"],
                yerr=line_data[self.plotted_column+"_se"],
                color=self.tools_format[tool]["col"],
                marker=self.tools_format[tool]['mark'],
                elinewidth=1,
                label=tool)
        ax1.legend(
            loc="upper left", numpoints=1, fontsize="small")

        ax1.set_xlabel("Length (Mb) for fixed sample of {}".format(self.fixed_sample_size))
        ax1.set_ylabel(self.y_label)

        df = self.data.query("length == @self.fixed_length")
        for tool in df.tool.unique():
            line_data = df.query("tool == @tool")
            ax2.errorbar(
                line_data.sample_size,
                line_data[self.plotted_column+"_mean"],
                yerr=line_data[self.plotted_column+"_se"],
                color=self.tools_format[tool]["col"],
                marker=self.tools_format[tool]['mark'],
                elinewidth=1,
                label=tool)

        ax2.set_xlabel("Sample size for fixed length of {} Mb".format(
            self.fixed_length))

        fig.tight_layout()

        self.save()


class FastargTsinferComparisonTimeFigure(FastargTsinferComparisonFigure):
    name = "fastarg_tsinfer_comparison_time"
    plotted_column = "cputime"
    y_label = "CPU time (hours)"
    def __init__(self):
        super().__init__()
        # Scale time to hours
        time_scale = 3600
        cpu_cols = [c for c in self.data.columns if c.startswith("cputime")]
        self.data[cpu_cols] /= time_scale



class FastargTsinferComparisonMemoryFigure(FastargTsinferComparisonFigure):
    name = "fastarg_tsinfer_comparison_memory"
    plotted_column = "memory"
    y_label = "Memory (GiB)"

    def __init__(self):
        super().__init__()
        # Scale memory to GiB
        mem_cols = [c for c in self.data.columns if c.startswith("memory")]
        self.data[mem_cols] /= 1024 * 1024 * 1024

class TsinferPerformanceLengthSamplesFigure(ToolsFigure):
    """
    Superclass for the performance metrics figures. Each of these figures
    has two panels; one for scaling by sequence length and the other
    for scaling by sample size. Different lines are given for each
    of the different combinations of tsinfer parameters
    """
    plotted_column = None
    y_axis_label = None

    def plot(self):
        df = self.dataset.data
        # Rescale the length to MB
        df.length /= 10**6
        # Set statistics to the ratio of observed over expected
        source_colour = "red"
        inferred_colour = self.tools[TSINFER]["col"]
        recombination_linestyles = [':', '-', '--']
        recombination_rates = df.recombination_rate.unique()
        mutation_rates = self.datasetClass.between_sim_params['mutation_rate']
        assert len(recombination_linestyles) >= len(recombination_rates)
        assert len(mutation_rates) ==1
        mu = mutation_rates[0]
        #inferred_markers =    {False:'s',True:'b'}
        fig, (ax1, ax2) = pyplot.subplots(1, 2, figsize=(12, 6), sharey=True)
        ax1.set_title("Fixed number of chromosomes ({})".format(self.datasetClass.fixed_sample_size))
        ax1.set_xlabel("Sequence length (MB)")
        ax1.set_ylabel(self.y_axis_label)
        for rho_index, rho in enumerate(recombination_rates):
            dfp = df[np.logical_and.reduce((
                df.sample_size == self.datasetClass.fixed_sample_size,
                df.recombination_rate == rho))]
            group = dfp.groupby(["length"])
            #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
            mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
            ax1.errorbar(
                [m['mu'] for m in mean_sem],
                [m['mean'][self.plotted_column] for m in mean_sem],
                yerr=[m['sem'][self.plotted_column] for m in mean_sem] if getattr(self, 'error_bars', False) else None,
                linestyle= recombination_linestyles[rho_index],
                color=inferred_colour,
                #elinewidth=1
                )


        ax2.set_title("Fixed sequence length ({:g} Mb)".format(self.datasetClass.fixed_length / 10**6))
        ax2.set_xlabel("Sample size")
        ax2.set_ylabel(self.y_axis_label)
        for rho_index, rho in enumerate(recombination_rates):
            dfp = df[np.logical_and.reduce((
                df.length == self.datasetClass.fixed_length / 10**6,
                df.recombination_rate == rho))]
            group = dfp.groupby(["sample_size"])
            #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
            mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
            ax2.errorbar(
                [m['mu'] for m in mean_sem],
                [m['mean'][self.plotted_column] for m in mean_sem],
                yerr=[m['sem'][self.plotted_column] for m in mean_sem] if getattr(self, 'error_bars', False) else None,
                linestyle= recombination_linestyles[rho_index],
                color=inferred_colour,
                #elinewidth=1
                )
        params = [
            pyplot.Line2D(
                (0,0),(0,0), color= inferred_colour,
                linestyle=recombination_linestyles[rho_index], linewidth=2)
            for rho_index, rho in enumerate(recombination_rates)]
        ax1.legend(
            params, [r"$\rho$ = {}".format("$\mu$" if rho==mu else r"{:g}$\mu$".format(rho/mu) if rho>mu else r"$\mu$/{:g}".format(mu/rho))
                for rho_index, rho in enumerate(recombination_rates)],
            loc="lower right", fontsize=10, title="Relative rate of\nrecombination")
    
        pyplot.suptitle(r'Tsinfer large dataset performance for $\mu$=${}$'.format(latex_float(mu)))
        self.savefig(fig)


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

    def single_metric_plot(self, df, x_variable, ax, av_method, 
        rho = None, markers = True):
        """
        A single plot on an ax. This requires plotting separate lines, e.g. for each tool
        If rho is give, plot an x=rho vertical line, assuming x is the mutation_rate
        """
        v_cols = ['sample_size', 'tool', 'polytomies']
        v_order = df[v_cols].drop_duplicates() #find unique combinations
        v_order = v_order.sort_values(v_cols, ascending=[True, True, False]) #sort for display
        ss_order = {v:k for k,v in enumerate(df.sample_size.unique())}
        for r in v_order.itertuples():
            query = []
            query.append("sample_size == @r.sample_size")
            query.append("tool == @r.tool")
            query.append("polytomies == @r.polytomies")
            line_data = df.query("(" + ") and (".join(query) + ")")
            if not line_data.empty:
                ax.errorbar(
                    line_data[x_variable],
                    line_data.treedist_mean,
                    yerr=line_data.treedist_se if self.error_bars else None,
                    linestyle=self.polytomy_and_averaging_format[r.polytomies][av_method]["linestyle"],
                    fillstyle=self.sample_size_format[ss_order[r.sample_size]]['fillstyle'],
                    color=self.tools_format[r.tool]["col"],
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
            len(error_params), figsize=(4*len(error_params), 15), sharey='row')
        for j, ((metric, root), rows) in enumerate(metric_and_rooting.items()):
            for k, error in enumerate(error_params):
                # we are in the j,k th subplot
                ax = axes[j][k] if len(error_params)>1 else axes[j]
                ax.set_xscale('log')
                display_order = self.single_metric_plot(
                    self.data.loc[rows].query("error_param == @error"), "mutation_rate",
                    ax, method, rho, markers = (len(sample_sizes)!=1))
                # Make the labels on each Y axis line up properly
                ax.get_yaxis().set_label_coords(-0.15,0.5)
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
            
        maintitle = "Average tree distance for neutral simulations"
        #if len(sample_sizes)>1:
        #    artists = [
        #        tuple([
        #            plt.Line2D((0,0),(0,0),
        #                color=self.tools_format[tool]['col'],
        #                marker=self.tools_format[tool]['mark'],
        #                fillstyle=fillstyle, 
        #                linestyle='None')
        #            for tool, poly in self.data.groupby(["tool", "polytomies"]).groups.keys()])
        #        for fillstyle in self.sample_size_format]
        #    axes[0][-1].legend(
        #        artists, ["Sample size = {}".format(n) for n in sample_sizes],
        #        loc="upper right", numpoints=1, 
        #        handler_map={tuple: matplotlib.legend_handler.HandlerTuple(ndivide=None, pad=2)})
        #else:
        maintitle += " of {} samples\n".format(sample_sizes[0])
        if any([isinstance(Ne, str) for Ne in eff_sizes]):
            maintitle += "(Model{}=".format('s' if len(eff_sizes) > 1 else "") \
                + ",".join(["“{}”".format(x.replace("."," ")) for x in eff_sizes])
        else:
            maintitle += "($N_e$=" + ",".join(["{}".format(x) for x in eff_sizes])
        maintitle += r"; $\rho="+ ",".join(["{}".format(latex_float(x)) for x in rhos])
        maintitle += "$; sequence length = " + " ,".join(["{:g}kb".format(x/1e3) for x in lengths])
        maintitle += ")"
        fig.suptitle(maintitle, fontsize=16)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        self.save()


class MetricsAllToolsAccuracyFigure(MetricsAllToolsFigure):
    """
    Show the metrics tending to 0 as mutation rate increases
    """
    name = "metrics_all_tools_accuracy"
    
class MetricsAllToolsAccuracyDemographyFigure(MetricsAllToolsFigure):
    """
    Simple figure that shows all the metrics at the same time for
    a genome under a more complex demographic model (the Gutenkunst 
    Out Of Africa model), as mutation rate increases to high values
    """
    name = "metrics_all_tools_accuracy_demography"


class MetricAllToolsFigure(TreeMetricsFigure):
    """
    Plot each metric in a different pdf file.
    For the publication: make symbols small and solid
    """
    name = "metric_all_tools"
    figsize = (9,4.5)
    y_axis_label="Average distance from true trees (mean $\pm$ s.e.)"
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

class MetricAllToolsAccuracyFigure(MetricAllToolsFigure):
    """
    Superclass of the metric all tools figure. 
    """
    name = "metric_all_tools_accuracy"    
    y_axis_label="Average distance from true trees (mean $\pm$ s.e.)"
    output_metrics = [("KC","rooted")] #can add extras in here if necessary


class MetricAllToolsAccuracySweepFigure(TreeMetricsFigure):
    """
    Figure for simulations with selection. 
    Each page should be a single figure for a particular metric, with error on the 
    """
    name = "metric_all_tools_accuracy_sweep"
    error_bars = True
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
                figsize=getattr(self,'figsize',(7*len(error_params), 3*len(output_freqs))),
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
            fig, axes = plt.subplots(len(lengths), len(error_params),
                figsize=(12, 12), sharey=True)
            for j, l in enumerate(lengths):
                for k, error in enumerate(error_params):
                    ax = axes[j][k]
                    display_order = self.single_metric_plot(
                        df.query("(error_param == @error) and (length == @l)"),
                        "subsample_size", ax, method, 
                        rho = None, markers = False)
                    if j == 0:
                        ax.set_title(self.error_label(error))
                    if k == 0:
                        ylab = getattr(self, 'y_axis_label', "Average " + \
                            self.metric_titles[metric] + " over {}kb".format(l//1000))
                        ax.set_ylabel(ylab)
                    if j == len(lengths) - 1:
                        ax.set_xlabel("Original sample size")
            if len(display_order)>1:
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
            fig.suptitle('Accuracy with increased sampling (trees subsampled down to {} tips)'.format(
                ", ".join(str(t) for t in tree_tips)))
            if len(self.output_metrics)==1:
                self.save()
            else:
                self.save("_".join([self.name, metric, rooting]))                

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
        "name", type=str, help="figure name", choices=list(name_map.keys()))
    
    args = parser.parse_args()
    fig = name_map[args.name]()
    fig.plot()


if __name__ == "__main__":
    main()
