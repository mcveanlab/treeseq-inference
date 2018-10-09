# Generates all the actual figures. Run like
# python3 src/plot.py PLOT_NAME

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


class MetricFigure(ToolsFigure):

    polytomy_and_averaging_format = collections.OrderedDict([
        ("breaking polytomies", {
            "per site":    {"linestyle":"--"},
            "per variant": {"linestyle":":"}}),
        ("leaving polytomies", {
            "per site":    {"linestyle":"-"},
            "per variant": {"linestyle":"-."}})
    ])

    fillstyles = ['full', 'none']
            
    def single_metric_plot(self, df, ax, sample_sizes, rho, averaging, markers = True):
        #now plot each line on this subplot
        for n, fillstyle in zip(sample_sizes, self.fillstyles):
            for tool in df.tool.unique():
                for poly in df.polytomy_treatment.unique():
                    query = []
                    query.append("sample_size == @n")
                    query.append("tool == @tool")
                    query.append("polytomy_treatment == @poly")
                    line_data = df.query("(" + ") and (".join(query) + ")")
                    if not line_data.empty:
                        ax.errorbar(
                            line_data.mutation_rate,
                            line_data.treedist_mean,
                            line_data.treedist_se,
                            linestyle=self.polytomy_and_averaging_format[poly][averaging]["linestyle"],
                            fillstyle=fillstyle,
                            color=self.tools_format[tool]["col"],
                            marker=self.tools_format[tool]['mark'] if markers else None,
                            elinewidth=1)
        ax.set_ylim(ymin=0)
        ax.axvline(x=rho, color = 'gray', zorder=-1, linestyle=":", linewidth=1)
        ax.text(rho, ax.get_ylim()[1]/40,  r'$\mu=\rho$',
            va="bottom",  ha="right", color='gray', rotation=90)
        ax.get_yaxis().set_label_coords(-0.15,0.5)


class MetricsAllToolsFigure(MetricFigure):
    """
    Simple figure that shows all the metrics at the same time.
    Assumes at most 2 sample sizes
    """
    name = "metrics_all_tools"

    def plot(self):
        averaging = self.data.averaging.unique()
        eff_sizes = self.data.Ne.unique()
        rhos = self.data.recombination_rate.unique()
        lengths = self.data.length.unique()
        assert len(averaging) == len(eff_sizes) == len(rhos) == 1 #these should not vary
        
        sample_sizes = self.data.sample_size.unique()
        # y-direction is different error rates
        error_params = self.data.error_param.unique()
        # x-direction is the permutations of metric + whether it is rooted
        metric_and_rooting = self.data.groupby(["metric", "rooting"]).groups
        # sort this (TO DO)
        fig, axes = plt.subplots(len(metric_and_rooting),
            len(error_params), figsize=(4*len(error_params), 15), sharey='row')
        for j, ((metric, root), rows) in enumerate(metric_and_rooting.items()):
            for k, error in enumerate(error_params):
                # we are in the j,k th subplot
                ax = axes[j][k] if len(error_params)>1 else axes[j]
                ax.set_xscale('log')
                if j == 0:
                    ax.set_title(self.error_label(error))
                if j == len(metric_and_rooting) - 1:
                    ax.set_xlabel("Mutation rate")
                if k == 0:
                    rooting_suffix = " (unrooted)" if root=="unrooted" else ""
                    if getattr(self, 'y_axis_label', None) is None:
                        ax.set_ylabel(metric + rooting_suffix)
                self.single_metric_plot(
                    self.data.loc[rows].query("error_param == @error"),
                    ax, sample_sizes, rhos[0], averaging[0],
                    markers = (len(sample_sizes)!=1))

        artists = [
            plt.Line2D((0,1),(0,0), linewidth=2,
                color=self.tools_format[tool]["col"],
                linestyle=self.polytomy_and_averaging_format[poly][averaging[0]]["linestyle"],
                marker = None if len(sample_sizes)==1 else self.tools_format[tool]['mark'])
            for tool, poly in self.data.groupby(["tool", "polytomy_treatment"]).groups.keys()]
        tool_labels = [tool + ("" if poly == "leaving polytomies" else (" " + poly)) 
            for tool, poly in self.data.groupby(["tool", "polytomy_treatment"]).groups.keys()]
        axes[0][0].legend(
            artists, tool_labels, numpoints=1, labelspacing=0.1)
            
        maintitle = "Average tree distance for neutral simulations"
        if len(sample_sizes)>1:
            artists = [
                tuple([
                    plt.Line2D((0,0),(0,0),
                        color=self.tools_format[tool]['col'],
                        marker=self.tools_format[tool]['mark'],
                        fillstyle=fillstyle, 
                        linestyle='None')
                    for tool, poly in self.data.groupby(["tool", "polytomy_treatment"]).groups.keys()])
                for fillstyle in self.fillstyles]
            axes[0][-1].legend(
                artists, ["Sample size = {}".format(n) for n in sample_sizes],
                loc="upper right", numpoints=1, 
                handler_map={tuple: matplotlib.legend_handler.HandlerTuple(ndivide=None, pad=2)})
        else:
            maintitle += " of {} samples".format(sample_sizes[0])
        maintitle += " spanning "  + ",".join(["{:g}kb\n".format(x/1e3) for x in lengths])
        if any([isinstance(Ne, str) for Ne in eff_sizes]):
            maintitle += "(Model{}=".format('s' if len(eff_sizes) > 1 else "") \
                + ",".join(["“{}”".format(x.replace("."," ")) for x in eff_sizes])
        else:
            maintitle += "($N_e$=" + ",".join(["{}".format(x) for x in eff_sizes])
        maintitle += r"; $\rho="+ ",".join(["{}".format(latex_float(x)) for x in rhos])
        maintitle += "$)"
        fig.suptitle(maintitle, fontsize=16)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        self.save()


class MetricsAllToolsAccuracyFigure(MetricsAllToolsFigure):
    """
    Show the metrics tending to 0 as mutation rate increases
    """
    name = "metrics_all_tools_accuracy"
    error_bars = True
    
class MetricsAllToolsAccuracyDemographyFigure(MetricsAllToolsFigure):
    """
    Simple figure that shows all the metrics at the same time for
    a genome under a more complex demographic model (the Gutenkunst 
    Out Of Africa model), as mutation rate increases to high values
    """
    name = "metrics_all_tools_accuracy_demography"


class MetricAllToolsFigure(MetricFigure):
    """
    Superclass of the metric all tools figure. Each subclass should be a
    single figure for a particular metric.
    """
    def plot(self):
        metrics = self.data.metric.unique()
        # y-direction is different error rates
        error_params = self.data.error_param.unique()

        averaging = self.data.averaging.unique()
        eff_sizes = self.data.Ne.unique()
        rhos = self.data.recombination_rate.unique()
        lengths = self.data.length.unique()
        assert len(averaging) == len(eff_sizes) == len(rhos) == len(metrics) == 1
        
        metric = metrics[0]
        sample_sizes = self.data.sample_size.unique()

        fig, axes = plt.subplots(1, len(error_params),
            figsize=getattr(self,'figsize',(12, 6)), sharey=True)
        lines = []
        for k, error in enumerate(error_params):
            ax = axes[k] if len(error_params)>1 else axes
            ax.set_title(self.error_label(error))
            ax.set_xlabel("Mutation rate")
            ax.set_xscale('log')
            if k == 0:
                ax.set_ylim(self.ylim)
                ax.set_ylabel(getattr(self, 'y_axis_label', metric))
            self.single_metric_plot(
                self.data.query("error_param == @error"),
                ax, sample_sizes, rhos[0], averaging[0])

        # Create legends from custom artists
        #artists = [
        #    pyplot.Line2D((0,1),(0,0), color= setting["col"], fillstyle=self.fillstyles[0],
        #        marker= setting["mark"], linestyle=setting["linestyle"])
        #    for tool,setting in self.tools_and_metrics_params.items()]
        #tool_labels = [(l.replace("_0", "").replace("_2"," breaking polytomies") if l.startswith(TSINFER) else l.replace("_0", "")) 
        #    for l in self.tools_and_metrics_params.keys()]
        #first_legend = axes[0].legend(
        #    artists, tool_labels, numpoints=1, labelspacing=0.1, loc="upper right")
        #    # bbox_to_anchor=(0.0, 0.1))
        ## ax = pyplot.gca().add_artist(first_legend)
        #if len(sample_sizes)>1:
        #    artists = [
        #        pyplot.Line2D(
        #            (0,0),(0,0), color="black", fillstyle=fillstyle, linewidth=2)
        #        for n, fillstyle in zip(sample_sizes, self.fillstyles)]
        #    axes[-1].legend(
        #        artists, ["Sample size = {}".format(n) for n, fillstyle in zip(sample_sizes, self.fillstyles)],
        #        loc="upper right")
        self.save()

class MetricAllToolsAccuracyFigure(MetricAllToolsFigure):
    """
    Superclass of the metric all tools figure. Each subclass should be a
    single figure for a particular metric.
    """
    pass


class KCAllToolsFigure(MetricAllToolsFigure):
    """For the publication: make symbols small and solid"""
    name = "kc_all_tools"
    ylim = (0, 24)
    fillstyles = ['full']
    figsize = (9,4.5)
    error_bars = True
    y_axis_label="Average distance from true trees (mean $\pm$ s.e.)"


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
