# Generates all the actual figures. Run like
# python3 src/plot.py PLOT_NAME

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

class Figure(object):
    """
    Superclass of figures for the paper. Each figure is a concrete subclass.
    """

    


class StoringEveryone(Figure):
    """
    Figure showing how tree sequences can store the entire human population
    worth of variation data.
    """
    name = "storing_everyone"

    def plot(self):
        GB = 1024**3

        # Fit the model for the observed data. Wek
        rho = 4 * Ne * recombination_rate * length

        def additive_model(n, a, b, c):
            # The model is based on the proof that we need O(n + rho log n) space
            # to store edge. We then need a further O(theta log n) space to
            # store the mutations. Because rho == mu, we don't need a separate term
            # for rho, and represent them both by z. Thus, a here is the per-edge
            # cost, b is the per site cost, and c is a constant for the overheads.
            z = rho * np.log(n)
            return a * (n + z) + b * z + c

        def mulplicative_model(n, a, b, c):
            # This is the number of sites.
            z = rho * np.log(n)
            return a * n * z + b * z + c

        df = pd.read_csv(datafile)
        df = df[df.sample_size > 100]
        projected_n = 10**np.arange(3, 11)

        msp_fit_params, _ = optimize.curve_fit(
            additive_model, df.sample_size[:-1], df.uncompressed[:-1])
        msp_fit = additive_model(projected_n, *msp_fit_params)
        mspz_fit_params, _ = optimize.curve_fit(
            additive_model, df.sample_size[:-1], df.compressed[:-1])
        mspz_fit = additive_model(projected_n, *mspz_fit_params)

        index = df.sample_size < 10**6
        vcf_fit_params, _ = optimize.curve_fit(
            mulplicative_model, df.sample_size[index], df.vcf[index])
        vcf_fit = mulplicative_model(projected_n, *vcf_fit_params)

        bcf_fit_params, _ = optimize.curve_fit(
            mulplicative_model, df.sample_size[index], df.bcf[index])
        bcf_fit = mulplicative_model(projected_n, *bcf_fit_params)
        bcf_fit_params, _ = optimize.curve_fit(
            mulplicative_model, df.sample_size[index], df.bcf[index])

        fig = pyplot.figure()
        ax1 = fig.add_subplot(111)
        xytext = (18, 0)

        index = df.vcf > 0
        line, = ax1.loglog(df.sample_size[index], df.vcf[index], "^", label="VCF")
        ax1.loglog(projected_n, vcf_fit, "--", color=line.get_color())
        ax1.annotate(
            humanize.naturalsize(vcf_fit[-1] * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(projected_n[-1], vcf_fit[-1]), xycoords="data")

        line, = ax1.loglog(df.sample_size[index], df.bcf[index], "s", label="BCF")
        ax1.loglog(projected_n, bcf_fit, "--", color=line.get_color())
        ax1.annotate(
            humanize.naturalsize(bcf_fit[-1] * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(projected_n[-1], bcf_fit[-1]), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.uncompressed, "o", label=".trees")
        ax1.loglog(projected_n, msp_fit, "--", color=line.get_color())
        ax1.annotate(
            humanize.naturalsize(msp_fit[-1] * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(projected_n[-1], msp_fit[-1]), xycoords="data")

        line, = ax1.loglog(
            df.sample_size, df.compressed, "*", label=".trees.gz")
        ax1.loglog(projected_n, mspz_fit, "--", color=line.get_color())
        ax1.annotate(
            humanize.naturalsize(mspz_fit[-1] * GB, binary=True, format="%d"),
            textcoords="offset points", xytext=xytext,
            xy=(projected_n[-1], mspz_fit[-1]), xycoords="data")

        ax1.set_xlabel("Number of chromosomes")
        ax1.set_ylabel("File size (GiB)")
        pyplot.legend()
        # pyplot.tight_layout()
        pyplot.savefig("figures/storing_everyone.pdf")



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
