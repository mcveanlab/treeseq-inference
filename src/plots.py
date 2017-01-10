"""
Code to run simulations, inference methods and generate all plots
in the paper.
"""

import argparse
import sys
import os.path

if sys.version_info[0] < 3:
    raise Exception("Python 3 only")

class Figure(object):
    """
    Superclass of all figures. Each figure depends on a dataset.
    """
    name = None
    """
    Each figure has a unique name. This is used as the identifier and the
    file name for the output plots.
    """

class Dataset(object):
    """
    A dataset is some collection of simulations which are run using
    the generate() method and stored in the raw_data_path directory.
    The process() method then processes the raw data and outputs
    the results into the data file.
    """
    name = None
    """
    Each dataset a unique name. This is used as the prefix for the data
    file and raw_data_dir directory.
    """

    data_dir = "data"

    def __init__(self):
        self.data_file = os.path.join(data_dir, "{}.csv".format(self.name))
        self.raw_data_dir = os.path.join(data_dir, "raw__NOBACKUP__", self.name)
        if not os.path.exists(self.raw_data_dir):
            print("making", self.raw_data_dir)

    def generate(self):
        pass

    def process(self):
        pass


def run_generate(cls, args):
    f = cls()
    print("generating")
    # f.generate()


def run_process(cls, args):
    f = cls()
    f.process()


def run_plot(cls, args):
    f = cls()
    f.plot()


def main():
    datasets = [
        Dataset,
    ]
    figures = [
        Figure,
    ]
    name_map = dict([(d.name, d) for d in datasets + figures])
    print(name_map)
    parser = argparse.ArgumentParser(
        description= "Generate datasets, process raw data and generate figures.")
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    generate_parser = subparsers.add_parser('generate')
    # TODO: something like this will be useful to run smaller runs for
    # testing purposes and to control the number of processes used.
    # generate_parser.add_argument(
    #     '-n', help="number of replicates", type=int, default=-1)
    # generate_parser.add_argument(
    #     "--processes", '-p', help="number of processes",
    #     type=int, default=1)
    generate_parser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    generate_parser.set_defaults(func=run_generate)

    process_parser = subparsers.add_parser('process')
    process_parser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the simulation identifier')
    process_parser.set_defaults(func=run_process)

    figure_parser = subparsers.add_parser('figure')
    figure_parser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the figure dentifier')
    figure_parser.set_defaults(func=run_plot)

    args = parser.parse_args()

    k = args.name[0]
    if k == "all":
        classes = datasets
        if args.func == run_plot:
            classes = figures
        for name, cls in name_map.items():
            if cls in classes:
                print("Running:", name)
                args.func(cls, args)
    else:
        cls = name_map[k]
        args.func(cls, args)

if __name__ == "__main__":
    main()
