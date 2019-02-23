# Evaluating tsinfer

This repository contains all code used for evaluating [tsinfer](https://tsinfer.readthedocs.io/en/latest/)
and producing the figures in the preprint:
[Inferring the ancestry of everyone](https://www.biorxiv.org/content/10.1101/458067v1). This 
includes code for comparing tsinfer with other tools using simulated data, as well as building 
tree sequences from human data. Because of the complexity of downloading and preparing 
real data, this code is kept isolated in the ``human-data`` directory. Except for the human 
data pipeline, everything is intended to be run from the repository root.

We assume from here on that you have cloned this github repository into a directory called e.g. 
`treeseq-inference`, and are running commands from within it, e.g. using

```
git clone https://github.com/mcveanlab/treeseq-inference
cd treeseq-inference
```

### General requirements

Code is primarily written in Python and requires Python >= 3.4. For benchmarking,
we  use an R package via [rpy2](https://rpy2.readthedocs.io/) and so a working 
R installation is also required. Some external software requires other libraries too
(e.g. cmake for SLiM, the GNU scientific library for msprime/tsinfer). 
These are detailed below.

#### Installing system prerequisites 
You will need to install install Python (version 3) with pip and the GNU scientific library (`gsl`).
For benchmarking, you will need R and cmake. Testing real-world data (e.g. for the
human analyses) uses the `cyvcf2` Python module, which requires pre-installation of 
cython and the curl libraries.

For example, to install all these on Ubuntu:

```
# Ubuntu-only: install pip, GNU scientific library, R, cmake for SLiM, cython & curl libs for cyvcf2
sudo apt-get install python3-pip libgsl-dev r-base-core cmake cython3 libssl-dev libcurl4-openssl-dev
```

#### Installing required python modules

The Python packages required are listed in the ``requirements.txt`` file. These can be 
installed with

```
$ python3 -m pip install -r requirements.txt
```

if you are using pip. Conda may also be used to install these dependencies.

## Human data

You will need the `cyvcf2` Python module to read VCF files. Once the requirements above have been installed you should simply be able to do:

```
$ python3 -m pip install cyvcf2 # only for human data analysis: needs to be installed *after* numpy
```

Please see the [README](human-data/README.md) in the ``human-data`` directory 
for further details on running the human data pipelines.

## Simulation benchmarks

### Requirements

For calculating ARG distance metrics, we require an `R` installation with certain packages, as well as
our own `ARGmetrics` package. To test other ARG inference software packages, we require them to be
available in the ``tools`` directory. Simple installation instructions for setting this up are below.

#### Installing necessary R packages

We require the `ape`, `phangorn`, and `Rcpp` packages. If you don't already have these installed
in your local R installation, you should be able to install them in the standard way. For example,
from within R you can issue the command 
`install.packages(c("ape", "phangorn", "Rcpp"), INSTALL_opts="--byte-compile")`, which may prompt
you for various bits of information, if they are not already known, e.g. your choice of CRAN repository.
If you have superuser (root) access to your machine, you can install packages without requiring any user interaction by

```
# Install latest required packages within R - this recompiles stuff so may take a few mins
sudo R -e 'install.packages(c("ape", "phangorn", "Rcpp"), repos="https://cran.r-project.org", INSTALL_opts="--byte-compile")'
```

You can then install our local `ARGmetrics` package, bundled in this github repository, by running `R CMD INSTALL` 
from within the github directory, as follows:

```
# Install ARGmetrics into R (you can omit `sudo` if installing locally)
sudo R CMD INSTALL ARGmetrics
```

#### Installing alternative ARG inference software

We compare our results against [ARGweaver](https://github.com/CshlSiepelLab/argweaver), [RentPlus](https://github.com/SajadMirzaei/RentPlus), and [fastARG](https://github.com/lh3/fastARG). 
We also use [SLiM](https://github.com/MesserLab/SLiM) to run forwards simulations. These stand-alone
software tools are kept in the ``tools`` directory and can be downloaded and built using 

```
# Download and compile other simulation and inference tools for testing
$ make -C tools
```

### Running evaluations

The code for running simulations is held in the ``src/evaluations.py`` file
which has several subcommands. For help, run

```
$ python3 src/evaluation.py --help
```

#### An example evaluation

All evaluations have 3 steps: setup, inference, and summarizing. The final summarised data can then be plotted using the script in `src/plots.py`. 

##### Setup
Run simulations to generate known ancestries, and sample haplotype files for ancestral 
inference. To see the various evaluation datasets that can be generated, run 
`python3 src/evaluation.py setup -h` (the special dataset "all" will generate all datasets,
see below). Here we show an example with the "all_tools" dataset, using the additional switches

* `-P` to show a progress monitor
* `-r 2` to only run 2 replicates (rather than hundreds), for speed
```
python3 src/evaluation.py setup -P -r 2 all_tools
```

##### Infer
Run inferences for various combinations of parameters / inference tools.
For the "all_tools" dataset, `tsinfer` plus 3 other inference tools are run.
While `tsinfer` takes only a few seconds or minutes to run, others (especially ARGweaver)
may may take a number of hours. The estimated time remaining is output as part of the
progress monitor.
```
python3 src/evaluation.py infer -P all_tools
```
In general, the inference step takes the most time. This script can be killed then rerun:
it will resume any not-yet-completed inference tasks.

##### Summarize
Data from inference is stored in a large csv file. Running the `summarize` command
takes this generated data and slims it down into single summary csv file corresponding
to a single evaluation plot.
```
python3 src/evaluation.py summarize metrics_all_tools
```

##### Plot using the csv file
```
python3 src/plot.py metrics_all_tools
```

The result should be an appropriately named pdf or png file in the `figures` directory 
(e.g. `figures/metrics_all_tools.pdf`)

#### Running all evaluations

To produce all the data in our paper, run the following, in order

```
python3 src/evaluation.py setup -P all # will take many hours/days
python3 src/evaluation.py infer -P all # will take many days/weeks
python3 src/evaluation.py summarize all #will take a few minutes
```

You can speed up the evaluations by using multiple processors, specified using the `-p` flag.
For instance, on a 64 core machine, using all cores:

```
python3 src/evaluation.py setup -p 64 all # will take a few hours
python3 src/evaluation.py infer -p 64 all # will take a few days (mostly to run ARGweaver)
python3 src/evaluation.py summarize all #will take a few minutes
```

The final figures can then be plotted using

```
python3 src/plot.py all
```
