# Evaluating tsinfer

This repository contains all code used for evaluating [tsinfer](https://tsinfer.readthedocs.io/en/latest/)
and producing the figures in the 
[Inferring the ancestry of everyone](https://www.biorxiv.org/content/10.1101/458067v1) preprint. This 
includes code for comparing tsinfer with other tools on simulated data, as well as building 
tree sequences from human data. Because of the complexity of downloading and preparing 
real data, this code is kept isolated in the ``human-data`` directory. Except for the human 
data pipeline, everything is intended to be run from the repository root.

We assume from here on that you have cloned this github repository into a directory called e.g. 
`treeseq-inference`, and are running commands from within it, e.g. using

```
git clone https://github.com/mcveanlab/treeseq-inference
cd treeseq-inference
```

### Requirements

Code is primarily written in Python and requires Python >= 3.4. For benchmarking,
we  use an R package via [rpy2](https://rpy2.readthedocs.io/) and so a working 
R installation is also required. Some external software requires other software too
(e.g. cmake for SLiM, the GNU scientific library for msprime/tsinfer). 
These are detailed below.

#### Installing system prerequisites 
You will need to install install python (3) with pip and the GNU scientific library (`gsl`).
For benchmarking, you will need R and cmake. For testing in real-world data such as the 
human analyses, you will need cython and the curl libraries.

For example, to install all these on Ubuntu:

```
# Install pip, GNU scientific library, R, cmake for SLiM, cython & curl libs for cyvcf2
sudo apt-get install python3-pip libgsl-dev r-base-core cmake cython3 libssl-dev libcurl4-openssl-dev
```

#### Installing required python modules

The Python packages required are listed in the ``requirements.txt`` file. These can be 
installed with

```
$ python3 -m pip install -r requirements.txt
$ python3 -m pip install cyvcf2 # only for human data analysis: needs to be installed *after* numpy
```

if your are using pip. Conda may also be used to install these dependencies.

## Human data

Please see the [README](human-data/README.md) in the ``human-data`` directory 
for details on running the human data pipelines.

## Benchmarking using simulations

### Requirements

For calculating ARG distance metrics, we require an `R` installation with certain packages, as well as
our own `ARGmetrics` package. To test other ARG inference software packages, we require them to be
available in the ``tools`` directory. Simple installation instructions for setting this up are below.

#### Installing R requirements

We require the `ape`, `phangorn`, and `Rcpp` packages. If you don't already have these installed
in your local R installation, you should be able to install them using:

```
# Install latest required packages within R - this recompiles stuff so may take a few mins
sudo R -e 'install.packages(c("ape", "phangorn", "Rcpp"), repos="https://cran.r-project.org", INSTALL_opts="--byte-compile")'
```

You can then install our local `ARGmetrics` package, bundled in this github repository, by running `R CMD INSTALL` 
from within the github directory, as follows:

```
# Install ARGmetrics into R
sudo R CMD INSTALL ARGmetrics
```

If can't run sudo, because you do not have superuser (root) access to your machine, you should be able to 
download pakages into a local R library folder by running the `install.packages(...)` command above 
from within an R session, which will prompt you to create a local folder. Any further installation can then be
done without `sudo`.

#### Installing alternative evaluation tools

The tools that we compare against are kept in the ``tools`` directory and can be 
downloaded and built using 

```
$ make -C tools
```



### Running evaluations

The code for running simulations is held in the ``src/evaluations.py`` file
which has several subcommands. For help, run

```
$ python3 src/evaluation.py --help
```

**TODO** Give desription of the various commands and give an example of 
running one of the evaluations from start to finish.

