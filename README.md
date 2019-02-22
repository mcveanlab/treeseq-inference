# Evaluating tsinfer

This repository contains all code used for evaluating [tsinfer](https://tsinfer.readthedocs.io/en/latest/)
and producing the figures in the 
[Inferring the ancestry of everyone](https://www.biorxiv.org/content/10.1101/458067v1) preprint. This 
includes code for comparing tsinfer with other tools on simulated data, as well as building 
tree sequences from human data. Because of the complexity of downloading and preparing 
real data, this code is kept isolated in the ``human-data`` directory. Except for the human 
data pipeline, everything is intended to be run from the repository root.

We assume from here on that you have cloned this github repository into a directory called 
`treeseq-inference`, e.g. using

```
git clone https://github.com/mcveanlab/treeseq-inference
```

## Human data

Please see the [README](human-data/README.md) in the ``human-data`` directory 
for details on running the human data pipelines.

## Benchmarking using simulations

### Requirements

The benchmarking code is primarily written in Python and requires Python >= 3.4. We
also use an R package via [rpy2](https://rpy2.readthedocs.io/) and so a working 
R installation is also required. 

#### Installing system prerequisites 
To install msprime & tsinfer you need to have the GNU scientific library (`gsl`) installed.
To calculate tree distance metrics, you will need to have `R` installed, with packages
as described below, and the `rpy2` python-to-R library. To install SLiM for simulating
selection you will need to install cmake, and to install the `cyvcf` library to read VCF
files, you will need  `curl` libraries too. These can all be installed e.g. on Ubuntu by:

```
# Install GNU scientific library, R and python2r interface, cmake for SLiM, cython & curl libs for cyvcf2
sudo apt-get install libgsl-dev r-base-core python3-pip cython3 cmake libssl-dev libcurl4-openssl-dev
```

#### Installing required python modules

The Python packages required are listed in the ``requirements.txt`` file. These can be 
installed with

```
$ python3 -m pip install -r treeseq-inference/requirements.txt
```

if your are using pip. Conda may also be used to install these dependencies.

#### Installing R requirements

For calculating ARG distance metrics, we require the `ape`, `phangorn`, and `Rcpp` packages.
To get the latest versions of these, if you don't have them installed already in your local
R installation, you should be able to do something like the following

```
# Install latest required packages within R - this recompiles stuff so may take a few mins
sudo R -e 'install.packages(c("ape", "phangorn", "Rcpp"), repos="https://cran.r-project.org", INSTALL_opts="--byte-compile")'
```

You can then install our local `ARGmetrics` package, bundled in this github repository.
Assuming this repository is in `treeseq-inference`, simply do

```
# Install ARGmetrics into R
sudo R CMD INSTALL treeseq-inference/ARGmetrics
```

If you don't have superuser (root) access to your machine, you should be able to [set a local R library folder]() 
local folder

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

