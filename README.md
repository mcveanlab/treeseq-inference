# Evaluating tsinfer

This repository contains all code used for evaluating [tsinfer](https://tsinfer.readthedocs.io/en/latest/)
and producing the figures in the 
[Inferring the ancestry of everyone](https://www.biorxiv.org/content/10.1101/458067v1) preprint. This 
includes code for comparing tsinfer with other tools on simulated data, as well as building 
tree sequences from human data. Because of the complexity of downloading and preparing 
real data, this code is kept isolated in the ``human-data`` directory. Except for the human 
data pipeline, everything is intended to be run from the repository root.

## Human data

Please see the [README](human-data/README.md) for details on running the 
human data pipelines.


## Benchmarking using simulations

### Requirements

The tools that we compare against are kept in the ``tools`` directory and can be 
downloaded and built using 

```
$ make -C tools
```

The benchmarking code is primarily written in Python and requires Python >= 3.4. We
also use an R package via [rpy2](https://rpy2.readthedocs.io/) and so a working 
R installation is also required. First, install the local R package by running:

```
$ R
> library(devtools)
> install("ARGmetrics")
```

We have assumed that the ``devtools`` package is installed---you will need to 
install this first if it is not present.

The Python packages required are listed in the ``requirements.txt`` file. These can be 
installed with

```
$ python3 -m pip install -r requirements.txt
```

if your are using pip. Conda may also be used to install these dependencies.

Once all the packages requirements have been installed, you can test it out 
by running 

```
$ python3 src/evaluation.py
```

### Running evaluations

The code for running simulations is held in the ``src/evaluations.py`` file
which has several subcommands. For help, run

```
$ python3 src/evaluation.py --help
```

**TODO** Give desription of the various commands and give an example of 
running one of the evaluations from start to finish.

