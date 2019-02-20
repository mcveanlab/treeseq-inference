# Human data pipelines

This directory contains the pipelines used to convert human data into 
``tsinfer`` input format and run tsinfer. It is based around a Makefile,
which breaks up the various steps and can avoid repeating parts of the
analysis.

We assume that all commands are run **within** this directory.

## Requirements

The pipelines require:

- bcftools
- tabix
- Python 3

The Python package requirements are listed in the ``requirements.txt`` file 
in the repository root.


## 1000 Genomes

To build a tree sequence for 1000 Genomes chromosome 20, run:

```
$ make 1kg_chr20.trees
```

This will download the data, and create the intermediate files and 
output the final tree sequence ``.trees`` file. Converting from 
VCF to the ``.samples`` file will take some time.

The number of threads used in steps of the pipeline that support 
threading can be controlled using the ``NUM_THREADS`` make 
variable, i.e.

```
$ make NUM_THREADS=20 1kg_chr20.trees
```

will tell ``tsinfer`` to use 20 threads where appropriate.

The pipeline for 1000 genomes is generic, and so any chromosome can be built
by running, e.g., ``make 1kg_chr1.trees``.


## SGDP

The pipeline for the SGDP chromosome 20 data can executed using

```
$ make sgdp_chr20.trees
```

## UKBB

The pipeline for UKBB can be run in a similar way, but because it is not a public
dataset the input data cannot be downloaded automatically. See the Makefile for 
the details of the required files and paths.
