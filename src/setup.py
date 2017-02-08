from __future__ import division
from __future__ import print_function

import subprocess
import numpy as np

# First, we try to use setuptools. If it's not available locally,
# we fall back on ez_setup.
try:
    from setuptools import setup, Extension
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, Extension

long_description = "TODO"

_tsinfer_module = Extension(
    '_tsinfer',
    sources=["_tsinfermodule.c", "ls.c"],
    # Enable asserts by default.
    undef_macros=["NDEBUG"],
    libraries=["m"],
    include_dirs = [np.get_include()],
)

setup(
    version = 0.2,
    name="tsinfer",
    description="Infer tree sequences from genetic data.",
    long_description=long_description,
    author="Jerome Kelleher",
    author_email="jerome.kelleher@well.ox.ac.uk",
    url="http://pypi.python.org/pypi/tsinfer",
    ext_modules=[_tsinfer_module],
)
