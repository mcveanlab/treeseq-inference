
from distutils.core import setup, Extension
import sys
import pathlib
import numpy as np

conda_base = pathlib.Path(sys.executable).parents[1]

module1 = Extension(
    'simplebgen', 
    libraries=["bgen", "zstd", "athr"],
    include_dirs=[str(conda_base / "include"), np.get_include()],
    library_dirs=[str(conda_base / "lib")],
    sources = ['simplebgenmodule.c'])

setup(
    name='simplebgen',
    version='0.0.1',
    description='Simple access to bgen data with numpy arrays.',
    ext_modules=[module1])
