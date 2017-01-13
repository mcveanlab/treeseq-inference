
deps:
	make -C fastARG
	# This will ultimately be redundant as we will pip install 
	# a released version.
	make -C msprime ext3
	# argweaver is Python 2 only, so we don't pip install it.
	cd argweaver && python2 setup.py build_ext --inplace

