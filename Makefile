
deps:
	make -C src
	make -C fastARG
	make -C argweaver
	# This will ultimately be redundant as we will pip install 
	# a released version.
	make -C msprime ext3

