help:
	echo WRITE SOME HELP

deps:
	make -C src
	# this should download and compile fastARG, ARGweaver, RentPlus, ftprime, etc
	make -C tools
