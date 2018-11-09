all: | fastARG argweaver RentPlus SLiM
	#the following two tools require stand-alone compilation
	make -C fastARG
	make -C argweaver
	make -C SLiM/build

fastARG: 
	git clone https://github.com/lh3/fastARG
	# checkout a specific revision (lastest as of Sept 2018)
	cd fastARG && git checkout 841b4971cdf33932f07044b1c6c2fb14000bdd94

argweaver: 
	git clone https://github.com/CshlSiepelLab/argweaver
	# checkout a specific revision (lastest as of Sept 2018)
	cd argweaver && git checkout 905d3280b0ce3d799767b9511862b39294b0c093
	
RentPlus:
	git clone https://github.com/SajadMirzaei/RentPlus
	# checkout a specific revision (lastest as of Sept 2018)
	cd RentPlus && git checkout 788287300e269921c37bd8e59d3c624f23806ace

SLiM:
	#this is used for the file selective_sweep.py
	git clone --branch v3.2 https://github.com/MesserLab/SLiM
	mkdir SLiM/build && cd SLiM/build && cmake ../

