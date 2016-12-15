#!/usr/bin/env python3
"""test an example file from ARGweaver

start=0	end=10000
name	event	age	pos	parents	children
1	coal	122.586947411	0	2	n2,n1
2	coal	1545.31450861	0	4	n3,1
3	coal	12061.8085146	0	7	5,4
4	recomb	8051.70236778	536	3,5	2
5	coal	12061.8085146	0	3	6,4
6	recomb	8051.70236778	1033	5,7	n0
7	coal	26970.5983226	0		3,6
n0	gene	0.0	0	6	
n1	gene	0.0	0	1	
n2	gene	0.0	0	1	
n3	gene	0.0	0	2	

(see http://mdrasmus.github.io/argweaver/doc/ for format)

msprime input (http://msprime.readthedocs.io/en/stable/api.html#msprime.load_txt) has 
[start  end  node, child1,child2  time population]


Larger simulations can be run e.g. by 

bin/arg-sim -k 8 -L 100000 -N 10000 -r 1.6e-8 -m 1.8e-8 -o ARGweaver_medtest.arg

"""
import sys
sys.path.insert(0,'../msprime/') # use the local copy of msprime in preference to the global one
from tempfile import NamedTemporaryFile
import msprime
msprime_filename="MStest_med.msprime"
ARGweaver_filename="AWtest_med.msprime"


from msprime_ARGweaver import ARGweaver_arg_to_msprime_txt
with open("../test_files/ARGweaver_medtest.arg", 'r+') as arg, \
     NamedTemporaryFile("w+") as msprime_in:
    conversion_table = ARGweaver_arg_to_msprime_txt(arg, msprime_in)
        
    ts = msprime.load_txt(msprime_in.name)

'''
print("ORIGINAL")
print(" Tree:")
for t in ts.trees():
    print(t.interval, t)
print(" Coalescense records:")
for r in ts.records():
    print("{}".format(r))


print("MINIMAL")
subset = ts.simplify()
print(" Tree:")
for t in subset.trees():
    print(t.interval, t)
'''
print(" ARGweaver coalescence records output to {}".format(ARGweaver_filename))
sim = ts.simplify()
with open(ARGweaver_filename, 'w+') as f:
    sim.write_records(f, precision=12)

print(" msprime coalescence records output to {}".format(msprime_filename))
ts = msprime.simulate(sample_size=8, Ne=10000, mutation_rate=1.8e-8, recombination_rate=1.6e-8,length=100000)
with open(msprime_filename, 'w+') as f:
    ts.write_records(f, precision=12)
