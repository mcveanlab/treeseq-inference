#!/usr/bin/env python3
"""test the example file from fastARG
the fastARG output from ex1.txt is 

E	11
N	4	3
R	4	2	0	3	0
C	5	4	0	2	0
C	5	3	0	3	0
C	6	0	0	3	1	0
C	6	4	2	3	0
C	7	1	0	3	0
C	7	5	0	3	1	1
C	8	6	0	3	0
C	8	7	0	3	1	2
S	8	000

(C & R nodes have [node, child, start, end, mut])

msprime input (http://msprime.readthedocs.io/en/stable/api.html#msprime.load_txt) has 
[start  end  node, child1,child2  time population]

so the corresponding input is

0 3 4 2    4 0
2 3 5 3    5 0
0 2 5 3,4  5 0
0 2 6 0    6 0
2 3 6 0,4  6 0
0 3 7 1,5  7 0
0 3 8 6,7  8 0

Where sequential records that don't span the same chunk of genome, like

C	5	4	0	2	0
C	5	3	0	3	0

have been broken up into pieces.

C	5	4	0	2	0
C	5	3	0	2	0
C	5	3	2	3	0
 

"""
import sys
sys.path.insert(0,'../msprime/') # use the local copy of msprime in preference to the global one
import tempfile
import msprime

records = """\
0 3 4 2    4 0
2 3 5 3    5 0
0 2 5 3,4  5 0
0 2 6 0    6 0
2 3 6 0,4  6 0
0 3 7 1,5  7 0
0 3 8 6,7  8 0
"""

with tempfile.NamedTemporaryFile("w+") as f:
    f.write(records)
    f.flush()
    ts = msprime.load_txt(f.name)
print("ORIGINAL")
print(" Tree:")
for t in ts.trees():
    print(t.interval, t)
print(" Coalescense records:")
for r in ts.records():
    print("{}".format(r))


print("MINIMAL")
subset = ts.subset([0,1,2,3])
print(" Tree:")
for t in subset.trees():
    print(t.interval, t)
print(" Coalescense records:")
for r in subset.records():
    print("{}".format(r))
