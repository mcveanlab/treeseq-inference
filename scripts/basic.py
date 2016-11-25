#!/usr/bin/env python3
import sys
sys.path[0:0] = '../msprime/' # puts the /foo directory at the start of your path
import tempfile
import msprime

records = """\
0   2   3   0,1     1   0
2   3   4   1,2     2   0
0   2   4   2       2   0
0   3   5   3,4     3   0
"""

with tempfile.NamedTemporaryFile("w+") as f:
    f.write(records)
    f.flush()
    ts = msprime.load_txt(f.name)

for t in ts.trees():
    print(t.interval, t)
