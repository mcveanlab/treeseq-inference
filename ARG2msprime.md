Correspondance between Ancestral Recombination graphs and coalescence records
------------



An ancestral recombination graph is classically represented as a series of coalescence nodes and recombination node, where a coalescence node has one parent and two children, and a recombination node has two parents, one child, and additional parameters which specify which parts of the genome are taken from which parent.

An example of such a format is the .arg format introduced by ARGweaver, documented at http://mdrasmus.github.io/argweaver/doc/. The example .arg file used

```
start=0	end=10000
nd	event	age (gens)			pos		parents	children
1	coal	122.586947411		0		2		n2,n1
2	coal	1545.31450861		0		4		n3,1
3	coal	12061.8085146		0		7		5,4
4	recomb	8051.70236778		536		3,5		2
5	coal	12061.8085146		0		3		6,4
6	recomb	8051.70236778		1033	5,7		n0
7	coal	26970.5983226		0		3,6
n0	gene	0.0	0	6	
n1	gene	0.0	0	1	
n2	gene	0.0	0	1	
n3	gene	0.0	0	2	
```

Where the 'pos' column is only relevant for recombination nodes, and marks a recombination breakpoint such that the first listed parent inherits the sequence to the left of that point (non-inclusive), and the second listed parent inherits sequence to the right (inclusive).  By ARGweaver convention the node numbers are usually allocated so that the node contributing the left hand side of the sequence is given a lower id number than that which contributes the right hand portion.

Also note in this case that node 7 has no parent, which marks it as the root node.

Such an ARG can be represented in a graph format, for instance as follows



But an alternative representation is to 