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

Note that having *both* parent and child nodes listed is redundant. Only one or the other is sufficient to reconstruct the ARG. Also that the node with no parent is the root node (node 7 in the example above). The general principle is that a coalescence node requires a time, an ID, and two children, whereas a recombination node requires an ID, a single child, and a recombination breakpoint. 

An ARG such as the above can be represented in a graph format, for instance using the following R code

```
library(visNetwork)

con <- file('/Users/yan/Documents/Research/Wellcome/treeseq-inference/test_files/ARGweaver_test.arg')
open(con)

line1 <- strsplit(readLines(con,1),'\t')
seq.lims <- data.frame(row.names=1,do.call(rbind,strsplit(line1[[1]],'=')), stringsAsFactors=FALSE)
seq.lims[,1] <- as.numeric(seq.lims[,1])

arg <- read.delim(con, colClasses=c(name='character', age='numeric', parents='character', pos='numeric', children='character'), fill=T)

rownames(arg) <- arg[['name']]

edge_output <- function(x) {
 #take a list with 'name'='1','parents='2,3', etc
 if(x[['parents']]!="") {
   parents<-strsplit(x[['parents']],",")[[1]]
   if (x[['event']]=='recomb') {
     split.indexes <- cbind(seq_along(parents), seq_along(parents)+1)
     splits=c(seq.lims['start',], x[['pos']], seq.lims['end',])
     lower <- splits[split.indexes[1,]]
     upper <- splits[split.indexes[2,]]
   } else {
     lower <- upper <- NA
   }
   label <- ifelse(is.na(lower) |is.na(upper), NA, paste(lower,upper,sep="-"))
   data.frame(from=parents, to=x[['name']], lower, upper, label)
  }
}

edges <- data.frame(do.call('rbind',by(arg, arg$name, edge_output)))

height = 100 #max height, pre-logging

nodes <- data.frame(id=arg$name, label=arg$name, level=log(arg$age/max(arg$age)*height+1), shape=ifelse(arg$event=='recomb','box','circle'), physics=ifelse(arg$parent=='',FALSE,FALSE), x=NA, stringsAsFactors = FALSE) 

nw <- visNetwork(nodes, edges, main="ARG representation") %>% 
  visEdges(arrows = "to", font=list(align='middle')) %>% 
  visNodes(shapeProperties=list(borderRadius=2)) %>% 
  visIgraphLayout(layout="layout_as_tree", flip.y=FALSE) %>%
  visHierarchicalLayout(direction="DU") %>%
  visInteraction(zoomView=TRUE, dragView=FALSE)

visSave(nw, ARG.html, FALSE)
```

But an alternative representation is to place lines