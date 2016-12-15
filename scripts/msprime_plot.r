#plot an msprime simulation from the coalescense record output format 



c.records <- read.delim("/Users/yan/Documents/Research/Wellcome/treeseq-inference/test_files/AWtest_med.msprime", stringsAsFactors=FALSE)

edge_output <- function(x) {
 #take a list with 'name'='1','parents='2,3', etc
   children<-as.numeric(strsplit(x[['children']],",")[[1]])
   lower <- x[['left']]
   upper <- x[['right']]
   label <- ifelse(is.na(lower) |is.na(upper), NA, paste(lower,upper,sep="-"))
   data.frame(from=x[['node']], to=children, lower, upper, label)
}

edges <- data.frame(do.call('rbind',by(c.records, rownames(c.records), edge_output)))

edges2 <- data.frame(do.call(rbind,by(edges, paste(edges$from, edges$to), function(x) (colMeans(x[,1:2])))))

height = 100 #max height, pre-logging

ids <- unique(c(edges2$from, edges2$to))
nodes <- data.frame(id=ids, label=ids, row.names=ids)
times <- by(c.records, c.records$node, function(x) mean(x$time))
nodes$level <- times[rownames(nodes)]/5000
#nodes$level <- log(c.records[rownames(nodes),'time']/max(c.records$time)*height+1)
nodes$level <- ifelse(is.na(nodes$level), 0, nodes$level) #set NA nodes (tips) to 0 - an msprime restriction for text outputs

nw <- visNetwork(nodes, edges2, main="Representation of ARGweaver 8 tip simulation as coalescence records") %>% 
  visNodes(shape = "circle") %>% 
  visEdges(arrows = "to", font=list(align='middle')) %>% 
  visIgraphLayout(layout="layout_as_tree", flip.y=FALSE) %>%
  visNodes(physics = TRUE) %>% 
  visPhysics(hierarchicalRepulsion=list(nodeDistance=1000)) %>% 
  visHierarchicalLayout(direction="DU") %>%
  visInteraction(zoomView=TRUE, dragView=FALSE)
  nw
visSave(nw, "ARGweaver_medium_size_msprime.html", FALSE)