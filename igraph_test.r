
#pipeline
library(magrittr)

library(igraph)
links <- read.table("./output/rrfile/P30443_df01_2.txt")
#rownames(links) <- colnames(links)
#a<- graph_from_data_frame(d=links,vertices=seq,directed=T)
links2 <- as.matrix(links)
net5 <- graph_from_adjacency_matrix(links2,mode = 'undirected',diag=F)

plot(net5,vertex.size=3,vertex.label=NA)
tkplot(net5,vertex.size=2)

rglplot(net5,vertex.size=2)


