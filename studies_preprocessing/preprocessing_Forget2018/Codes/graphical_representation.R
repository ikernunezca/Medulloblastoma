 source("~/MultilayerNetworks/required_packages.R", echo=FALSE)
 load("~/MultilayerNetworks/data/MultilayerNetwork.RData")
 rm(list=setdiff(ls(), c("exp_files", "set.samples", "multi.layer.net.with.ggi",
				"GGI", "set.common.genes")))
 cytoscapePing()
 
 group.attributes <- names(exp_files)
 
 nodes <- data.frame(id = set.samples, 
                     group = "samples" ) 
 nodes.att.data <- mclapply(group.attributes, function(x){
    sel.id <- grep(x, multi.layer.net.with.ggi$to)
    nodes.att <- data.frame(id = as.character(unique(multi.layer.net.with.ggi$to[sel.id])), 
                            group= x)
    return(nodes.att)
 })
 nodes.att.data <- do.call(rbind, nodes.att.data)
 
 nodes.att.ggi <- data.frame(id= as.character(GGI$nodes$PARTICIPANT), 
                             group= "GGI")
 
 nodes<- rbind(nodes, nodes.att.data, nodes.att.ggi)
 nodes<- data.frame(nodes, 
               score=1 , 
               stringsAsFactors= FALSE)

 
 edges <- data.frame(source=multi.layer.net.with.ggi$from,
                     target=multi.layer.net.with.ggi$to,
                     group= multi.layer.net.with.ggi$layer,
                     weight=1, # numeric
                     stringsAsFactors=FALSE)

 #createNetworkFromDataFrames(nodes, edges) ## Takes a lot of time
 ## Too many nodes and edges.. 
 ig <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
 igraph::vcount(ig)
 igraph::ecount(ig)
 
 #create a subgraph  THIS IS JUST AN ILLUSTRATION 
 ## Select some nodes
 sel.nodes <- c(which(nodes$group=="samples")[1:3], which(nodes$group=="Methylation")[4:6], 
                which(nodes$group=="Phospho")[1:3], which(nodes$group=="Protein")[c(1:3,100)], 
                which(nodes$group=="RNAseq")[c(1:3,44)], 
                which(nodes$id %in% c("SFRI", "DBNL", "SLC10A3", "WAS", "FAP", "GFAP", "RPL4", "RPL41")) )
 sig <- subgraph(graph=ig, v=sel.nodes)
 createNetworkFromIgraph(sig)

 
 column <- 'group'
 values.node <- unique(nodes$group)
 col <- c('#FF00FF', '#00FFFF', '#FF329F', '#FFFF00', '#00FF00', '#FF0000')
 names(col) <- values.node
 shapes <- getNodeShapes()[1:length(values)]
 setNodeShapeMapping (column, values.node, shapes)
 setNodeColorMapping(column, table.column.values = values.node , colors= col, 
                     mapping.type = 'd')

 values.edge <- unique(edges$group)
 col.edg <- c(col[-1], '#000000') 
 setEdgeColorMapping(column, table.column.values = values.edge, colors =col.edg, 
                     mapping.type = 'd')
 setEdgeLineStyleDefault("DOT")
 
 