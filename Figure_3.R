# Load dependencies
library(CmmD)
library(dendroextras)
library(dendextend)
library(circlize)

# Create a vector with the paths where Molti's Output files are saved. 
structures_12 <- paste0("data/Molti_Output/",seq(0.5,12,0.5),".csv")
genes_of_interest <- unique(as.character(read.table(file = "data/gene_mentions.lst")[,1])) # Genes of interest come from text-mining activities of the paper.
curie_to_12_full <- CmmD_from_community_structures(nodelist = NULL, community_structures = structures_12, resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = "hamming",threads = 7)
curie_to_12_full$hamming_distance_matrix = as.matrix(curie_to_12_full$distance_matrix) * 24 # This transformation is needed because parallel dist is weighted.
# 24 = length(seq(0.5,12,0.5)) -> number of resolution values analyzed


# Dendrogram of Molti's Clustering. Final plot (Figure 3) is a representation on how communities dissagregate during the process.
# Index in rows of distance matrix (d_matrix) is the same as curie_to_12$constant
d_matrix <- matrix(NA,ncol=length(curie_to_12$l_constant),nrow=length(curie_to_12$l_constant))
for(i in 1:nrow(d_matrix)){
  a <- strsplit(names(curie_to_12$l_constant[i]),"_")[[1]]
  for(j in 1:ncol(d_matrix)){
    b <- strsplit(names(curie_to_12$l_constant[j]),"_")[[1]]
    d_matrix[i,j] <- hamming.distance(a,b)
    d_matrix[j,i] <- d_matrix[i,j]
  }
  print(i)
}
hc <- hclust(as.dist(d_matrix),"ward.D2")

set.seed(2019) # Set the seeds to generate the same dendrogram

# Generate the dendrogram object
dend3 <- as.dendrogram(hc)
# Generate a vector of colors (as) for the names of the trajectories. Except for the ones we are looking for, the rest are set to white in order to hide them.
as <- rep("white",length(curie_to_12$l_constant))
# Get the trajectory of each gene of interest
wnt <- which(hc$order== grep("\\b1499\\b",curie_to_12$l_constant))
shh <- which(hc$order==grep("\\b6469\\b",curie_to_12$l_constant))
mycn <- which(hc$order==grep("\\b4613\\b",curie_to_12$l_constant))
myc <- which(hc$order==grep("\\b4609\\b",curie_to_12$l_constant))
cdk6 <- which(hc$order==grep("\\b1021\\b",curie_to_12$l_constant))
erbb4 <- which(hc$order==grep("\\b2066\\b",curie_to_12$l_constant))
src <- which(hc$order==grep("\\b6714\\b",curie_to_12$l_constant))
# Set red color for trajectories of interest in color vector (as)
as[wnt] <- "red"
as[shh] <- "red"
as[mycn] <- "red"
as[myc] <- "red"
as[cdk6] <- "red"
as[erbb4] <- "red"
as[src] <- "red"
# Generate a letter size vector for the names of the trajectories, the same actions as with the color.
letter_size <- rep(0.1,888)
letter_size[wnt] <- 0.8
letter_size[shh] <- 0.8
letter_size[mycn] <- 0.8
letter_size[myc] <- 0.8
letter_size[cdk6] <- 0.8
letter_size[erbb4] <- 0.8
letter_size[src] <- 0.6
# Keep original labels
old_labels <- labels(dend3)
old_labels[wnt] <- "CTNBB1"
old_labels[shh] <- "SHH"
old_labels[mycn] <- "MYCN"
old_labels[myc] <- "MYC"
old_labels[cdk6] <- "CDK6"
old_labels[erbb4] <- "ERBB4"
old_labels[src] <- "CTNNB1 - SRC"
dend4 <- dend3 %>%
  set("labels", seq_len(nleaves(dend3))) %>% # Generate new trajectory labels in order to use the next set statements.
  set("by_labels_branches_col",value=c(wnt,shh,mycn,myc,cdk6,erbb4,src),type="any") %>% # Set red colors for th etrajectories
  set("by_labels_branches_lwd",value=1:length(curie_to_12),type="any",TF_values = 0.3) %>%
  set("by_labels_branches_lwd",value=c(wnt,shh,mycn,myc,cdk6,erbb4,src),type="any",TF_values = 1) %>%
  set("labels_col",as) %>% # Set labels color
  set("labels_cex", letter_size) %>% # Set labels names
  set("hang") %>%
  set("labels",old_labels) # Re-set original labels.

xy <- dend4 %>% get_nodes_xy()
is_internal_node <- is.na(dend4 %>% get_nodes_attr("leaf"))
is_internal_node[which.max(xy[,2])] <- FALSE
xy <- xy[is_internal_node,]
par(mar=c(0,0,0,0))
# Plot
pdf("data/Plots/Figure_3.pdf")
plot(circlize_dendrogram(dend4))
# plot(dend4)
# text(xy[,1]+.2, xy[,2]+.2, labels=format(xy[,2], digits=2), col="red")
dev.off()
