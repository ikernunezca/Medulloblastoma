ascuas <- commandArgs()[6] # Patient number(1 to 38) to obtain its enrichments

# Load depencencies

library(neat)
library(igraph)

# Load input data and put it into a unique data.frame
a <- read.table("data_for_enrichments/genes_per_patient_final.csv",sep="\t")
b <- read.table("data_for_enrichments/patient_clusters_final.csv",sep= ",",header=T)
colnames(a)[1] <- "names"
colnames(b)[1] <- "names"
genes_per_patient_and_clusters <- merge(b,a,by="names",all = T)


# Load layers
biogrid <- read.table("biogrid_interactions.csv",sep = "\t",header = F, stringsAsFactors=F)
class(biogrid[,1]) <- "character"
class(biogrid[,2]) <- "character"
class(biogrid[,3]) <- "character"
biogrid_graph <- graph_from_data_frame(biogrid,directed = FALSE)
E(biogrid_graph)$interaction_type <- biogrid[,3]

#############################################

# Explore enrichments in protein-interaction Data

# Obtain all vertices of biogrid ppi layer
proteins_available <- names(V(biogrid_graph)) # Actually it is the genes available in the network (nodes)

# genes per patient appearing in the ppi network
genes_lists <- list()
for(i in 1:nrow(genes_per_patient_and_clusters)){
  genes <- strsplit(as.character(genes_per_patient_and_clusters[i,4]),"_")[[1]]
  genes_lists[[i]] <- unique(genes[which(genes %in% proteins_available)])
}
names(genes_lists) <- genes_per_patient_and_clusters[,1] # List of genes associated to each patient appearing in the drug network

# For each protein, obtain adjacent proteins

genes_per_protein <- lapply(adjacent_vertices(graph = biogrid_graph,v = names(V(biogrid_graph))),function(x)
names(x))

message(paste0("Starting protein interaction enrichment for patient ",names(genes_lists)[as.numeric(ascuas)]))
Sys.time()
protein_enrichments <- neat(genes_lists[as.numeric(ascuas)],blist = genes_per_protein,network = biogrid_graph,nettype = 'undirected',nodes=proteins_available,mtc.type = 'fdr',alpha = 0.05)
message(paste0("Finished protein interaction enrichment for patient ",names(genes_lists)[as.numeric(ascuas)]))
Sys.time()

# Save enrichment results
write.table(protein_enrichments,file=paste0("data/Protein_enrichments/protein_enrichments_",names(genes_lists)[as.numeric(ascuas)],".csv"),sep="\t",col.names=F,row.names=F)
