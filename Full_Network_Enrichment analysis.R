library(neat)
library(igraph)

# Load input data and put it into a unique data.frame
  a <- read.table("data/genes_per_patient.csv",sep="\t")
  b <- read.table("data/patient_clusters.csv",sep= ",",header=T)
  colnames(a)[1] <- "names"
  colnames(b)[1] <- "names"
  genes_per_patient_and_clusters <- merge(b,a,by="names",all = T)
  write.table(genes_per_patient_and_clusters,file="Supplementary_Tables/Supplementary_Table_4.csv",sep="\t",col.names=T,row.names=F)

# Load all layers -Except for biogrid interactions, See Supplementary_Table_6.R-
  kegg <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/cirillodavide/gene_multilayer_network/master/networks/KEGG_drugs.19-10-2019.gr",sep = " ",header = F, stringsAsFactors=F))
  class(kegg[,1]) <- "character"
  class(kegg[,2]) <- "character"
  class(kegg[,3]) <- "character"
  kegg_graph <- graph_from_data_frame(kegg[,1:2],directed = FALSE)
  E(kegg_graph)$drug <- kegg[,3]
  mondo <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/cirillodavide/gene_multilayer_network/master/networks/MoNDO_diseases.19-10-2019.gr",sep = " ",header = F, stringsAsFactors=F))
  class(mondo[,1]) <- "character"
  class(mondo[,2]) <- "character"
  class(mondo[,3]) < "character"
  mondo_graph <- graph_from_data_frame(mondo,directed = FALSE)
  E(mondo_graph)$mondo_ID <- mondo[,3]
  reactome <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/cirillodavide/gene_multilayer_network/master/networks/Reactome_pathways.19-10-2019.gr",sep = " ",header = F, stringsAsFactors=F))
  class(reactome[,1]) <- "character"
  class(reactome[,2]) <- "character"
  class(reactome[,3]) <- "character"
  reactome_graph <- graph_from_data_frame(reactome,directed = FALSE)
  E(reactome_graph)$reactome_ID <- reactome[,3]
  recon3d <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/cirillodavide/gene_multilayer_network/master/networks/Recon3D_metabolites.19-10-2019.gr",sep = " ",header = F, stringsAsFactors=F))
  class(recon3d[,1]) <- "character"
  class(recon3d[,2]) <- "character"
  class(recon3d[,3]) <- "character"
  recon3d_graph <- graph_from_data_frame(recon3d,directed = FALSE)
  colnames(kegg)[3] <- "kegg"
  colnames(mondo)[3] <- "mondo"
  colnames(reactome)[3] <- "reactome"
  colnames(recon3d)[3] <- "recon3d"

############################################# 
# 1. Explore enrichments in Drug Data 
    # Obtain all vertices of kegg drugs layer
    drugs_available <- names(V(kegg_graph))
    
    # genes per patient appearing in the drug network
    genes_lists <- list()
    for(i in 1:nrow(genes_per_patient_and_clusters)){ 
    genes <- strsplit(as.character(genes_per_patient_and_clusters[i,4]),"_")[[1]]
    genes_lists[[i]] <- unique(genes[which(genes %in% drugs_available)])
    }
    names(genes_lists) <- genes_per_patient_and_clusters[,1] # List of genes associated to each patient appearing in the drug network
    
    # For each drug, obtain genes associated because of it
    genes_per_drug <- list()
    drugs <- unique(kegg[,3])
    for(i in 1:length(drugs)){
      droga <- drugs[i]
      genes_per_drug[[i]] <- unique(unlist(kegg[which(kegg[,3]==droga),][,1:2]))
    }
    names(genes_per_drug) <- drugs # List of drugs, and the genes associated to each drug.
    
    message('Starting drug enrichments. Please note this process may take a while')
    drug_enrichments <- neat(genes_lists,blist = genes_per_drug,network = kegg_graph,nettype = 'undirected',nodes=drugs_available,mtc.type = 'fdr',alpha = 0.05)
    
    # Save enrichment results
    write.table(drug_enrichments,file="Supplementary_Tables/Drug_Enrichments.csv",sep="\t",col.names=T,row.names=F)
    rm(drug_enrichments)

############################################# 
# 2. Explore enrichments in Disease Data
    
    # Obtain all vertices of mondo disease layer
    diseases_available <- names(V(mondo_graph))
    
    # genes per patient appearing in the drug network
    genes_lists <- list()
    for(i in 1:nrow(genes_per_patient_and_clusters)){ 
      genes <- strsplit(as.character(genes_per_patient_and_clusters[i,4]),"_")[[1]]
      genes_lists[[i]] <- unique(genes[which(genes %in% diseases_available)])
    }
    names(genes_lists) <- genes_per_patient_and_clusters[,1] # List of genes associated to each patient appearing in the disease network
    
    # For each disease, obtain genes associated because of it
    genes_per_disease <- list()
    diseases <- unique(mondo[,3])
    for(i in 1:length(diseases)){
      disease <- diseases[i]
      genes_per_disease[[i]] <- unique(unlist(mondo[which(mondo[,3]==disease),][,1:2]))
    }
    names(genes_per_disease) <- diseases # List of disease, and the genes associated to each disease.
    
    message('Starting disease enrichments. Please note this process may take a while')
    disease_enrichments <- neat(genes_lists,blist = genes_per_disease,network = mondo_graph,nettype = 'undirected',nodes=diseases_available,mtc.type = 'fdr',alpha = 0.05)
    
    # Save enrichment results
    write.table(disease_enrichments,file="Supplementary_Tables/Disease_Enrichments.csv",sep="\t",col.names=T,row.names=F)
    rm(disease_enrichments)
   
 #############################################  
# 3. Reactome Data
    
    # Obtain all vertices of reactome net layer
    reactome_available <- names(V(reactome_graph))
    
    # genes per patient appearing in the pathway network
    genes_lists <- list()
    for(i in 1:nrow(genes_per_patient_and_clusters)){ 
      genes <- strsplit(as.character(genes_per_patient_and_clusters[i,4]),"_")[[1]]
      genes_lists[[i]] <- unique(genes[which(genes %in% reactome_available)])
    }
    names(genes_lists) <- genes_per_patient_and_clusters[,1] # List of genes associated to each patient appearing in the reactome network
    
    # For each pathway, obtain genes associated because of it
    genes_per_pathway <- list()
    pathways <- unique(reactome[,3])
    for(i in 1:length(pathways)){
      pathway <- pathways[i]
      genes_per_pathway[[i]] <- unique(unlist(reactome[which(reactome[,3]==pathway),][,1:2]))
    }
    names(genes_per_pathway) <- pathways # List of pathways, and the genes associated to each.
    
    message('Starting pathway enrichments. Please note this process may take a while')
    pathway_enrichments <- neat(genes_lists,blist = genes_per_pathway,network = reactome_graph,nettype = 'undirected',nodes=reactome_available,mtc.type = 'fdr',alpha = 0.05) #Expected run time= 17 min. Up to 16 GB of ram occupied.
    # Save enrichment results
    write.table(pathway_enrichments,file="Supplementary_Tables/Pathway_Enrichments.csv",sep="\t",col.names=T,row.names=F)
    rm(pathway_enrichments)

#############################################  
# 4. Explore enrichments in metabolic Data
    
    # Obtain all vertices of metabolome layer
    metabolites_available <- names(V(recon3d_graph))
    
    # genes per patient appearing in the metabolome network
    genes_lists <- list()
    for(i in 1:nrow(genes_per_patient_and_clusters)){ 
      genes <- strsplit(as.character(genes_per_patient_and_clusters[i,4]),"_")[[1]]
      genes_lists[[i]] <- unique(genes[which(genes %in% metabolites_available)])
    }
    names(genes_lists) <- genes_per_patient_and_clusters[,1] # List of genes associated to each patient appearing in the metabolome network
    
    # For each metabolite, obtain genes associated because of it
    genes_per_metabolite <- list()
    metabolites <- unique(recon3d[,3])
    for(i in 1:length(metabolites)){
      metabolite <- metabolites[i]
      genes_per_metabolite[[i]] <- unique(unlist(recon3d[which(recon3d[,3]==metabolite),][,1:2]))
    }
    names(genes_per_metabolite) <- metabolites # List of metabolites, and the genes associated to each metabolite.
    
    metabolite_enrichments <- neat(genes_lists,blist = genes_per_metabolite,network = recon3d_graph,nettype = 'undirected',nodes=metabolites_available,mtc.type = 'fdr',alpha = 0.05)
    
    # Save enrichment results
    write.table(metabolite_enrichments,file="Supplementary_Tables/Metabolite_enrichments.csv",sep="\t",col.names=T,row.names=F)

