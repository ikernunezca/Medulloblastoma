# Define jaccard index function

jaccard_ind <- function(x){
  res <- matrix(data=NA,nrow=ncol(x),ncol=ncol(x))
  rownames(res) <- colnames(x)
  colnames(res) <- colnames(x)
  for(i in 1:ncol(x)){
    uno <- rownames(x)[which(x[,i]==1)]
    for(j in 1:ncol(x)){
      if(is.na(res[i,j])){
        dos <- rownames(x)[which(x[,j]==1)]
        numerador <- length(intersect(uno,dos))
        denom <- length(union(uno,dos))
        out <- numerador/denom
        res[i,j] <- out
        res[j,i] <- res[i,j]
      }
    }
  }
  attr(res,"method") <- "jaccard_ind" 
  return(as.dist(1-res)) ###Also: We have to be calculating a dissimilarity matrix with the function, so 1-x.
}

# Load dependencies
library(sigclust2)
library(pvclust)
library(fpc)
library(jaccard)
library(CmmD)
library(readr)
library(knitr)
library(data.table)

# Create a vector with the paths where Molti's Output files are saved. 
structures_12 <- paste0("data/Molti_Output/",seq(0.5,12,0.5),".csv")
# Detect community trajectories and tree distances between each gene. 
curie_to_12_full <- CmmD_from_community_structures(nodelist = NULL, community_structures = structures_12, resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = "hamming",threads = 7)
curie_to_12_full$hamming_distance_matrix = curie_to_12_full$distance_matrix * 24 # This transformation is needed because parallel dist is weighted.

# Load genes associated to each patient from CURIE data
tata <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/iPC-project-H2020/ipcrg/master/scripts/CURIE2gr/multi.layer.net.gr",sep = "\t",header = F, stringsAsFactors=F))
splited_patients <- split(tata[,2],tata[1])

# Generate ground truth table from CURIE data
ground_truth_patients <- c("MB30","MB31","MB34",
                           "MB04","MB05","MB06","MB24","MB40","MB25","MB43","MB46","MB49","MB55",
                           "MB01","MB02","MB03","MB14","MB19","MB47","MB50","MB51","MB52","MB53",
                           "MB07","MB08","MB09","MB13","MB15","MB16","MB17","MB20","MB22","MB39","MB48","MB54")
ground_truth <- matrix(nrow= length(ground_truth_patients),ncol=2)
ground_truth[,1] <- ground_truth_patients
ground_truth[,2] <- c(rep("WNT",3),rep("SHH",10),rep("G3",10),rep("G4",12))
colnames(ground_truth) <- c("Patient","Real_class")
rownames(ground_truth) <- ground_truth[,1]


message("Performing filtering based on tetha (0 to 10)")
genes_per_patient_list <- list()
for(k in 0:10){
  genes_per_patient <- list()
  for(i in 1:length(splited_patients)){
    genes_interesantes <- names(which(table(splited_patients[[i]])>=1))
    
    # Filter CURIE genes to those included in the multilayer network.
    genes_interesantes <- genes_interesantes[genes_interesantes %in% rownames(curie_to_12_full$hamming_distance_matrix)]
    matricin <- curie_to_12_full$hamming_distance_matrix[genes_interesantes,genes_interesantes]
    
    # For each value of theta between 0 and 10, we create a vector (suc), where for each gene of the multilayer we calculate the
    # number of other multilayer genes that are below a maximum k value of theta -distance in the tree-.
    names_matricin <- colnames(matricin)
    suc <- vector("numeric",length = length(names_matricin))
    names(suc) <- names_matricin
    for(j in 1:ncol(matricin)){
      leng_matricin <- length(matricin[,j][matricin[,j]<=k]) ### k (in our case, a representation of theta) is maximum distance allowed in the clustering plot. (Sauron's eye)
      suc[j] <- leng_matricin
    }
    suc <- suc[suc>1] # Filter genes that have no patient associated partners
    genes_per_patient[[i]] <- suc - 1 # As a value of 2 mean that only there is one more gene with the gene being analyzed, we substract 1.
  }
  names(genes_per_patient) <- names(splited_patients)
  genes_per_patient_list[[k+1]] <- genes_per_patient
  names(genes_per_patient_list)[[k+1]] <- as.character(k)
  message(paste0("theta= ",k))
}

message("Tetha based filtering finished. Calculating clustering accuracies for tetha 0 to 10 and lambda 1 to 20")

# Start a 11 x 20 matrix to be filled with the hierarchical clustering accuracy values.
final_accuracy_matrix <- matrix(0, ncol= 20, nrow= 11)
final_kk_used <- matrix(0, ncol= 20, nrow= 11)
rownames(final_accuracy_matrix) <- as.character(0:10)
rownames(final_kk_used) <- as.character(0:10)

u <- 1 # theta = 0
val <- 6 #lambda
    preserve_genes_per_patient <- genes_per_patient_list[[u]] # u==k == tetha + 1
    genes_per_patient <- preserve_genes_per_patient
    genes_per_patient <- lapply(genes_per_patient,function(x) x[x<=val]) #val = lambda. We filter the genes that are over the lambda value tested
    
    genes_per_patient_names <- lapply(genes_per_patient,function(x) names(x))
    all_genes_possible <- unique(unlist(genes_per_patient_names,use.names=F))
    
    # Generate a 0-1 patient x genes matrix that acts as input for pamk, jaccard_ind and hclust.
    n_genes_p_patients <- matrix(data= 0, nrow= 38,ncol= length(all_genes_possible))
    colnames(n_genes_p_patients) <- all_genes_possible
    rownames(n_genes_p_patients) <- names(genes_per_patient)
    
    for(rowi in 1:nrow(n_genes_p_patients)){
      n_genes_p_patients[rowi,] <- as.integer(colnames(n_genes_p_patients) %in% genes_per_patient_names[[rowi]])
    }
    
    # Save Table 1
    # values <- as.matrix(jaccard_ind(t(patient_matrix)))
    # table_2 <- 1-values[,c("MB10","MB21","MB33")]
    # write.table(table_2,file="Tables/Table_1.csv", row.names=T, col.names=T, sep="\t")
    
    # Exclude validation set
    WHATEVER <- c("MB10","MB21","MB33")
    patient_matrix <- n_genes_p_patients
    
    patient_matrix2 <- patient_matrix[- which(rownames(patient_matrix) %in% WHATEVER),] # Exclude patients with missing data from clustering
    
    #Obtain optimal clusters
    pamk.best <- pamk(patient_matrix2)
    kk <- pamk.best$nc ## kk is the optimal number of clusters for the particular optimization.
    
    # Perform hierarchical clustering with the suggested number of clusters
    patient_matrix3 <- t(patient_matrix2)
    set.seed(2020)
    
    # Pvclust Bootstrap
    res_pvclust <- pvclust(patient_matrix3, method.dist=jaccard_ind,method.hclust = "ward.D2", nboot=10000,store = TRUE,parallel = TRUE) 
    # Add 2 decimal positions for plotting
    res_pvclust_100 <- res_pvclust
    res_pvclust_100$edges <- res_pvclust$edges * 100
    # shc Bootstrap
    vfun <- function(x, y) {1 - jaccard(x, y)}
    mfun <- function(x) {
      as.dist(outer(split(x, f = row(x)), split(x, f = row(x)),
                    Vectorize(vfun)))
    }
    
    res_shc <- shc(patient_matrix2, matmet=mfun, linkage="ward.D2", n_sim = 10000)
    
    # Original hclust output
    res_hclust <- hclust(jaccard_ind(patient_matrix3),"ward.D2")
    
# Plotting

pdf("data/Plots/clustering_pv_shc_hclust.pdf",height= 14, width= 20)
plot(res_hclust,main= "Ward Hierarchical clustering of Medulloblastoma patients")
plot(res_pvclust_100,main= "Ward Hierarchical Clustering Bootstrap AU p−values (%) (package pvclust)",print.pv= "au")
# Plot empirical values
plot(res_shc,alpha=0.5,ci_emp=T) # we change alpha in order to plot the values. With predefined alpha not all partition values would be plotted
plot(res_shc,alpha=0.5,ci_emp=F) # Plot Gaussian aproximate p-values
dev.off()

