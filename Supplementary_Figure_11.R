
# Load dependencies
library(sigclust2)
library(pvclust)
library(fpc)
library(jaccard)
library(CmmD)
library(readr)
library(knitr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
# Mini test of proteome data from Archer et al 2018



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

setwd('~/Escritorio/Medulloblastoma/')
# Create a vector with the paths where Molti's Output files are saved. 
structures_12 <- paste0("data/Molti_Output/",seq(0.5,12,0.5),".csv")
# Detect community trajectories and tree distances between each gene. 
curie_to_12_full <- CmmD_from_community_structures(nodelist = NULL, community_structures = structures_12, resolution_start = 0.5,resolution_end = 12,interval = 0.5)
curie_to_12_full$hamming_distance_matrix = curie_to_12_full$distance_matrix * 24 # This transformation is needed because parallel dist is weighted.
# 24 = length(seq(0.5,12,0.5)) -> number of resolution values analyzed


# Load patient data from Archer et al. 2018

tata <- as.data.frame(read.csv(file = "data/multilayer_archer.gr",sep = "\t",header = F, stringsAsFactors=F))
splited_patients <- split(tata[,2],tata[1])

# Generate ground truth table from Archer et al 2018. They suggest 6 groups.
ground_truth_patients <- c("x1M6","x5M15","Burdenko_1360",
                           "MB037","MB244","MB268","MB274","MB275","MB102","MB104","MB234","MB239","MB284",
                           "MB136","MB266","MB287","MB088","MB206",
                           "MB018","MB095","MB106","MB170","MB226","MB247","MB248","MB260",
                           "MB277","MB164","MB271","MB166","MB278","MB288",
                           "MB177","MB199","MB281","MB174","MB227","MB265","MB270","MB264","MB282","MB099","MB091","MB118","MB269")
ground_truth <- matrix(nrow= length(ground_truth_patients),ncol=2)
ground_truth[,1] <- ground_truth_patients
ground_truth[,2] <- c(rep("WNT",3),rep("SHHa",10),rep("SHHb",5),rep("G3a",8),rep("G3b",6),rep("G4",13))
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
    suc <- suc[suc>1] ###Filetr genes that have no patient associated partners
    genes_per_patient[[i]] <- suc - 1 ### As a value of 2 mean that only there is one more gene with the gene being analyzed, we substract 1.
  }
  names(genes_per_patient) <- names(splited_patients)
  genes_per_patient_list[[k+1]] <- genes_per_patient
  names(genes_per_patient_list)[[k+1]] <- as.character(k)
  message(paste0("theta= ",k))
}

message("Tetha based filtering finished. Calculating clustering accuracies for tetha 0 to 10 and lambda 1 to 20")

# Start a 11 x 20 matrix to be filled with the hierarchical clustering accuracy values.
final_accuracy_matrix <- matrix(0, ncol= 20, nrow= 11)
final_matthews_matrix <- matrix(0, ncol= 20, nrow= 11)
final_kk_used <- matrix(0, ncol= 20, nrow= 11)
mean_per_pair <- matrix(0, ncol= 20, nrow= 11)
rownames(final_accuracy_matrix) <- as.character(0:10)
rownames(final_kk_used) <- as.character(0:10)

for(u in 1:11){ # 1 to 11 because of non 0-based language: If we want, for example, theta to be 0, we set u=1
  for(val in 1:20){ #lambda
    preserve_genes_per_patient <- genes_per_patient_list[[u]] # u==k == tetha + 1
    genes_per_patient <- preserve_genes_per_patient
    genes_per_patient <- lapply(genes_per_patient,function(x) x[x<=val]) #val = lambda. We filter the genes that are over the lambda value tested
    
    genes_per_patient_names <- lapply(genes_per_patient,function(x) unique(names(x)))
    all_genes_possible <- unique(unlist(genes_per_patient_names,use.names=F))
    
    # Generate a 0-1 patient x genes matrix that acts as input for pamk, jaccard_ind and hclust.
    n_genes_p_patients <- matrix(data= 0, nrow= 45,ncol= length(all_genes_possible))
    colnames(n_genes_p_patients) <- all_genes_possible
    rownames(n_genes_p_patients) <- names(genes_per_patient)
    
    for(rowi in 1:nrow(n_genes_p_patients)){
      n_genes_p_patients[rowi,] <- as.integer(colnames(n_genes_p_patients) %in% genes_per_patient_names[[rowi]])
    }
    
    # WHATEVER <- c("MB10","MB21","MB33")
    patient_matrix <- n_genes_p_patients
    
    patient_matrix2 <- patient_matrix # Exclude patients with missing data from clustering
    
    #Get mean gene length
    all_lengths <- unlist(lapply(genes_per_patient_names,function(x) length(x)))
    media <- mean(all_lengths)
    
    #Obtain optimal clusters
    pamk.best <- pamk(patient_matrix2)
    kk <- pamk.best$nc ## kk is the optimal number of clusters for the particular optimization.
    
    # Perform hierarchical clustering with the suggested number of clusters
    patient_matrix3 <- t(patient_matrix2)
    set.seed(2020)
    res_hclust <- hclust(jaccard_ind(patient_matrix3),"ward.D2")
    
    # Save Best fitted pair
      if(u == 2){
        if(val==5){
          pdf("data/Plots/Supplementary_Figure_11.pdf",height= 14, width= 20)
          plot(res_hclust,main= "Ward Hierarchical clustering of Archer et al. 2018 Medulloblastoma patients")
          dev.off()
          }}
    
    # Calculate two 0-1 matrices in order to compare our clustering with the ground truth.
    arbol <- cutree(res_hclust,kk)
    arbol_splited <- split(names(arbol),arbol)
    splited_ground_truth <- split(ground_truth[,1],ground_truth[,2])
    arbol_splited_mat <- matrix(0,ncol= nrow(ground_truth),nrow= nrow(ground_truth))
    ground_truth_mat <- matrix(0,ncol= nrow(ground_truth),nrow= nrow(ground_truth))
    dimnames(arbol_splited_mat) <- list(rownames(ground_truth),rownames(ground_truth))
    dimnames(ground_truth_mat) <- list(rownames(ground_truth),rownames(ground_truth))
    for(f in 1:nrow(arbol_splited_mat)){
      current_patient_row <- rownames(ground_truth_mat)[f]
      for(g in 1:ncol(arbol_splited_mat)){
        current_patient_col <- colnames(ground_truth_mat)[g]
        cluster_pat_row_ground_truth <- grep(current_patient_row,splited_ground_truth)
        cluster_pat_col_ground_truth <- grep(current_patient_col,splited_ground_truth)
        cluster_pat_row_arbol_splited <- grep(current_patient_row,arbol_splited)
        cluster_pat_col_arbol_splited <- grep(current_patient_col,arbol_splited)
        # print(paste('f=',f,'g=',g,'length of row_gt=',length(cluster_pat_row_ground_truth),'length of col_gt=',length(cluster_pat_col_ground_truth),'length of row_ab=',length(cluster_pat_row_arbol_splited),'length of col_ab=',length(cluster_pat_col_arbol_splited)))
        if(cluster_pat_row_ground_truth==cluster_pat_col_ground_truth){
          ground_truth_mat[f,g] <- 1
          ground_truth_mat[g,f] <- ground_truth_mat[f,g]
        }
        if(cluster_pat_row_arbol_splited==cluster_pat_col_arbol_splited){
          arbol_splited_mat[f,g] <- 1
          arbol_splited_mat[g,f] <- arbol_splited_mat[f,g]
        }
      }
    }
    
    sum_matrix <- arbol_splited_mat + ground_truth_mat
    tab_sum_matrix <- table(sum_matrix)
    zeros <- tab_sum_matrix["0"] # True Negatives
    if(is.na(zeros)){
      zeros <- 0
    }
    twos <- tab_sum_matrix["2"] # True Positives
    if(is.na(twos)){
      twos <- 0
    }
    # Accuracy
    accuracy <- (zeros+twos)/sum(tab_sum_matrix,na.rm = T)
    final_accuracy_matrix[u,val] <- accuracy
    # Matthew's Coefficient
    sum_matrix[which(sum_matrix==1,arr.ind=T)] <- 6 #Set all false to value 6.
    dif_sum_mat <- sum_matrix - arbol_splited_mat
    dif_tab_sum_matrix <- table(dif_sum_mat)
    sixes <- dif_tab_sum_matrix["6"] # False Positives
    if(is.na(sixes)){
      sixes <- 0
    }
    fives <- dif_tab_sum_matrix["5"] # False Negatives
    if(is.na(fives)){
      fives <- 0
    }
    tn <- as.double(unname(zeros))
    tp <- as.double(unname(twos))
    fp <- as.double(unname(fives))
    fn <- as.double(unname(sixes))
    numerador <- (tp * tn) - (fp * fn)
    den_1 <- tp + fp
    den_2 <- tp + fn
    den_3 <- tn + fp
    den_4 <- tn + fn
    pro_den <- den_1 * den_2 * den_3 * den_4
    denominador <- sqrt(pro_den)
    matthews <- numerador/denominador
    final_matthews_matrix[u,val] <- matthews
    final_kk_used[u,val] <- kk
    mean_per_pair[u,val] <- media
  }
  message(paste0("Acuracies for tetha=",u-1," calculated."))
}

write.table(final_accuracy_matrix,row.names=T,col.names=F,file="Supplementary_Tables/Supplementary_Table_10.tsv") # Supplementary 10
write.table(final_kk_used,row.names=T,col.names=F,file="Supplementary_Tables/Supplementary_Table_11.tsv") # Supplementary 11
write.table(final_matthews_matrix,row.names=T,col.names=F,file="Supplementary_Tables/Supplementary_Table_12.tsv") # Supplementary 12
