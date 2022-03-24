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
        if(is.na(out)){
          out <- 0
        }
        res[i,j] <- out
        res[j,i] <- res[i,j]
      }
    }
  }
  attr(res,"method") <- "jaccard_ind" 
  return(as.dist(1-res)) # In order to calculate a dissimilarity matrix.
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
# 24 = length(seq(0.5,12,0.5)) -> number of resolution values analyzed

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

# Define variables to fill at each iteration
counts <- 0
genes_per_patient_X <- vector(mode = "list", length = length(splited_patients)) # start empty gene list
thetas <- c() # vector for best tetha values from each iteration (value corresponds to theta + 1)
lambdas <- c() # vector for lambda values from each iteration
accuracies <- c() # vector of best accuracy values from each iteration
accuracy_matrices <- list() # list with the accuracy matrices from each iteration.
n_erased_genes <- c() # mean number of genes to be erased in the next iteration. n_erased_genes[i] = X genes to be erased in the iteration i+1. 
best_genes_per_iteration <- list() # A list with the genes per patient that generated the best accuracy

repeat{

  counts <- counts + 1

  message(paste("Starting iteration",counts,"at",Sys.time()))

  message("Performing filtering based on tetha (0 to 10)")

  genes_per_patient_list <- list()

  for(k in 0:10){ # k ---> Theta
    genes_per_patient <- list()
    for(i in 1:length(splited_patients)){
     genes_interesantes <- names(which(table(splited_patients[[i]])>=1))
    
      # Filter CURIE genes to those included in the multilayer network.
      genes_interesantes <- genes_interesantes[genes_interesantes %in% rownames(curie_to_12_full$hamming_distance_matrix)]  
    
      # Filter genes that yielded the best accuracy at the previous iteration. If conditional is needed in order not to erase all genes.
      if(length(which(genes_interesantes %in% genes_per_patient_X[[i]])>0)){
       genes_interesantes <- genes_interesantes[-which(genes_interesantes %in% genes_per_patient_X[[i]])]
      }
      matricin <- curie_to_12_full$hamming_distance_matrix[genes_interesantes,genes_interesantes]
    
      # For each value of theta between 0 and 10, we create a vector (suc), where for each gene of the multilayer we calculate the
      # number of other multilayer genes that are below a maximum k value of theta -distance in the tree-.
      names_matricin <- colnames(matricin)
      suc <- vector("numeric",length = length(names_matricin))
      names(suc) <- names_matricin
      for(j in 1:ncol(matricin)){
        leng_matricin <- length(matricin[,j][matricin[,j]<=k])
        suc[j] <- leng_matricin
      }
      suc <- suc[suc>1] # Filter genes that have no patient associated partners.
      genes_per_patient[[i]] <- suc - 1 #  As a value of 2 mean that only there is one more gene with the gene being analyzed, we substract 1.
    }
    names(genes_per_patient) <- names(splited_patients)
    genes_per_patient_list[[k+1]] <- genes_per_patient
    names(genes_per_patient_list)[[k+1]] <- as.character(k)
   message(paste0("tetha=",k))
  }

  message("Tetha based filtering finished. Calculating clustering accuracies for tetha 0 to 10 and lambda 1 to 20")

  # Start a 11 x 20 matrix to be filled with the hierarchical clustering accuracy values.
  final_accuracy_matrix <- matrix(0, ncol= 20, nrow= 11)
  final_kk_used <- matrix(0, ncol= 20, nrow= 11)
  rownames(final_accuracy_matrix) <- as.character(0:10)
  rownames(final_kk_used) <- as.character(0:10)

  for(u in 1:11){ # 1 to 11 because of non 0-based language: If we want, for example, theta to be 0, we set u=1
    for(val in 1:20){ # lambda 
      preserve_genes_per_patient <- genes_per_patient_list[[u]] # u==k == tetha + 1
    
      genes_per_patient <- preserve_genes_per_patient
      genes_per_patient <- lapply(genes_per_patient,function(x) x[x<=val]) #val = lambda. We filter the genes that are over the lambda value tested
    
      # We calculate the accuracy if at least one patient has genes associated. At high iterations we may have excluded too much genes and therefore our accuracy is set to 0 by default as there is genes to compare.
      if(length(unlist(genes_per_patient))>0){ 
        genes_per_patient_names <- lapply(genes_per_patient,function(x) names(x))
        all_genes_possible <- unique(unlist(genes_per_patient_names,use.names=F))
      
        # Generate a 0-1 patient x genes matrix that acts as input for pamk, jaccard_ind and hclust.
        n_genes_p_patients <- matrix(data= 0, nrow= 38,ncol= length(all_genes_possible))
        colnames(n_genes_p_patients) <- all_genes_possible
        rownames(n_genes_p_patients) <- names(genes_per_patient)
      
        for(rowi in 1:nrow(n_genes_p_patients)){
          n_genes_p_patients[rowi,] <- as.integer(colnames(n_genes_p_patients) %in% genes_per_patient_names[[rowi]])
        }
      
        WHATEVER <- c("MB10","MB21","MB33") 
        patient_matrix <- n_genes_p_patients
      
        patient_matrix2 <- patient_matrix[- which(rownames(patient_matrix) %in% WHATEVER),] # Exclude patients with missing data from clustering
      
        # Obtain optimal clusters
        pamk.best <- pamk(patient_matrix2)
        kk <- pamk.best$nc ## kk is the optimal number of clusters for the teration.
      
        # Perform hierarchical clustering with the suggested number of clusters
        patient_matrix3 <- t(patient_matrix2)
        set.seed(2020)
        res_hclust <- hclust(jaccard_ind(patient_matrix3),"ward.D2")
      
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
      
        sum_matrix <- arbol_splited_mat + ground_truth_mat #sum both matrices
        tab_sum_matrix <- table(sum_matrix)
        zeros <- tab_sum_matrix["0"] # True Negatives
        if(is.na(zeros)){
          zeros <- 0
        }
        twos <- tab_sum_matrix["2"] #True positives
        if(is.na(twos)){
          twos <- 0
        }
        accuracy <- (zeros+twos)/sum(tab_sum_matrix,na.rm = T)
        final_accuracy_matrix[u,val] <- accuracy
        final_kk_used[u,val] <- kk
      }
    }
    message(paste0("Acuracies for tetha=",u-1," calculated."))
  }

  message("All acuracies calculated. Obtaining best accuracy and theta-lambda values.")

  # Get best accuracy and the corresponding theta-lambda pair
  maximized <- which(final_accuracy_matrix==max(final_accuracy_matrix),arr.ind = T)
  if(length(maximized)==2){
    best_theta <- maximized[1]
    best_lambda <- maximized[2]
  }
  if(length(maximized)>2){ #If various cases present the best accuracy, we choose the lowest value if theta possible
    maxi <- maximized[which(maximized[,1]==min(maximized[,1]))[1],]
    best_theta <- maxi[1]
    best_lambda <- maxi[2]
  }

  u <- unname(best_theta)
  val <- unname(best_lambda)

  # Get genes to filter in the next iteration
  preserve_genes_per_patient_Y <- genes_per_patient_list[[u]] ### u==k == tetha + 1

  genes_per_patient_Y <- preserve_genes_per_patient_Y
  genes_per_patient_Y <- lapply(genes_per_patient_Y,function(x) x[x<=val])
  genes_per_patient_Y <- lapply(genes_per_patient_Y,function(x) names(x)) 

  message("Calculated. Saving data of iteration and joining genes to exclude in next iteration...")

  # Fill output variables
  best_genes_per_iteration <- c(best_genes_per_iteration,genes_per_patient_Y)
  accuracies <- c(accuracies,max(final_accuracy_matrix))
  accuracy_matrices <- c(accuracy_matrices,final_accuracy_matrix)
  lambdas <- c(lambdas,best_lambda)
  thetas <- c(thetas, best_theta)


  genes_per_patient_X <- mapply(c, genes_per_patient_X, genes_per_patient_Y) #join genes from previous iterations for the be filtered in the next one
  genes_per_patient_X <- lapply(genes_per_patient_X, function(x) unique(x))

  mean_tamanos <- mean(unlist(lapply(genes_per_patient_X,function(x) length(x))))
  n_erased_genes <- c(n_erased_genes,mean_tamanos)

  if(max(final_accuracy_matrix) == 0){
    message(paste("Process stopped at iteration",counts))
    break()
  }

  finish_time <- Sys.time()

  message(paste("Iteration",counts,"finished at",finish_time,". Accuracy=",max(final_accuracy_matrix),", Theta= ",best_theta-1," Lambda= ",best_lambda))

} #repeat{}


# Generate data matrix for the plot (Supplementary_Figure_3)

thetas <- thetas - 1
n_erased_genes <- n_erased_genes[-(length(n_erased_genes))]
n_erased_genes <- c(0,n_erased_genes)
matrix_for_plot <- matrix(ncol= length(lambdas),nrow=4) 
matrix_for_plot[1,] <- accuracies
matrix_for_plot[2,] <- thetas
matrix_for_plot[3,] <- lambdas
matrix_for_plot[4,] <- n_erased_genes

# Plotting

ngenes_p_patient_allgenes <- matrix_for_plot

data <- ngenes_p_patient_allgenes[1:3,]
data[1,] <- data[1,]*100

data <- as.matrix(data)
colnames(data) <- as.character(round(as.matrix(ngenes_p_patient_allgenes[4,]),digits = 2))
colnames(data)[1] <- "0"

barplot(data,beside = T)
png("data/Plots/Supp_Figure_7.png",width= 3250,height= 850)
xx <- barplot(data,beside = T,col = c("turquoise","violet","red"),xlab = "Mean removed genes per patient",ylab="value",ylim=c(0,120))
text(x = xx, y = data, label = round(data,2), pos = 3, cex = 0.8, col = "red")
dev.off()

