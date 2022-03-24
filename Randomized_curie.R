# The following script is intended to be run at clusters running SLURM and Greasy. It is programed to be run from bash in the following way:
# RScript Randomized_curie.R 1 50 
# With such command, the script will perform randomizations 1 to 50 and save 2 files, the accuracy matrix and the number of clusters suggested for each theta-lambda pair.

set.seed(NULL) # Restart all seeds (Pure randomizing)

# Load Dependencies
library(pvclust)
library(fpc)
library(readr)
library(knitr)
library("AnnotationDbi")
library("igraph")
library("stringr")
library("e1071")

# Define Jaccard index function to use with hclust
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
  return(as.dist(1-res))
}

# Other useful function
unsplit_to_data_frame <- function(x){
  unlisted <- unlist(x,use.names=T)
  col_1 <- unlist(x,use.names=F)
  col_2 <- substr(names(unlisted),start=1,stop=4)
  tabla <- as.data.frame(matrix(nrow=length(unlisted),ncol=2))
  tabla[,1] <- col_1
  tabla[,2] <- col_2
  return(tabla)
}

# Define CmmD Package functions (Not installed at our supercomputer)
CmmD <- function (nodelist = NULL, input_layers, resolution_start, resolution_end,
          interval, destfile_community_analysis)
{
  require("AnnotationDbi")
  require("igraph")
  require("stringr")
  require("e1071")
  if (length(input_layers) < 2) {
    stop("ERROR: Input_layers argument must be a list of at least 2 network files")
  }
  if (class(resolution_end) != "numeric") {
    stop("ERROR: Resolution parameter must be a number")
  }
  if (class(resolution_start) != "numeric") {
    stop("ERROR: Resolution parameter must be a number")
  }
  if (class(interval) != "numeric") {
    stop("ERROR: Interval value must be a number")
  }
  if (class(destfile_community_analysis) != "character") {
    stop("ERROR: destfile_community_analysis expects a character string")
  }
  layers <- paste0(input_layers, collapse = " ")
  message(paste0("Resolution parameter starts at: ", resolution_start))
  message(paste0("Resolution parameter ends at: ", resolution_end))
  resolution_interval <- seq(from = resolution_start, to = resolution_end,
                             by = interval)
  desfile_vector <- paste0(destfile_community_analysis, resolution_interval,
                           ".csv")
  message(paste0("Starting community analysis."))
  start_time <- Sys.time()
  for (i in 1:length(resolution_interval)) {
    current_resolution <- resolution_interval[i]
    current_destfile <- desfile_vector[i]
    current_layers <- layers
    message(paste0("Resolution parameter: ", current_resolution))
    message(Sys.time())
    system_order <- paste("molti-console", "-o", current_destfile,
                          "-p", current_resolution, layers)
    system(system_order)
  }
  message(paste0("Reading MolTi output files. Calculating Gene/Community matrix"))
  output_files <- list.files(destfile_community_analysis)
  to_be_forgotten <- grep("_", output_files)
  output_files <- output_files[-to_be_forgotten]
  alllists <- list()
  for (i in 1:length(output_files)) {
    red <- readLines(paste0(destfile_community_analysis,
                            output_files[i]))
    cluster_ids <- grep("Cluster", red)
    lista <- list()
    for (j in 1:length(cluster_ids)) {
      st <- cluster_ids[j]
      if (j == length(cluster_ids)) {
        en <- length(red)
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-length(current_cluster)]
      }
      else {
        en <- cluster_ids[j + 1]
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-c(length(current_cluster),
                                               length(current_cluster) - 1)]
      }
      lista[[j]] <- current_cluster2[2:length(current_cluster2)]
      names(lista)[j] <- paste0("Cluster_", j)
    }
    kaz <- output_files[i]
    assign(paste0("com_", kaz), value = lista)
    alllists[[i]] <- lista
  }
  names(alllists) <- output_files
  tamano_alllists <- length(alllists)
  allgenes <- unique(unlist(alllists))
  if (length(nodelist) > 0) {
    inter_nodes <- intersect(allgenes, nodelist)
    allgenes <- inter_nodes
  }
  print(paste0("Files red. Calculating Gene/Community matrix"))
  res_matrix <- matrix(ncol = tamano_alllists + 1, nrow = length(allgenes))
  rownames(res_matrix) <- allgenes
  colnames(res_matrix) <- c(output_files, "Pattern")
  for (i in 1:length(allgenes)) {
    gen <- rownames(res_matrix)[i]
    for (j in 1:tamano_alllists) {
      searched <- unlist(lapply(alllists[[j]], function(x) gen %in%
                                  x))
      comunidad <- unname(which(searched == TRUE))
      res_matrix[i, j] <- comunidad
    }
    res_matrix[i, "Pattern"] <- paste0(res_matrix[i, 1:(ncol(res_matrix) -
                                                          1)], collapse = "_")
    percentage <- round((i/length(allgenes)), digits = 4) *
      100
    porcentajes <- seq(from = 0, to = 100, by = 5)
    progresos <- paste0("Progress: ", porcentajes, "%")
    porcentajes[1] <- 1
    if ((percentage %in% porcentajes) == TRUE) {
      cual_percentage <- which(porcentajes == percentage)
      to_post <- paste0("Progress: ", percentage, "%")
      message(to_post)
    }
  }
  message(paste0("Gene/Community matrix calculated, calculating Hamming distances for all gene pairs. This process may take a while: It takes about 14 min with an Intel Xeon E-2124 processor"))
  genes_same_communities <- split(rownames(res_matrix), res_matrix[,
                                                                   "Pattern"])
  final_res_matrix_length <- ncol(res_matrix) - 1
  distance_matrix <- hamming.distance(res_matrix[, 1:final_res_matrix_length])
  final_output <- list(res_matrix[, 1:final_res_matrix_length],
                       genes_same_communities, distance_matrix)
  names(final_output) <- c("gene_community_matrix", "l_constant",
                           "hamming_distance_matrix")
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  message(paste0("Run Time: ", diff_time))
  return(final_output)
}

CmmD_from_community_structures <- function (nodelist = NULL, community_structures, resolution_start,
          resolution_end, interval)
{
  start_time <- Sys.time()
  require("AnnotationDbi")
  require("igraph")
  require("stringr")
  require("e1071")
  if (class(resolution_end) != "numeric") {
    stop("ERROR: Resolution parameter must be a number")
  }
  if (class(resolution_start) != "numeric") {
    stop("ERROR: Resolution parameter must be a number")
  }
  if (class(interval) != "numeric") {
    stop("ERROR: Interval value must be a number")
  }
  message(paste0("Resolution parameter starts at: ", resolution_start))
  message(paste0("Resolution parameter ends at: ", resolution_end))
  resolution_interval <- seq(from = resolution_start, to = resolution_end,
                             by = interval)
  message(paste0("Reading MolTi output files. Calculating Gene/Community matrix"))
  output_files <- community_structures
  alllists <- list()
  for (i in 1:length(output_files)) {
    red <- readLines(output_files[i])
    cluster_ids <- grep("Cluster", red)
    lista <- list()
    for (j in 1:length(cluster_ids)) {
      st <- cluster_ids[j]
      if (j == length(cluster_ids)) {
        en <- length(red)
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-length(current_cluster)]
      }
      else {
        en <- cluster_ids[j + 1]
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-c(length(current_cluster),
                                               length(current_cluster) - 1)]
      }
      lista[[j]] <- current_cluster2[2:length(current_cluster2)]
      names(lista)[j] <- paste0("Cluster_", j)
    }
    kaz <- output_files[i]
    assign(paste0("com_", kaz), value = lista)
    alllists[[i]] <- lista
  }
  names(alllists) <- output_files
  tamano_alllists <- length(alllists)
  allgenes <- unique(unlist(alllists))
  if (length(nodelist) > 0) {
    inter_nodes <- intersect(allgenes, nodelist)
    allgenes <- inter_nodes
  }
  message(paste0("Files red. Calculating Gene/Community matrix"))
  res_matrix <- matrix(ncol = tamano_alllists + 1, nrow = length(allgenes))
  rownames(res_matrix) <- allgenes
  colnames(res_matrix) <- c(output_files, "Pattern")
  for (i in 1:length(allgenes)) {
    gen <- rownames(res_matrix)[i]
    for (j in 1:tamano_alllists) {
      searched <- unlist(lapply(alllists[[j]], function(x) gen %in%
                                  x))
      comunidad <- unname(which(searched == TRUE))
      res_matrix[i, j] <- comunidad
    }
    res_matrix[i, "Pattern"] <- paste0(res_matrix[i, 1:(ncol(res_matrix) -
                                                          1)], collapse = "_")
    percentage <- round((i/length(allgenes)), digits = 4) *
      100
    porcentajes <- seq(from = 0, to = 100, by = 5)
    progresos <- paste0("Progress: ", porcentajes, "%")
    porcentajes[1] <- 1
    if ((percentage %in% porcentajes) == TRUE) {
      cual_percentage <- which(porcentajes == percentage)
      to_post <- paste0("Progress: ", percentage, "%")
      message(to_post)
    }
  }
  message(paste0("Gene/Community matrix calculated, calculating Hamming distances for all gene pairs. This process may take a while: It takes about 14 min with an Intel Xeon E-2124 processor"))
  genes_same_communities <- split(rownames(res_matrix), res_matrix[,
                                                                   "Pattern"])
  final_res_matrix_length <- ncol(res_matrix) - 1
  distance_matrix <- hamming.distance(res_matrix[, 1:final_res_matrix_length])
  final_output <- list(res_matrix[, 1:final_res_matrix_length],
                       genes_same_communities, distance_matrix)
  names(final_output) <- c("gene_community_matrix", "l_constant",
                           "hamming_distance_matrix")
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  message(paste0("Run Time: ", diff_time))
  return(final_output)
}

# Define two variables to set the interval of randomizations to do. It only affects the naming of the final files. 
as <- commandArgs()[6]
bs <- commandArgs()[7]

en_cual_guardo <- as:bs

# Create a vector with the paths where Molti's Output files are saved.
structures_12 <- paste0("data/Molti_Output/",seq(0.5,12,0.5),".csv")
# Detect community trajectories and tree distances between each gene. 
curie_to_12_full <- CmmD_from_community_structures(nodelist = NULL, community_structures = structures_12, resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = "hamming",threads = 7)
curie_to_12_full$hamming_distance_matrix = curie_to_12_full$distance_matrix * 24 # This transformation is needed because parallel dist is weighted.
# 24 = length(seq(0.5,12,0.5)) -> number of resolution values analyzed

# Load genes associated to each patient from CURIE data
tata <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/iPC-project-H2020/ipcrg/master/scripts/CURIE2gr/multi.layer.net.gr",sep = "\t",header = F, stringsAsFactors=F))
splited_patients <- split(tata[,2],tata[1])
all_genes <- unique(tata[,2])
splited_patients_original <- split(tata[,2],tata[1])

for(p in 1:50){ # Perform 50 randomizations
  splited_patients <- list()
  for(i in 1:length(splited_patients_original)){
    splited_patients[[i]] <- sample(all_genes,size = length(splited_patients_original[[i]]),replace = T ) # Get a random sample of genes that is of the same size as the number of originally associated genes
  }
  
  # Generate ground truth table fro CURIE data
  names(splited_patients) <- names(splited_patients_original)
  ground_truth_patients <- c("MB30","MB31","MB34",
                             "MB04","MB05","MB06","MB24","MB40","MB25","MB43","MB46","MB49","MB55",
                             "MB01","MB02","MB03","MB14","MB19","MB47","MB50","MB51","MB52","MB53",
                             "MB07","MB08","MB09","MB13","MB15","MB16","MB17","MB20","MB22","MB39","MB48","MB54")
  ground_truth <- matrix(nrow= length(ground_truth_patients),ncol=2)
  ground_truth[,1] <- ground_truth_patients
  ground_truth[,2] <- c(rep("WNT",3),rep("SHH",10),rep("G3",10),rep("G4",12))
  colnames(ground_truth) <- c("Patient","Real_class")
  rownames(ground_truth) <- ground_truth[,1]


  genes_per_patient_list <- list()
  for(k in 0:10){ # k ---> Theta
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
      suc <- suc[suc>1] # Filter genes that have no patient associated partners.
      genes_per_patient[[i]] <- suc - 1 #  As a value of 2 mean that only there is one more gene with the gene being analyzed, we substract 1.
    }
    names(genes_per_patient) <- names(splited_patients)
    genes_per_patient_list[[k+1]] <- genes_per_patient
    names(genes_per_patient_list)[[k+1]] <- as.character(k)
    message(paste0("tetha=",k))
  }

  # Start a 11 x 20 matrix to be filled with the hierarchical clustering accuracy values.
  final_accuracy_matrix <- matrix(0, ncol= 20, nrow= 11)
  final_kk_used <- matrix(0, ncol= 20, nrow= 11)
  rownames(final_accuracy_matrix) <- as.character(0:10)
  rownames(final_kk_used) <- as.character(0:10)

  for(u in 1:11){ # 1 to 11 because of non 0-based language: If we want, for example, theta to be 0, we set u=1
    for(val in 1:20){ # lambda
      preserve_genes_per_patient <- genes_per_patient_list[[u]] # u==k == tetha + 1
      genes_per_patient <- preserve_genes_per_patient
      genes_per_patient <- lapply(genes_per_patient,function(x) x[x<=val]) # val = lambda. We filter the genes that are over the lambda value tested

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

      patient_matrix2 <- patient_matrix[- which(rownames(patient_matrix) %in% WHATEVER),]  #Exclude patients with missing data from clustering

      # Obtain optimal clusters

      pamk.best <- pamk(patient_matrix2)
      kk <- pamk.best$nc ## kk is the optimal number of clusters for the particular optimization.

      # Perform hierarchical clustering with the suggested number of clusters
      patient_matrix3 <- t(patient_matrix2)
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

      sum_matrix <- arbol_splited_mat + ground_truth_mat
      tab_sum_matrix <- table(sum_matrix)
      zeros <- tab_sum_matrix["0"]
      if(is.na(zeros)){
        zeros <- 0
      }
      twos <- tab_sum_matrix["2"]
      if(is.na(twos)){
        twos <- 0
      }
      accuracy <- (zeros+twos)/sum(tab_sum_matrix,na.rm = T)
      final_accuracy_matrix[u,val] <- accuracy
      final_kk_used[u,val] <- kk
    }
    message(paste0("Acuracies for tetha=",u-1," calculated."))
  }
  write.table(final_accuracy_matrix,file=paste("data/Randomizations/Accuracymatrix_",en_cual_guardo[p],".csv"),row.names=T, col.names = F)
  write.table(final_kk_used,file=paste("data/Randomizations/Clustersused_",en_cual_guardo[p],".csv"),row.names=T, col.names = F)

  message(paste0("Random_",en_cual_guardo[p]," done"))
}
