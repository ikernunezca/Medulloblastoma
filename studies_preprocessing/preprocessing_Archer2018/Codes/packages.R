if(!"BiocManager" %in% installed.packages()){
  install.packages("BiocManager")
}
library(BiocManager)  ## Package installation
if(!"RCy3" %in% installed.packages()){
  BiocManager::install("RCy3")
}
library(RCy3) ## Network visualization 
if(!"paxtoolsr" %in% installed.packages()){
  BiocManager::install("paxtoolsr")
}
library(paxtoolsr)  ## Interact with bioPAX (download of PPI)
if(!"igraph" %in% installed.packages()){
  install.packages("igraph")
}
library(igraph)  ## Network visualization
if(!"plyr" %in% installed.packages()){
  install.packages("plyr")
}
library(plyr)
if(!"biomaRt" %in% installed.packages()){
  install.packages("biomaRt")
}
library(biomaRt)  ## Retrieve annotated HSapiens genes 
if(!"parallel" %in% installed.packages()){
  install.packages("parallel")
}
library(parallel)  ## Parallelisation
