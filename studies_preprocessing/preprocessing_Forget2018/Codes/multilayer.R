source("studies_preprocessing/preprocessing_Forget2018/Codes/required_packages.R", echo=FALSE) ## (internet connection is assummed to be on)
source("studies_preprocessing/preprocessing_Forget2018/Codes/required_functions.R", echo=FALSE)

## CONSTRUCTION OF MULTILAYER NETWORK
## from gene-gene interaction network and 
## experimenal data at different layers of observation. 

### DOWNLOAD A GGI NERTWORK FROM PATHWAY COMMONS 
#GGI <- downloadPc2(version=8, verbose = TRUE) ## Load for instance 103. Note that the format of the choosen file is ".hgnc.txt.gz".
## A new version can be similarly download and upload by setting version = 12

## LOAD EXPERIMENTAL DATA
## Assume the experimental data is all in one folder
# setwd("~/MultilayerNetworks/data") ## Set the working directory to this folder
## THE DATA SHOULD BE ORGANIZED AS FOLLOWS: 
## Samples as columns and genes as rows
## Sample/gene names should be provided as data matrix column/row names, respectively. 
## In addition, sample ordering must be the same through datasets (not necessary for the genes). 
## In case of missing values or missed observed sample layer, add a zero value/column, respectively. 
## Note: Name the files wrt the biol. layer 
names_exp_files <- list.files(path= 'studies_preprocessing/preprocessing_Forget2018/OmicsData',pattern = ".csv") 
names_exp_files <- paste0('studies_preprocessing/preprocessing_Forget2018/OmicsData/',names_exp_files)

exp_files  <- list()
for(i in 1:length(names_exp_files)){
  exp_files[[i]] <- read.csv(names_exp_files[i], sep=",", header=TRUE, row.names=1)
}
names(exp_files)<- unlist(lapply(names_exp_files, function(x) substr(x, start=1, stop=nchar(x)-4)))  ## Data matrices are now stored in "exp_files" list
names(exp_files) <- c("Methylation","Phospho","Protein","RNA_seq")

## DEFINE THE SET OF GENES TO BE USED IN THE CONSTRUCTION OF THE MULTILAYER NETWORK 
#set.selected.genes <- GGI$nodes$PARTICIPANT

## Refinements:
## REDUCE THE GGI NETWORK TO THE SET OF SELECTED GENES
## REDUCE THE DATA MATRICES TO THE SET OF SELECTED GENES

## CONSTRUCT BINARY DATA MATRICES IF CONTINUOUS DATA MATRICES ARE PROVIDED
binary.data.from.exp <- from.exp.data.to.bin.data(exp_files, perc.abnormals = 0.3)


## CONSTRUCT MULTILAYER NETWORK (without GGI)
set.samples <- colnames(exp_files$Methylation)
multi.layer.net.no.ggi <- mclapply(set.samples, per.sample.network, list.bin.data.matrices = binary.data.from.exp, name.sample=TRUE)
multi.layer.net.no.ggi <- do.call(rbind, multi.layer.net.no.ggi)

write.table(multi.layer.net.no.ggi,'studies_preprocessing/preprocessing_Forget2018/Networks/multi.layer.no.ggi.net',sep='\t',row.names=F,quote=F)

print("Multi-omic set generated.")

# ## CONSTRUCTION THE GGI NETWORK 
# ## Here we roughly consider all from-to edges but the code could be refine considering the type of interaction 
# ## described in GGI$edges$INTERACTION_TYPE
# to.del <- which((GGI$edges$PARTICIPANT_A=="")|(GGI$edges$PARTICIPANT_B=="")|is.na(GGI$edges$PARTICIPANT_A)|is.na(GGI$edges$PARTICIPANT_B))
# GGI$edges <- GGI$edges[-to.del, ]
# ggi.net <- data.frame(from=GGI$edges$PARTICIPANT_A, to=GGI$edges$PARTICIPANT_B, layer="GGI") 
# 
# ## CONSTRUCT THE WHOLE NETWORK 
# # There can be links only between genes present on both the mln and the ggi 
# set.common.genes <- set.selected.genes[which(set.selected.genes %in% c(rownames(exp_files$Phospho), 
#                                       rownames(exp_files$Protein), rownames(exp_files$RNAseq), rownames(exp_files$Methylation)))]
# multi.layer.net.with.ggi <- connect.ggi.with.multilayer.network(multi.layer.net.no.ggi, ggi.net= ggi.net, set.selected.genes =  set.common.genes)


# save.image("MultiLayerNetwork.RData")



