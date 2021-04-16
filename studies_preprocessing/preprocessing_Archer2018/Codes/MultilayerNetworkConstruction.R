### Archer 2018 data
names_exp_files <- list.files(path= 'studies_preprocessing/preprocessing_Archer2018/OmicsData',pattern = ".csv") 
names_exp_files <- paste0('studies_preprocessing/preprocessing_Archer2018/OmicsData/',names_exp_files)
exp_files  <- list()
for(i in 1:length(names_exp_files)){
  exp_files[[i]] <- read.csv(names_exp_files[i], sep=";", header=TRUE,
                             dec=",", stringsAsFactors=FALSE)
}
names(exp_files)<- unlist(lapply(names_exp_files, function(x) substr(x, start=1, stop=nchar(x)-4)))  ## Data matrices are now stored in "exp_files" list
try1 <- lapply(exp_files, function(x){
  gene.names <- unique(x$GeneSymbol)
  median.per.gene <- lapply(gene.names, function(gn){
    sel.rows <- which(x$GeneSymbol==gn)
    if (length(sel.rows)==1){
      median.val <- x[, -c(1:2)][sel.rows,]
    }else{
      sel.associated.genes <- x[, -c(1:2)][sel.rows,]
      ### Here I consider the median
      median.val <- apply(sel.associated.genes, 2, median, na.rm=TRUE)
      median.val <- unlist(median.val)
    }
    return(median.val)
  })
  median.per.gene <- do.call(rbind, median.per.gene)
  rownames(median.per.gene) <- gene.names
  return(median.per.gene)
})
exp_files <- try1
names(exp_files) <- c("Acetyl","phospho","phospho_psty","proteomics","RNA_seq")

## CONSTRUCT BINARY DATA MATRICES IF CONTINUOUS DATA MATRICES ARE PROVIDED
binary.data.from.exp <- from.exp.data.to.bin.data(exp_files, perc.abnormals = 0.3)


## CONSTRUCT MULTILAYER NETWORK (without GGI)
## Here select the sample you want to consider:
set.samples <- unique(c(unique(colnames(exp_files$Acetyl)), 
                      unique(colnames(exp_files$phospho)), 
                      unique(colnames(exp_files$phospho_psty)),
                      unique(colnames(exp_files$proteomics)), 
                      unique(colnames(exp_files$RNA_seq))))
binary.data.from.all.samples <- lapply(binary.data.from.exp, function(x){
  ordered.expr <- try(x[, set.samples], silent = TRUE)
  if (inherits(ordered.expr, "try-error")){
    ordered.expr <- matrix(data=0, nrow = nrow(x), ncol = length(set.samples),)
    colnames(ordered.expr) <- set.samples
    rownames(ordered.expr) <- rownames(x)
    ordered.expr[,colnames(x)] <- as.matrix(x)
  }
  return(ordered.expr)
})

## CONSTRUCT MULTILAYER NETWORK (without GGI)
multi.layer.net.no.ggi <- mclapply(set.samples, per.sample.network, list.bin.data.matrices = binary.data.from.all.samples, name.sample=TRUE)
multi.layer.net.no.ggi <- do.call(rbind, multi.layer.net.no.ggi)

write.table(multi.layer.net.no.ggi,'studies_preprocessing/preprocessing_Archer2018/Networks/multi.layer.no.ggi.net',sep='\t',row.names=F,quote=F)


# ### DOWNLOAD A GGI NERTWORK FROM PATHWAY COMMONS 
# GGI <- downloadPc2(version=11, verbose = TRUE) ## Load for instance 48. Note that the format of the choosen file is ".hgnc.txt.gz".
# ## A new version can be similarly download and upload by setting version = 12
# 
# ## DEFINE THE SET OF GENES TO BE USED IN THE CONSTRUCTION OF THE MULTILAYER NETWORK 
# set.selected.genes <- GGI$nodes$PARTICIPANT
# 
# ## CONSTRUCTION THE GGI NETWORK 
# ## Here we roughly consider all from-to edges but the code could be refine considering the type of interaction 
# ## described in GGI$edges$INTERACTION_TYPE
# to.del <- which((GGI$edges$PARTICIPANT_A=="")|(GGI$edges$PARTICIPANT_B=="")|is.na(GGI$edges$PARTICIPANT_A)|is.na(GGI$edges$PARTICIPANT_B))
# GGI$edges <- GGI$edges[-to.del, ]
# ggi.net <- data.frame(from=GGI$edges$PARTICIPANT_A, to=GGI$edges$PARTICIPANT_B, layer="GGI") 
# 
# ## CONSTRUCT THE WHOLE NETWORK 
# # There can be links only between genes present on both the mln and the ggi 
# set.common.genes <- set.selected.genes[which(set.selected.genes %in% c(rownames(exp_files$Acetyl), 
#                                                                        rownames(exp_files$phospho), 
#                                                                        rownames(exp_files$phospho_psty),
#                                                                        rownames(exp_files$proteomics), 
#                                                                        rownames(exp_files$RNA_seq)))]
# multi.layer.net.with.ggi <- connect.ggi.with.multilayer.network(multi.layer.net.no.ggi, ggi.net= ggi.net, set.selected.genes =  set.common.genes)
# write.table(multi.layer.net.with.ggi,'studies_preprocessing/preprocessing_Archer2018/Networks/multi.layer.net',sep='\t',row.names=F,quote=F)

