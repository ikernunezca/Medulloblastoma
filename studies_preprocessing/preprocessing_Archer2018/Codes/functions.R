binary.transformation <- function(sel.not.bin, ## index of the non binary dataset (integer)
                                  exp.matrix, ## list containing all datasets 
                                  perc.abnormals= 0.2 ## construction of the threshold
){
  exp.M <- exp.matrix[[sel.not.bin]]
  threshold.given.perc <- stats::quantile(as.vector(exp.M), prob=1-perc.abnormals, na.rm = TRUE)
  if (length(which(is.na(exp.M)))!=0) exp.M[is.na(exp.M)] <- 0
  exp.M.bin <- exp.M
  exp.M.bin[exp.M<threshold.given.perc] <- rep(0, length(which(exp.M<threshold.given.perc)))
  exp.M.bin[!(exp.M<threshold.given.perc)] <- rep(1, length(which(!(exp.M<threshold.given.perc))))
  return(exp.M.bin)
}

from.exp.data.to.bin.data <- function(exp_files, ## List with gene information at multiple layer of observation 
                                      perc.abnormals = 0.3 ## A priori threshold of abnormality
                                      ){
  is.binary <- lapply(exp_files, function(x) return(sum(na.omit(x)==0) + sum(na.omit(x)==1) == length(na.omit(x)))) #Check whether the matrix is binary
  is.binary <- unlist(is.binary)
  sel.not.bin <- which(!is.binary) # Indexes of the continuous matrices (vector)
  if (all(is.binary)){
    bin.exp.data <- exp_files
  }else{
    bin.exp.data <- mclapply(sel.not.bin, binary.transformation, exp.matrix = exp_files, 
                             perc.abnormals = perc.abnormals)
    if (any(is.binary)){
      bin.exp.data <- c(bin.exp.data, exp_files[which(is.binary)])
    }
  }
  return(bin.exp.data)
}


per.sample.network <- function(sel.sample,  ## Either an integer or sample name
                               list.bin.data.matrices, ## List with the binary data matrices
                               name.sample =FALSE  ## Logical indicating if sample name is provided as sel.sample
                               ){
  if (!name.sample) {sample.name <- paste("S", sel.sample, sep="") }
  else { sample.name <- sel.sample}
  sample.exp.data <- lapply(list.bin.data.matrices, function(x) {
    data.vec <- x[ ,sel.sample]
    names(data.vec) <- rownames(x)
    return(data.vec)
    })
  list.per.layer <- mclapply(1:length(sample.exp.data), per.layer.network, exp.data.per.sample= sample.exp.data, sample.name=sample.name )
  df.per.sample.network <- do.call(rbind, list.per.layer)
  return(df.per.sample.network)
} 

per.layer.network <- function(sel.layer, exp.data.per.sample, sample.name){
  layer.exp.data.per.sample <- exp.data.per.sample[[sel.layer]]
  from <- sample.name 
  to <- paste(names(layer.exp.data.per.sample)[which(layer.exp.data.per.sample == 1)], names(exp.data.per.sample)[sel.layer], sep=".")
  per.layer.net.data.frame <- data.frame(from=from, to=to, layer=names(exp.data.per.sample)[sel.layer])
  return(per.layer.net.data.frame)
}

connect.ggi.with.multilayer.network <- function(multi.layer.net.no.ggi, set.selected.genes, ggi.net){
  ggi.to.mln.network <- mclapply(set.selected.genes, from.gene.to.gene.layers, multi.layer.net = multi.layer.net.no.ggi)
  cl.ggi.to.mln.network <- ggi.to.mln.network[-which(sapply(ggi.to.mln.network, is.null))]
  ggi.to.mln.network <- do.call(rbind, cl.ggi.to.mln.network)
  full.ml.network <- rbind(multi.layer.net.no.ggi, ggi.net, ggi.to.mln.network)
  return(full.ml.network)
}

from.gene.to.gene.layers <- function(sel.gene, multi.layer.net){
  pos.link.gene.gene.layer <- grep(sel.gene, multi.layer.net$to)
  if (length(pos.link.gene.gene.layer)==0) {
    df.per.gene.ggi.to.mln <- NULL
  } else {
    df.per.gene.ggi.to.mln <- data.frame(from=sel.gene, to= unique(multi.layer.net[pos.link.gene.gene.layer, ]$to), layer= "ggi.mln" )
  }
  return(df.per.gene.ggi.to.mln)
}












