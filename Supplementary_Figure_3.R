# Load dependencies
library(mgcv)
library(MASS)
library(CmmD)

# First :Define functions needed.
# Function from: gavinsimpson/derivSimulCI.R (https://gist.github.com/gavinsimpson/ca18c9c789ef5237dbc6)
`derivSimulCI` <- function(mod, n = 200, eps = 1e-7, newdata, term,
                           samples = 10000) {
  stopifnot(require("MASS"))
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x) - (2*eps), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  newDF <- data.frame(newD) ## needs to be a data frame for predict
  X0 <- predict(mod, newDF, type = "lpmatrix")
  newDF <- newDF + eps
  X1 <- predict(mod, newDF, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  ## sample draws from the posterior distribution of model coefficients
  Rbeta <- t(mvrnorm(n = samples, coef(mod), vcov(mod)))
  ## loop over the terms
  for(i in seq_len(nt)) {
    want <- grep(t.labs[i], colnames(X1))
    lD[[i]] <- list(deriv = Xp[, want] %*% coef(mod)[want],
                    simulations = Xp[, want] %*% Rbeta[want, ])
  }
  class(lD) <- "derivSimulCI"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}

plot.derivSimulCI <- function(x, alpha = 0.05, polygon = TRUE,
                              sizer = FALSE, term,
                              eval = 0, lwd = 3,
                              col = "lightgrey", border = col,
                              ylab, xlab, main, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else {
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  ciFUN <- function(x, alpha) {
    ahalf <- alpha / 2
    apply(x$simulations, 1, quantile, probs = c(ahalf, 1 - ahalf))
  }
  CI <- lapply(x[seq_len(l)], ciFUN, alpha = alpha)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  on.exit(layout(1))
  for(i in term) {
    lwr <- CI[[i]][1,]
    upr <- CI[[i]][2,]
    ylim <- range(upr, lwr)
    plot(x$eval[,i], x[[i]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,i], rev(x$eval[,i])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,i], upr, lty = "dashed")
      lines(x$eval[,i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
    }
  }
  invisible(x)
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

# Create a vector with the paths where Molti's Output files are saved.
structures_12 <- paste0("data/Molti_Output/",seq(0.5,30,0.5),".csv")
# Detect community trajectories and tree distances between each gene. 
curie_to_12_full <- CmmD_from_community_structures(nodelist = NULL, community_structures = structures_12, resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = "hamming",threads = 7)
curie_to_12_full$hamming_distance_matrix = curie_to_12_full$distance_matrix * 24 # This transformation is needed because parallel dist is weighted.
curie = curie_to_12_full
# 24 = length(seq(0.5,12,0.5)) -> number of resolution values analyzed
class(curie$gene_community_matrix) <- "numeric"

# Plotting Mean gene number per community & Number of communities
maxis <- c()

for(i in 1:ncol(curie$gene_community_matrix)){
  maxim_cur <- max(curie$gene_community_matrix[,i])
  maxis <- c(maxis,maxim_cur)
}
ngenes <- 18948
mean_gene_per_community <- ngenes/maxis
pdf("data/Plots/curie_mean_gene_per_comm_vs_n_communities_log.pdf")
plot(log(maxis),
     axes=FALSE,
     main="Curie Multilayer Mean gene size per community & N communities",
     ylab="",
     col="blue",xlab="MolTi's resolution parameter",ylim=c(0,8))
axis(2)
axis(1,at = 0:60, labels=seq(0,30,0.5))
points(log(mean_gene_per_community), col="red")
legend("bottomright",
       c("log(N communities)","log(Mean gene size per community)"),
       fill=c("blue","red")
       ,box.lty = 0)
dev.off()

# Code adapted from: https://rpubs.com/hrlai/gam_inflection to perform curve fitting

m <- gam(log(mean_gene_per_community) ~ s(maxis)) # Generate GAM model
logYpred <- predict(m)

# Supplementary figure 3
pdf("data/Plots/predicted_log_mean_gene_per_community_vs_smoothed_n_communities.pdf")
plot(maxis, log(mean_gene_per_community),main= "F(x): \n Predicted log Mean gene size per community \n vs Smoothed n communities",ylab="Predicted log Mean gene size per community" ,xlab= "s(N communities)")
lines(maxis, logYpred,col="blue")
dev.off()

fd <- derivSimulCI(m)

pdf("data/Plots/Derivative_x.pdf")
plot(fd, sizer = TRUE,main= "",xlab="")
title(main="F'(x)",
      xlab="N communities")
dev.off()

CI <- lapply(fd[1], function(x) t(apply(x$simulations, 1, quantile, probs = c(0.025, 0.975))))
first.zero.slope.index <- min(which(sign(CI$maxis[, "2.5%"]) != sign(CI$maxis[, "97.5%"])))
fd$eval[first.zero.slope.index] # Number of communities at the inflection point
