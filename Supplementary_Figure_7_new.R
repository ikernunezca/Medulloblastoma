acc <- c(94.94,94.94,94.94,94.94,88.57,93.31,79.27,75.67,79.76,77.96,76.98,77.31,74.20,74.20,73.88,72.90,73.39,72.73,72.41,72.73,72.57,71.75,71.92,73.39,72.08,70.45,69.96,72.24,74.20,57.06,53.14,47.59,43.02,31.43,31.43,31.43,0.0,0.0,0.0,0.0)
thet <- c(0,3,7,10,2,10,9,0,7,0,10,5,8,5,9,7,10,0,9,7,1,4,4,7,10,7,5,5,10,1,3,8,6,1,5,8,0,0,0,0)
lamb <- c(5,5,7,5,20,9,3,15,4,18,7,4,6,2,6,5,8,20,2,5,2,1,7,1,4,4,2,14,20,4,2,1,2,2,2,2,1,1,1,1)
remgen <- c(0,1743.63,2534.45,3329.66,3850.84,5138.58,5408.61,5438.61,5565.97,5594.34,5678.76,5697.95,5725.16,5747.26,5767.24,5772.55,5774.11,5778.03,5849.03,5851.74,5860.66,5881.47,5889.21,5900.05,5901.45,5905.5,5907.11,5908.45,5921.45,5944.13,5947.87,5949.11,5949.47,5950.34,5950.47,5950.61,5950.63,5950.63,5950.63,5950.63)




# tab <- read.csv("tab.txt",header=F)

par(mar=c(5, 4, 4, 6) + 0.1)

pdf('data/Plots/Supp_Figure_7_paper_version.pdf',width = 13, height=10)
plot(1:40, acc, type='o', col='#3D9F88', lwd=2, xlab = 'Optimization cycles', ylab="", ylim=c(0,100))
points(1:40, thet, type='o', col='magenta', lwd=2)
points(1:40, lamb, type='o', col='red', lwd=2)
mtext(expression(paste("Accuracy and optimal ", theta," and ",lambda," values",sep="")), side=2, line=2.5)
box()

par(new = TRUE)
plot(1:40, remgen, type='o', col='grey60', lwd=2, axes=FALSE, ylab="", xlab="")
axis(side=4, at = pretty(range(remgen)))
mtext("Average removed genes per patient (cumulative sum)", side=4, line=2.5)

legend(33,5000, legend=c("Accuracy", expression(paste("Optimal ",theta,sep="")),
	expression(paste("Optimal ",lambda,sep="")),"Removed genes"),col=c('#3D9F88','magenta','red','grey60'), lty=rep(1,4), cex=0.8)

dev.off()
