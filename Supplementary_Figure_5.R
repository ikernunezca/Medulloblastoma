# Load dependencies
library(data.table)
library(ggplot2)

# Define variables to fill with the data coming from each randomization
all_accuracies <- c()
all_clusters <- c()
all_maximums <- c()

for(i in 1:10000){
  suc <- fread(paste0("data/Randomizations/Accuracymatrix_ ",i," .csv"))
  clus <- fread(paste0("data/Randomizations/Clustersused_ ",i," .csv"))
  # theta = 0 , lambda = 6
  all_accuracies <- c(all_accuracies,unname(unlist(suc[1,6]))) # Vector with accuracy value at [0,6] from each randomization
  all_clusters <- c(all_clusters,unname(unlist(clus[1,6]))) # Vector with number of suggested clusters at [0,6] from each randomization
  all_maximums <- c(all_maximums,unname(unlist(max(suc)))) # Vector with the maximum value of accuracy of each randomization. 
  print(i)
}

p_value_5_clusters <- unname(table(all_clusters==5)["TRUE"]) / 10000 # Is suggesting 5 clusters at [0,6] significant?
max_accuracy <- max(all_accuracies,na.rm = T) 
min_accuracy <- min(all_accuracies,na.rm = T)
mean_accuracy <- mean(all_accuracies,na.rm = T)
sd_accuracy <- sd(all_accuracies,na.rm = T)

# Generate data for plotting
data_for_plot <- as.data.frame(matrix(data=NA, ncol= 3, nrow=2))
colnames(data_for_plot) <- c("X","Accuracy","sd")
data_for_plot[1,] <- c("Original gene associations",0.9493,0) # 0.9493 is the original accuracy
data_for_plot[2,] <- c("Mean Randomized",mean_accuracy,sd_accuracy)
data <- data_for_plot
class(data[,2]) <- "numeric"
class(data[,3]) <- "numeric"

# Create ggplot with error bars indicating the standard deviation
pdf("data/Plots/Supplementary_Figure_5.pdf")
ggplot(data) +
  geom_bar( aes(x=X, y=Accuracy), stat="identity", fill="skyblue", alpha=0.7,orientation = "x") +
  geom_errorbar( aes(x=X, ymin=Accuracy-sd, ymax=Accuracy+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) 
dev.off()
