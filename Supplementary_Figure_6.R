# Get randomization output files
files <- list.files("data/Randomizations/")
acc_files <- files[grep("Accuracy",files)]

example_mat <- read.table(acc_files[1])
sum_mat <- matrix(0,nrow=nrow(example_mat),ncol=ncol(example_mat))

# Perform the sum for all the accuracy matrices...

for(i in 1:length(acc_files)){
  sup <- read.table(acc_files[i])
  sum_mat <- sup + sum_mat
  print(i)
}

# ... and then obtain the mean of each theta-lamba pair by dividing by the number of files.
mean_mat <- sum_mat / length(acc_files)


## Plotting
to_plot <- mean_mat[,-1] # Take out theta value column
to_plot <- unlist(to_plot)

real_acc <- read.table("data/Supplementary_Tables/Supplementary_Table_2.csv")  # Real obtained accuracies
real_acc <- real_acc[,-1]
real_acc <- unlist(real_acc)

# Formatting
dat_real_acc <- data.frame(real_acc)
colnames(dat_real_acc) <- 'Accuracy'
dat_rand_acc <- data.frame(to_plot)
colnames(dat_rand_acc ) <- 'Accuracy'
dat_real_acc$legend <- 'Observed accuracies'
dat_rand_acc$legend <- 'Mean randomized accuracies'
densities_what <- rbind(dat_real_acc, dat_rand_acc)

library(ggplot2)

png(file= "data/Plots/Supplementary_Figure_6.png",width= 1920, height=1080)
ggplot(densities_what, aes(Accuracy,fill = legend)) + geom_density(alpha = 0.2)
dev.off()
