setwd("C:/Users/adi44/Desktop/ATM/Spring2018/STAT646/EXAM2")

set.seed(101)
Y <- read.csv("exam2_data2.csv", row.names = 1) #Load the data
m <- nrow(Y) # Get the rows
Y_train <- Y[, 1:225] # Training set
Y_test <- Y[, 226:300] # Testing set
GRP_train <- factor(rep(c("A", "B", "C"), each = 75)) # Create training labels
GRP_test <- factor(rep(c("A", "B", "C"), each = 25)) # Create testing labels
Y_train_df <- data.frame("GRP" = GRP_train, t(Y_train)) # Create the training dataset
Y_test_df <-data.frame("GRP" = GRP_test, t(Y_test)) # Create Testing dataset
Y_df <- rbind(Y_train_df, Y_test_df)# Whole data set


library(class) # Contains KNN

pred_KNN <- rep(NA, 10) # array for containing precition accuracy
for(k in 1:10) { # Compute KNN accuracy rate for K= 1 to 10
  pred_KNN_k <- knn(Y_train_df[, -1], Y_test_df[, -1], Y_train_df$GRP, k = k)
  pred_KNN[k] <- mean(pred_KNN_k == Y_test_df$GRP)
}
error=1-pred_KNN # Compute the error
plot(1:10, pred_KNN, xlab = "K", ylab = "Accuracy", xaxt = "n", type = "l", lwd = 2)# Plot accuracy vs K
axis(1, at = 1:10)

# Find the K for which the test accuracy is the highest
which(pred_KNN==max(pred_KNN))
print(sessionInfo())
