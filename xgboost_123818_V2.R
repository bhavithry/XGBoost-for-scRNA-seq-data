# Clear environment
########
rm(list=ls())
gc()
########

# Import libraries
########
library(dplyr)
library(ggplot2)
library(ROCR)
library(xgboost)
########

# Import data
########
setwd("C:/Users/bhavi/Desktop/Research")

# Import data
setwd("C:/Users/bhavi/OneDrive/Desktop/Research/Data")

# Import data
class1 <- read.csv("Root_single_cell_shr_datamatrix.csv")
class2 <- read.csv("Root_single_cell_wt_datamatrix.csv")
gene <- class1[,1]
# Transpose and change to data frame
class1 <- data.frame(t(data.matrix(class1[,2:length(class1)])))
class2 <- data.frame(t(data.matrix(class2[,2:length(class2)])))
colnames(class1) <- gene
colnames(class2) <- gene
########

# Data Pre-processing
########
# Shuffle rows
set.seed(100)
class1<- class1[sample(nrow(class1)),]
set.seed(100)
class2<- class2[sample(nrow(class1)),]

class1['y'] <- 1
class2['y'] <- 0
l1 <- nrow(class1)
l2 <- nrow(class2)

class_all <- rbind(class1,class2)

# Remove genes not expressed or equally expressed in all cells
removeZeroVar <- function(df){
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}
data <- removeZeroVar(class_all)

data['y'] <- class_all['y']

gc()
########

# Pre-allocate variables
########
df_genes <- data.frame(genes=colnames(data[,-ncol(data)]), beta = numeric(ncol(data)-1))
df_genes$beta <- 0
row.names(df_genes) <- df_genes$genes

k=1
j=l1+1
nfolds = 2
l1fold = floor(l1/nfolds)
l2fold = floor(l2/nfolds)

xgboost.val <- rep(0,nfolds)

# create parameter list
params <- list(
  eta = .1,
  max_depth = 5,
  min_child_weight = 2,
  subsample = .8,
  colsample_bytree = .9
)
#########

# K-Fold Cross validation starts here
########
for (i in 1:nfolds){
  
# Train Test Split
#######
  
  if(i==nfolds){
    test_rows <- c(k:l1,j:(l1+l2))
  } else {
    test_rows <- c(k:(k+l1fold-1),j:(j+l2fold-1))
  }
  
  x_test <- as.matrix(data[test_rows,1:length(data)-1])
  y_test <- as.matrix(data[test_rows, length(data)])
  x_train <- as.matrix(data[-test_rows,1:length(data)-1])
  y_train <- as.matrix(data[-test_rows, length(data)])
  k=k+l1fold
  j=j+l2fold
#######


# XGBoost
#######
#xgboost.fit <- xgb.cv( params = params
#                       , data = x_train
#                      , label = y_train
#                       , nrounds = 100
#                      , nfold = 10
#                      , objective = "binary:logistic"  # for classification models
#                       , verbose = 0   # silent
#                       , early_stopping_rounds = 10
#                      )

start.time <- Sys.time()
set.seed(100)
xgboost.fit <- xgboost( params = params
                       , data = x_train
                       , label = y_train
                       , nrounds = 100
                       , objective = "binary:logistic"  # for classification models
                       , verbose = 0   # silent
                       , early_stopping_rounds = 10
)
end.time <- Sys.time()
comp.time <- end.time - start.time
#######

# AUC
#######
xgboost.probabilities <- predict(xgboost.fit, newdata = x_test)
xgboost.y_pred <- ifelse(xgboost.probabilities > 0.5, 1, 0)
xgboost.pred <- prediction(xgboost.y_pred, y_test)
xgboost.auc.perf <- performance(xgboost.pred, measure = "auc")
xgboost.val[i] <- xgboost.auc.perf@y.values
#######

}
#######

# Get top important genes
#######
importance_matrix <- xgb.importance(model = xgboost.fit)
xgb.plot.importance(importance_matrix, top_n = 20, measure = "Gain")
#######

# Get number of trees that minimize error
########
xgboost.fit$evaluation_log %>%
  dplyr::summarise(
    ntrees.train = which(train_logloss_mean == min(train_logloss_mean))[1],
    rmse.train   = min(train_logloss_mean),
    ntrees.test  = which(test_logloss_mean == min(test_logloss_mean))[1],
    rmse.test   = min(test_logloss_mean),
  )
#######

# plot error vs number trees
########
ggplot(xgboost.fit$evaluation_log) +
  geom_line(aes(iter, train_logloss_mean), color = "red") +
  geom_line(aes(iter, test_logloss_mean), color = "blue")
########