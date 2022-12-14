# Clear environment
##############
rm(list=ls())
gc()
########

# Import libraries
##############
library(dplyr)
library(ggplot2)
library(xgboost)
########

# Import data
##############
setwd("C:/Users/bhavi/Desktop/Research")

# Import data
class1 <- read.csv("GSE81861_CRC_tumor_epithelial_cells_FPKM.csv")
class2 <- read.csv("GSE81861_CRC_NM_epithelial_cells_FPKM.csv")
gene <- class1$X

# Transpose and change to data frame
class1 <- data.frame(t(data.matrix(class1[,-1])))
class2 <- data.frame(t(data.matrix(class2[,-1])))

colnames(class1) <- gene
colnames(class2) <- gene
########

# Data Pre-processing
##############
# Shuffle rows
set.seed(100)
class1<- class1[sample(nrow(class1)),]
set.seed(100)
class2<- class2[sample(nrow(class2)),]

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
##############
df_genes <- data.frame(genes=colnames(data[,-ncol(data)]), beta = numeric(ncol(data)-1))
df_genes$beta <- 0
row.names(df_genes) <- df_genes$genes

k=1
j=l1+1
nfolds = 1
l1fold = floor(l1/nfolds)
l2fold = floor(l2/nfolds)


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
##############
for (i in 1:nfolds){
  
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
  
  
# XGBoost
##############
  set.seed(100)
  xgboost.fit <- xgb.cv( params = params
                         , data = x_train
                         , label = y_train
                         , nrounds = 100
                         , nfold = 10
                         , objective = "binary:logistic"  # for classification models
                         , verbose = 0   # silent
                         , early_stopping_rounds = 10
  )
  
##################
}

# get number of trees that minimize error
#######
xgboost.fit$evaluation_log %>%
  dplyr::summarise(
    ntrees.train = which(train_logloss_mean == min(train_logloss_mean))[1],
    rmse.train   = min(train_logloss_mean),
    ntrees.test  = which(test_logloss_mean == min(test_logloss_mean))[1],
    rmse.test   = min(test_logloss_mean),
  )
#######

###### plot error vs number trees
########
ggplot(xgboost.fit$evaluation_log) +
  geom_line(aes(iter, train_logloss_mean), color = "red") +
  geom_line(aes(iter, test_logloss_mean), color = "blue")
########

# Apparently 61 trees is optimal


# create hyperparameter grid
###########
hyper_grid <- expand.grid(
  eta = c(.01, .05, .1, .3),
  max_depth = c(1, 3, 5, 7),
  min_child_weight = c(1, 3, 5, 7),
  subsample = c(.65, .8, 1), 
  colsample_bytree = c(.8, .9, 1),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)
#############

# grid search
###########
for(i in 1:nrow(hyper_grid)) {
  
  # create parameter list
  params <- list(
    eta = hyper_grid$eta[i],
    max_depth = hyper_grid$max_depth[i],
    min_child_weight = hyper_grid$min_child_weight[i],
    subsample = hyper_grid$subsample[i],
    colsample_bytree = hyper_grid$colsample_bytree[i]
  )
  
  # reproducibility
  set.seed(100)
  
  # train model
  xgb.tune <- xgb.cv(
    params = params,
    data = x_train,
    label = y_train,
    nrounds = 1000,
    nfold = 10,
    objective = "binary:logistic",  # for regression models
    verbose = 0,               # silent,
    early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_logloss_mean)
  hyper_grid$min_RMSE[i] <- min(xgb.tune$evaluation_log$test_logloss_mean)
}

hyper_grid %>%
  dplyr::arrange(min_RMSE) %>%
  head(10)
##########