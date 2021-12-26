library(survival)
library(AdaptiveLasso)
library(MultiRNG)
library(CPE)
library(parallel)
library(foreach)
library(doParallel)

dat <- read.csv("~/Desktop/LungReadyStudy.csv", row.names = 1)
dat <- dat[, which(apply(dat, 2, sum)/nrow(dat) > 0.03)]

set.seed(1)
train <- dat
test <- dat

time <- train$time2
delta <- train$status
zz <- as.matrix(train[,-c(1:3)])

results <-  alasso_cox_tuned(NN = 15, time, delta, zz)


####################################
######### ENSEMBLE LEARNING ########
####################################
dat <- read.csv("~/Desktop/LungReadyStudy.csv", row.names = 1)
dat <- dat[, which(apply(dat, 2, sum)/nrow(dat) > 0.03)]

runs <- 10 # 30
cores <- 3

cl <- parallel::makeCluster(cores,
                            setup_strategy = "sequential",
                            outfile = paste0("Adaptive_LASSO_log.txt"))
registerDoParallel(cl)

final.lasso <- list()

LASSO <- foreach(run=1:runs) %dopar% {
  library(survival)
  library(AdaptiveLasso)
  set.seed(run+210793)
  cat(paste("Run : ", run,"\n",sep=""),
      file=paste0("Run_Adaptive_LASSO.txt"), append=TRUE)

  rm.samples <- sample(1:nrow(dat), ceiling(nrow(dat)*1/3),replace = FALSE)
  train <- dat[-rm.samples,]
  test <- dat[rm.samples,]

  time <- train$time2
  delta <- train$status
  zz <- as.matrix(train[,-c(1:3)])

  results <- alasso_cox_tuned(NN = 10, time, delta, zz)

  return(results)
}
save(LASSO, file = "results/Lung_ensemble.Rdata")
