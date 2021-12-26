library(survival)
library(AdaptiveLasso)
library(MultiRNG)
library(CPE)

#### Generate Data ####
simulWeib <- function(N, lambda, rho, beta, rateC,X)
{
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda * exp(X %*% beta)))^(1 / rho)
  C <- rexp(n=N, rate=rateC)
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  # data set
  return(data.frame(#id=1:N,
    time=time,
    delta=status,
    X))
}

sample.size <- c(150, 300)
cor.set <- 0
covs <- 30
datasets <- 25

for(i in 1:length(sample.size)){
  n = sample.size[i]

  print(paste("Current Sample Size :" ,n))

  MSEs <- as.data.frame(matrix(nrow = datasets, ncol = 3))
  MSEs_var <- as.data.frame(matrix(nrow = datasets, ncol = 3))
  Corrects <- as.data.frame(matrix(nrow = datasets, ncol = 3))
  Incorrects <- as.data.frame(matrix(nrow = datasets, ncol = 3))
  CPEs <- as.data.frame(matrix(nrow = datasets, ncol = 3))

  for(j in 1:datasets){
    set.seed(j+210793)
    print(paste0("Dataset: ",j))

    # cov.mat <- matrix(rep(cor.set,covs^2),nrow=covs,ncol=covs)
    # diag(cov.mat) <- 1
    # X <- draw.d.variate.uniform(no.row=n,d=covs,cov.mat)
    X <- matrix(data = rnorm(n = n*covs), nrow = n, ncol = covs)
    # X <- matrix(data = rbinom(n = n*covs, size = 1, prob = 1/2), nrow = n, ncol = covs)
    colnames(X) <- paste0("Cov",1:covs)
    rownames(X) <- paste0("Patient",1:n)

    # generate risk and time variables #
    beta <- c(rep(1,5),rep(-1,5),rep(0,covs-10))
    dat <- simulWeib(N = n, lambda = 0.6, rho = 1, beta = beta, rateC = 0.01, X = X)
    time <- dat$time
    delta <- dat$delta
    sum(delta)/length(delta)
    zz <- as.matrix(X)

    results <- alasso_cox_tuned(NN = 15, time, delta, zz)

    # MSE #
    MSEs[j,] <- c(t(as.matrix(results$beta.raw - beta)) %*% cov(zz) %*% as.matrix(results$beta.raw - beta),
                  t(as.matrix(results$beta.lasso - beta)) %*% cov(zz) %*% as.matrix(results$beta.lasso - beta),
                  t(as.matrix(results$beta.alasso - beta)) %*% cov(zz) %*% as.matrix(results$beta.alasso - beta)
    )

    MSEs_var[j,] <- c(t(as.matrix(results$beta.raw - beta)) %*% diag(summary(coxph(Surv(time,delta)~zz))$coefficients[,3]) %*% as.matrix(results$beta.raw - beta),
                      t(as.matrix(results$beta.lasso - beta)) %*% diag(results$beta.sd.lasso) %*% as.matrix(results$beta.lasso - beta),
                      t(as.matrix(results$beta.alasso - beta)) %*% diag(results$beta.sd.alasso) %*% as.matrix(results$beta.alasso - beta)
    )

    # correct 0's #
    Corrects[j,] <- c(sum(results$beta.raw[-c(1:10)] == 0)/20,
                      sum(results$beta.lasso[-c(1:10)] == 0)/20,
                      sum(results$beta.alasso[-c(1:10)] == 0)/20
    )
    # incorrect 0's #
    Incorrects[j,] <- c(sum(results$beta.raw[c(1:10)] == 0)/10,
                        sum(results$beta.lasso[c(1:10)] == 0)/10,
                        sum(results$beta.alasso[c(1:10)] == 0)/10
    )
    # CPE #
    # new random dataset #
    X_test <- matrix(data = rnorm(n = n*covs), nrow = n, ncol = covs)
    # X <- matrix(data = rbinom(n = n*covs, size = 1, prob = 1/2), nrow = n, ncol = covs)
    colnames(X_test) <- paste0("Cov",1:covs)
    rownames(X_test) <- paste0("Patient",1:n)

    # generate risk and time variables #
    beta <- c(rep(0.5,5),rep(-0.5,5),rep(0,covs-10))
    dat_test <- simulWeib(N = n, lambda = 0.2, rho = 1, beta = beta, rateC = 0.01, X = X_test)
    time_test <- dat_test$time
    delta_test <- dat_test$delta

    # fit new data with trained parameters #
    fit_raw <- coxph(Surv(time,delta) ~ ., data = dat_test, init = results$beta.raw,iter=0)
    fit_lasso <- coxph(Surv(time,delta) ~ ., data = dat_test, init = results$beta.lasso,iter=0)
    fit_alasso <- coxph(Surv(time,delta) ~ ., data = dat_test, init = results$beta.alasso,iter=0)

    # Get CPE #
    CPEs[j,] <- c(
      as.numeric(phcpe(coxph(Surv(dat_test$time, dat_test$delta) ~predict(fit_raw, newdata=as.data.frame(dat_test))))),
      as.numeric(phcpe(coxph(Surv(dat_test$time, dat_test$delta) ~predict(fit_lasso, newdata=as.data.frame(dat_test))))),
      as.numeric(phcpe(coxph(Surv(dat_test$time, dat_test$delta) ~predict(fit_alasso, newdata=as.data.frame(dat_test)))))
    )

  }
  out_results <- list("MSEs" = MSEs, "MSEs_var" = MSEs_var,"Corrects" = Corrects,
                      "Incorrects" = Incorrects, "CPEs" = CPEs)
  save(out_results, file = paste0("results/StrongModel_", n,"_", "10cens","R.data"))

}
