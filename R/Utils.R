#' Calculate partial log likelihood
#'
#' @export

loglik = function(n,delta,z,beta){
  L = 0.0
  for(i in 1:n){
    temp = 0.0
    for(j in i:n) temp = temp + exp(sum(beta*z[j,]))
    L = L - delta[i]*(sum(beta*z[i,])-log(temp))
  }
  return(L)
}


#' Calculate derivative of partial log likelihood (score function)
#'
#' @export

dloglik = function(n,delta,z,beta){
  p = length(beta)
  L = numeric(p)
  for(i in 1:n){
    temp = 0.0
    temp1 = numeric(p)
    for(j in i:n){
      temp = temp + exp(sum(beta*z[j,]))
      temp1 = temp1 + z[j,]*exp(sum(beta*z[j,]))
    }
    L = L - delta[i]*(z[i,] - temp1/temp)
  }
  return(L)
}

#' Calculate double derivative of partial log likelihood (Hessian)
#'
#' @export

ddloglik = function(n,delta,z,beta){
  p = length(beta)
  L = matrix(0,nrow=p,ncol=p)
  for(i in 1:n){
    temp = 0.0
    temp1 = numeric(p)
    temp2 = matrix(0,nrow=p,ncol=p)
    for(j in i:n){
      temp = temp + exp(sum(beta*z[j,]))
      temp1 = temp1 + z[j,]*exp(sum(beta*z[j,]))
      temp2 = temp2 + z[j,]%*%t(z[j,])*exp(sum(beta*z[j,]))
    }
    if(delta[i] == 0.0) L = L + 0.0
    else L = L + delta[i]*(temp2/temp - temp1%*%t(temp1)/(temp*temp))
  }
  return(L)
}


#' Shooting algorithm
#'
#' @export

wshoot <- function(n,p,x,y,init,weight,lambda,maxiter,tol)
{
  Q = t(x)%*%x
  B = t(x)%*%y
  i=0
  status = 0

  lams =lambda*weight
  oldbeta <- init
  tmpbeta <- oldbeta

  while (i<maxiter && status==0){
    for (j in 1:p){
      s<-ss2(j,tmpbeta,Q,B)
      if (s > lams[j])
        tmpbeta[j]<-(lams[j]-s)/(2*Q[j,j])
      else if (s < (-lams[j]) )
        tmpbeta[j]<-(-lams[j]-s)/(2*Q[j,j])
      else
        tmpbeta[j]<- 0.0
    }
    dx<-max(abs(tmpbeta-oldbeta))
    oldbeta <- tmpbeta
    if (dx<=tol)
      status <- 1
    i <- i+1
  }
  tmpbeta
}


#' Calculate partial log likelihood
#'
#' @export

ss2 <- function(j,tmpb,Q,B)
{
  a <- sum(tmpb*Q[,j])-tmpb[j]*Q[j,j]
  s <- 2*(a-B[j])
  return(s)
}

#' Calculate partial log likelihood
#'
#' @export

ginv = function(X, tol = sqrt(.Machine$double.eps)){
  s = svd(X)
  nz = s$d > tol * s$d[1]
  if(any(nz)) s$v[,nz] %*% (t(s$u[,nz])/s$d[nz])
  else X*0
}

#' Calculate partial log likelihood
#'
#' @export

normalize = function(x){
  y = (x-mean(x))/sqrt(sum((x-mean(x))^2))
  return(y)
}

#' Compute least square error
#'
#' @export

lse <- function(x,y)
{ Q <- t(x)%*%x
P <- t(x)%*%y
beta <- solve(Q,P)
beta
}

#' Compute adaptive lasso estimates
#'
#' @export

alasso_cox <- function(NN = 10, time, delta, z, lambda)
{
  iter = 100
  tol = 1.0e-10
  n = length(time)
  p = length(z[1,])
  true.sd = sqrt(apply(z,2,var)*(n-1))
  delta = delta[order(time)]
  ordz = z[order(time),]
  z = apply(ordz,2,normalize)
  sd = sqrt(apply(z,2,var)*(n-1))
  #computing initial estimates
  ii = 0
  beta = numeric(p)
  while(ii < NN){
    fn=loglik(n,delta,z,beta)
    G=dloglik(n,delta,z,beta)
    H=ddloglik(n,delta,z,beta)

    X = chol(H)
    vecY = forwardsolve(t(X),H%*%beta-G)
    lsbeta = lm(vecY~-1+X)$coef
    beta1 = lsbeta
    dx = max(abs(beta1-beta))
    ii = ii + 1
    istop = ii
    if(dx <= 1.0e-5) ii = NN
    beta = beta1
  }
  inibeta = beta
  #computing adaptive-Lasso solutions
  sd = sd/abs(inibeta)
  ii = 0
  beta = numeric(p)
  while(ii < NN){
    fn=loglik(n,delta,z,beta)
    G=dloglik(n,delta,z,beta)
    H=ddloglik(n,delta,z,beta)

    X = chol(H)
    vecY = forwardsolve(t(X),H%*%beta-G)
    lsbeta = lm(vecY~-1+X)$coef
    beta1 = wshoot(p,p,X,vecY,init=lsbeta,sd,lambda,iter,tol)
    dx = max(abs(beta1-beta))
    ii = ii + 1
    istop = ii
    if(dx <= 1.0e-5) ii = NN
    beta = beta1
  }
  beta = beta/true.sd
  inibeta = inibeta/true.sd

  w = diag(2*abs(beta))
  ginvw = ginv(w)
  A = H + lambda*sd*ginvw
  ps = sum(diag(solve(A)%*%H)) - sum(beta == 0)
  GCV = fn/(n*(1-ps/n)^2)
  return(rbind(t(c(inibeta,GCV)),t(c(beta,GCV))))
}


#########################################

#########################################

#' Compute adaptive tuned lasso estimates
#'
#' @export

alasso_cox_tuned <- function(NN = 10, time, delta, z, gd = 10)
{

  # set up #
  ps = numeric(gd)
  GCV = numeric(gd)
  newGCV = numeric(gd)
  n = length(time)
  p = length(z[1,])
  iter = 100
  tol = 1.0e-10
  beta.s = matrix(0,nrow=gd,ncol=p)
  newbeta.s = matrix(0,nrow=gd,ncol=p)
  beta.sd = matrix(0,nrow=gd,ncol=p)
  newbeta.sd = matrix(0,nrow=gd,ncol=p)

  # get initial values #
  true.sd = sqrt(apply(z,2,var)*(n-1))
  delta = delta[order(time)]
  ordz = z[order(time),]
  beta.raw = coxph(Surv(time,delta)~z)$coef
  z = apply(ordz,2,normalize)
  sd = sqrt(apply(z,2,var)*(n-1))
  beta.ini = coxph(Surv(time,delta)~z)$coef

  ###Computing Lasso estimates
  for(j in 1 : gd){
    lam = 2^((j-1)/3) - 1.0
    beta = beta.ini
    ii = 0
    while(ii < NN){
      fn=loglik(n,delta,z,beta)
      G=dloglik(n,delta,z,beta)
      H=ddloglik(n,delta,z,beta)

      X = chol(H)
      vecY = forwardsolve(t(X),H%*%beta-G)
      lsbeta = lm(vecY~-1+X)$coef

      beta1=wshoot(p,p,X,vecY,init=lsbeta,sd,lam,iter,tol)
      dx = max(abs(beta1-beta))
      ii = ii + 1
      istop = ii
      if(dx <= 1.0e-5) ii = NN
      beta = beta1
    }
    beta.s[j,] = beta
    fn = loglik(n,delta,z,beta)
    cat(istop,dx,fn,"\n")
    w = diag(2*abs(beta))
    ginvw = ginv(w)
    A = H + lam*sd*ginvw
    ps[j] = sum(diag(solve(A)%*%H)) - sum(beta.s[j,] == 0)
    GCV[j] = fn/(n*(1-ps[j]/n)^2)
    #variance computation
    if(j == 1) beta.sd[j,]=sqrt(diag(solve(H)))/true.sd
    else {
      w = numeric(p)
      w[abs(beta) > 0] = 1/abs(beta[abs(beta) > 0])
      w[beta == 0] = 1.0e10
      A = H + diag(0.5*lam*w)
      invA = solve(A)
      beta.sd[j,]=sqrt(diag(invA%*%H%*%invA))/true.sd
    }
  }

  ###Computing Adaptive Lasso estimates
  newbeta.s[1,] = beta.s[1,]
  newbeta.sd[1,] = beta.sd[1,]
  newGCV[1] = GCV[1]

  beta = newbeta.s[1,]
  sd = sd/abs(beta)
  for(j in 2 : gd){
    lam = 3^((j-1)/2) - 1.0
    beta = beta.ini
    ii = 0
    while(ii < NN){
      fn=loglik(n,delta,z,beta)
      G=dloglik(n,delta,z,beta)
      H=ddloglik(n,delta,z,beta)

      X = chol(H)
      vecY = forwardsolve(t(X),H%*%beta-G)
      lsbeta = lm(vecY~-1+X)$coef
      vecX = as.vector(X)

      beta1=wshoot(p,p,X,vecY,init=lsbeta,sd,lam,iter,tol)
      dx = max(abs(beta1-beta))
      ii = ii + 1
      istop = ii
      if(dx <= 1.0e-5) ii = NN
      beta = beta1
    }
    newbeta.s[j,] = beta
    fn = loglik(n,delta,z,beta)
    cat(istop,dx,fn,"\n")
    w = diag(2*abs(beta))
    ginvw = ginv(w)
    A = H + lam*sd*ginvw
    ps[j] = sum(diag(solve(A)%*%H)) - sum(newbeta.s[j,] == 0)
    newGCV[j] = fn/(n*(1-ps[j]/n)^2)
    #variance computation
    w = numeric(p)
    w[abs(beta) > 0] = 1/abs((beta[abs(beta) > 0])^2)
    w[beta == 0] = 1.0e10
    A = H + diag(0.5*lam*w)
    invA = solve(A)
    v = numeric(p)
    v[abs(beta) > 0] = 1/abs((beta[abs(beta) > 0])^2)
    v[beta == 0] = 0.0
    B = H + diag(0.5*lam*v)
    invH = solve(H)
    newbeta.sd[j,]=sqrt(diag(invA%*%B%*%invH%*%B%*%invA))/true.sd
  }

  for(j in 1:gd){
    beta.s[j,] = beta.s[j,]/true.sd
    newbeta.s[j,] = newbeta.s[j,]/true.sd
  }

  beta.sd = round(beta.sd*1000)/1000
  newbeta.sd = round(newbeta.sd*1000)/1000

  ###solutions
  beta.lasso = beta.s[rank(GCV)==1,]
  beta.alasso = newbeta.s[rank(newGCV)==1,]

  return(list("beta.lasso" = beta.lasso,
              "beta.alasso" = beta.alasso,
              "beta.ini" = beta.ini,
              "beta.raw" = beta.raw,
              "beta.sd.lasso" = beta.sd[rank(GCV)==1,],
              "beta.sd.alasso" = newbeta.sd[rank(newGCV)==1,]
              ))
}



#########################################

#########################################

#' Modified shooting algo for 2 lambdas
#'
#' @export


wshoot_enet <- function(n,p,x,y,init,weight,lambda,maxiter,tol)
{
  Q = t(x)%*%x
  B = t(x)%*%y
  i=0
  status = 0

  lams = apply(lambda*weight,1,sum)
  oldbeta <- init
  tmpbeta <- oldbeta

  while (i<maxiter && status==0){
    for (j in 1:p){
      s<-ss2(j,tmpbeta,Q,B)
      if (s > lams[j])
        tmpbeta[j]<-(lams[j]-s)/(2*Q[j,j])
      else if (s < (-lams[j]) )
        tmpbeta[j]<-(-lams[j]-s)/(2*Q[j,j])
      else
        tmpbeta[j]<- 0.0
    }
    dx<-max(abs(tmpbeta-oldbeta))
    oldbeta <- tmpbeta
    if (dx<=tol)
      status <- 1
    i <- i+1
  }
  tmpbeta
}


#########################################

#########################################

#' Compute adaptive tuned ENET estimates
#'
#' @export

aenet_cox_tuned <- function(NN = 10, time, delta, z, gd = 10)
{

  # set up #
  ps = numeric(gd)
  GCV = numeric(gd)
  newGCV = numeric(gd)
  n = length(time)
  p = length(z[1,])
  iter = 100
  tol = 1.0e-10
  beta.s = matrix(0,nrow=gd,ncol=p)
  newbeta.s = matrix(0,nrow=gd,ncol=p)
  beta.sd = matrix(0,nrow=gd,ncol=p)
  newbeta.sd = matrix(0,nrow=gd,ncol=p)

  # get initial values #
  true.sd = sqrt(apply(z,2,var)*(n-1))
  delta = delta[order(time)]
  ordz = z[order(time),]
  beta.raw = coxph(Surv(time,delta)~z)$coef
  z = apply(ordz,2,normalize)
  sd = sqrt(apply(z,2,var)*(n-1))
  beta.ini = coxph(Surv(time,delta)~z)$coef

  ###Computing Lasso estimates
  for(j in 1 : gd){
    lam = 2^((j-1)/3) - 1.0
    beta = beta.ini
    ii = 0
    while(ii < NN){
      fn=loglik(n,delta,z,beta)
      G=dloglik(n,delta,z,beta)
      H=ddloglik(n,delta,z,beta)

      X = chol(H)
      vecY = forwardsolve(t(X),H%*%beta-G)
      lsbeta = lm(vecY~-1+X)$coef

      beta1=wshoot(p,p,X,vecY,init=lsbeta,sd,lam,iter,tol)
      dx = max(abs(beta1-beta))
      ii = ii + 1
      istop = ii
      if(dx <= 1.0e-5) ii = NN
      beta = beta1
    }
    beta.s[j,] = beta
    fn = loglik(n,delta,z,beta)
    cat(istop,dx,fn,"\n")
    w = diag(2*abs(beta))
    ginvw = ginv(w)
    A = H + lam*sd*ginvw
    ps[j] = sum(diag(solve(A)%*%H)) - sum(beta.s[j,] == 0)
    GCV[j] = fn/(n*(1-ps[j]/n)^2)
    #variance computation
    if(j == 1) beta.sd[j,]=sqrt(diag(solve(H)))/true.sd
    else {
      w = numeric(p)
      w[abs(beta) > 0] = 1/abs(beta[abs(beta) > 0])
      w[beta == 0] = 1.0e10
      A = H + diag(0.5*lam*w)
      invA = solve(A)
      beta.sd[j,]=sqrt(diag(invA%*%H%*%invA))/true.sd
    }
  }

  ###Computing Adaptive Lasso estimates
  newbeta.s[1,] = beta.s[1,]
  newbeta.sd[1,] = beta.sd[1,]
  newGCV[1] = GCV[1]

  beta = newbeta.s[1,]

  # new sd #
  sd = cbind(sd/abs(beta) , sd^2/beta^2)
  for(j in 2 : gd){
    lam = c(3^((j-1)/2) - 1.0,3^((j-1)/2) - 1.0)
    beta = beta.ini
    ii = 0
    while(ii < NN){
      fn=loglik(n,delta,z,beta)
      G=dloglik(n,delta,z,beta)
      H=ddloglik(n,delta,z,beta)

      X = chol(H)
      vecY = forwardsolve(t(X),H%*%beta-G)
      lsbeta = lm(vecY~-1+X)$coef
      vecX = as.vector(X)

      beta1=wshoot_enet(p,p,X,vecY,init=lsbeta,sd,lam,iter,tol)
      dx = max(abs(beta1-beta))
      ii = ii + 1
      istop = ii
      if(dx <= 1.0e-5) ii = NN
      beta = beta1
    }
    newbeta.s[j,] = beta
    fn = loglik(n,delta,z,beta)
    cat(istop,dx,fn,"\n")
    w = diag(2*abs(beta))
    ginvw = ginv(w)
    A = H + apply(lam*sd,1,sum)*ginvw
    ps[j] = sum(diag(solve(A)%*%H)) - sum(newbeta.s[j,] == 0)
    newGCV[j] = fn/(n*(1-ps[j]/n)^2)
    #variance computation
    w = numeric(p)
    w[abs(beta) > 0] = 1/abs((beta[abs(beta) > 0])^2)
    w[beta == 0] = 1.0e10
    A = H + diag(0.5*lam*w)
    invA = solve(A)
    v = numeric(p)
    v[abs(beta) > 0] = 1/abs((beta[abs(beta) > 0])^2)
    v[beta == 0] = 0.0
    B = H + diag(0.5*lam*v)
    invH = solve(H)
    newbeta.sd[j,]=sqrt(diag(invA%*%B%*%invH%*%B%*%invA))/true.sd
  }

  for(j in 1:gd){
    beta.s[j,] = beta.s[j,]/true.sd
    newbeta.s[j,] = newbeta.s[j,]/true.sd
  }

  beta.sd = round(beta.sd*1000)/1000
  newbeta.sd = round(newbeta.sd*1000)/1000

  ###solutions
  beta.lasso = beta.s[rank(GCV)==1,]
  beta.alasso = newbeta.s[rank(newGCV)==1,]

  return(list("beta.lasso" = beta.lasso,
              "beta.aenet" = beta.alasso,
              "beta.ini" = beta.ini,
              "beta.raw" = beta.raw,
              "beta.sd.lasso" = beta.sd[rank(GCV)==1,],
              "beta.sd.aenet" = newbeta.sd[rank(newGCV)==1,]
  ))
}

