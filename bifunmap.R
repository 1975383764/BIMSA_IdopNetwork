library(MASS)
library(mvtnorm)
library(deSolve)
library(ggplot2)
library(orthopolynom)
pheno_df <- read.csv("F_2015-HT-DIA.csv",header = F)
row.names(pheno_df) <- pheno_df[,1]
pheno_df <- pheno_df[,-1]
geno_df <- read.csv("geno_df.csv",header = T)
geno_df <- geno_df[,-1:-3]
geno_df <- geno_df[,as.numeric(row.names(pheno_df))]

############################################
logistic <- function(t, miu_par) {
  miu_par[1] / (1 + miu_par[2]*exp(-miu_par[3]*t))
}

get_SAD1_covmatrix <- function(par,n){
  phi <- par[1]; gamma <- par[2];
  sigma <- array(dim=c(n,n))
  diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
  sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^seq(1:(n-c))*diag(sigma)[c]))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  return(gamma^2*sigma)
}

get_biSAD1 <- function(par,n){
  sig1 <- get_SAD1_covmatrix(par[1:2],n)
  sig2 <- get_SAD1_covmatrix(par[3:4],n)
  sig12 <- array(0, dim=c(n,n))
  sigma1 <- cbind(sig1,sig12)
  sigma2 <- cbind(sig12,sig2)
  sigma <- rbind(sigma1,sigma2)
  return(sigma)
}

L0 <- function(par){
  miu = c(logistic(t,par[1:3]),logistic(t,par[4:6]))
  SAD1 = get_biSAD1(par[7:10],length(t))
  L0 = -sum(dmvnorm(y_all,miu,SAD1,log = T))
  L0
  if (is.infinite(L0)) {
    return(1e10)
  }
  return(L0)
}
L1 <- function(par){
  SAD1 <- get_biSAD1(par[13:16],length(t))
  L_0 <- -sum(dmvnorm(pheno_0,c(logistic(t,par[1:3]),logistic(t,par[4:6])),SAD1,log = T))
  L_1 <- -sum(dmvnorm(pheno_1,c(logistic(t,par[7:9]),logistic(t,par[10:12])),SAD1,log = T))
  LL <- L_0+L_1
  if (is.infinite(LL)) {
    return(1e10)
  }
  return(LL)
}
L2 <- function(par){
  SAD1 <- get_biSAD1(par[19:22],length(t))
  L_0 <- -sum(dmvnorm(pheno_0,c(logistic(t,par[1:3]),logistic(t,par[4:6])),SAD1,log = T))
  L_1 <- -sum(dmvnorm(pheno_1,c(logistic(t,par[7:9]),logistic(t,par[10:12])),SAD1,log = T))
  L_2 <- -sum(dmvnorm(pheno_2,c(logistic(t,par[13:15]),logistic(t,par[16:18])),SAD1,log = T))
  LL <- L_0+L_1+L_2
  if (is.infinite(LL)) {
    return(1e10) 
  }
  return(LL)
}

  genotsable <- geno_df[i,]
  marker=genotsable
  pheno_0 <- pheno_df[which(marker==0),]
  pheno_1 <- pheno_df[which(marker==1),]
  pheno_2 <- pheno_df[which(marker==2),]
  pheno_9 <- pheno_df[which(marker==9),]
  y_all <- rbind(pheno_0,pheno_1,pheno_2)
  
  t <- c(1:11)
  
  get_init_pars <- function(par){
    y <- as.numeric(colMeans(y_all))
    y1 <- c(logistic(t,par[1:3]),logistic(t,par[4:6]))
    ssr <- sum((y1-y)^2)
  }
  init_curve_par <- optim(c(mean(y_all[,11]), 1.6,  0.6,mean(y_all[,22]),1.7,0.38)
                          ,get_init_pars,method = "BFGS")$par
  
  init_par <- c(init_curve_par,1.05,3.6,1.05,0.2)
  lower <- c(-Inf, -Inf,-Inf,-Inf, -Inf,-Inf, 0,0,0,0)
  upper <- c(Inf, Inf,Inf,Inf, Inf,Inf,Inf,Inf, Inf,Inf)
  NH_0 <- optim(init_par,L0,method = "L-BFGS-B",upper=upper,lower=lower,control = list(maxit = 10000, factr = 1e1, pgtol = 1e-10))
  NH_0
  
  if (nrow(pheno_2)==0) {
    h1_pars <- c(NH_0$par[1:6],NH_0$par[1:6],NH_0$par[7:10])
    NH_1 <- optim(h1_pars,L1,method="L-BFGS-B",upper=upper,lower=lower,control = list(maxit = 10000, factr = 1e1, pgtol = 1e-10))
  }else{
    h2_pars <- c(NH_0$par[1:6],NH_0$par[1:6],NH_0$par[1:6],NH_0$par[7:10])
    NH_1 <- optim(h2_pars,L2,method="L-BFGS-B",upper=upper,lower=lower,control = list(maxit = 10000, factr = 1e1, pgtol = 1e-10))
  }
  LR <- 2*(NH_0$value - NH_1$value)
  

