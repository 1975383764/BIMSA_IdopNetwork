library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)
pheno_df <- read.csv("F_2015-HT-DIA.csv",header = F)
row.names(pheno_df) <- pheno_df[,1]
pheno_df <- pheno_df[,-1]
geno_df <- read.csv('geno_df.csv')
data_9<- c(132,173,197,224,236,286,307,315,343,393,496,687,704,775,927,1116,1165,1171,1239,1378,269,460,588,845,861,936,1008,1058,1366)
geno_df <- geno_df[,-1:-3]
geno_df <- geno_df[,as.numeric(row.names(pheno_df))]
geno_df <- geno_df[-data_9,]
data_big <- read.csv("data_big.csv")
data_big <- data_big[,-1]
data_small <- read.csv("data_small.csv")
data_small <- data_small[,-1]

X <- as.matrix(data_big)
t <- c(1:11)
start_params1 <- c(111, 1.62405,  0.6)
start_params2 <- c(9, 1.7,  0.35)


logistic <- function(t, miu_par) {
  miu_par[1] / (1 + miu_par[2]*exp(-miu_par[3]*t))
}

loss <- function(params, t, y) {
  y_pred <- c(logistic(t,params))
  sum((y - y_pred)^2)
}
get_init_par <- function(data,k){
  
  get_logistic_par <- function(y,t,start_params){
    logistic_par <- optim(start_params, fn=loss, t = t, y=y, method = "BFGS")$par
    return(logistic_par)
  }
  
  #get initial pars based on k-means
  init_cluster <- kmeans(data,centers = k,iter.max = 1000)
  cuM <- init_cluster$centers #??ֵ????
  init_curve_par <- cbind(t(sapply(1:k,function(c)get_logistic_par(y=as.numeric(cuM[c,1:11]),t=t,start_params=start_params1))),
                          t(sapply(1:k,function(c)get_logistic_par(y=as.numeric(cuM[c,12:22]),t=t,start_params=start_params2))))
  init_SAD_par <- c(0.8,0.35,0.8,IQR(diag(cov(X[,12:22]))))
  init_pro <- table(init_cluster$cluster)/nrow(data)
  return_object <- list(init_SAD_par,init_curve_par,init_pro)
  names(return_object)<-c("init_SAD_par","init_curve_par","init_pro")
  return(return_object)
}

get_cluster <- function(data,k){
  
  input <- biFunClu_intial_pars
  Delta <- 100; iter <- 1; itermax <- 100;
  get_biSAD1 <- function(par){
    n=ncol(data)/2
    get_SAD1_covmatrix <- function(par){
      phi <- par[1]; gamma <- par[2]; 
      sigma <- array(dim=c(n,n))
      diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
      sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^seq(1:(n-c))*diag(sigma)[c]))
      sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
      return(gamma^2*sigma)
    }
    
    sig1 <- get_SAD1_covmatrix(par[1:2])
    sig2 <- get_SAD1_covmatrix(par[3:4])
    sig12 <- array(0, dim=c(n,n))
    sigma1 <- cbind(sig1,sig12)
    sigma2 <- cbind(sig12,sig2)
    sigma <- rbind(sigma1,sigma2)
    return(sigma)
  }
  mle <- function(par,data,prob){
    par1 <- par[1:4]
    par2 <- matrix(par[-c(1:4)],nrow = k,ncol = 6)
    miu <- t( sapply(1:k, function(c)c(logistic(t,par2[c,1:3]),
                                       logistic(t,par2[c,(4:6)]))))
    temp_S <- sapply(1:k,function(c)dmvnorm(data,miu[c,],get_biSAD1(par1))*prob[c] )
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  cat(paste0("Start biFunClu Calculation ","\n","Cluster_number=",k))
  while ( Delta > 1 && iter <= itermax ) {
    # initiation
    if(iter == 1){
      init_SAD_par <- input[[1]]
      init_curve_par <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_SAD_par,init_curve_par)
    LL_mem <- mle(old_par,data,pro)
    miu <- t( sapply(1:k, function(c)c(logistic(t,init_curve_par[c,1:3]),
                                       logistic(t,init_curve_par[c,(4:6)]))))
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             miu[c,],
                                             get_biSAD1(init_SAD_par))*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_SAD_par <- new_par$par[1:4]
    init_curve_par <- matrix(new_par$par[-c(1:4)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    if (Delta > 20000)
      break
    cat("iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
}
  k=k
  biFunClu_intial_pars <- get_init_par(data=X,k=k)
  biFunClu_results <- get_cluster(data=X,k=k)










