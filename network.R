library(mvtnorm)
library(pbapply)
library(parallel)
library(deSolve)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)
setwd("D:/data/H_D/figure/网络73")
load("B_FunClu_results 5 .RData")###
times <- seq(1,11,length=30)
k=5###
logistic <- function(t, miu_par) {
  miu_par[1] / (1 + miu_par[2]*exp(-miu_par[3]*t))
}
get_interaction <- function(data,col){
  n <- nrow(data)
  clean_data <- data
  gene_list <- list()
  m <- clean_data[,col]
  M <- clean_data[,-col]
  x_matrix <- M
  x_matrix <- as.matrix(x_matrix)
  vec <- sapply(1:length(M[1,]),function(c)cor(m,M[,c]))
  x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
  x_matrix <- as.matrix(x_matrix)
  name <- colnames(clean_data)
  ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 6,alpha = 0)
  best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
  
  fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 6,alpha = 0.99,
                       penalty.factor = abs(1/best_ridge_coef),
                       keep = TRUE)
  best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
  
  gene_list_one <- list()
  gene_list_one[[1]] <- name[col]
  gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
  gene_list_one[[3]] <- best_alasso_coef1@x[-1]
  gene_list[[col]] <- gene_list_one
  
  return(gene_list_one)
}
cluster_mean <- t(sapply(1:k, function(c)logistic(times,biFunClu_results$curve_par[c,4:6])))
rownames(cluster_mean) <- 1:k
module_relationship <- pblapply(1:k,function(c) get_interaction(t(cluster_mean),c))

#----------------------

get_module_result <- function(k,times,cluster_result){
  get_effect <- function(pars, effect, times, y0) {
    
    Lorenz <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        dE <- a1 * E * (1 - (E / k1)^a1)
        list(c(dE))
      })
    }
    
    parameters <- list(a1 = pars[1], k1 = pars[2])
    rk4 <- function(f, state, t, dt, parameters) {
      k1 <- dt * f(t, state, parameters)[[1]]
      k2 <- dt * f(t + dt / 2, state + k1 / 2, parameters)[[1]]
      k3 <- dt * f(t + dt / 2, state + k2 / 2, parameters)[[1]]
      k4 <- dt * f(t + dt, state + k3, parameters)[[1]]
      state + dt / 6 * (k1 + 2 * (1 - 1 / sqrt(2)) * k2 + 2 * (1 + 1 / sqrt(2)) * k3 + k4)
    }
    dt <- times[2] - times[1]
    y <- numeric(length(times))
    y[1] <- y0
    for (i in 2:length(times)) {
      y[i] <- rk4(Lorenz, c(E = y[i - 1]), times[i - 1], dt, parameters)
    }
    dy_fit <- effect * c(0, diff(y))
    return(cumsum(dy_fit))
  }
  legendre <- function(pheno,t,par){
    min_value <- min(pheno)
    max_value <- max(pheno)
    xn <- 2 * (pheno - min_value) / (max_value - min_value) - 1
    y <- t*(par[1]+par[2]*xn+par[3]*(0.5*(3*xn^2-1)))+par[4]
    y[1] <- 0
    return(y)
  }
  fn <- function (par, pheno_ind,pheno_dep,state,t){
    ind_pars <- matrix(par[1:2])
    dep_pars <- matrix(par[-1:-2],ncol=4)
    if ( nrow(dep_pars)==1 ) {
      y <- c(get_effect(ind_pars,pheno_ind,times,state)+state+
               legendre(pheno_dep,times,dep_pars))
    }else{
      ind_y <- c(get_effect(ind_pars,pheno_ind,times,state)+state);
      dep_y <- sapply(1:nrow(dep_pars), function(c)
        legendre(pheno_dep[c,],times,dep_pars[c,]));
      y=ind_y+rowSums(dep_y)
    }
    #add penalty
    alpha=5e-4
    ridge <- sum((pheno_ind-y)^2+alpha*(sum(ind_pars^2)+sum(dep_pars^2)))
    return(ridge)}
  
  get_par <- function(times,par) {
    init_index <- c(1:k)
    get_value <- function(effect,data,times,c){
      #input
      ind <- data[[1]]
      dep <- data[[2]]
      ind_no <- as.numeric(which(colnames(effect)==ind))
      dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effect)==dep[c])))
      sample_m <- as.numeric(effect[,ind_no])
      sample_m1 <- matrix(as.numeric(t(effect[,dep_no])),nrow=length(dep_no))
      H0_state <- c(E=as.numeric(sample_m[1]))
      if(length(dep_no)!=0){
        init_paras_list <- list(
          c(0.808017175318673,4.48837019409984,-0.268279149197042,-0.000200216192752123,0.153506491798908,0.567714696284384,-0.268279149197042,-0.000200216192752123,0.153506491798908,0.567714696284384),
          c(0.700877960724756,4.95456727733836,0.909578449092805,0.336549715138972,-0.320503944065422,0.291405486408621,0.909578449092805,0.336549715138972,-0.320503944065422,0.291405486408621),
          c(0.383851149189286,4.35051672067493,-0.435603402089328,-0.204167664051056,0.091941628139466,0.481124949641526,-0.435603402089328,-0.204167664051056,0.091941628139466,0.481124949641526),
          c(0.451851344690658,6.55534715438262,-0.597769915591925,0.163165043108165,-0.674537206534296,-0.808200033847243),
          c(0.536085225921124,4.19288107077591,-0.113792103715241,0.438245038036257,0.2798729124479,-0.01774957543239,-0.113792103715241,0.438245038036257,0.2798729124479,-0.01774957543239))}
      
      LS_init_pars <- init_paras_list[[init_index[c]]]
      saveRDS(LS_init_pars, file = paste0("LS_init_pars_", c, ".rds"))
      LS1  = optim(LS_init_pars, fn = fn, pheno_ind = sample_m,pheno_dep = sample_m1, 
                   state = H0_state,t=times, method = "Nelder-Mead", control = list(maxit = 1000))#Nelder-Mead
      effect_par <- LS1$par
      
      return(effect_par)
    }
    lop_par <- pblapply(1:nrow(cluster_mean),function(c)get_value(t(cluster_mean),
                                                                  module_relationship[[c]],times,c))
    return(list(lop_par,module_relationship))
  }
  all_lop_par <- get_par(times=times, par=cluster_result$curve_par)
  #output for result-
  get_output <- function(relationship,par,effect,times){
    output <- list()
    output[[1]] <- relationship[[1]]  
    output[[2]] <- relationship[[2]]  
    output[[3]] <- par[1:2]
    output[[4]] <- matrix(par[-1:-2],ncol=4)
    ind_no <- as.numeric(which(colnames(effect)==output[[1]]))
    dep_no <- as.numeric(sapply(1:length(output[[2]]), 
                                function(c) which(colnames(effect)==output[[2]][c])))
    inital_value <- c(E=effect[,ind_no][1])
    sample_m <- effect[,ind_no]
    sample_m1 <- matrix(as.numeric(t(effect[,dep_no])),nrow=length(dep_no))
    ind_effect <- get_effect(output[[3]],sample_m,times,y0=inital_value)+inital_value
    
    if (length(dep_no)==1) {
      dep_effect <- legendre(effect[,dep_no],times,output[[4]])
    }else{
      dep_effect <- sapply(1:length(dep_no), function(c)
        legendre(sample_m1[c,],times,output[[4]][c,]))
      colnames(dep_effect) <- dep_no
    }
    #------------
    all_effect <- cbind(ind_effect,dep_effect)
    effect_mean <- apply(all_effect,2,mean)
    output[[5]] <- effect_mean
    output[[6]] <- all_effect
    output[[7]] <- sample_m
    return(output)
  }

  module_relationship <- all_lop_par[[2]]
  all_net <- pblapply(1:k,function(c)get_output(module_relationship[[c]],all_lop_par[[1]][[c]],t(cluster_mean),times=times))
  
  
  get_net_output <- function(j){
    get_after <- function(i){
      temp <- matrix(NA,nrow = length(i[[2]]),ncol=3)
      temp[,1] <- i[[2]]
      temp[,2] <- i[[1]]
      temp[,3] <- i[[5]][2:(length(i[[2]])+1)]
      colnames(temp) <- c('from','to','dep_effect')
      temp <- data.frame(temp)
      temp[,3] <- as.numeric(as.character(temp[,3]))
      return(temp)
    }
    links <- do.call(rbind,lapply(j, get_after))
    
    get_link_color <- function(i){
      tmp <- links$dep_effect[i]
      if (tmp >= 0 ) {
        tmp2 <- '+1'
      } else {
        tmp2 <- '-1'
      }
      return(tmp2)
    }
    links$cor <- sapply(1:nrow(links),function(c)get_link_color(c))
    
    get_ind <- function(i){
      temp <- i[[5]][1]
      return(temp)
    }
    nodes <- data.frame(unique(links[,2]),paste0('M',1:k),sapply(j,get_ind))
    colnames(nodes) <- c("id","name","ind_effect")
    nodes$size <- table(cluster_result$clustered_height$cluster)
    nodes$influence <- aggregate(dep_effect ~ to, data = links, sum)[order(
      as.numeric(aggregate(dep_effect ~ to, data = links, sum)[,1])),2]
    
    links$dep_effect <- abs(links$dep_effect)
    links$from <- paste0('M',links$from)
    links$to <- paste0('M',links$to)
    
    return(list(links,nodes[,c(2,4,5)]))
  }
  ODE_result <- get_net_output(all_net)
  return(list(ODE_result,all_net))
}

net_result <- get_module_result(k=k,times,cluster_result = biFunClu_results)

