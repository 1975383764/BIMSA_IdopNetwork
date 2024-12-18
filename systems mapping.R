library(MASS)
library(mvtnorm)
library(deSolve)
library(ggplot2)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(reshape2)
library(patchwork)
pheno_df <- read.csv("new30.F_2015-HT-DIA.csv",header = T)
height <- pheno_df[,1:30]
diameter <- pheno_df[,31:60]
times <- seq(1,11,len = 30)

Lorenz <- function(t,state, parameters)
{
  with( as.list(c(state, parameters)),
        {
          dE <- a1*E*(1-(E/k1)^a1)
          
          list(c(dE))
        }
  )
}

legendre <- function(xn,t,par)
{
  y <- t*(par[1]+par[2]*xn+par[3]*(0.5*(3*xn^2-1))+par[4]*(0.5*(5*xn^3-3*xn)))+par[5]
  return(y)
}

get_miu = function(miu_par,state,xn,t)
{
  parameters <- c(a1 = miu_par[1],
                  k1 = miu_par[2]
  )
  par <- c( miu_par[3],
            miu_par[4],
            miu_par[5],
            miu_par[6],
            miu_par[7]
  )
  
  miu <- as.numeric(ode(y = state, times = t, func = Lorenz, parms = parameters, method = 'rk4')[,2])
  miu1 <- as.numeric(legendre(xn=xn,t=t,par=par))
  return (miu+miu1);
}
fn <- function (par, pheno,state,xn,t){
  y <- c(get_miu(par,xn, state = state,t = t));
  return(sum((pheno-y)^2))}

get_value = function(miu_par,state,xn,t)
{
  parameters <- c(a1 = miu_par[1],
                  k1 = miu_par[2]
  )
  par <- c( miu_par[3],
            miu_par[4],
            miu_par[5],
            miu_par[6],
            miu_par[7]
  )
  
  miu <- as.numeric(ode(y = state, times = t, func = Lorenz, parms = parameters, method = 'rk4')[,2])
  miu1 <- as.numeric(legendre(xn=xn,t=t,par=par))
  return (miu+miu1);
}

i=73
t <- c(seq(1,11,len = 30))
pheno <- pheno_df[i,]
x <- pheno[31:60]
u <- -1;
v <- 1;
xmin <- min(x);
xmax <- max(x);
xn  <- u + ((v - u) * (x - xmin)) / (xmax - xmin);
xn <- as.numeric(xn)

x1 <- pheno[1:30]
u <- -1;
v <- 1;
xmin <- min(x1);
xmax <- max(x1);
xn1  <- u + ((v - u) * (x1 - xmin)) / (xmax - xmin);
xn1 <- as.numeric(xn1)

sample_m <- as.numeric(pheno[1:30])
sample_m1 <- as.numeric(pheno[31:60])
H0_state <- c(E=as.numeric(pheno[1]))
H0_state1 <- c(E=as.numeric(pheno[31]))
lower <- c(0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf) 
upper <- c(Inf, Inf, Inf, Inf, Inf, Inf, Inf)  
#######73
LS_init_par <- c(1.1,  100,0.077859573,	-0.013871414,	-0.014060134,	0.000677726,	-0.078544293)#0.8730265,  max(pheno[1:11]),   7.6368376, -12.9871982,   3.2353543,   0.6795738, -24.1510614
LS_init_par1 <- c(0.3,  5, -0.40210699, -0.09103915,  0.16905234,  0.04101178,  0.20989313)
LS  = optim(LS_init_par, fn = fn, pheno = sample_m, state = H0_state,xn=xn,t=t, method = "BFGS")
LS1  = optim(LS_init_par1, fn = fn, pheno = sample_m1, state = H0_state1, xn=xn1,t=t, method = "Nelder-Mead")
#######60
LS_init_par <- c(0.4108854, 166.1587706,   0.3191281 , -0.1629520,  -0.2384959 , -0.3425897 , -0.2871979)
LS_init_par1 <- c(0.20279362, 13.57443357,  0.05399541, -0.99686189,  0.18596807, -0.96544587, -0.87030851)
LS  = optim(LS_init_par, fn = fn, pheno = sample_m, state = H0_state,xn=xn,t=t, method = "Nelder-Mead")
LS1  = optim(LS_init_par1, fn = fn, pheno = sample_m1, state = H0_state1, xn=xn1,t=t, method = "Nelder-Mead")




