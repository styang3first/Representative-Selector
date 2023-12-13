path_root = dirname(rstudioapi::getActiveDocumentContext()$path); setwd(path_root)

library(lars)
library(MASS)
library(magrittr)
library(doParallel)
library(foreach)
library(doRNG)
library(iterators)
num_cluster = 10
if(!exists('cl')){cl = makeCluster(num_cluster); registerDoParallel(cl)}
# stopCluster(cl)
source('functions_ver28.R')

### Noiseless
cor=0.8; sigma=0; m=10
# recovery_plot(cor=cor, m=m, sigma=sigma, runs=20, lam1s = 1e-5)

###########
## Noisy ##
###########
par(mfrow=c(2, 2))
sigma=1; m=10
## Lasso
lasso = lapply( ((1:4)*2-1) / 10, function(cor) recovery_plot(cor=cor, m=m, sigma=sigma, runs=1:50, model='lasso', plot_only=F))
CARS = lapply( ((1:4)*2-1) / 10, function(cor) recovery_plot(cor=cor, m=m, sigma=sigma, runs=1:50, model='CARS', plot_only=F))

# lapply( 0.1, function(cor) recovery_plot(cor=cor, m=m, sigma=sigma, runs=1:50, model='lasso', plot_only=F))



##############
# source('functions_ver28.R')
# p=512; scale=3; cor=0.8; sigma=1; m=10
# datat = datagenerate11(p=p, scale=scale, cor=cor, m=m, sigma=sigma)
# 
# lam1s=NULL; nlam1=100
# # lasso = cv.vs(datat=datat, lam1s=lam1s, nlam1=nlam1, show=F)
# # CARS = RepreSelect(datat=datat, lam1=lam1s, nlam1=nlam1, lam2=0.1, nlam2=1, g_given=datat$g, K=max(datat$g), rho.min=0.5)
# opt_checker(X=datat$x, Y=datat$y, S=datat$S)
# opt_checker(X=datat$x, Y=datat$y, S=datat$S, g=datat$g)
