datagenerate11 = function(p, n=NULL, scale=2, sigma=0, cor=0.8, m=1){
  # n: sample size
  # p: dimensionality
  s = ifelse(sigma==0,ceiling(p/10), ceiling(sqrt(p)))
  if(is.null(n)) n = ifelse(sigma==0, round(scale*(s*log(p/s))), round(scale*(s*log(p-s))))
  
  K = round(s/2); Sigma = matrix(cor, nrow=m, ncol=m); diag(Sigma)=1
  x = matrix(rnorm(n*(p-K*m)), nrow=n)
  x = cbind(do.call(cbind, lapply(1:K, function(i) mvrnorm(n, rep(0,m), Sigma=Sigma) )), x)
  x = stand(x)$X
  
  idx = cumsum(c(1, rep(m,K), rep(1,s-K-1)))
  b = sample(c(-3,3), s, replace=T)
  w = rnorm(n); w=w-mean(w)
  y = x[,idx]%*%b + sigma*w
  g = rep(0,p); g[1:(K*m)]=rep(1:K, each=m)
  return(list(y=y, x=x, g=g, b=b, S=idx, s=s, w=w))
}
data_scaling = function( datat, datav=NULL, xnames=NULL, idx_valid=NULL, reorder=FALSE, compressed=F){
  # datav: specific testing data
  # idx_valid: specify index of datat to create testing data
  if(!is.null(idx_valid)){
    datav = list(x=datat$x[ idx_valid,,drop=F], y=datat$y[ idx_valid])
    datat = list(x=datat$x[-idx_valid,,drop=F], y=datat$y[-idx_valid])
  }
  if(is.null(xnames)) xnames=names(datat$x)
  stan=stand(datat$x)
  x_center=stan$center; x_scale=stan$scale; y_center=mean(datat$y)
  datat = list(x=stan$X, y=datat$y-y_center)
  
  if(!is.null(datav)){
    datav = list(x=sapply(1:ncol(datav$x), function(i) (datav$x[,i]-x_center[i])/x_scale[i] ), y=datav$y-y_center)
    if(!is.matrix(datav$x)) datav$x=matrix(datav$x,nrow=1)
  }
  rm = which(datat$x[1,] %>%is.nan)
  if(length(rm)>0){
    datat$x=datat$x[,-rm]
    datav$x=datav$x[,-rm]
  }
  if(reorder){
    od = order(abs(cor(datat$x, datat$y)), decreasing=T)
    xnames = xnames[od]
    datat$x = datat$x[,od]
    if(!is.null(datav)){datav$x = datav$x[,od]}
    od_inv = sapply(1:ncol(datav$x), function(i) which(od==i))
    return(list(datat=datat, datav=datav, xnames=xnames, x_center=x_center, x_scale=x_scale, y_center=y_center, idx_valid=idx_valid, od=od, od_inv=od_inv))
  }
  return(list(datat=datat, datav=datav, xnames=xnames, x_center=x_center, x_scale=x_scale, y_center=y_center, idx_valid=idx_valid))
}
lam1lam2 = function(lam=NULL, lam1=NULL, lam1min=NULL, lam1max=NULL, lam2=NULL, nlam1=10, nlam2=10, l2_fac=1, cutend=T){
  if(is.null(lam1min)) lam1min=0.001
  if(!is.null(lam)){
    return(lam)
  }
  # else
  if(is.null(lam1)){ 
    if(lam1max<=lam1min){
      lam1s=lam1min; nlam1=1; nlam2=1
    }else{
      lam1s = exp(seq(from=log(lam1min), to=log(lam1max), length.out=nlam1))
      if(cutend) lam1s[nlam1] = (lam1s[nlam1]+lam1s[nlam1-1])/2
    }
  }else{ # if lam1s are given as a sequence of lambda_1
    nlam1 = length(lam1); lam1s = lam1
  }
  if(is.null(lam2)){
    lam2s = exp((-1):(nlam2-2))*l2_fac; lam2s[1]=0
    lam_ls = lapply(lam1s, function(lam1) cbind(lam1,lam1*lam2s) )
    lam = NULL
  }else{
    nlam2 = length(lam2); lam2s = lam2
    lam_ls = lapply(lam1s, function(lam1) cbind(lam1,lam2s) )
    lam = NULL
  }
  lam=do.call(rbind, lam_ls)
  return(lam)
}
aweights = function(obj){
  b = obj$b; g = obj$g; med = obj$medoids
  w = sapply(1:length(b), function(i){ if(g[i]==0){1/abs(b[i])}else{1/abs(b[med[g[i]]])}})
  return(w)
}
complete = function(self,dist=NULL){
  name_lst = names(self); for(name in name_lst){assign(name, self[[name]])}
  if(is.null(w)){w=rep(1,length(b))}
  if(!is.null(g)){
    if(is.null(self$objective)) self$objective=Objective(self=self,dist=dist)$obj
    c_ind = g[b!=0 & g!=0] # index of important clusters
    b_ind = which(b!=0); names(b_ind)=NULL
    b_val = b[b!=0]
  }
  output = list();  for(name in name_lst) output[[name]] = get(name)
  return(output)
}
stand = function(x, scale=0){
  n = nrow(x); p = ncol(x)
  n = n - scale
  center = apply(x, 2, mean)
  X = sapply(1:p, function(i) (x[,i]-center[i]))
  scale = sqrt(apply(X^2, 2, sum))/sqrt(n)
  X = sapply(1:p, function(i) X[,i]/scale[i])
  output = list(X=X, center=center, scale=scale)
  return(output)
}
mydist = function(X, Y=NULL, n=NULL){
  if(is.null(n)){n=nrow(X)};if(is.null(Y)){Y=X}
  rho = round(abs(t(X)%*%Y)/n,10)
  return(1-rho)
}
Soft = function(a, b){return(sign(a)*pmax(abs(a)-b,0))}
checkrepresent = function(self){
  c1 = all(self$b[self$medoids[self$g[self$b!=0 & self$g!=0]]]!=0) # important clusters' medoids are nonzero
  c2 = all(table(self$g[self$g!=0 & self$b!=0])==1) # representatives are selected
  return(c1&c2)
}

# objective: lst sq/2 + lam1*(l1-norm + c*pairwise) + lam2*PAM
# type: objective value of 1. "RS" = representative selection, 2. "SPPAM" = SPPAM clustering
Objective = function(self, dist=NULL, type="all", sum=T, XY=T){
  if(!(type %in% c("all", "RS", "SPPAM"))){print("wrong type");return(0)}
  name_lst = names(self); for(name in name_lst){assign(name, self[[name]])}
  
  ######### function starts here
  if(!sum){type="SPPAM"}
  
  ## SPPAM ignores loss
  ## loss = 0.5||Y-Xb||_2^2
  loss = ifelse(type=="SPPAM", 0, 0.5*sum((Y-X%*%b)^2)) 
  
  ## SPPAM ignores l1_penalty
  ## L1-norm = ||b||_1
  L1_penalty = ifelse(type=="SPPAM", 0, lam1*sum(abs(w[w<Inf]*b[w<Inf])))
  
  ## both SPPAM and CD need this term
  ## sum_{i \neq j} |b[i]b[j]| = 2*sum_{i < j} |b[i]b[j]| = ||b||_1^2 - ||b^2||_1
  pairwise_penalty = lam1*cc*sapply(1:max(g),function(k) sum(abs(b[g==k]))^2- sum(b[g==k]^2)) 
  
  ## RS ignores within cluster variation
  PAM = if(type=="RS"){0}else{
    wc_var = if(is.null(dist)){
      sapply(1:max(g),function(k) sum( mydist(X[,medoids[k],drop=F],X[,g==k,drop=F])))
    }else{
      sapply(1:max(g),function(k) sum( dist[medoids[k],g==k]) )
    }
    c(lam2*(1-rho.min)*sum(g==0),lam2*wc_var )}
  # 
  objective = if(sum){
    loss+L1_penalty+sum(pairwise_penalty)+sum(PAM)
  }else{
    # C0, C1,...,CK
    c(0,pairwise_penalty)+PAM
  }
  ######### function ends here
  output = list();  for(name in name_lst){output[[name]] = get(name)}
  if(!XY) output$X=NULL; output$Y=NULL
  return(output)
}
sort_cluster=function(self){
  ## sort clusters according to the order of medoids
  od = order(self$medoids)
  rk = rank(self$medoids)
  self$g[self$g!=0] = rk[self$g]
  self$medoids = self$medoids[od]
  return(self)
}

#################################
#####  spare estimation  ########
#################################
## Local quadratic approximation (LQA)
LQA_grad=function(self, dist, ct.max=10, A=NULL, tol=1e-7, debugging=F, learning_rate = 0.8){
  # Y, X, b, w, g, lam1, cc
  ######### function starts here
  ct = 0
  obj = Objective(self=self,dist=dist,type="RS")$obj
  lam1 = self$lam1
  Y = self$Y
  cc = self$cc
  repeat{
    # A2 movable, A1 don't move
    obj0 = obj
    self0 = self
    if(is.null(A)){A = abs(self$b)>1e-7}
    
    b=self$b[A]; w=self$w[A]; X=self$X[,A, drop=F]; g=self$g[A]; 
    b[b==0] = tol
    
    if(ncol(X)==1){break;}
    temp2=lam1*w
    if(max(self$g)>0){
      for(k in setdiff(unique(g),0)){temp2[g==k] = (w[g==k]*lam1+2*lam1*cc*(sum(abs(b[g==k]))-abs(b[g==k])))}
    }
    grad2 = -t(X)%*%(Y-X%*%b) + temp2*sign(b)
    hessian2 = t(X)%*%X + diag(as.vector(temp2/abs(b)))
    step2 = try(-solve(hessian2, grad2), silent=T )
    
    if(!is.matrix(step2)){break;}
    self$b[A] = self$b[A] + step2*learning_rate
    obj = Objective(self=self,dist=dist,type="RS")$obj
    if(ct>ct.max){break;}
    if(obj-obj0>=0){break;}
    if(obj0/obj-1<tol){break}
    ct=ct+1
  }
  return(self0)
}
UpdateBetas = function(self, dist, alpha=0.1, tol=1e-7, ct.max=10, boost=c(1,2,3), iter.max=10,pt=F){
  # boost1: LQA fit
  # boost2: lars fit
  # boost3: switch representatives
  # Y, X, b, w, g, lam1, lam2
  
  ######### function starts here
  exclude = which(self$w==Inf)
  obj0 = Objective(self=self,dist=dist)$obj
  obj_min = Inf
  self$converge = F
  
  ### proposed loop: coordinate-descend + boost1 + boost2
  ct=0
  repeat{
    ct=ct+1
    ### Do LQA
    if(1 %in% boost) self = LQA_grad(self=self,dist=dist,ct.max=ct.max,tol=tol)
    
    # join active dimensions back
    if(ct==0){
      join = setdiff(1:length(self$b), exclude)
    }else{
      temp2=self$lam1*self$w
      Xb = A=which(self$b!=0); Xb = self$X[,A,drop=F]%*%self$b[A]
      if(max(self$g)>0){for(k in setdiff(unique(self$g),0)){temp2[self$g==k] = (self$w[self$g==k]*self$lam1+2*self$lam1*self$cc*(sum(abs(self$b[self$g==k]))-abs(self$b[self$g==k])))}}
      a = temp2; b = abs(t(self$X)%*%(self$Y-Xb)); join = union(which(a<b-tol & self$b==0),self$medoids)
    }
    Ck = intersect(which((self$g==0)),join)
    if(length(Ck)>0){
      Xw = sapply(Ck, function(j) self$X[,j]/self$w[j])
      lars = lars(x=Xw,y=self$Y-self$X[,-Ck,drop=F]%*%self$b[-Ck],type="lasso",use.Gram = F,normalize=F )
      self$b[Ck] = predict(lars, s = self$lam1, type = "coefficients", mode="lambda")$coefficients/self$w[Ck]
    }
    
    Xb = A=which(self$b!=0); Xb = self$X[,A,drop=F]%*%self$b[A]
    for(k in 1:max(self$g)){
      A = which(self$b!=0)
      Ck = intersect(which((self$g==k)), union(A,join))
      rCk = Ck[sample(length(Ck), length(Ck))]
      for(j in rCk){
        Xj <- self$X[,j, drop=F]; s = sum(Xj^2)
        yj <- self$Y - (Xb - Xj*self$b[j]) ## partial residual
        p2 = sum(abs(self$b[Ck]))-abs(self$b[j])
        bj.old = self$b[j]
        penalty = (self$w[j]*self$lam1+2*self$cc*self$lam1*p2); # if(k!=0){penalty=penalty/sum(self$b[Ck]!=0)}
        step = as.numeric(Soft( sum(Xj*yj), penalty ))/s - bj.old
        self$b[j] = bj.old + sign(step)*min(abs(step), alpha); Xb = Xb + Xj*(self$b[j]-bj.old)
      }
      if(sum(self$b[Ck]!=0)==1){ # this step resets obj0 if medoid_k is switched
        candi = which(self$g==k)
        dist_candi = if(is.null(dist)){mydist(self$X[,candi,drop=F])}else{dist[candi, candi, drop=F]}
        wc_var = apply(dist_candi, 2, sum)
        Xk=self$X[,candi,drop=F];wk=self$w[candi]
        Yk=self$Y-(Xb - self$X[,Ck,drop=F]%*%self$b[Ck])
        bk=Soft(t(Xk)%*%Yk, self$lam1*wk)/apply(Xk^2, 2, sum)
        obj_mix = self$lam2*wc_var + sapply(1:sum(self$g==k), function(j) 0.5*sum((Yk-Xk[,j]*bk[j])^2))+ ifelse(bk==0,0,self$lam1*wk*abs(bk))
        opt = obj_mix %>% which.min()
        if(self$medoids[k] != candi[opt]){self$medoids[k]=candi[opt];obj0=Inf}
        j.old = which(self$b[candi]!=0); bj.old=self$b[candi[j.old]]
        self$b[candi[j.old]]=0; Xb = Xb + Xk[,j.old]*(0-bj.old)
        self$b[candi[opt]]=bk[opt]; Xb = Xb + Xk[,opt]*(bk[opt]-0)
      };# print(Objective(self=self,dist=dist)$obj)
    }
    
    if(2 %in% boost){ ## weighted exact lasso
      A_join = which(self$b!=0)
      if(length(A_join)>1){
        # w = a[A_join]
        w = sapply(A_join, function(j){output=if(self$g[j]!=0){self$w[j]*self$lam1+2*self$cc*self$lam1*(sum(abs(self$b[self$g==self$g[j]]))-abs(self$b[j]))}else{self$w[j]*self$lam1}});names(w)=NULL
        Xw = sapply(1:length(A_join), function(j) self$X[,A_join[j]]/w[j])
        lars = lars(x=Xw,y=self$Y,type="lasso",use.Gram = F,normalize=F )
        self$b[A_join] = predict(lars, s = 1, type = "coefficients", mode="lambda")$coefficients/w
      }
    }
    obj = Objective(self=self,dist=dist)$obj
    if(obj<obj_min){
      obj_min=obj
      solution = list(b=self$b,med=self$medoids)
    }
    # because opt_mix, objective value could increase
    if( obj0/obj-1<tol ){self$converge=T;break}
    obj0 = obj; 
    if(ct>iter.max){break}
  }
  self$b=solution$b; self$medoids=solution$med
  return(self) 
}

################
## clustering ##
################
myPPAM = function(self, dist=NULL, package=T, PPAM=NULL, nstart=1){
  p = ncol(self$X)
  if(is.null(PPAM)){
    if(is.null(dist)){dist=mydist(self$X)}
    # BUILD phase. medoid p+1 is fixed all the time
    # Using supervised distance
    dist = rbind(cbind(dist,1-self$rho.min),1-self$rho.min); dist[p+1,p+1] = 0
    
    val = Inf
    for(i in 1:nstart){
      while(T){
        PPAM_try = fastpam(as.dist(dist),n=p+1,k=self$K+1,seed=i) # build-in PAM algorithm
        if(p%in%PPAM_try@medoids){
          if(PPAM_try@cost<val){PPAM = PPAM_try;val=PPAM@cost}
          break
        }
      }
    }
  }
  
  if(class(PPAM)=='list'){
    return(PPAM)
  }else if(class(PPAM)=="KmedoidsResult"){
    self$medoids = PPAM@medoids+1
    self$g = PPAM@assignment
    
    k_ex = self$g[p+1]
    self = sort_cluster(self)
    
    self$g[self$g==(self$K+1)]=0
    self$g=self$g[1:p]
    self$medoids=self$medoids[1:self$K]
    return(self)
  }
}

mySPPAM = function(self, dist, nstart=1, iter.max=20){
  # fixed parameters
  K = self$K
  rho.min = self$rho.min
  b = self$b
  lam1cclam2 = self$lam1*self$cc/self$lam2
  self$converge = F
  
  ct=0
  A = which(b!=0); Ac = which(b==0)
  while(ct<iter.max){
    ct = ct+1;# print(ct)
    # iterative parameters
    self0=self
    g = self$g
    med = self$medoids
    
    ## update cluster label:
    dist_med = if(is.null(dist)){ mydist(self$X[,med,drop=F],self$X) }else{dist[med,,drop=F]}
    dist0 = rbind(1-rho.min, dist_med) # distance matrix with C0
    # simultaneously update unimportant variables' cluster labels. Because these updates do not affect cluster penalties
    g[Ac] = apply(dist0[, Ac], 2, which.min)-1
    sumb = sapply(1:K, function(k) sum(abs(b[g==k])))
    
    
    # sequentially update important variables' cluster labels. (Because these updates "DO" affect cluster penalties)
    # the closer j to medoid, the earlier j joins cluster
    A_order = A[ order(apply(dist_med[,A,drop=F], 2, min)) ]
    # A_order = sample(A)
    for(j in rev(A_order)){ 
      k = g[j] # old cluster label
      dist_penalty = sumb*abs(b[j]); if(k!=0){dist_penalty[k] = dist_penalty[k]-b[j]^2}
      dist_penalized = dist0[,j] + c(0, ifelse(dist_penalty==0, 0, dist_penalty*2*lam1cclam2)) # prevent lam1cclam2=Inf
      g[j] = which.min(dist_penalized)-1
      
      # update sumb
      if(k!=0) sumb[k] = sumb[k]-abs(b[j])
      if(g[j]!=0) sumb[g[j]] = sumb[g[j]]+abs(b[j])
    }
    self$medoids=med;self$g=g
    ## update medoid: need to watchout empty cluster (members might be dropped in previous steps)
    for(k in 1:K){ 
      if(sum(g==k)<=1){ # if empty or one-member cluster 
        # remove single-member-cluster
        idx_single = which(g==k); g[idx_single] = 0 
        # extract a new cluster from C0
        # the new cluster could include nonzero beta. Therefore, need to update its medoid right after
        candi = setdiff(which(g==0),idx_single)
        dist_candi = if(is.null(dist)){mydist(self$X[,candi,drop=F])}else{dist[candi, candi, drop=F]}
        wc_var = apply(dist_candi, 2, function(d){
          new_cluster = (d<1-self$rho.min)
          return(sum(d[new_cluster]) + sum(!new_cluster)*(1-self$rho.min))
        })
        opt = wc_var %>% which.min()
        med[k] = candi[opt]
        g[candi[dist_candi[opt,]<1-self$rho.min]] = k
      }
   
      candi = which(g==k)
      dist_candi = if(is.null(dist)){mydist(self$X[,candi,drop=F])}else{dist[candi, candi, drop=F]}
      wc_var = apply(dist_candi, 2, sum)
      if(all(self$b[candi]!=0)){
        med[k] = candi[wc_var %>% which.min()]
      }else{ # boost 3
        Xk=self$X[,candi,drop=F];wk=self$w[candi]
        Yk=self$Y-self$X[,-candi,drop=F]%*%self$b[-candi]
        bk=Soft(t(Xk)%*%Yk, self$lam1*wk)/apply(Xk^2, 2, sum)
        obj_mix = self$lam2*wc_var + sapply(1:sum(g==k), function(j) 0.5*sum((Yk-Xk[,j]*bk[j])^2)) + ifelse(bk==0,0,self$lam1*wk*abs(bk))
        opt = obj_mix %>% which.min()
        med[k] = candi[opt]
        self$b[candi[self$b[candi]!=0]]=0
        self$b[candi[opt]]=bk[opt]
      }
    }
    if(all(self$g==self0$g) & all(self$med==self0$med)){self$converge=T;return(self)}
    self$medoids=med;self$g=g
  }
  return(self)
}

################ algorithm 1 in section 4.5
algo = function(self, dist=NULL, ct=0, cc_init=1e-5, nstart=1, ct.max=100, alpha=0.1, tol=1e-7, bhat=NULL, g_given = NULL, boost=c(1,2,3), nonstop=T, seed){
  cat('solving lam1=', round(self$lam1,6), ', lam2=', round(self$lam2,6),': ', sep='')
  n = nrow(self$X); p = ncol(self$X); if(self$cc==0){self$cc=cc_init*exp(ct)}
  b_ini=self$b; g_ini=self$g; med_ini=self$medoids
  set.seed(seed,sample.kind = "Rejection")
  if(is.null(g_given)) self=mySPPAM(self, dist=dist, nstart=nstart)
  if(length(g_given)==p){ self$g=g_given; self$medoids = sapply(1:max(self$g), function(k) which(self$g==k)[1]) }
  repeat{
    if(all(self$w==Inf) | self$lam2==0){self=Objective(self=self,dist=dist,XY=F);t_g=t_b=rep(0,3);break}
    # step 1: increase c
    ct = ct+1; cat(ct,', ', sep='')
    self$cc = self$cc*exp(1)
    repeat{
      b0 = self$b; g0=self$g
      
      # step 2: update beta
      tic=proc.time(); self=UpdateBetas(self=self,dist=dist,alpha=alpha,tol=tol,ct.max=ct.max,boost=boost); toc=proc.time()
      converge_beta=self$converge
      if(exists("t_b")){t_b = t_b+toc-tic}else{t_b=toc-tic}
      
      # step 3: update g
      if(is.null(g_given)){ 
        tic=proc.time(); self=mySPPAM(self, dist=dist, nstart=nstart); toc=proc.time()
        converge_cluster=self$converge
      }else{tic=proc.time(); toc=proc.time(); converge_cluster=T}
      if(exists("t_g")){t_g = t_g+toc-tic}else{t_g=toc-tic}
      
      if( all(self$g==g0) & max(abs(self$b - b0))<tol ) break
      if( nonstop ) break
    }
    medoids_unimp = self$medoids[setdiff(1:self$K, unique(self$g[self$b!=0]))]
    c1 = if(length(medoids_unimp)==0){T}else{all(mydist(self$X[,self$b!=0 & self$g==0], self$X[,medoids_unimp])>1-self$rho.min)}
    if(length(g_given)>0) c1=T
    # if( checkrepresent(self) & c1 & converge_beta & converge_cluster){
    if( (checkrepresent(self) & c1)){
      A = which(self$b!=0 | self$g==0) # representatives and independents
      Xw = sapply(A, function(j)self$X[,j]/self$w[j])
      lars = lars(x=Xw,y=self$Y,type="lasso",use.Gram = F,normalize=F )
      self$b[A] = predict(lars, s = self$lam1, type = "coefficients", mode="lambda")$coefficients/self$w[A]
      cyclic = if(exists('self_hist')){apply(self_hist,2,function(o){all(o$g==self$g)&(max(abs(o$b-self$b)))<tol})}else{F}
      self_hist=if(exists('self_hist')){cbind(self_hist,Objective(self=self,dist=dist,XY=F))}else{cbind(Objective(self=self,dist=dist,XY=F))}
      if(any(cyclic)){  # check cyclic loop
        min = which.min(apply(self_hist,2, function(o) o$obj))
        self = self_hist[,min]
        break
      }
    }
    if(ct>50){stop('ct exceed limit')}
  }
  cat('\n user: g=', round(t_g[1],3),', b=',round(t_b[1],3), ', clock: g=', round(t_g[3],3),', b=',round(t_b[3],3), '\n', sep='')
  output=complete(self=self,dist=dist); output$iter=ct; output$timecost = list(t_g=t_g, t_b=t_b)
  return(output)
}
################################
############# script ###########
################################

RepreSelect = function(datatv=NULL, datat=NULL, datav=NULL, rho.min, K, w=NULL, path=NULL, 
                       lam=NULL, lam1=NULL, lam1min=NULL, lam1max=NULL, lam2=NULL, nlam1=10, nlam2=10, l2_fac=1, break_lam2=T, full_lam2=F, lam_start=1, lam_end=NULL,
                       dist=T, PPAM=NULL, lars=NULL, # preprocesses
                       subsample=NULL, nfold=NULL, # for cross-validation
                       g_given=NULL, nstart=1,  # clustering parameters: nstart: clustering, 
                       ct.max=10, alpha=0.1, tol=1e-5, once=T, boost=c(1,2,3), # RS parameters: ct.max=iter_num_LQA, alpha=max_step
                       parallel=0, path_function=NULL, filename=NULL){ # parallel computing parameters
  bhat=NULL # if bhat is given, adaptive weights variate through iteration
  ## append: True -> append #clusters when Update_g
  ## updatebeta: True -> if not append #clusters, mergecluster then Update_beta
  ## data transformation
  # obj atttributes:
  # X: standard X
  # Y: standard Y
  # dist: input; TRUE=compute and save; FALSE=compute everytime
  # lam: fixed grids
  # lam1: para
  # lam2: para
  # b: regression coefficients
  # g: cluster label
  # w: adaptive weights
  # medoids: cluster medoids
  # objective: objective value
  # important: important cluster index ( effective variables)
  # unimportant: unimportant cluster index ( effecitve variables)
  # b_important: value of important b (representatives)
  # b_unimportant: value of unimportant b (independent b)
  # cluster_cor: correlation between representatives
  if(parallel & is.null(path_function)){stop('need to give function library path')}
  
  #### preprocesses start ####
  ## Training / Validation(Testing) data
  if(!is.null(datatv)){datat=datatv$datat;datav=datatv$datav}
  if(is.null(datav)){datav = datat}
  datatv = data_scaling(datat=datat,datav=datav, reorder=F)
  X = datatv$datat$x; Y = datatv$datat$y
  n = nrow(X); p = ncol(X)
  
  ## Parameter setup
  if(is.null(w)){w=rep(1,p)}; w[w==0]=1e-10
  
  # tuning paramemter assignments
  if(is.null(lam1max)){lam1max=max(abs(t(X)%*%Y));cutend=T}else{cutend=F} 
  if(is.null(lam1min)){lam1min=0.001*n}
  lam = lam1lam2(lam=lam, lam1=lam1, lam1min=lam1min, lam1max=lam1max, lam2=lam2, nlam1=nlam1, nlam2=nlam2,l2_fac=l2_fac,cutend=cutend)
  if(is.null(lam_end)) lam_end = nrow(lam)
  
  # distance matrix. This could be very RAM comsuming if p is large
  if(isTRUE(dist)){dist = mydist(X)}else if(is.matrix(dist)){if( abs(dist[1,2]-mydist(X[,1:2])[1,2])>1e-10){stop('dist does not match')}}
  # Penalized Partition Around Medoids (PPAM)
  if(is.null(PPAM) & length(g_given)<=1 ){ # large computation
    self = list(X=X, K=K, rho.min=rho.min)
    PPAM = myPPAM(self=self, dist=dist, nstart=1)
  }
  
  # After initializing PPAM, update adaptive weights
  if(!is.null(bhat)){
    w_ini = sapply(1:length(b), function(i){if(PPAM$g[i]==0){1/abs(bhat[i])}else{min(1/abs(bhat[PPAM$g==PPAM$g[i]]))}})
    if(!is.null(w)){if( !all(w==w_ini) ){stop('the given w does not match (PPAM + bhat)')}}
    w = w_ini
  }
  
  # exact adaptive lasso solution
  if(is.null(lars)){
    #glmnet and gglasso use decreasing lambda
    Xw = sapply(1:ncol(X), function(j) X[,j]/w[j]) ## ||Y-Xb||^2_2 + ||wb||_1 = ||Y-(X/w)wb||^2_2 + ||wb||_1
    lars = lars(x=Xw,y=Y,type="lasso",use.Gram = F,normalize=F)
    # lam1s = sort(unique(unlist(lam[,1])),decreasing=T)
    # solpath = glmnet(x=Xw,y=Y,lambda=lam1s/n,alpha=0.1,intercept=F,standardize=F,thresh=1e-12)$beta %>% apply(2, function(b)b/w)
    
    # g_temp = PPAM$g; g_temp[g_temp==0]= (K+1:sum(g_temp==0))
    # solpath = gglasso(x=Xw,y=Y,group=g_temp,lambda=lam1s/n,intercept=F,eps=1e-16)$beta %>% apply(2, function(b)b/w)
  }else{
    if( any(w!=lars$w)){stop('lars input does not match w')}
    lars = lars$lars
  }
  #### preprocesses end ####
  
  run = function(lamj,seed=1){
    lam1 = lam[lamj, 1]; lam2 = lam[lamj,2]
    b=if(all(w==Inf)){rep(0,length(w))}else{b = predict(lars, s = lam1, type = "coefficients", mode="lambda")$coefficients/w}
    self = list(X=X, Y=Y,
                lam1 = lam1, cc = 0, lam2 = lam2,
                K = K, rho.min = rho.min,
                w = w, b = as.numeric(b), bmax=as.numeric(b),
                g = PPAM$g, medoids = PPAM$medoids,
                objective = NULL,
                c_ind = NULL, b_ind = NULL, b_val = NULL)
    while(T){
      self_try = try({
        algo(self=self,dist=dist,nstart=nstart,ct.max=ct.max,alpha=alpha,
             tol=tol,bhat=bhat,g_given=g_given,boost=boost,seed=seed)
      },silent=F)
      if(class(self_try)!='try-error'){self=self_try;break}
      print('algo has error')
      seed=seed+1
    }
    return(self)
  }
  ## grid parameter solutions
  if(!is.null(path)){
    lam = t(matrix(unlist(path[c(1,3),]),nrow=2))
  }else if(parallel){
    lam1s = sort(unique(lam[,1]))
    idx_lam1 = sapply(lam1s, function(lam1) min(which( abs(lam[,1]-lam1)<1e-10)) )
    path = foreach(start=idx_lam1, .packages=c('MASS', 'lars', 'cluster', 'magrittr'), .combine=cbind, .verbose = T, .inorder=F )%dopar%{
    # path_lst = lapply(idx_lam1, function(start){
      source(path_function)
      path=NULL; next_lam2=T; o=list(lam1=0); lamj=start
      for(lamj in lam_start:lam_end){ 
        if( abs(lam[lamj,1] - lam[start,1])>1e-10){break}
        print(lamj)
        if(next_lam2 | o$lam1<lam[lamj,1] ){
          path = cbind(path, o<-run(lamj))
        }else if(full_lam2){
          o$lam2=lam[lamj,2];o$X=X;o$Y=Y; o=Objective(self=o,dist=dist,XY=F)
          path=cbind(path, o)
        }
        if(break_lam2){ next_lam2 = (any(o$medoids!=PPAM$medoids)|any(o$g!=PPAM$g)) }
      }
      return(path)
    }
    print('finish')
  }else{
    path=NULL; next_lam2=T; o=list(lam1=0); count=lam_start-1
    for(lamj in lam_start:lam_end){
      load=F; print(lamj)
      if(!next_lam2 & abs(o$lam1-lam[lamj,1])<1e-10 & !full_lam2){next}
      count=count+1; 
      if(!is.null(filename)){ 
        path_file=paste0(filename, '_', count,'.Rdata')
        if(file.exists(path_file)){
          o=get(load(path_file)); 
          if(abs(o$lam1-lam[lamj,1])>1e-10 | abs(o$lam2-lam[lamj,2])>1e-10){stop('saved file inconsistent. wrong full_lam2?')}
          load=T
        }
      }
      if(!load){
        if(!next_lam2 & abs(o$lam1-lam[lamj,1])<1e-10){
          o$lam2=lam[lamj,2];o$X=X;o$Y=Y; o=Objective(self=o,dist=dist,XY=F)    
        }else{
          o = run(lamj)
        }
        if(!is.null(filename)) save(o, file=paste0(filename, '_', count,'.Rdata'))
      }
      path=cbind(path, o)
      if(break_lam2){ next_lam2 = (any(o$medoids!=PPAM$medoids)|any(o$g!=PPAM$g)) }
    }
  }
  colnames(path)=NULL
  MSEt0 = sum(datatv$datat$y^2)
  MSEts = apply(path, 2, function(o) sum((datatv$datat$y-datatv$datat$x[,o$b!=0,drop=F]%*%o$b[o$b!=0])^2) ) # within-sample MSE
  
  MSEv0 = sum(datatv$datav$y^2)
  MSEvs = apply(path, 2, function(o) sum((datatv$datav$y-datatv$datav$x[,o$b!=0,drop=F]%*%o$b[o$b!=0])^2) ) # out-of=sample MSE
  badtune = apply(path, 2, function(o){
    medoids_imp = o$medoids[unique(o$g[o$b!=0])]
    return(any(mydist(datatv$datat$x[,medoids_imp,drop=F], datatv$datat$x[,o$g==0 & o$b!=0,drop=F])<1-o$rho.min))
  })
  badtune2 = apply(path, 2, function(o){
    return(any(mydist(datatv$datat$x[,o$g==0 & o$b!=0,drop=F], datatv$datat$x[,o$g==0 & o$b!=0,drop=F])<1-o$rho.min))
  })
  MSEt3s = MSEt2s = MSEts; MSEt2s[badtune]=Inf; MSEt3s[badtune2]=Inf
  
  output = list(X=X,Y=Y,Xv=datatv$datav$x,Yv=datatv$datav$y,path=path, MSEt0=MSEt0, MSEts=MSEts, MSEt2s=MSEt2s, MSEt3s=MSEt3s, MSEv0=MSEv0, MSEvs=MSEvs)
  if(!is.null(subsample)) nfold=length(subsample)
  if(!is.null(nfold)){
    if(is.null(subsample)) subsample = fold_sample(n=n, fold=nfold, seed=1)
    lam = t(matrix(unlist(path[c(1,3),]),nrow=2))
    CARS_CV = sapply(1:nfold, function(k){
      print(paste('fold', k))
      filename_CV = if(!is.null(filename)){paste0(filename,'_cv',k)}else{NULL}
      datatv_CV=data_scaling(datat=datat, idx_valid=subsample[[k]] )
      dist=ifelse(is.matrix(dist),T,NULL) 
      CARS_k = RepreSelect(datatv=datatv_CV, rho.min=rho.min, K=K, w=w, lam=lam,
                         full_lam2=T, dist=dist, g_given=g_given, nstart=nstart,
                         ct.max=ct.max, alpha=alpha, tol=tol, once=once, boost=boost,
                         parallel=parallel, path_function=path_function, filename=filename_CV)
      return(CARS_k)
    })
    output$MSEcvs = apply(apply(CARS_CV, 2, function(CARS_k) CARS_k$MSEvs), 1, sum)
    output$MSEcv0 = apply(CARS_CV, 2, function(CARS_k) sum(CARS_k$Yv^2)) %>% sum
    output$CrossValidate = CARS_CV
  }
  return(output)
}

path_output = function(SIS, K, rho.min, model, job=NULL){ 
  output = 'HTCresult/'
  char_job = if(is.null(job)){NULL}else{paste0('_job',job)}
  return(paste0(output, model, '_SIS',SIS,'_K',K,'_r',rho.min,char_job,'_',data_name))
}

mypath = function(prefix, par_lst, suffix='.Rdata', job=NULL){ 
  # prefix is for folder/object
  # par_lst has to be a list
  output = NULL
  if(!is.null(prefix)) output = paste0(output, prefix)
  for(name in names(par_lst)) output = paste0(output, '_',name,par_lst[[name]])
  if(!is.null(job)) output = paste0(output, '_job',job)
  if(!is.null(suffix)) output = paste0(output, '_',suffix)
  return(output)
}


cv.vs = function(datatv=NULL, datat=NULL, path=NULL, datav=NULL, nlam1=100, w=NULL, lam1s=NULL, lam1max=NULL, njob=NULL, # parameters
                 tune=T, subsample=NULL, fold=NULL, model='lasso', dftype=1, tol=1e-7, # tuning info
                 g=NULL, clustering='PAM', K=NULL, rho.min=NULL,  # clustering
                 MSEt=1, # only applied if clustering exists
                 parallel=0, path_function=NULL, pathonly=F, nopath=F, show=T){ # summary
  # model: 1. lasso 2. glasso 3. crl 4. SCAD 5. MCP
  # model: lasso2 use ncvreg (cd algorithm) to approach lasso
  # method: 1. CV 2. BIC 3. AIC 4. GCV
  # clustering: 1. PAM 2. PPAM 3. CARS-PAM
  if(model%in%c('glasso','crl')){
    if( is.null(K) & is.null(g)){
      print('need to specify K or g');return(0)
    }else if(clustering=='PPAM' & is.null(rho.min)){
      print('need to specify rho.min');return(0)
    }
  }else if((model=='CARS' | model=='aCARS') & is.null(njob) & is.null(path)){stop('need njob')}
  # check parallel
  if(parallel & is.null(path_function)){print('need to give function library path'); return(0)}
  
  ## scale datat and datav
  if(!is.null(datatv)){datat=datatv$datat;datav=datatv$datav}
  if(is.null(datav)){datav = datat}
  datatv = data_scaling(datat=datat,datav=datav, reorder=F)
  X = datatv$datat$x; Y = datatv$datat$y
  n = nrow(X); p = ncol(X)
  
  ## assign tuning parameter
  if(is.null(lam1max)) lam1max = max(abs(t(Y)%*%X))
  if(is.null(lam1s)){lam1s=lam1lam2(lam1min=0.001*n, lam1max=lam1max, nlam1=nlam1, nlam2=1, l2_fac=1, cutend=F)[,1]}
  
  
  # training: compute path anyway
  if(show) cat('computing path...')
  if(is.null(path)){
    if(model%in%c('alasso', 'alasso2')){
      package=if(model=='alasso'){'lasso'}else{'lasso2'}
      lassos = vs(datat=datatv$datat,lam1s=lam1s, g=g,mode=package,tol=tol)
      path = NULL
      for(i in 1:length(lam1s)){
        print(i)
        o = lassos[,i]
        w = 1/abs(o$b)
        if(all(w==Inf)){
          pathi = sapply(lam1s, function(lam1) list(lam1=lam1,g=o$g,b=o$b))
        }else{
          pathi = vs(datat=datatv$datat,lam1s=lam1s, g=g,mode=package,w=w,tol=tol)    
        }
        pathi = sapply(1:ncol(pathi), function(j) append(list(lam1_0=o$lam1, w=w, b_0=o$b), pathi[,j])  )
        path = cbind(path, pathi)
      }
    }else{
      path = vs(datat=datatv$datat,lam1s=lam1s, g=g,w=w,mode=model,K=K,rho.min=rho.min,clustering=clustering,tol=tol)
    }
  }
  if(pathonly) return(path)
  
  if(show) cat('computing MSEt...')
  MSEt0 = sum(datatv$datat$y^2)
  MSEts = apply(path, 2, function(o) sum((datatv$datat$y-datatv$datat$x[,o$b!=0,drop=F]%*%o$b[o$b!=0])^2) )
  
  badtune = apply(path, 2, function(o){
    medoids_imp = o$medoids[unique(o$g[o$b!=0])]
    return(any(mydist(datatv$datat$x[,medoids_imp,drop=F], datatv$datat$x[,o$g==0 & o$b!=0,drop=F])<1-o$rho.min))
  })
  badtune2 = apply(path, 2, function(o){
    return(any(mydist(datatv$datat$x[,o$g!=0 & o$b!=0,drop=F], datatv$datat$x[,o$g!=0 & o$b!=0,drop=F])<1-o$rho.min))
  })
  MSEt3s = MSEt2s = MSEts; MSEt2s[badtune]=Inf; MSEt3s[badtune2]=Inf
  
  
  if(!is.null(fold) | !is.null(subsample)){
    if(is.null(subsample)) subsample = fold_sample(n=nrow(X), fold=fold)
    cat(fold,'-fold cross-validation: ',sep='')
    MSEcv0=0
    for(k in 1:fold){
      cat(k,' ')
      temp = data_scaling(datat, idx_valid = subsample[[k]])
      result = cv.vs(datatv=temp, lam1s=lam1s, g=g,w=w,model=model,K=K,rho.min=rho.min,clustering=clustering,tol=tol,show=show,tune=F)
      MSEcv0 = MSEcv0+result$MSEv0
      if(exists('MSEcvs')){MSEcvs=MSEcvs+result$MSEvs}else{MSEcvs=result$MSEvs}
    }
    cat('\n'); rm(result)
  }else{
    MSEcvs=NULL; MSEcv0=NULL
  }
  training = list(X=datatv$datat$x, Y=datatv$datat$y, path=path, MSEts=MSEts, MSEt2s=MSEt2s, MSEt3s=MSEt3s, MSEt0=MSEt0, MSEcvs=MSEcvs, MSEcv0=MSEcv0)
  
  ## testing data
  if(show) cat('computing MSEv...')
  if(!is.null(datav)){
    MSEv0 = sum(datatv$datav$y^2)
    MSEvs = apply(path, 2, function(o) sum((datatv$datav$y-datatv$datav$x[,o$b!=0,drop=F]%*%o$b[o$b!=0])^2) )
    testing = list(Xv=datatv$datav$x, Yv=datatv$datav$y, MSEvs=MSEvs, MSEv0=MSEv0)
  }else{testing = NULL}
  output = append(training, testing)
  rm(path, training, testing)
  
  ## tuning of hyper parameters: 
  if(MSEt==2 & !grepl( 'CARS', model, fixed = TRUE)){stop('MSEt=2 is only for CARS/aCARS')}
  if(tune){output = summary.vs(output, MSEt=MSEt, show=show)}
  if(nopath) output$path=NULL
  return(output)
}

summary.vs = function(output, n=NULL, dftype=1, p=NULL, fac_HBIC2=2, Rsq_min=0, gamma1=1, gamma2=1, MSEt=1, show=T, hist=F, MSEplot=F, main=NULL, path=NULL){
  if(gamma1>1) stop('gamma1 in (0,1)')
  if(gamma2<1) stop('gamma2 > 1')
  if(is.null(n)) n=nrow(output$X)
  if(is.null(p)) p=if(is.null(output$path)){ncol(output$X)}else{length(output$path[,1]$g)}
  if(is.null(output$df1s)){
    if(show) print('recompute df1s...')
    output$df1s = apply(output$path, 2, function(o){
      if(is.null(o$w)) o$w=rep(1,length(o$b))
      o$b = round(o$b,7)
      A = (o$b!=0)
      if(sum(A)==0){return(0)}
      if(sum(A)==1){return( n/ (n+o$lam1*o$w[A]/abs(o$b[A])) )}
      XtX_A = t(output$X[,A,drop=F])%*%output$X[,A,drop=F]
      ( solve( XtX_A + diag(o$lam1*o$w[A]/abs(o$b[A]))) %*% XtX_A ) %>% diag %>% sum
      # check (solve( XtX_A + diag(o$lam1*o$w[A]/abs(o$b[A]))) %*% t(X[,A]) %*% Y) - o$b[A]
    })
  }
  if(is.null(output$df2s)){
    if(show) print('recompute df2s...')
    output$df2s = apply(output$path, 2, function(o) sum(o$b!=0) )
  }
  
  dfs = if(dftype==1){output$df1s}else{output$df2s}
  MSEts = if(MSEt==1){output$MSEts}else if(MSEt==2){output$MSEt2s}else if(MSEt==3){output$MSEt3s}
  # var = MSEts/(n-dfs-1)
  var = 1
  # AIC  = n*log(MSEts/n)+2*dfs*var
  AIC  = n*log(MSEts) + dfs*var*2
  BIC  = n*log(MSEts) + dfs*var*log(n)
  EBIC1 = MSEts + dfs*log(n)+2*gamma1*log(p) # reviewed in 2011
  EBIC2 = log(MSEts/n)*n + dfs*(log(n)+2*gamma1*log(p)) # reviewed in 2011
  HBIC1 = log(MSEts/n)*n + dfs*2*gamma2*log(p) # Tao Wang and Lixing Zhu 2011
  HBIC2 = log(MSEts/n)*n + dfs*fac_HBIC2*log(log(n))*log(p) # Lan Wang, Yongdai Kim, and Runze Li 
  
  # HBIC1[1-MSEts/output$MSEt0<Rsq_min]=Inf
  # HBIC2[1-MSEts/output$MSEt0<Rsq_min]=Inf
  
  OPT = output$MSEvs; OPT[MSEts==Inf] = Inf
  if(!is.null(output$MSEcvs)){CV = output$MSEcvs; MSEcv0=output$MSEcv0}else{MSEcv0=NULL}
  
  for(method in c('CV', 'AIC', 'BIC', 'HBIC1', 'HBIC2', 'OPT')){
    if(!exists(method)){next}
    loss = get(method)
    opt = which.min(loss)
    MSEt = MSEts[opt]
    MSEv = output$MSEvs[opt]
    df = dfs[opt]
    if(!is.null(path)){
      obj = append(path[, opt], list(loss=loss, opt=opt, MSEt=MSEt, MSEv=MSEv, df=df, MSEcv0=MSEcv0))
    }else if(!is.null(output$path)){
      obj = append(output$path[, opt], list(loss=loss, opt=opt, MSEt=MSEt, MSEv=MSEv, df=df, MSEcv0=MSEcv0))
    }else{
      if(opt==output[[method]]$opt){
        obj = output[[method]]
        obj$loss=loss; obj$MSEt=MSEt; obj$MSEv=MSEv; obj$MSEcv0=MSEcv0
      }else{
        obj = list(loss=loss, opt=opt, MSEt=MSEt, MSEv=MSEv, df=df, MSEcv0=MSEcv0)
        print('obj has not estimates')
      }
    }
    output[[method]]=obj
    
    if(show){
      cat('Selection by ', method,": ", sep='')
      if(!is.null(obj$lam1_0)) cat("lam1_0=", round(obj$lam1_0,3),', ', sep='')
      if(!is.null(obj$lam2_0)) cat("lam2_0=", round(obj$lam2_0,3),', ', sep='')
      if(!is.null(obj$lam1))   cat("lam1=",   round(obj$lam1,3),', ', sep='')
      if(!is.null(obj$lam2))   cat("lam2=",   round(obj$lam2,3),', ', sep='')
      cat('Rsq=',       round(1-obj$MSEt/output$MSEt0,3),
          ', Rsq_oos=', round(1-obj$MSEv/output$MSEv0,3),sep='')
      if(method=='CV') cat(', Rsq_CV=', round(1-min(loss)/MSEcv0,3),sep='')
      cat('\n')
    }
    if(method %in% hist){
      abs(cor(datat$x[,obj$b!=0])) %>% hist(ylim=c(0,20), xlim=c(0.5,1), 
              main=paste0('Hist. of corr of selected variables by using ', method))
    }
  }
  if(MSEplot){
    vss = vs(datat=datat,lam1s=lam1s,mode=model,K=K)
    MSEts = apply(vss, 2, function(vs) sum((Y - X%*%vs$b)^2))
    yaxis = c(MSEts, MSE_CV)
    plot(lam1s, MSEts, 'l',ylim=c(0, max(yaxis)), xlab='lambdba', ylab='MSE', main=main )
    points(lam1s, MSE_CV, 'l', col=2)
    text(output$lam1, min(MSE_CV), 'optimal lambda')
    legend('bottomright', c('MSE train', "MSE CV"), col=c(1,2), lty=1)
  }
  return(output)
}


fold_sample = function(n, fold=NULL, seed=1){
  if(is.null(fold)){fold=n}
  if(!is.null(seed)){set.seed(seed,sample.kind = "Rejection")}
  if(fold!=n){
    ord = sample(1:n, n)
    mod = n%/%fold
    remainder = n%%fold
    subsample = list()
    for(k in 1:fold){subsample[[k]] = ord[1:mod + (k-1)*mod]}
    if(remainder>0){for(r in 1:remainder){subsample[[r]] = c(subsample[[r]], ord[mod*k+r])}}
    for(k in 1:fold){subsample[[k]] = sort(subsample[[k]])}
  }else{
    subsample = lapply(1:n, function(i) return(i))
  }
  return(subsample)
}

vs = function(datat, lam1s, w=NULL, g=NULL, PPAM=NULL, model=NULL, K=NULL, rho.min=NULL, clustering=NULL, tol=1e-7){
  stan = stand(datat$x)
  X=stan$X
  Y=datat$y-mean(datat$y)
  if(model %in% c('lasso', 'MCP', "SCAD", 'lasso2')){ # noncluster-based method
    g = rep(0, ncol(X))
    g_ori=g
  }else if( model %in% c('glasso', 'crl') ){ # cluster-based methods
    if(clustering == 'PAM'){
      if(is.null(g)){
        tic = proc.time()
        kmedoids <- cluster::pam(1-abs(cor(X)), K, diss=TRUE) # create k-medoids clustering with K clusters
        t_g=proc.time()-tic
        g <- kmedoids$cluster
      }
      g_ori=g
    }else if(clustering =='PPAM'){
      self = list(X=X, Y=Y, lam1 = 0, cc = 0, lam2 = 1, K = K, rho.min = rho.min, w = 0, b = 0, bmax=0, g = 0, medoids = NULL, objective = 0, c_ind = NULL, b_ind = NULL, b_val = NULL)
      tic = proc.time()
      self = myPPAM(self=self, PPAM=PPAM, nstart=5)
      t_g=proc.time()-tic
      g = self$g
      C0 = which(g==0)
      g_ori=g; g[C0] = 1:length(C0)+K
    }
  }
  
  ## special X manipulation for crl
  if(model=='crl'){ # cluster representative lasso
    # create representatives by taking mean of each group
    X0=X; g0=g
    X = sapply(1:max(g), function(k){ 
      if(k>K){
        return( X[,g==k])
      }else{
        return( apply(X[,g==k,drop=F],1,mean) )
      }
    })
    # X = sapply(1:max(g), function(k) apply(X[,g==k,drop=F],1,mean))
    g = 1:max(g)
  }
  
  if(model=='glasso'){
    b = gglasso(x=X, y=Y, group=g, lambda=lam1s/nrow(X))$beta
    b = b[,length(lam1s):1, drop=F] # make it p*nlam
    # stop('check this modification')
    b = apply(b, 2, function(bb){
      for(k in 1:max(g)){
        if(any(bb[g==k]!=0)){
          bb[g==k & bb==0]=1e-10
        }
      }
      return(bb)
    })
  }else if(model=='lasso' | model=='crl'){
    if(model=='lasso' & any(w!=1)){
      Xw = sapply(1:ncol(X), function(j) X[,j]/w[j])
      lars = lars(x=Xw,y=Y,type="lasso",use.Gram = F,normalize=F)
    }else{
      lars = lars(x=X,y=Y,type="lasso",use.Gram = F,normalize=F)
    }
    b = predict(lars, s = lam1s, type = "coefficients", mode="lambda")$coefficients
    if(length(lam1s)==1){b = cbind(b)}else{b=t(b)}
    if(!is.null(w)) b = sapply(1:ncol(b), function(j) b[,j]/w)
    if( model=='crl'){
      b = apply(b, 2, function(blam1){
        bb = rep(0, ncol(X0))  
        for(k in 1:max(g0)){ bb[g0==k] = blam1[k]/sum(g0==k)}
        return(bb)
      })
    }
  }else if(model %in% c('SCAD','MCP','lasso2')){
    if(model=='lasso2') model='lasso'
    b = matrix(0, nrow=ncol(X), ncol=length(lam1s))
    b[,length(lam1s)] = ncvfit(X=X, y=Y, penalty=model, lambda=lam1s[length(lam1s)]/nrow(X), eps=tol, max.iter=as.integer(1/1e-7) )$beta
    if(length(lam1s)>1){
      for(i in (length(lam1s)-1):1){
        b[,i] = ncvfit(X=X, y=Y, penalty=model, lambda=lam1s[i]/nrow(X), eps=tol, max.iter=as.integer(1/1e-7) , init=b[,i+1] )$beta
      }
    }
  }
  path = sapply(1:length(lam1s), function(i)return(list(lam1=lam1s[i], g=g_ori, b=b[,i])))
  return(path)
}

###################################
########## summary tools ##########
###################################
upcor = function(X, idx){cor= round(cor(X[,idx]),3);return(cor[upper.tri(cor)])}
table_evaluation = function(simu, g=NULL, A=NULL, R=NULL, xvar=NULL, b=NULL, nv=NULL, eva_func = evaluation4, table_type=NULL, tune='OPT'){
  # table_type: 1. Lasso 2. CARS 3. without 4. with
  rerun = 1:ncol(simu)
  if(is.null(nv)) nv = nrow(simu[,1]$Xv)
  for(i in rerun){
    if(!is.null(simu[,i]$b0)){b=simu[,i]$b0; g=simu[,i]$g0; A=simu[,i]$A0; xvar=simu[,i]$xvar}
    if(i==1){table = eva_func(b=b, g=g, A=A, xvar=xvar, vs = simu[,i][[tune]], nv=nv); next}
    table = table + eva_func(b=b, g=g, A=A, xvar=xvar, vs = simu[,i][[tune]], nv=nv)
  }
  table = table/length(rerun)
  table[,1:2] = round(table[,1:2]*100, 2)
  table[,3:ncol(table)] = round(table[,3:ncol(table)], 3)
  
  tab_name = colnames(table)
  # if('PS' %in% tab_name) tab_name[which(tab_name=='PS')] = '|A_hat|'
  # if('TPS' %in% tab_name) tab_name[which(tab_name=='TPS')] = '|A_hat & A|'
  # colnames(table)=tab_name
  col_idx = tab_name %in% c('BIAS','RMSE')
  table[2:nrow(table),col_idx]=NA
  if(table_type=='lasso'){
    n_row = nrow(table); n_col = ncol(table)
    table=apply(table, 2, sum)
    table=table[-2:0+n_col]
    table[3]=table[3]/n_row
  }else{
    table[nrow(table), 3:5]=NA
  }
  print(table)
  return(table)
}
evaluation4 = function(b, g, A, xvar, vs, nv, method=c("TPC", "FPC", "ERS", "WRS", "SRS", "PS", "TPS", "BIAS", "RMSE")){
  if( !all(c('b','g') %in% names(vs))){stop('obj needs b and g in list')}
  MSEv = vs$MSEv
  ghat = vs$g; K = max(ghat)
  bhat = vs$b; Ahat = which(bhat!=0)
  
  if(!is.null(b)){
    R = A
    if(is.null(xvar)) stop('R and xvar are needed')
    V = sapply(R, function(j){
      if(g[j]==0){return(j)}
      idx = which.min(abs(b[j]/sqrt(xvar[g==g[j]])-bhat[g==g[j]]))
      return(which(g==g[j])[idx])
    })
    
    b2 = rep(0,length(b)); b2[V] = b[R]/sqrt(xvar[V])
    BIAS = sqrt(sum((b2-bhat)^2))
  }
  
  
  ## deciding corresponding cluster labels
  index = sapply(1:max(g), function(k1){
    Ck = which(g==k1)
    gamma=sapply(1:max(ghat), function(k2){ 
      Ckhat = which(ghat==k2)
      return(length(intersect(Ckhat, Ck))/length(Ck))
    })
    return(which.max(gamma))
  })
  method_cluster = c('TPC', 'FPC', 'ERS', 'WRS', 'SRS')
  # cbind(1:length(g), bhat, ghat)[bhat!=0, ]
  output = sapply(sort(unique(g)), function(k){
    Ck = which(g==k)
    Ckc = which(g!=k)
    Ckhat = if(k==0){which(ghat==k)}else{which(ghat==index[k])}
    Ckchat = if(k==0){which(ghat!=k)}else{which(ghat!=index[k])}
    
    TPC = length(intersect(Ck, Ckhat))/length(Ck)
    FPC = if(length(Ckhat)==0){0}else{length(setdiff(Ckhat, Ck)) / length(Ckhat)}
    ERS = length(Reduce(intersect, list(Ahat, Ckhat, Ck)))
    WRS = length(Reduce(intersect, list(Ahat, Ckhat, Ckc)))
    SRS = length(Reduce(intersect, list(Ahat, Ckchat, Ck)))
    
    #PS = length(Reduce(intersect, list(Ahat, Ckhat, Ck)))
    #TPS = length(Reduce(intersect, list(Ahat, A, Ckhat, Ck)))
    PS = length(Reduce(intersect, list(Ahat, Ck)))
    TPS = length(Reduce(intersect, list(Ahat, A, Ck)))
    
    PSALL = length(Ahat)
    TPSALL = length(intersect(Ahat, A))
    RMSE = sqrt(MSEv/nv)
    output = c()
    for(m in method){output=c(output, get(m))}
    return(output)
  })
  if(length(method)>1){
    output = cbind(output[,-1], output[,1])
    colnames(output) = c(paste0('C', 1:max(g)), 'C0')
    rownames(output)=method
    output = t(output)
  }else{
    output = c(output[-1], output[1])
  }
  return(output)
}

#######################
## summary functions ##
#######################
tsneplot2 = function(CARS, cv_modified = F){
  if(cv_modified){MSE = CARS$MSEv2}else{MSE=CARS$MSEv}
  o = CARS$path[, which.min(MSE)]  
  g = o$g
  dist = 1-abs(cor(CARS$Xv))
  obj_min = o$rho.min*sum(g==0)+sum(sapply(1:max(g),function(k) sum(dist[medoids[k],g==k]))) %>% round(2)
  MSE_min = sqrt(min(MSE)/nrow(CARS$Xv)) %>% round(2)
  tsne = tsneplot(tsne=tsne, g=o$g+1, b=o$b, corner='topright', main=paste0("K=", o$K,', rho_min=', o$rho.min, ', MSE=', MSE_min, ', SWCV=', obj_min))
}


recovery_plot = function(cor, m, sigma, p_lst=c(128,256,512), scale_lst=(1:12)/4, runs=1:20, 
                         model='lasso', lam1s=NULL, nlam1=100, plot_only=T){
  if(!is.null(lam1s)) nlam1=length(lam1s)
  grids = expand.grid(p_lst, scale_lst, runs)
  rec = foreach(i=1:nrow(grids), .packages=c('MASS', 'lars', 'magrittr'), .combine=cbind)%dopar%{
    source('functions_ver28.R')
    p=grids[i,1]; scale=grids[i,2]; run=grids[i,3]
    set.seed(run)
    datat = datagenerate11(p=p, scale=scale, cor=cor, m=m, sigma=sigma)
    g = if(model=='lasso'){NULL}else if(model=='CARS'){g=datat$g}
    exist = opt_checker(X=datat$x, Y=datat$y, S=datat$S, g=g, lam1s=lam1s, nlam1=nlam1)
    return(c(p, scale, exist))
  }
  rec = apply(expand.grid(p_lst, scale_lst), 1, function(grids){
    p=grids[1]; scale=grids[2]; prob=rec[3, rec[1,]==p & rec[2,]==scale]
    return(c(p,scale, prob))
  })
  plot(NULL, ylim=c(0,1), xlim=c(0,3), xlab='Rescaled sample size', main=paste0(model, ': m=',m,', r=',cor,', sigma=',sigma))  
  for(j in 1:length(p_lst) ){rec_sub = rec[,rec[1,]==p_lst[j],drop=F]; points( rec_sub[2,], apply(rec_sub[-(1:2),,drop=F], 2, mean), col=j, 'l' )}
  if(!plot_only) return(rec)
}

opt_checker = function(X, Y, S, g=NULL, lam1s=NULL, nlam1=100){
  if(is.null(g)) g = rep(0, ncol(X))
  if(!is.null(lam1s)){nlam1=length(lam1s)}else{
    lam1max = max(abs(t(X)%*%Y))
    lam1s = exp(seq(from=log(nrow(X)/1e3), to=log(lam1max), length.out=nlam1))
  }
  
  # stage 1: check C0
  lars = lars(x=X[,S],y=Y,type="lasso",use.Gram = F,normalize=F )
  # check any false positive in C0
  C0_sparsity = sapply(lam1s, function(lam1){
    b = predict(lars, s = lam1, type = "coefficients", mode="lambda")$coefficients
    obj = sum((Y - X[,S]%*%b)^2)/2 + lam1*sum(abs(b))
    opt = all(abs(t(X[,g==0])%*%(Y - X[,S]%*%b))<lam1+1e-5) & all(b!=0)
    return(c(obj, opt))
  })
  
  if(!any(C0_sparsity[2,]==1)){ return(FALSE) }
  if(max(g)==0) return(TRUE)
  
  # stage 2: check exchangeability
  lam1s = lam1s[C0_sparsity[2,]==1]; C0_sparsity = C0_sparsity[,C0_sparsity[2,]==1,drop=F]
  for(k in 1:max(g)){
    Ck = which(g==k); if(length(Ck)==1){next}
    for(j in Ck){
      S2 = sort(c(j, setdiff(S, Ck))); if(setequal(S, S2)){next}
      lars = lars(x=X[,S2],y=Y,type="lasso",use.Gram = F,normalize=F )
      # check any false positive in C0 after exchanging
      C0_sparsity2 = sapply(lam1s, function(lam1){
        b = predict(lars, s = lam1, type = "coefficients", mode="lambda")$coefficients
        obj = sum((Y - X[,S2]%*%b)^2)/2 + lam1*sum(abs(b))
        opt = all(abs(t(X[,g==0])%*%(Y - X[,S2]%*%b))<lam1+1e-5) & all(b!=0)
        return(c(obj, opt))
      })
      RS = C0_sparsity2[2,]==1
      if( all( !RS ) ){next}
      if( all( C0_sparsity[1,RS] > C0_sparsity2[1,RS]) ){return(FALSE)}
      # if( all( C0_sparsity[1,RS]>C0_sparsity2[1,RS]) ){cat('k=',k, ', j=',j, '\n');return(FALSE)}
    }
  }
  return(TRUE)
}
