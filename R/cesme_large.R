#' R package: cesme.
#'
#' High-dimensional Cluster Analysis with Latent Semiparametric Mixture Models.
#'
#' @param y Response.
#' @param K Number of clusters.
#' @param thre Threshold for estimating transformation. If thre is NULL, use thre=1/(4*n^0.25*sqrt(pi*log(n))).
#' @param N Number of iterates. 
#' @param C0 A constant for penalty parameter of SCAD. 
#' @param C2 A constant for penalty parameter of SCAD. 
#' @param kap A constant for penalty parameter of SCAD. 
#' @param z.initial Initial clustering assignment. If z.initial is NULL, use spectral clustering to obtain an initialization.
#' @param mu.initial Initial clustering assignment. If mu.initial is NULL, use sample cluster means with z.initial to obtain an initialization.
#' @param maxit_msda Maximum number of iterates in SCAD solver.
#' @param eps_msda Precision of the penalized regression solver.
#' @param ncvx Whether to use the nonconvex penalty (SCAD).
#' @param sanity.check For debug use.
#' @param lam_max_print For testing the code: print the smallest lambda such that all coefficients are zero.
#' @param g.method Which method to use to estimate the monotone function.
#' @param kmnstart Number of initial starts in spectral clustering.
#' 
#' @details
#' An algorithm for high-dimensional model-based clustering via alternating maximizations. Ues bigmemory to store sample covariance matrix of data. 
#' @return
#' A list consisting of
#'   \item{z.iter}{estimate of z}
#'   \item{mu.iter}{estimate of mu}
#'   \item{sigma.iter}{estimate of sigma}
#'   \item{g.iter}{estimate of g}
#'   \item{bic1.iter}{bic values}
#'   \item{bic2.iter}{bic values}
#'   \item{time.iter}{run time}
#' @keywords clustering
#' @importFrom stats cov kmeans model.matrix qnorm sd
#' @import Matrix
#' @useDynLib cesme, .registration=TRUE
#' @export
#'
#-----------------------------------------------------------------------------

##cesme with SCAD
#BIC twice
npn.clust.bic.large = function(y, K, thre=NULL, N=100, C0=1, C2=1, kap=0.8, z.initial=NULL, mu.initial=NULL,
                         maxit_msda=as.integer(1e+3), eps_msda=1e-4, ncvx=TRUE,
                         sanity.check=FALSE, lam_max_print=FALSE, g.method=1, kmnstart=100) {
  
  # g.method = 1 empirical cdf
  # g.method = 2 shaped restricted regression
  # g.method = 3 hybrid estimation
  # g.method = 4 trivial identical mapping
  
  ## set parameters
  N = as.integer(N)
  K = as.integer(K)
  p = dim(y)[2]
  n.smp = dim(y)[1]
  n.inv = 1 / n.smp
  ## ten-step ECM
  S = 10L
  #twenty rho's
  #Crho = c(seq(0.2,1.1,0.3),seq(1.4,2,0.2),seq(2.1,3,0.1))
  Crho = seq(0.1,2,0.1)
  L = length(Crho)
  
  ## setting for msda
  nk = K - 1L
  vnames = paste("V", seq(p), sep = "")
  ## maxit_msda = as.integer(1e+06)
  maxit_msda = as.integer(maxit_msda)
  ## pf = rep(1, p)
  ## eps_msda = as.double(1e-04)
  eps_msda = as.double(eps_msda)
  dfmax = as.integer(n.smp)
  pmax = as.integer(min(dfmax * 2 + 20, p))
  
  ## ## start points
  z.iter = array(0L, dim=c(N, n.smp, L))
  z2.iter = array(0L, dim=c(S, n.smp, L))
  mu.iter = array(0, dim=c(K, p, S+1))
  sigma.iter = list()
  for(i in 1:(S+1)){ 
    sigma.iter = c(sigma.iter, big.matrix(p, p, type = "short", init = 0))
  }
  g.iter = array(0, dim=c(n.smp, p, 2))
  pro.iter = matrix(0, ncol=K, nrow=S+1)
  
  bic1.iter = matrix(Inf, ncol=S, nrow=N)
  bic2.iter = matrix(Inf, ncol=L, nrow=N)
  
  
  time.iter = matrix(0, 2, N)
  
  ## sparse
  gamma.iter = array(0, dim=c(p, K-1, N*L))
  gamma2.iter = array(0, dim=c(p, K-1, L*S))
  rho.iter = c(C0 + C2 * sqrt(log(p)*n.inv), rep(0,N))
  
  if(is.null(z.initial)){
    ## initialized by SC
    z.iter[1,,1] = as.integer(initial.sc(y, K, rho.iter[1], nstart=kmnstart))
  }else{
    z.iter[1,,1] = z.initial
  }
  
  if ((length(unique(z.iter[1,,1])) < K) || (min(table(z.iter[1,,1])) < 2L)) {
    miss_z = setdiff(seq(K), unique(z.iter[1,,1]))
    large_z = which(z.iter[1,,1]==which.max(table(z.iter[1,,1])))
    for (k in miss_z) {
      z.iter[1,large_z[k],1] = k
      z.iter[1,large_z[k+K],1] = k
    }
  }
  
  pro.iter[1,] = as.numeric(table(z.iter[1,,1]) * n.inv)
  
  if (is.null(thre)) {
    thre <- 1 / (4 * n.smp^0.25 * sqrt(pi * log(n.smp)))
  }
  delta0 = thre
  ## 1/(4*(n.smp)^0.25 * sqrt(pi * log(n.smp)))
  
  alpha.temp = array(0, dim=c(L+1, n.smp, K))
  
  if(g.method == 1){
    ## estimate monotone function empirically
    tt1 = Sys.time()
    g.iter[,,1] = g.est.Rversion(y=y, z=z.iter[1,,1], n=n.smp, p=p, K=K, 
                        pro0=pro.iter[1,], delta0=delta0)
    tt2 = Sys.time()
    if (sanity.check) print("new_g.iter_cdf")
    if (sanity.check) print(ttt <- tt2 - tt1)
    if (sanity.check){ if( ttt > 1) { 
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    }}
  }else if(g.method ==4){
    g.iter[,,1] = y
  }
  
  
  
  
  ############################################################################
  # if (sanity.check) {
  #   ## print("check")
  #   tt1 = Sys.time()
  #   tmp = sapply(c(1:p), g.est.old, y=y, z=z.iter[1,,1], 
  #     K=K, pro0=pro.iter[1,], delta0=delta0)
  #   tt2 = Sys.time()
  #   print("old_g.iter")
  #   print(tt2 - tt1)
  #   print("------------------------------------")
  #   if(max(abs(tmp - g.iter[,,1])) > 1.0e-5) stop("g.est")
  # }
  ############################################################################
  ### scale it
  g.iter[,,1] = t((t(g.iter[,,1])-apply(t(g.iter[,,1]),1,mean))/(apply(t(g.iter[,,1]),1,sd)+1e-4))
  
  ## ## starting value for mu
  center2.g = g.iter[,,1]
  if(is.null((mu.initial))){
    for (k in seq(K)) {
      if (sum(z.iter[1,,1]==k) == 1L) {
        mu.iter[k,,1] = g.iter[z.iter[1,,1] == k,,1]
      } else {
        mu.iter[k,,1] = colMeans(g.iter[z.iter[1,,1] == k,,1])
        center2.g[z.iter[1,,1] == k, ] = apply(center2.g[z.iter[1,,1] == k, ], 
                                               2, scale, scale=FALSE)
      }
    }
  }else{
    mu.iter[,,1] = mu.initial
  }
  
  ## starting value for variance
  sigma.iter[[1]] = sigma.est.fn.large(mu=mu.iter[seq(K),,1], 
                                       g=g.iter[,,1], alpha=matrix(1/K,n.smp,K), 
                                       n=n.smp, p=p, K=K, n.inv=n.inv)
  ## apply(simplify2array(sapply(c(1:n.smp), 
  ##   function(index){return((center2.g[index,]%*%t(center2.g[index,]))/n.smp)}, 
  ##   simplify = F)), c(1:2), sum)
  
  
  ### start iteration
  for (iter in seq(N - 1L)) {
    ## if (iter %% 5 == 0)
    
    if (sanity.check) { 
      print("iter")    
      print(iter)  
    }  
    l = 1
    #ten-step ECM
    
    i = 1
    
    bic3.iter = matrix(Inf, ncol=L, nrow=S)
    l.iter = rep(0,S)
    
    for(i in seq(S)){
      ### E-step
      
      ############################################################################
      ## reduce the number of computation of ginv
      ############################################################################
      
      if (sanity.check) {
        print("i")
        print(i)
      }
      
      time1 = Sys.time()
      
      #det_sigma = determinant(sigma.iter[,,i],logarithm=TRUE)$modulus[1]
      #eig_sigma = eigen(sigma.iter[,,i])
      #eig_sigma_value = eig_sigma$values
      #eig_sigma_vector = eig_sigma$vector
      #det_sigma = sum(log(eig_sigma_value[which(eig_sigma_value>1e-6)]))
      #ginv_sigma = eig_sigma$vector[,which(eig_sigma_value>1e-6)]%*%
      #  diag(1/eig_sigma_value[which(eig_sigma_value>1e-6)])%*%
      #  t(eig_sigma$vector[,which(eig_sigma_value>1e-6)])
      
      
      
      #########################################################################
      ## apply SCAD to get posterior
      
      sep_p = split(seq(p),rep(seq(ceiling(p/1000)),each = 1000)[seq(p)])
      llh0 = 0
      
      for(j in 1:length(sep_p)){
        
        p0 = length(sep_p[[j]])
      
      sigma = sigma.iter[[i]][sep_p[[j]],sep_p[[j]]]
      dSigma = diag(sigma)+1e-2
      delta = crossprod(mu.iter[,sep_p[[j]],i], rbind(-1,diag(1,K-1)))
      #scale_delta = sqrt(apply(delta^2,2,mean)) + 1e-4
      #delta = t(t(delta)/scale_delta)
      lambda = rho.iter[iter] * Crho
      lambda.factor = ifelse((n.smp - K) <= p0, 0.2, 0.001)
      ulam = as.double(rev(sort(lambda)))
      pfmat = matrix(1, ncol=L, nrow=p0)
      
      ### (for test) print the smallest lambda such that all gamma is zero
      if (lam_max_print) {
        foo = .Fortran("msda_ncx", obj=double(3L), nk=nk, p=p0, 
                       sigma=as.double(sigma), delta=as.double(delta), 
                       pfmat=as.double(matrix(1, ncol=3, nrow=p0)), 
                       dfmax=dfmax, pmax=as.double(pmax), nlam=3L, flmin=lambda.factor, 
                       ulam=as.double(rep(1, 3)), eps=eps_msda, 
                       maxit=maxit_msda, sml=as.double(1e-06), verbose=as.integer(FALSE), 
                       nalam=integer(1), theta=double(pmax * nk * 3), itheta=integer(pmax), 
                       ntheta=integer(3), alam=double(3), npass=integer(1), 
                       jerr=integer(1))
        print("print lambda max.")
        print(foo$alam[1])
      }
      ############################################################################
      ### local linear algorithm
      ## lasso-initial
      tt1 = Sys.time()
      
      fit = .Fortran("msda_ncx", obj=double(L), nk=nk, p=p0, 
                     sigma=as.double(sigma), delta=as.double(delta), pfmat=pfmat, 
                     dfmax=dfmax, pmax=pmax, nlam=L, flmin=1, ulam=ulam, eps=eps_msda, 
                     maxit=maxit_msda, sml=as.double(1e-06), verbose=as.integer(FALSE), 
                     nalam=integer(1), theta=double(pmax * nk * L), itheta=integer(pmax), 
                     ntheta=integer(L), alam=double(L), npass=integer(1), 
                     jerr=integer(1))
      outlist0 = formatoutput(fit, maxit_msda, pmax, p0, vnames[sep_p[[j]]], nk)
      
      tt2 = Sys.time()
      if (sanity.check) print("lasso_initial")
      if (sanity.check) print(ttt <- tt2 - tt1)
      if (sanity.check){ if( ttt > 1) { 
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      }}
      
      
      if (ncvx) {
        ## compute the weights according to the local linear algorithm
        L_ncx = fit$nalam
        for (l in seq(L_ncx)) {
          btnm = apply(outlist0$theta[[l]], 1, function(x) sqrt(sum(x * x)) )
          pfmat[ ,l] = as.vector(derivative.pen(abs(btnm), ulam[l]))
        }
        ulam_ncvx = rep(1, L_ncx)
        
        tt1 = Sys.time()
        ## ncvx penalization
        fit = .Fortran("msda_ncx", obj=double(L_ncx), nk=nk, p=p0, 
                       sigma=as.double(sigma), delta=as.double(delta), pfmat=pfmat[, seq(L_ncx)], 
                       dfmax=dfmax, pmax=pmax, nlam=L_ncx, flmin=1, ulam=ulam_ncvx, eps=eps_msda, 
                       maxit=maxit_msda, sml=as.double(1e-06), verbose=as.integer(FALSE), 
                       nalam=integer(1), theta=double(pmax * nk * L_ncx), itheta=integer(pmax), 
                       ntheta=integer(L_ncx), alam=double(L_ncx), npass=integer(1), 
                       jerr=integer(1))
        tt2 = Sys.time()
        if (sanity.check) print("SCAD")
        if (sanity.check) print(ttt <- tt2 - tt1)
        if (sanity.check){ if( ttt > 1) { 
          print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        }}
        
        outlist = formatoutput(fit, maxit_msda, pmax, p0, vnames[sep_p[[j]]], nk)
      }
      
      l0 = 1
      l1 = min(fit$nalam,L)
      
      for(l in seq(fit$nalam)){
        if(outlist$df[fit$nalam+1-l]<2){
          l0 = max(l0,l)
        }
        gamma2.iter[sep_p[[j]],,(i-1)*L+l] = as.matrix(outlist$theta[[fit$nalam+1-l]])
      }
      
      llh0 = llh0 - sum((g.iter[,sep_p[[j]],1]-mu.iter[1,sep_p[[j]],i])^2/dSigma)/2
      
      }
      
      
      #llh0 = 0
      
      ############################################################################
      for(l in seq(l1)){
        
        tt1 = Sys.time()
        logw.tmp = alpha.fn.Rversion(g.tmp=g.iter[,,1], mu.tmp=mu.iter[,,i],
                            pro.tmp=pro.iter[i,], gamma.tmp=t(gamma2.iter[,,(i-1)*L+l]),
                            n=n.smp, p=p, K=K)
        w.tmp <- exp(logw.tmp)
        w.tmp[is.na(w.tmp)] <- 0
        w.tmp[is.infinite(w.tmp)] <- 10^6
        alpha.temp[l,,] = w.tmp/rowSums(w.tmp)
        pro.tmp = colMeans(alpha.temp[l,,])
        
        if((l<=l0)&&(g.method==1)){
          z2.iter[i,,l] = z.iter[iter,,1]
          alpha.temp[l,,] =  rep(1,n.smp)%*%t(pro.iter[i,])
        }else{
          z2.iter[i,,l] = apply(alpha.temp[l,,],1,which.min)
        }
        
        
        if((length(unique(z2.iter[i,,l])) < K) || (min(table(z2.iter[i,,l])) < 2L)){
          #miss_z = setdiff(seq(K), unique(z2.iter[i,,l]))
          #large_z = which(z2.iter[i,,l]==which.max(table(z2.iter[i,,l])))
          #for (k in miss_z) {
          #  z2.iter[i,large_z[k],l] = k
          #  z2.iter[i,large_z[k+K],l] = k
          #}
          #z.tmp = z2.iter[i,,l]
          #z.tmp[which(z.tmp==which.max(table(z.tmp)))] = 
          #  rep(1:K,each = ceiling(max(table(z.tmp))/K))[seq(max(table(z.tmp)))]
          #z2.iter[i,,l] = z.tmp
          z2.iter[i,,l] = z.iter[iter,,1]
        }
        
        tt2 = Sys.time()
        if (sanity.check) print("new_alpha.fn")
        if (sanity.check) print(ttt <- tt2 - tt1)
        if (sanity.check){ if( ttt > 1) { 
          print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        }}
        
        #BIC
        
        
        llh = llh0
        #llh = sum(table(z2.iter[i,,l])*log(pro.tmp))
        if (sanity.check) print("new_dmvnorm_log")
        if (sanity.check) print(ttt <- tt2 - tt1)
        if (sanity.check){ if( ttt > 1) { 
          print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        }}
        
        for(k in 1:K){
          yid = z2.iter[i,,l]==k
          ## ytmp = g.iter[yid,,iter]
          llh = llh + sum(logw.tmp[yid,k])
        }  
        
        bic3.iter[i,l] = -2 * llh + log(n.smp) * outlist$df[fit$nalam+1-l] + (log(n.smp)+log(p)) * (K-1)
      }
      
      time2 = Sys.time()
      
      #best l
      l = max(which.min(bic3.iter[i,]), 1)
      l.iter[i] = l
      
      ## M-step  part 1: update cluster proportions p_k
      pro.iter[i+1,] = colMeans(alpha.temp[l,,])
      
      ## M-step  part 2: update cluster means mu_j
      ## mu.iter[,,(iter-1)*S+i+1] = t(sapply(1:K, mu.est.fn, y = g.iter[,,iter], 
      ##   alpha = alpha.temp[(iter-1)*S+i+1,,]))
      mu.iter[,,i+1] = crossprod(alpha.temp[l,,], 
                                 g.iter[,,1]) / colSums(alpha.temp[l,,])
      
      ## M-step  part 3: update cluster covariance sigma
      # if (sanity.check) {
      #   ## print("check")        
      #   tmp = sigma.iter[,,i+1] 
      # } 
      tt1 = Sys.time()
      sigma.iter[[i+1]] = 
        sigma.est.fn.large(mu=mu.iter[seq(K),,i+1], 
                     g=g.iter[,,1], alpha=alpha.temp[l,,seq(K)], 
                     n=n.smp, p=p, K=K, n.inv=n.inv)
      tt2 = Sys.time()
      if (sanity.check) print("new_sigma.est.fn")
      if (sanity.check) print(ttt <- tt2 - tt1)
      if (sanity.check){ if( ttt > 1) { 
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      }}
      ############################################################################
      # if (sanity.check) {
      #   ## print("check")  
      #   tt1 = Sys.time()
      #   for(k in 1:K){
      #     tmp = tmp + 
      #       apply(simplify2array(sapply(1:n.smp, sigma.est.fn.old, 
      #         mu=mu.iter[k,,i+1], 
      #         g=g.iter[,,1], alpha = alpha.temp[(iter-1)*S+i+1,,k], 
      #         n.smp=n.smp, simplify = F)), c(1,2), sum)
      #   }
      #   tt2 = Sys.time()
      #   print("old_sigma.est.fn")
      #   print(tt2 - tt1)
      #   print("------------------------------------")
      #   if(max(abs(tmp - sigma.iter[,,i+1])) > 1.0e-5) stop("sigma.iter")
      # }
      ############################################################################
      ## BIC
      
      bic1.iter[iter,i] = bic3.iter[i,l]
      
      time3 = Sys.time()
      
      ### stop inner ECM iteration if BIC increases
      if (i > 1 && bic1.iter[iter, i] >= bic1.iter[iter, i-1]){
        if(i<S){
          bic1.iter[iter, (i+1):S] = bic1.iter[iter,i] + 1
        }
        break
      }
      
      
    }
    ############################################################################
    #best i
    i = max(which.min(bic1.iter[iter, ]), 1)
    
    #let gamma and z be best among inner iteration
    z.iter[iter,,] = z2.iter[i,,]
    z.iter[iter+1,,1] = z2.iter[i,,l.iter[i]]
    gamma.iter[,,((iter-1)*L+1):(iter*L)] = gamma2.iter[,,((i-1)*L+1):(i*L)]
    
    #store bic for best i
    bic2.iter[iter,] = bic3.iter[i,]
    
    
    
    ## rho
    rho.iter[iter+1] = kap * rho.iter[iter] + C2 * sqrt(log(p) * n.inv)
    
    if(g.method==1){
      ## estimate monotone function empirically
      tt1 = Sys.time()
      g.iter[,,2] = g.est.Rversion(y=y, z=z.iter[iter+1,,1], n=n.smp, p=p, K=K, 
                          pro0 = pro.iter[i+1,], delta0=delta0)
      tt2 = Sys.time()
      if (sanity.check) print("new_g.est_cdf")
      if (sanity.check) print(ttt <- tt2 - tt1)
      if (sanity.check){ if( ttt > 1) { 
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      }}
      ############################################################################  
      # if (sanity.check) {
      #   ## print("check")
      #   tt1 = Sys.time()
      #   tmp = sapply(c(1:p), g.est.old, y=y, z=z.iter[iter+1,,1], K=K, 
      #                pro0 = pro.iter[i+1,], delta0=delta0)
      #   tt2 = Sys.time()
      #   print("old_g.est")
      #   print(tt2 - tt1)
      #   print("------------------------------------")
      #   if(max(abs(tmp - g.iter[,,2])) > 1.0e-5) stop("g.est")
      # }
      ############################################################################  
    }else if(g.method==2){
      ## estimate monotone function thourgh shape restricted regression
      tt1 = Sys.time()
      
      g.iter[,,2] = g.est.alt2(y=y, z=z.iter[iter+1,,1], mu=mu.iter[,,i+1], n=n.smp, p=p,
                               sanity.check=sanity.check)
      tt2 = Sys.time()
      if (sanity.check) print("new_g.iter_srr")
      if (sanity.check) print(ttt <- tt2 - tt1)
      if (sanity.check){ if( ttt > 1) { 
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      }}
    }else if(g.method==3){
      ## estimate monotone function empirically for unselected entries
      ## estimate monotone function thourgh srr for selected entries
      tt1 = Sys.time()
      if(K==2){
        selected = which(gamma.iter[,,(iter-1)*L+l.iter[i]]!=0)
      }else{
        selected = which(apply((gamma.iter[,,(iter-1)*L+l.iter[i]]!=0),1,sum)!=0)
      }
      
      
      
      unselected = setdiff(1:p,selected)
      
      
      if (length(selected) > 0) {
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(c("selected ",length(selected)))
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        
      }  
      
      
      
      if(length(unselected)!=0){
        g.iter[,unselected,2] = g.est.Rversion(y=as.matrix(y[,unselected]), z=z.iter[iter+1,,1], n=n.smp, p=length(unselected), K=K, 
                                      pro0 = pro.iter[i+1,], delta0=delta0)
        if (sanity.check) print(c("unselected ",length(unselected)))
      }
      
      if(length(selected)!=0){
        g.iter[,selected,2] = g.est.alt2(y=as.matrix(y[,selected]), z=z.iter[iter+1,,1], mu=mu.iter[,selected,i+1], 
                                         n=n.smp, p=length(selected), sanity.check = sanity.check)
        if (sanity.check) print(c("selected ",length(selected)))
      }
      
      
      
      tt2 = Sys.time()
      if (sanity.check) print("new_g.iter_srr")
      if (sanity.check) print(ttt <- tt2 - tt1)
      if (sanity.check){ if( ttt > 1) { 
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      }}
    }else if(g.method ==4){
      g.iter[,,2] = y
    }
    
    ### scale it
    g.iter[,,2] = t((t(g.iter[,,2])-apply(t(g.iter[,,2]),1,mean))/(apply(t(g.iter[,,2]),1,sd)+1e-4))
    
    
    time.iter[1,iter] = time2 - time1
    time.iter[2,iter] = time3 - time2
    ### use relative difference as stopping criterion
    ### add code to check stopping criterion
    if(iter>2){
      #if(min(bic2.iter[iter,])>min(bic2.iter[iter-1,])){
      #  iter=iter-1
      #  break
      #}
      if(sum(abs(g.iter[,,2] - g.iter[,,1]))/sum(abs(g.iter[,,1]))+
         sum(abs(mu.iter[,,i+1] - mu.iter[,,1]))/
         (sum(abs(mu.iter[,,1])) + 1)+ 
         sum(abs(pro.iter[i+1,] - pro.iter[1,]))/
         sum(pro.iter[1,]) < 1e-2) {
        break
      }
    }
    
    
    ### let the update with smallest BIC be the starting point of next iteration
    mu.iter[,,1] = mu.iter[,,i+1]
    pro.iter[1, ] = pro.iter[i+1, ]
    sigma.iter[[1]] = sigma.iter[[i+1]]
    g.iter[,,1] = g.iter[,,2]
    
    
  }
  return(list(z=z.iter, mu=mu.iter, sigma=sigma.iter, g=g.iter, gamma=gamma.iter,
              bic1=bic1.iter, bic2= bic2.iter, time.iter=time.iter, iter=iter))
}





