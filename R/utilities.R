### gtilde
g.est.old = function(j, y, z, K, pro0, delta0){
  n.smp = dim(y)[1]
  zmat = matrix(0, ncol=K, nrow = n.smp)
  zmat[,c(1:K) %in% names(table(z))] = model.matrix(~as.factor(z)-1)
  prob1 = t(diag(1/colSums(zmat)))%*% t(zmat) %*% outer(y[,j], y[,j], FUN = '<=')
  prob1[is.na(prob1)]=0
  winsor1 = delta0*(prob1 <= delta0) + prob1*(prob1 > delta0 & prob1 <= 1-delta0) + (1-delta0)*(prob1 > 1-delta0)
  return(pro0%*%qnorm(winsor1))
}

# Boxiang: g.est.Rversion is an R function of g.est
g.est.Rversion = function(y, z, n, p, K, pro0, delta0){
  one_minus_delta0 = 1 - delta0
  tz = table(z)
  prob1 = array(NA, c(K, n, p))
  y_rank = z_sort = matrix(0L, n, p)
  for (j in seq(p)) {
    y_rank[,j] = rank(y[,j])
    z_sort[,j] = factor(z[order(y_rank[,j])], levels=seq.int(K))
  }
  tz_inv = 1.0 / tz
  ct = matrix(0L, n, K)
  for (j in seq(p)) {
    ct[1, ] = 0L
    ct[1, z_sort[1, j]] = 1L
    for (i in seq(2, n)) {
      ct[i, ] = ct[i-1, ]
      ct[i, z_sort[i, j]] = ct[i, z_sort[i, j]] + 1
    }
    for (i in seq(n)) {
      prob1[, i, j] = ct[y_rank[i, j], ] * tz_inv
    }
  }
  prob1[is.na(prob1)]=0
  winsor1_qnorm = qnorm(ifelse(prob1 <= delta0, delta0,
                               ifelse(prob1 > one_minus_delta0, one_minus_delta0, prob1)))
  return(apply(winsor1_qnorm, 3, function(x) pro0 %*% x))
}

# Boxiang: g.est replaces g.est.old. Now g.est returns prob for all columns.
g.est = function(y, z, n, p, K, pro0, delta0){
  tz = as.integer(table(z))
  y_rank = z_sort = matrix(0L, n, p)
  for (j in seq(p)) {
    y_rank[,j] = rank(y[,j])
    z_sort[,j] = factor(z[order(y_rank[,j])], levels=seq.int(K))
  }
  matrix(.Fortran("g_est", as.integer(y_rank), z_sort, tz, n, p, K, 
    delta0, pro0, ret=double(n * p))$ret, nrow=n, ncol=p)
}

# Lyuou: an alternative estimate of g based on monreg
ksmooth.gcv <- function(x, y){
  nobs <- length(y)
  xs <- sort(x, index.return = TRUE)
  x <- xs$x
  y <- y[xs$ix]
  xdif <- outer(x, x, FUN = "-")
  tune.ksmooth <- function(h){
    xden <- dnorm(xdif / h)
    xden <- xden / rowSums(xden)
    df <- sum(diag(xden))
    fit <- xden %*% y
    mean((fit - y)^2) / (1 - df/nobs)^2
  }
  xrng <- diff(range(x))
  oh <- optimize(tune.ksmooth, interval = c(xrng/nobs, xrng))$minimum
  if(any(oh == c(xrng/nobs, xrng)))
    warning("Minimum found on boundary of search range.\nYou should retune model with expanded range.")
  return(h = min(oh,0.1))
}

g.est.alt.single = function(j, y, z, mu, n, sanity.check){
  sca = sqrt(sum(mu[,j]^2))
  mu0 = mu[z,j]/sca
  if (sanity.check) {
    tt1 = Sys.time()
  }
  hr = ksmooth.gcv(x=mu0, y=y[,j])
  if (sanity.check) {
    tt2 = Sys.time()
  }
  gest = monreg(x=mu0, y=y[,j], t=n, hd=hr^2, hr=hr, inverse=1)
  if (sanity.check) {
    tt3 = Sys.time()
    print("ksmooth.gcv")
    print(tt2 - tt1)
    print("monreg")
    print(tt3 - tt2)
  }
  x = (y[,j]-min(gest$t))/(max(gest$t)-min(gest$t))*n
  x[x<=1] = 1
  x[x>=n-1] = n-1
  xfloor = floor(x)
  gest0 = gest$estimation[xfloor]+(x-xfloor)*gest$estimation[xfloor+1]
  gest0[gest0==NA] = gest$estimation[n/2]
  return(gest0*sca)
}

g.est.alt = function(y, z, mu, n, p, sanity.check){
  gest1 = sapply(1:p, g.est.alt.single, y=y, z=z, mu=mu, n=n, sanity.check=sanity.check)
  return(gest1)
}



#Lyuou: shaped restrcited regression
g.est.alt2.single = function(j, y, z, mu, n, sanity.check){
  #set minimal value for scale
  sca = sqrt(sum(mu[,j]^2)) + 1e-4
  
  #avoid response to be the same
  mu0 = mu[z,j]/sca+rnorm(n,0,1e-4)
  
  if (sanity.check) {
    tt1 = Sys.time()
  }
  
  fit = shapereg(mu0 ~ incr(y[,j]))
  
  if (sanity.check) {
    tt2 = Sys.time()
  }
  
  #avoid giving the same estimation at all points
  gest = (fit$constr.fit+rnorm(n,0,1e-4)) * sca
  
  if (sanity.check) {
    tt3 = Sys.time()
    print("shapereg")
    print(tt2 - tt1)
    print("rescale")
    print(tt3 - tt2)
  }
  return(gest)
}

g.est.alt2 = function(y, z, mu, n, p, sanity.check){
  gest1 = sapply(1:p, g.est.alt2.single, y=y, z=z, mu=mu, n=n, sanity.check=sanity.check)
  return(gest1)
}

### gtilde function
gtilde=function(x,j, y, z, K, pro0, delta0){
  g=0
  for(k in 1:K){
    prob1=sum((y[,j]<=x)*(z==k))/sum(z==k)
    prob1[is.na(prob1)]=0
    winsor1 = delta0*(prob1 <= delta0) + prob1*(prob1 > delta0 & prob1 <= 1-delta0) + (1-delta0)*(prob1 > 1-delta0)
    g=g+pro0[k]*qnorm(winsor1)
  }
  return(g)
}




###ghat function
ghat=function(x,j, y, z, K, pro0, delta0){
  gt=sapply(y[,j],gtilde,j=j, y=y, z=z, K=K, pro0=pro0, delta0=delta0)
  mg=mean(gt)
  sg=sd(gt)
  g=gtilde(x,j, y, z, K, pro0, delta0)
  return((g-mg)/sg)
}


#loglikelihood
dmvnorm_log.old <- function(index, mu, sigma, y) {
  loglik=-diag(t(t(y)-mu[index,])%*%ginv(sigma)%*%t(t(t(y)-mu[index,])))/2-
    log(2*pi)*dim(y)[2]/2-determinant(sigma,logarithm=TRUE)$modulus[1]
  #loglik=mvtnorm::dmvnorm(x=y, mean=mu[index,], sigma=sigma, log=TRUE)
  #loglik[which(loglik==-Inf)]=-10^(3)
  return(loglik)
}

# Boxiang: dmvnorm_log replaces dmvnorm_log.old 
#   Now dmvnorm_log returns loglik for all clusters.
#loglikelihood
dmvnorm_log = function(y, mu, ginv_sigma, det_sigma, n, p, K) {
  matrix(.Fortran("loglik", y, n, p, K, mu, ginv_sigma, det_sigma, 
    ret=double(n*K))$ret, nrow=n, ncol=K)
}

#w_k
alpha.fn.old = function(k, log.dens = temp.normal, pro.temp = pro.iter[t,]){
  if(pro.temp[k] == 0){
    out.alpha = rep(0,dim(log.dens)[1])
  } else {
    log.mix1 = sweep(log.dens,2, log(pro.temp), '+')
    log.ratio1 = sweep(log.mix1, 1, log.mix1[,k], '-')
    out.alpha = 1/rowSums(exp(log.ratio1))
    out.alpha[which(rowSums(exp(log.ratio1)) == 0)] = 0
  }
  return(out.alpha)
}

# Lyuou: new alpha.fn use gamma from SCAD
# return f_k(X_i)
alpha.fn.Rversion <- function(g.tmp, mu.tmp, pro.tmp, gamma.tmp, n, p, K, seq_K = 2:K){
  ret <- matrix(0,n,K)
  g.tmp = t(g.tmp)
  mu.tmp = t(mu.tmp)
  for(k in 2:K){
    ret[,k] <- -gamma.tmp[k-1,]%*%(g.tmp-(mu.tmp[,1]+mu.tmp[,k])/2)
  }
  exp.ret= exp(ret + t(matrix(log(pro.tmp),K,n)))
  #ret <- log(exp.ret/apply(exp.ret,1,sum))
  ret <- log(exp.ret)
  return(ret)}

# Boxiang: alpha.fn replaces alpha.fn.old 
#   Now alpha.fn returns alpha for all clusters.
#alpha.fn.previous = function(log.dens, pro.temp, n, K){
#  .Fortran("alpha_fn", log.dens, pro.temp, n, K, ret=double(n * K))$ret
#}


# Lyuou: new alpha.fn use gamma from SCAD
alpha.fn = function(g.tmp, mu.tmp, pro.tmp, gamma.tmp, n, p, K){
  ret = matrix(.Fortran("z_fn", g.tmp, mu.tmp, pro.tmp, gamma.tmp, 
                        n, p, K, ret=double(n*K))$ret, nrow=n, ncol=K)
  #ret = t(-(t(ret)-log(pro.tmp))+log(pro.tmp))
  return(ret)
}



#objective function
opt.fn=function(gamma,mu.tmp,sigma.tmp,rho.tmp){
  return(t(gamma)%*%sigma.tmp%*%gamma/2+abs(t(gamma)%*%(mu.tmp[2,]-mu.tmp[1,]))+sum(SCAD(gamma,rho.tmp)))
}

#derivative of objective function
drv.fn=function(gamma,mu.tmp,sigma.tmp,rho.tmp){
  return(sigma.tmp%*%gamma+abs(mu.tmp[2,]-mu.tmp[1,])+SCAD.derivative(gamma,rho.tmp))
}

# Boxiang: slightly change mu.est.fn
#muhat
mu.est.fn.old = function(k.index, y, alpha){
  mu.k1 = alpha[,k.index] %*% y/sum(alpha[,k.index])
  return(mu.k1)
}

mu.est.fn = function(k.index, y, alpha){
  alpha[,k.index] %*% y / sum(alpha[,k.index])
}

## Boxiang: sigma.est.fn replaces sigma.est.fn.old 
#   Now sigma.est.fn returns sigma estimates for all observations.
#sigmahat
sigma.est.fn.old = function(index, mu, g, alpha ,n.smp=n.smp){
  return(alpha[index]*(1/n.smp)*((g[index,] - mu)%*%t(g[index,] - mu)))
}

#sigmahat
sigma.est.fn = function(mu, g, alpha ,n, p, K, n.inv){
  matrix(.Fortran("sigma_est", mu, g, alpha, n, p, K, n.inv, 
    ret=double(p * p))$ret, p, p)
  # return(alpha[index]*(1/n.smp)*((g[index,] - mu)%*%t(g[index,] - mu)))
}

#sigmahat 
#for large matrix
sigma.est.fn.large = function(mu, g, alpha ,n, p, K, n.inv){
  sep_p = split(seq(p),rep(seq(ceiling(p/1000)),each = 1000)[seq(p)])
  sigma = big.matrix(p, p, type = "short", init = 0)
  for(i in 1:length(sep_p)){
    sigma[sep_p[[i]],sep_p[[i]]] = matrix(.Fortran("sigma_est", mu[,sep_p[[i]]], g[,sep_p[[i]]], alpha, n, 
      length(sep_p[[i]]), K, n.inv, ret=double(length(sep_p[[i]]) * length(sep_p[[i]])))$ret, 
      length(sep_p[[i]]), length(sep_p[[i]]))
  }
  return(sigma)
  # return(alpha[index]*(1/n.smp)*((g[index,] - mu)%*%t(g[index,] - mu)))
}

## Boxiang: may have a bug in the original function z.fn.bug
####################bug fixed####################
#zhat
z.fn.bug = function(i,g.tmp,mu.tmp,pro.tmp,gamma.tmp){
  s.tmp=c(log(pro.tmp[1]),(diag(t(g.tmp[i,]-t(mu.tmp[-1,]+mu.tmp[1,])/2)%*%t(gamma.tmp)))+log(pro.tmp[-1]))
  z.tmp=which.max(s.tmp)
  return(z.tmp)
}

z.fn.old = function(i,g.tmp,mu.tmp,pro.tmp,gamma.tmp){
  s.tmp=c(log(pro.tmp[1]),(diag(t(g.tmp[i,]-((t(mu.tmp[-1,])+mu.tmp[1,])/2))%*%t(gamma.tmp)))+log(pro.tmp[-1]))
  z.tmp=which.max(s.tmp)
  return(z.tmp)
}

#zhat
z.fn = function(g.tmp, mu.tmp, pro.tmp, gamma.tmp, n, p, K){
  ret = matrix(.Fortran("z_fn", g.tmp, mu.tmp, pro.tmp, gamma.tmp, 
    n, p, K, ret=double(n*K))$ret, nrow=n, ncol=K)
  apply(ret, 1, which.min)
}


### mis-clustering function, kendall's tau
miss.clus = function(a, true.label){
  pair.label.true = outer(true.label, true.label, '==')
  pair.label.est = outer(a, a, '==')
  return(sum(abs(pair.label.est[upper.tri(pair.label.est)] - pair.label.true[upper.tri(pair.label.true)]))/sum(upper.tri(pair.label.true)))
}

### mis-clustering function, hamming distance
error_clus=function(clus,true_clus){
  m=max(c(clus,true_clus))
  Ptmp=gtools:::permutations(m,m)
  etmp=rep(0,factorial(m))
  for(i in 1:factorial(m)){
    etmp[i]=mean((as.numeric(factor(clus,level=Ptmp[i,]))-true_clus!=0))
  }
  return(min(etmp))
}

zeromat=function(nvars,lmu,vnames,stepnames){
  ca=rep(0,lmu)
  ia=seq(lmu+1)
  ja=rep(1,lmu)
  dd=c(nvars,lmu)
  new("dgCMatrix", Dim = dd,
      Dimnames = list(vnames,stepnames),
      x = as.vector(ca),
      p = as.integer(ia - 1), i = as.integer(ja - 1))
}

# extract fortran outputs and format it into sparse matries
formatoutput <- function(fit, maxit, pmax, nvars, vnames, nk) {
  nalam <- fit$nalam
  ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta)
  lam <- fit$alam[seq(nalam)]
  obj <- fit$obj[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  resnames <- paste("delta", seq(nk), sep = "")
  
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
  
  dd <- c(nvars, nk)
  df <- rep(0, nalam)
  if (nthetamax > 0) {
    ja <- fit$itheta[seq(nthetamax)]
    oja <- order(ja)
    ja <- rep(ja[oja], nk)
    itheta <- cumsum(c(1, rep(nthetamax, nk)))
    pos <- rep(1:nalam, each = nk * pmax)
    theta <- split(fit$theta[seq(nk * pmax * nalam)], pos)
    for (l in 1:nalam) {
      theta[[l]] <- matrix(theta[[l]], pmax, nk, byrow = TRUE)[seq(nthetamax), 
                                                               , drop = FALSE]
      df[l] <- sum(rowSums(abs(theta[[l]])) != 0)
      theta[[l]] <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames, 
                                                               resnames), x = as.vector(theta[[l]][oja, ]), p = as.integer(itheta - 
                                                                                                                             1), i = as.integer(ja - 1))
    }
  } else {
    theta <- list()
    for (l in 1:nalam) {
      theta[[l]] <- zeromat(nvars, nk, vnames, resnames)
    }
    df <- rep(0, nalam)
  }
  list(theta = theta, df = df, dim = dd, lambda = lam, obj = obj)
}

# error messages from fortran
err <- function(n, maxit, pmax) {
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    # fatal error
    if (n < 7777) 
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    # non fatal error
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                   maxit, " iterations; solutions for larger lambdas returned.\n", 
                   sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", 
                   pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", 
                   sep = "")
    if (n < -20000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", 
                   pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", 
                   sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}



derivative.pen = function(u, lambda, a = 3.7) {
  u = abs(u) # u must be nonnegative
  lambda * (u <= lambda) + (a * lambda - u) / (a - 1) *
    (u > lambda) * (u <= a * lambda)
}

initial.kmeans = function(y,K,rho,nstart=100){
  p = dim(y)[2]
  n.smp = dim(y)[1]
  Crho = seq(0.1,2,0.1)
  L = length(Crho)
  
  nk = K - 1L
  vnames = paste("V", seq(p), sep = "")
  maxit_msda = as.integer(1e+3)
  eps_msda = as.double(1e-4)
  dfmax = as.integer(n.smp)
  pmax = as.integer(min(dfmax * 2 + 20, p))
  
  z.initial = matrix(0,nstart,n.smp)
  bic1.iter = rep(0,nstart)
  
  for(i in 1:nstart){
    centers <- y[sample.int(n.smp, K), ]
    kms1 = kmeans(y, centers=centers, iter.max=100)
    z.initial[i,] = as.integer(kms1$cluster)
    if ((length(unique(z.initial[i,])) < K) || (min(table(z.initial[i,])) < 2L)) {
      for (k in seq.int(K)) {
        z.initial[i,k] = k
        z.initial[i,k+K] = k
      }
    }
    pro.tmp = as.numeric(table(z.initial[i,]) /n.smp)
    delta0 =  1/(4*(n.smp)^0.25 * sqrt(pi * log(n.smp)))
    
    g.tmp = g.est(y=y, z=z.initial[i,], n=n.smp, p=p, K=K, 
                        pro0=pro.tmp, delta0=delta0)
    g.tmp = apply(g.tmp, 2, scale)
    
    center2.g = g.tmp
    mu.tmp = matrix(0,K,p)
    for (k in seq(K)) {
      mu.tmp[k,] = colMeans(g.tmp[z.initial[i,] == k,])
      center2.g[z.initial[i,] == k, ] = apply(center2.g[z.initial[i,] == k, ], 
                                               2, scale, scale=FALSE)
    }
    
    sigma = (n.smp-1L) * cov(center2.g) /n.smp
    
    delta = crossprod(mu.tmp, rbind(-1,diag(1,K-1)))
    lambda = rho * Crho
    lambda.factor = ifelse((n.smp - K) <= p, 0.2, 0.001)
    ulam = as.double(rev(sort(lambda)))
    pfmat = matrix(1, ncol=L, nrow=p)
    
    fit = .Fortran("msda_ncx", obj=double(L), nk=nk, p=p, 
                   sigma=as.double(sigma), delta=as.double(delta), pfmat=pfmat, 
                   dfmax=dfmax, pmax=pmax, nlam=L, flmin=1, ulam=ulam, eps=eps_msda, 
                   maxit=maxit_msda, sml=as.double(1e-06), verbose=as.integer(FALSE), 
                   nalam=integer(1), theta=double(pmax * nk * L), itheta=integer(pmax), 
                   ntheta=integer(L), alam=double(L), npass=integer(1), 
                   jerr=integer(1))
    outlist = formatoutput(fit, maxit_msda, pmax, p, vnames, nk)
    
    L_ncx = fit$nalam
    for (l in seq(L_ncx)) {
      btnm = apply(outlist$theta[[l]], 1, function(x) sqrt(sum(x * x)) )
      pfmat[ ,l] = as.vector(derivative.pen(abs(btnm), ulam[l]))
    }
    ulam_ncvx = rep(1, L_ncx)
    
    ## ncvx penalization
    fit = .Fortran("msda_ncx", obj=double(L_ncx), nk=nk, p=p, 
                   sigma=as.double(sigma), delta=as.double(delta), pfmat=pfmat[, seq(L_ncx)], 
                   dfmax=dfmax, pmax=pmax, nlam=L_ncx, flmin=1, ulam=ulam_ncvx, eps=eps_msda, 
                   maxit=maxit_msda, sml=as.double(1e-06), verbose=as.integer(FALSE), 
                   nalam=integer(1), theta=double(pmax * nk * L_ncx), itheta=integer(pmax), 
                   ntheta=integer(L_ncx), alam=double(L_ncx), npass=integer(1), 
                   jerr=integer(1))
    outlist = formatoutput(fit, maxit_msda, pmax, p, vnames, nk)
    
    bic2.iter = rep(0,fit$nalam)
    
    for(l in seq(fit$nalam)){
      for(k in seq(K-1)){
        gamma2.tmp = t(outlist$theta[[fit$nalam+1-l]][,k])
      } 
      tt1 = Sys.time()
      logw.tmp = -abs(alpha.fn(g.tmp=g.tmp, mu.tmp=mu.tmp,
                               pro.tmp=pro.tmp, gamma.tmp=t(gamma2.tmp),
                               n=n.smp, p=p, K=K))
      #BIC
      llh = 0
      
      for(k in 1:K){
        yid = z.initial[i,]==k
        llh = llh + sum(logw.tmp[yid,k])
      }  
      
      bic2.iter[l] = -2 * llh + log(n.smp) * outlist$df[fit$nalam+1-l]
    }
    
    bic1.iter[i] = min(bic2.iter)
  }
  
  i = max(which.min(bic1.iter), 1)
  
  return(z.initial[i,])
  
}

initial.sc = function(y,K,rho,nstart=100){
  p = dim(y)[2]
  n.smp = dim(y)[1]
  Crho = seq(0.1,2,0.1)
  L = length(Crho)
  
  nk = K - 1L
  vnames = paste("V", seq(p), sep = "")
  maxit_msda = as.integer(1e+3)
  eps_msda = as.double(1e-4)
  dfmax = as.integer(n.smp)
  pmax = as.integer(min(dfmax * 2 + 20, p))
  
  z.initial = matrix(0,nstart,n.smp)
  bic1.iter = rep(0,nstart)
  
  for(i in 1:nstart){
    z.initial[i,] = as.numeric(specc(y, centers=K))
    if ((length(unique(z.initial[i,])) < K) || (min(table(z.initial[i,])) < 2L)) {
      for (k in seq.int(K)) {
        z.initial[i,k] = k
        z.initial[i,k+K] = k
      }
    }
    pro.tmp = as.numeric(table(z.initial[i,]) /n.smp)
    delta0 =  1/(4*(n.smp)^0.25 * sqrt(pi * log(n.smp)))
    
    g.tmp = g.est(y=y, z=z.initial[i,], n=n.smp, p=p, K=K, 
                  pro0=pro.tmp, delta0=delta0)
    g.tmp = apply(g.tmp, 2, scale)
    
    center2.g = g.tmp
    mu.tmp = matrix(0,K,p)
    for (k in seq(K)) {
      mu.tmp[k,] = colMeans(g.tmp[z.initial[i,] == k,])
      center2.g[z.initial[i,] == k, ] = apply(center2.g[z.initial[i,] == k, ], 
                                              2, scale, scale=FALSE)
    }
    
    sigma = (n.smp-1L) * cov(center2.g) /n.smp
    
    delta = crossprod(mu.tmp, rbind(-1,diag(1,K-1)))
    lambda = rho * Crho
    lambda.factor = ifelse((n.smp - K) <= p, 0.2, 0.001)
    ulam = as.double(rev(sort(lambda)))
    pfmat = matrix(1, ncol=L, nrow=p)
    
    fit = .Fortran("msda_ncx", obj=double(L), nk=nk, p=p, 
                   sigma=as.double(sigma), delta=as.double(delta), pfmat=pfmat, 
                   dfmax=dfmax, pmax=pmax, nlam=L, flmin=1, ulam=ulam, eps=eps_msda, 
                   maxit=maxit_msda, sml=as.double(1e-06), verbose=as.integer(FALSE), 
                   nalam=integer(1), theta=double(pmax * nk * L), itheta=integer(pmax), 
                   ntheta=integer(L), alam=double(L), npass=integer(1), 
                   jerr=integer(1))
    outlist = formatoutput(fit, maxit_msda, pmax, p, vnames, nk)
    
    L_ncx = fit$nalam
    for (l in seq(L_ncx)) {
      btnm = apply(outlist$theta[[l]], 1, function(x) sqrt(sum(x * x)) )
      pfmat[ ,l] = as.vector(derivative.pen(abs(btnm), ulam[l]))
    }
    ulam_ncvx = rep(1, L_ncx)
    
    ## ncvx penalization
    fit = .Fortran("msda_ncx", obj=double(L_ncx), nk=nk, p=p, 
                   sigma=as.double(sigma), delta=as.double(delta), pfmat=pfmat[, seq(L_ncx)], 
                   dfmax=dfmax, pmax=pmax, nlam=L_ncx, flmin=1, ulam=ulam_ncvx, eps=eps_msda, 
                   maxit=maxit_msda, sml=as.double(1e-06), verbose=as.integer(FALSE), 
                   nalam=integer(1), theta=double(pmax * nk * L_ncx), itheta=integer(pmax), 
                   ntheta=integer(L_ncx), alam=double(L_ncx), npass=integer(1), 
                   jerr=integer(1))
    outlist = formatoutput(fit, maxit_msda, pmax, p, vnames, nk)
    
    bic2.iter = rep(0,fit$nalam)
    
    for(l in seq(fit$nalam)){
      for(k in seq(K-1)){
        gamma2.tmp = t(outlist$theta[[fit$nalam+1-l]][,k])
      } 
      tt1 = Sys.time()
      logw.tmp = -abs(alpha.fn(g.tmp=g.tmp, mu.tmp=mu.tmp,
                               pro.tmp=pro.tmp, gamma.tmp=t(gamma2.tmp),
                               n=n.smp, p=p, K=K))
      #BIC
      llh = 0
      
      for(k in 1:K){
        yid = z.initial[i,]==k
        llh = llh + sum(logw.tmp[yid,k])
      }  
      
      bic2.iter[l] = -2 * llh + log(n.smp) * outlist$df[fit$nalam+1-l]
    }
    
    bic1.iter[i] = min(bic2.iter)
  }
  
  i = max(which.min(bic1.iter), 1)
  
  return(z.initial[i,])
  
}


initial.random = function(y,K,C0){
  p = dim(y)[2]
  n.smp = dim(y)[1]
  
  dir = r_unif_sph(K,p)[,,1]
  
  mu0 = C0*(p*log(n.smp)/n.smp)^(1/4)*dir
  
  z.initial = kmeans(y, centers = mu0)$cluster
  
  return(list(z.initial = z.initial, mu.initial = mu0))
  
}
