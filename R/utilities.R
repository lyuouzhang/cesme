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

# Lyuou: an alternative estimate of g based on monreg
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

# Lyuou: an alternative estimate of g based on monreg
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

#Lyuou: shaped restrcited regression
g.est.alt2 = function(y, z, mu, n, p, sanity.check){
  gest1 = sapply(1:p, g.est.alt2.single, y=y, z=z, mu=mu, n=n, sanity.check=sanity.check)
  return(gest1)
}

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


#zhat
z.fn = function(g.tmp, mu.tmp, pro.tmp, gamma.tmp, n, p, K){
  ret = matrix(.Fortran("z_fn", g.tmp, mu.tmp, pro.tmp, gamma.tmp, 
    n, p, K, ret=double(n*K))$ret, nrow=n, ncol=K)
  apply(ret, 1, which.min)
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
    # do not report non fatal errors
    n <- 0
  }
  list(n = n, msg = msg)
}



derivative.pen = function(u, lambda, a = 3.7) {
  u = abs(u) # u must be nonnegative
  lambda * (u <= lambda) + (a * lambda - u) / (a - 1) *
    (u > lambda) * (u <= a * lambda)
}

initial.sc = function(y,K,nstart=10){
  A = t(y)%*%y
  eig.A = eigen(A)
  U0 = eig.A$vectors[,1:K]
  km.result = kmeans(y%*%U0,centers = K,nstart = nstart)
  sc.select = km.result$cluster
  sc.mean = t(U0%*%t(km.result$centers))
  
  return(list(z.initial=sc.select,mu.initial=sc.mean))
  
}


