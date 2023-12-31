% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cesme_large.R
\name{npn.clust.bic.large}
\alias{npn.clust.bic.large}
\title{R package: cesme.}
\usage{
npn.clust.bic.large(
  y,
  K,
  thre = NULL,
  N = 100,
  C0 = 1,
  C2 = 1,
  kap = 0.8,
  z.initial = NULL,
  mu.initial = NULL,
  maxit_msda = as.integer(1000),
  eps_msda = 1e-04,
  ncvx = TRUE,
  sanity.check = FALSE,
  lam_max_print = FALSE,
  g.method = 1,
  kmnstart = 100
)
}
\arguments{
\item{y}{Response.}

\item{K}{Number of clusters.}

\item{thre}{Threshold for estimating transformation. If thre is NULL, use thre=1/(4*n^0.25*sqrt(pi*log(n))).}

\item{N}{Number of iterates.}

\item{C0}{A constant for penalty parameter of SCAD.}

\item{C2}{A constant for penalty parameter of SCAD.}

\item{kap}{A constant for penalty parameter of SCAD.}

\item{z.initial}{Initial clustering assignment. If z.initial is NULL, use spectral clustering to obtain an initialization.}

\item{mu.initial}{Initial clustering assignment. If mu.initial is NULL, use sample cluster means with z.initial to obtain an initialization.}

\item{maxit_msda}{Maximum number of iterates in SCAD solver.}

\item{eps_msda}{Precision of the penalized regression solver.}

\item{ncvx}{Whether to use the nonconvex penalty (SCAD).}

\item{sanity.check}{For debug use.}

\item{lam_max_print}{For testing the code: print the smallest lambda such that all coefficients are zero.}

\item{g.method}{Which method to use to estimate the monotone function.}

\item{kmnstart}{Number of initial starts in spectral clustering.}
}
\value{
A list consisting of
  \item{z.iter}{estimate of z}
  \item{mu.iter}{estimate of mu}
  \item{sigma.iter}{estimate of sigma}
  \item{g.iter}{estimate of g}
  \item{bic1.iter}{bic values}
  \item{bic2.iter}{bic values}
  \item{time.iter}{run time}
}
\description{
High-dimensional Cluster Analysis with Latent Semiparametric Mixture Models.
}
\details{
An algorithm for high-dimensional model-based clustering via alternating maximizations. Ues bigmemory to store sample covariance matrix of data.
}
\keyword{clustering}
