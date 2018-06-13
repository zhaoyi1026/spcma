\name{spcma}
\alias{spcma}
\title{ Sparse principal component based mediation analysis
}
\description{ This function performs sparse principal component based mediation analysis with multiple mediators under the linear structural equation modeling framework. The new mediators, which are causally independent, are linear projections of the original ones with sparse loadings.
}
\usage{
spcma(X, M, Y, adaptive = FALSE, var.per = 0.8, n.pc = NULL, D = NULL, 
  gamma = 0, eps = 1e-04, maxsteps = 2000, per.jump = 0.7, boot = TRUE, 
  sims = 1000, boot.ci.type = c("bca", "perc"),   conf.level = 0.95, 
  p.adj.method=c("BH","bonferroni","BY"), PC.run = TRUE)
}
\arguments{
  \item{X}{ a vector of length \eqn{n}, the randomized treatment assignment.
}
  \item{M}{ a \eqn{n\times p} data matrix, the mediators.
}
  \item{Y}{ a vector of length \eqn{n}, the outcome of interest.
}
  \item{adaptive}{ a logic variable, if \code{adaptive = TRUE}, the number of principal components (PCs) is determined based on the percentage of variance explained. Default is \code{FALSE}.
}
  \item{var.per}{ a numeric variable, the designated percentage of variance explained, which determines the number of PCs. Will be ignored if \code{adaptive = FALSE}.
}
  \item{n.pc}{ an integer variable, the number of PCs considered. Will be ignored if \code{adpative = TRUE}. Default is \code{NULL}. If \code{n.pc = NULL}, will set \code{adaptive} to \code{TRUE}.
}
  \item{D}{ the penalty matrix prespecified for structured regularization in the generalized lasso. Default is \code{NULL}. When \code{D = NULL}, it is set to be a diagonal matrix, that is the lasso penalty is imposed.
}
  \item{gamma}{ a numeric variable greater than or equal to 0, indicating the ratio of the two tuning parameters, one for the fusion penalty, and the other for the pure \eqn{\ell_{1}} penalty. Default is 0. See \code{\link[genlasso]{fusedlasso}}.
}
  \item{eps}{ a numeric variable indicating the multiplier for the ridge penalty, in the case that \code{M} is wide (more columns than rows). If numeric problems occur, make \code{eps} larger. Default is \eqn{1e-4}. See \code{\link[genlasso]{fusedlasso}}.
}
  \item{maxsteps}{ an integer specifying the maximum number of steps for the algorithm to take before termination. Default is 2000. See \code{\link[genlasso]{fusedlasso}}.
}
  \item{per.jump}{ a numeric value, used for tuning parameter selection, the quantile of cut-off for the total variance change under different tuning parameter values. Default is \eqn{0.7}. The higher the value, the sparser the loading matrix.
}
  \item{boot}{ a logic variable, if \code{FALSE} a quasi-Bayesian approximation is used for confidence intervals; if \code{TRUE} nonparametric bootstrap will be used. Default is \code{TRUE}.
}
  \item{sims}{ a numeric value, the number of bootstrap replications. Default is 1000.
}
  \item{boot.ci.type}{ a character of the way of calculating bootstrap confidence interval. If \code{boot.ci.type = "bca"}, the bias corrected confidence interval is returned; if \code{boot.ci.type = "perc"}, the percentile confidence interval is returned.
}
  \item{conf.level}{ a numeric value, the designated significance level. Default is \eqn{0.05}, i.e., return \eqn{95\%} confidence interval.
}
  \item{p.adj.method}{ a character of the method implemented for multiple testing adjustment. See \code{\link[stats]{p.adjust}}.
}
  \item{PC.run}{ a logic value, whether perform principal component analysis based mediation analysis. Default is \code{TRUE}.
}
}
\details{ Consider the case of multiple mediators, and there is a linear projection of the original mediators \eqn{\tilde{M}^{(j)}=M\phi_{j}}, where \eqn{\Phi=(\phi_{1},\dots,\phi_{p})} is an orthonormal matrix, such that for \eqn{j\neq k},
  \deqn{\tilde{M}^{(j)}~\perp~\tilde{M}^{(k)}~|~X.}
We can perform a series of marginal mediation analysis on the new projected mediators, i.e.,
  \deqn{\tilde{M}_{i}^{(j)}=\alpha_{0j}+\alpha_{j}X_{i}+\xi_{ij},}
  \deqn{Y_{i}=\beta_{0j}+\gamma X_{i}+\beta_{j}\tilde{M}_{i}^{(j)}+\eta_{ij},}
where \eqn{\xi_{ij}} and \eqn{\eta_{ij}} are model errors with mean zero and independent of each other.

In order to yield intuitive interpretations, the loading matrix is sparsified based on a structured lasso penalty. This package considers a special case as the fused lasso. The structure is given in the matrix \eqn{D} (also see \code{\link[genlasso]{fusedlasso}}).
}
\value{ When \code{PC.run = TRUE}, two output lists, one for PCA-based analysis (\code{PCA}) and one for sparse PCA-based analysis (\code{SPCA}). The value of \code{PCA} is the same as the outcome in \code{\link[spcma]{mcma_PCA}}.

  For \code{SPCA}
  \item{IE}{an outcome matrix, the estimate, \eqn{p}-value, lower and upper bound of the confidence interval, and the adjusted \eqn{p}-value of the indirect effect for each mediator.}

  \item{DE}{an outcome matrix, for each mediator, this provides the estimate of the direct effect, as well as the \eqn{p}-value, lower and upper bound of the confidence interval.}

  \item{alpha}{an outcome matrix, the estimate of the treatment effect on each mediator, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

  \item{beta}{an outcome matrix, the estimate of the mediator effect on the outcome with a one-unit increment, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

  \item{gamma}{an outcome matrix, th estimate, \eqn{p}-value, lower and upper bound of the confidence interval, and the adjusted \eqn{p}-value of \eqn{\gamma_{j}} in the model. For each mediator, this is the effect not explained by the mediator.}

  \item{IE.total}{an outcome matrix, the estimate, \eqn{p}-value, lower and upper bound of the total indirect effect.}

  \item{W}{the sparsified loading matrix.}

  \item{var.per}{the cumulative proportion of variance explained by the PCs.}
}
\references{
Huang and Pan (2016). ``Hypothesis test of mediation effect in causal mediation model with high-dimensional continuous mediators''. \emph{Biometrics}, 72(2):402-413.

Zou, Hastie and Tibshirani (2006). ``Sparse principal component analysis''. \emph{Journal of computational and graphical statistics}, 15(2):265-286.
}
\author{
Yi Zhao, Johns Hopkins University, \email{zhaoyi1026@gmail.com};

Martin A. Lindquist, Johns Hopkins University, \email{mal2053@gmail.com};

Brian S. Caffo, Johns Hopkins University, \email{bcaffo@gmail.com}.
}
\seealso{
\code{\link[genlasso]{fusedlasso}}.
}
\examples{
#############################################
data(env.example)
X<-get("X",env.example)
M<-get("M",env.example)
Y<-get("Y",env.example)

# sparse principal component based mediation analysis
spcma(X,M,Y,adaptive=TRUE,var.per=0.75,boot=FALSE,PC.run=FALSE)
#############################################
}
\keyword{ models }
