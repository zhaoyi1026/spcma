\name{mcma_PCA}
\alias{mcma_PCA}
\title{ Principal component analysis based multiple mediation analysis
}
\description{ This function performs marginal causal mediation analysis under the linear structural equation modeling framework. The new mediators are linear projections of the original ones, such that they are causally independent.
}
\usage{
mcma_PCA(X, M, Y, adaptive = FALSE, var.per = 0.8, n.pc = NULL, boot = TRUE, 
  sims = 1000, boot.ci.type = c("bca", "perc"), conf.level = 0.95, 
  p.adj.method = c("BH", "bonferroni", "BY"))
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
}
\details{ Consider the case of multiple mediators, and there is a linear projection of the original mediators \eqn{\tilde{M}^{(j)}=M\phi_{j}}, where \eqn{\Phi=(\phi_{1},\dots,\phi_{p})} is an orthonormal matrix, such that for \eqn{j\neq k},
  \deqn{M^{(j)}~\perp~M^{(k)}~|~X}.
We can perform a series of marginal mediation analysis on the new projected mediators, i.e.,
  \deqn{\tilde{M}_{i}^{(j)}=\alpha_{0j}+\alpha_{j}X_{i}+\xi_{ij},}
  \deqn{Y_{i}=\beta_{0j}+\gamma X_{i}+\beta_{j}\tilde{M}_{i}^{(j)}+\eta_{ij},}
where \eqn{\xi_{ij}} and \eqn{\eta_{ij}} are model errors with mean zero and independent of each other.
}
\value{
  \item{IE}{an outcome matrix, the estimate, \eqn{p}-value, lower and upper bound of the confidence interval, and the adjusted \eqn{p}-value of the indirect effect for each mediator.}

  \item{DE}{an outcome matrix, for each mediator, this provides the estimate of the direct effect, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

  \item{alpha}{an outcome matrix, the estimate of the treatment effect on each mediator, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

  \item{beta}{an outcome matrix, the estimate of the mediator effect on the outcome with a one-unit increment, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

  \item{gamma}{an outcome matrix, th estimate, \eqn{p}-value, lower and upper bound of the confidence interval, and the adjusted \eqn{p}-value of \eqn{\gamma_{j}} in the model. For each mediator, this is the effect not explained by the mediator.}

  \item{IE.total}{an outcome matrix, the estimate, \eqn{p}-value, lower and upper bound of the total indirect effect.}

  \item{U}{the estimated loading matrix.}

  \item{var.per}{the cumulative proportion of variance explained by the PCs.}
}
\references{ Huang and Pan (2016). \emph{Hypothesis test of mediation effect in causal mediation model with high-dimensional continuous mediators.} Biometrics, 72(2):402-413.
}
\author{
Yi Zhao, Johns Hopkins University, \email{zhaoyi1026@gmail.com};

Martin A. Lindquist, Johns Hopkins University, \email{mal2053@gmail.com};

Brian S. Caffo, Johns Hopkins University, \email{bcaffo@gmail.com}.
}
\examples{
#############################################
data(env.example)
X<-get("X",env.example)
M<-get("M",env.example)
Y<-get("Y",env.example)

# principal component based mediation analysis
mcma_PCA(X,M,Y,adaptive=TRUE,var.per=0.75,boot=FALSE)
#############################################
}
\keyword{ models }
