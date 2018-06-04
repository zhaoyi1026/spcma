\name{mcma_BK}
\alias{mcma_BK}
\title{ Multiple/marginal causal mediation analysis
}
\description{ This function performs marginal causal mediation analysis under the linear structural equation modeling framework. It assumes the mediators are causally independent.
}
\usage{
mcma_BK(Z, M, Y, sims = 1000, boot = TRUE, boot.ci.type = c("bca", "perc"), 
    conf.level = 0.95, p.adj.method = c("BH", "bonferroni", "BY"))
}
\arguments{
  \item{Z}{ a vector of length \eqn{n}, the randomized treatment assignment.
}
  \item{M}{ a \eqn{n\times p} data matrix, the mediators.
}
  \item{Y}{ a vector of length \eqn{n}, the outcome of interest.
}
  \item{sims}{ a numeric value, the number of bootstrap replications. Default is 1000.
}
  \item{boot}{ a logic variable, if \code{FALSE} a quasi-Bayesian approximation is used for confidence intervals; if \code{TRUE} nonparametric bootstrap will be used. Default is \code{TRUE}. See \code{\link[mediation]{mediate}}.
}
  \item{boot.ci.type}{ a character of the way of calculating bootstrap confidence interval. If \code{boot.ci.type = "bca"}, the bias corrected confidence interval is returned; if \code{boot.ci.type = "perc"}, the percentile confidence interval is returned.
}
  \item{conf.level}{ a numeric value, the designated significance level. Default is \eqn{0.05}, i.e., return \eqn{95\%} confidence interval.
}
  \item{p.adj.method}{ a character of the method implemented for multiple testing adjustment. See \code{\link[stats]{p.adjust}}.
}
}
\details{ Consider the case of multiple mediators, when the mediators are causally independent, i.e., given the treatment assignment, the potential outcome of different mediators are conditionally independent, it is equivalent to conduct a seres of marginal mediation analysis on each individual mediators.
    \deqn{M_{i}^{(j)}=\alpha_{0j}+\alpha_{j}Z_{i}+\xi_{ij},}
    \deqn{Y_{i}=\beta_{0j}+\gamma_{j} Z_{i}+\beta_{j}M_{i}^{(j)}+\eta_{ij},}
where \eqn{\xi_{ij}} and \eqn{\eta_{ij}} are model errors with mean zero and independent of each other. See \code{\link[mediation]{mediate}} for details about single mediator analysis.
}
\value{
    \item{IE}{an outcome matrix, the estimate, \eqn{p}-value, lower and upper bound of the confidence interval, and the adjusted \eqn{p}-value of the indirect effect for each mediator.}

    \item{DE}{an outcome matrix, for each mediator, this provides the estimate of the effects not explained by the mediator, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

    \item{alpha}{an outcome matrix, the estimate of the treatment effect on each mediator, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

    \item{beta}{an outcome matrix, the estimate of the mediator effect on the outcome with a one-unit increment, as well as the \eqn{p}-value, lower and upper bound of the confidence interval and the adjusted \eqn{p}-value.}

    \item{gamma}{an outcome matrix, th estimate, \eqn{p}-value, lower and upper bound of the confidence interval, and the adjusted \eqn{p}-value of \eqn{\gamma_{j}} in the model.}

    \item{IE.total}{an outcome matrix, the estimate, \eqn{p}-value, lower and upper bound of the total indirect effect.}
}
\author{
Yi Zhao, Johns Hopkins University, \email{zhaoyi1026@gmail.com};

Martin A. Lindquist, Johns Hopkins University, \email{mal2053@gmail.com};

Brian S. Caffo, Johns Hopkins University, \email{bcaffo@gmail.com}.
}
\seealso{
\code{\link[mediation]{mediate}}
}
\examples{
#############################################
data(env.example)
X<-get("X",env.example)
M<-get("M",env.example)
Y<-get("Y",env.example)
Phi<-get("Phi",env.example)

# marginal mediation analysis on causally independent mediators
M.tilde<-M\%*\%Phi
mcma_BK(X,M.tilde,Y,boot=FALSE)
#############################################
}
\keyword{ models }
