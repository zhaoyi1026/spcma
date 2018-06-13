\name{plot_spcma}
\alias{plot_spcma}
\title{ Plot the outcome
}
\description{ This function realizes the visualization of the result of multiple mediation analysis.
}
\usage{
plot_spcma(object, plot.coef = c("alpha", "beta", "IE"), 
  cex.lab = 1, cex.axis = 1, pt.cex = 1, ...)
}
\arguments{
  \item{object}{ \code{mcma_PCA} or \code{mcma_BK} object.
}
  \item{plot.coef}{ a character indicating the parameter to be plotted.
}
  \item{cex.lab}{ the magnification to be used for \code{x} and \code{y} labels relative to the current setting of \code{cex}. See \code{\link[graphics]{par}}.
}
  \item{cex.axis}{ the magnification to be used for axis annotation relative to the current setting of \code{cex}. See \code{\link[graphics]{par}}.
}
  \item{pt.cex}{ a numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. See \code{\link[graphics]{par}}.
}
  \item{\dots}{ additional argument to be passed.
}
}
\details{ Visualization of the parameter estimates in the multiple mediation analysis.
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
Phi<-get("Phi",env.example)

# marginal mediation analysis on causally independent mediators
M.tilde<-M\%*\%Phi
re.BK<-mcma_BK(X,M.tilde,Y,boot=FALSE)
plot_spcma(re.BK,plot.coef="IE")

# principal component based mediation analysis
re.PCA<-mcma_PCA(X,M,Y,adaptive=TRUE,var.per=0.75,boot=FALSE)
plot_spcma(re.PCA,plot.coef="IE")

# sparse principal component based mediation analysis
re.SPCA<-spcma(X,M,Y,adaptive=TRUE,var.per=0.75,boot=FALSE,PC.run=FALSE)
plot_spcma(re.SPCA$SPCA,plot.coef="IE")
#############################################
}
\keyword{ aplot }
