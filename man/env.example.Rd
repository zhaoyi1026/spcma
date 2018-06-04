\name{env.example}
\alias{env.example}
\docType{data}
\title{ Simulated data
}
\description{ "env.example" is an R environment containing the data generated from the proposed high-dimensional mediation model with \eqn{p=50}.
}
\usage{data("env.example")}
\format{ An R environment
	\describe{
		\item{\code{X}}{a \eqn{n\times 1} vector, the binary treatment assignment of \eqn{n} subjects.}
		\item{\code{M}}{a \eqn{n\times p} matrix, the \eqn{p} mediators.}
		\item{\code{Y}}{a \eqn{n\times 1} vector, the outcome of interest.}
		\item{\code{Phi}}{a \eqn{p\times p} matrix, the true loading/projection matrix.}
		\item{\code{alpha}}{a \eqn{p\times 1} vector, the true \eqn{\alpha} value in the model.}
		\item{\code{beta}}{a \eqn{p\times 1} vector, the true \eqn{\beta} value in the model.}
		\item{\code{gamma}}{a scalar, the true \eqn{\gamma} value in the model.}
	}
}
\details{ The data was generated from the following mediation models:
	\deqn{\tilde{M}_{i}^{(j)}=\alpha_{j}X_{i}+\xi_{ij},}
	\deqn{Y_{i}=\gamma X_{i}+\beta_{1}\tilde{M}_{i}^{(1)}+\cdots+\beta_{p}\tilde{M}_{i}^{(p)}+\eta_{i},}
	where \eqn{\tilde{M}_{i}^{(j)}=M_{i}'\phi_{j}}, \eqn{\Phi=(\phi_{1},\dots,\phi_{p})}; \eqn{\xi_{ij}} and \eqn{\eta_{i}} are model errors normally distributed with mean zero.
}
\examples{
data(env.example)
X<-get("X",env.example)
M<-get("M",env.example)
Y<-get("Y",env.example)
}
\keyword{datasets}
