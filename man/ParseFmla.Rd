\name{ParseFmla}
\alias{ParseFmla}
\title{Create Matrices used by CreateMatrix}
\description{Computes Matrices with the experience plane and the model}
\usage{ParseFmla(fmla,data)}
\arguments{
 \item{fmla}{Formula with fixed and random effects

  \code{y ~ 1 + x + (1|a) + I( (1|b) + (1|c) )}

  where I() used as in formulas.}
  \item{data}{Experience plane, a dataframe, with fixed and random effects}
}
\value{
 \item{X}{fixed effects application matrix}
 \item{namesX}{Names of the fixeds effect which are columns of the fixed effects application matrix}
 \item{A}{forms of random effects covariance matrices}
 \item{namesA}{Names of random effects}
}

\author{Jean-Benoist Leger <jleger@agroparistech.fr>}

\examples{
data(dental)

model <- distance ~ 1 + Sex + (1|age) + (1|Subject)
XA<-ParseFmla(model,dental)

}
