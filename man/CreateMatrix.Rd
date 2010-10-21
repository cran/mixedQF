\name{CreateMatrix}
\alias{CreateMatrix}
\title{Create Matrices used by VarEstim}
\description{Computes MAtrices with the experience plane and the model}
\usage{CreateMatrix(XA)}
\arguments{
 \item{XA}{List which contains X, A, namesX, namesA, as the result of \code{ParseFmla} function}
}
\value{
 \item{X}{fixed effects application matrix}
 \item{namesX}{Names of the fixeds effect which are columns of the fixed effects application matrix}
 \item{A}{forms of random effects covariance matrices}
 \item{namesA}{Names of random effects}
 \item{iG}{invert of gram matrix}
 \item{piX}{projector on the subspaces of columns of fixxed effects matrix}
 \item{piXo}{orthogonal projector of the subspaces of columns of fixxed effects matrix}
}

\author{Jean-Benoist Leger <jleger@agroparistech.fr>}

\examples{
data(dental)

model <- distance ~ 1 + Sex + (1|age) + (1|Subject)
XA<-ParseFmla(model,dental)
G<-CreateMatrix(XA)

}
