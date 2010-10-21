\name{VarEstim}
\alias{VarEstim}
\title{Estimate Variance, fixed effects components, and do test}
\description{Estimate Variance, fixed effects components, and do test}
\usage{VarEstim(Y,Mats,contrast)}
\arguments{
 \item{Y}{Observed data}
  \item{Mats}{Matrices of the model created by CreateMatrix}
 \item{contrast}{Contrast of the fixed effects for doing test}
}
\value{
 \item{sigmas}{estimate of components of the variance}
 \item{Sigma}{Estimate of the variance covariance matrix of Y}
 \item{theta}{fixeds effects}
 \item{T}{Value of Student statistic}
 \item{pval}{p-value of null hypothesis of the contrast}
 \item{cli}{estimate of the contrast}
 \item{ksi}{estimate of variance of the contrast}
}

\author{J.Benoist Leger <jleger@agroparistech.fr>}

\examples{
data(dental)

model <- distance ~ 1 + Sex + (1|age) + (1|Subject)
G<-CreateMatrix(model,dental)

# G$namesX contains the names of columns to build contrast
# contrast is difference between Male and Female
contrast <- matrix(c(0,1,-1),1,3);

estim<-VarEstim(dental$distance,G,contrast);

# estim$pval contains the pvalue to test the null hypothesis on the contrast
# (null hypothesis is "there is no difference between mean of male an mean of
# female)

}
