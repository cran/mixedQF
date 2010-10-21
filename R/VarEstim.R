VarEstim<-function(Y, Mats, contrast)
{
  # Donnees
  # Mats : matrices précalculées
  # Matrice de projection de theta pour obtenir la contrast

  tr <- function(M) { sum(diag(M)) }

  n <- Mats$n;
  r <- Mats$r;

  YR <- t(Mats$PXoR) %*% Y;

  FQ <- matrix(0,length(Mats$AR),1);
  for(i in 1:length(Mats$AR))
  {
    FQ[i,1] <- crossprod(YR,Mats$AR[[i]]%*%YR);
  }

  ## Puis on estime les parametres de variance

  varestim <- Mats$iG %*% FQ;

  # On seuille
  varestim <- (varestim>0)*varestim;

  # On estime Sigma
  Sigma <- matrix(0,n,n);

  for(i in 1:length(Mats$AR))
  {
    Sigma <- Sigma + varestim[i] * Mats$A[[i]];
  }


  # On calcule la matrice d'estimation
  if(rcond(Sigma)>1e-10)
  {
    iSigma <- solve(Sigma);
    M <- MASS::ginv(t(Mats$X) %*% iSigma %*% Mats$X) %*% t(Mats$X) %*% iSigma;
    mcg <- TRUE;
  }
  else
  {
    M <- MASS::ginv(t(Mats$X) %*% Mats$X) %*% t(Mats$X);
    mcg <- FALSE;
  }

  alpha <- matrix(0,length(Mats$AR),1);

  for(i in 1:length(Mats$AR))
  {
    alpha[[i]] <- contrast %*% M %*% Mats$A[[i]] %*% t(M) %*% t(contrast);
  }

  ksic <- t(alpha) %*% varestim;

  # construisons le vecteur gamma
  gamma <- as.vector(c(varestim) %o% c(varestim));

  # construisons l'estimé de la matrice de covar
  v<-Mats$C%*%solve((diag(1,length(Mats$AR)^2,length(Mats$AR)^2)+Mats$C),gamma);
  Mcovar <- matrix(v,length(Mats$AR),length(Mats$AR));

  # les coefs interveants

  varcksic <- t(alpha) %*% Mcovar %*% alpha;

  etac <- 2*ksic^2/varcksic;

  # theta
  theta <- M %*% Y;
  
  T <- (contrast%*%theta)/sqrt(ksic);

  pval <- 2*pt(abs(T),etac,lower.tail=FALSE);

  list(sigmas = varestim, Sigma = Sigma, eta = etac, theta = theta, T = T, pval = pval, ksi=ksic, cli=contrast%*%theta,mcg=mcg);
}

