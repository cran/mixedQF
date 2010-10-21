CreateMatrix <- function(XA)
{
  X <- XA$X;
  nomColX <- XA$nomColX;
  aleaMats <- XA$A;
  nomAlea <- XA$nomAlea;

  tr <- function(M) { sum(diag(M)) }

  ## On crée les projecteurs
  n <- dim(X)[1];
  r <- dim(X)[2];
  piX <- X %*% MASS::ginv( t(X) %*% X) %*% t(X);
  piXo <- diag(1,n) - piX;

  ## On diagonalise et on calcule le PXoR
  dpiXo <- eigen(piXo,symmetric=TRUE);
  PXoR <- dpiXo$vectors[,1:(n-r)];

  ## on calcule AR
  AR=list()
  for(i in 1:length(aleaMats))
  {
    AR <- c(AR, list( t(PXoR) %*% aleaMats[[i]] %*% PXoR));
  }

  G <- matrix(0, length(aleaMats), length(aleaMats));
  ## Que le triangle sup
  for(i in 1:length(aleaMats))
  {
    for(j in i:length(aleaMats))
    {
      G[i,j] <- tr(AR[[i]] %*% AR[[j]])
    }
  }

  ## On complete par symetrie
  G <- G+t(G)-diag(diag(G));

  iG <- solve(G);

  ## On crée la matrice C
  C <- matrix(0,length(AR)^2,length(AR)^2);

  for(p in 1:length(AR))
  {
    for(q in 1:length(AR))
    {
      k<-(q-1)*length(AR)+p;
      # colone par colone
      Prov <- matrix(0,length(AR),length(AR));

      for(i in 1:length(aleaMats))
      {
        for(j in i:length(aleaMats))
        {
          Prov[i,j] <- tr(AR[[i]] %*% AR[[p]] %*% AR[[j]] %*% AR[[q]])
        }
      }

      Prov <- Prov + t(Prov) - diag(diag(Prov));

      C[,k] <- as.vector(2 * iG %*% Prov %*% iG);
    }
  }

  list(X = X, namesX = nomColX, A = aleaMats, AR=AR, namesA = nomAlea, iG = iG, piX = piX, piXo = piXo, PXoR = PXoR, C = C, n = n, r = r);
} 
