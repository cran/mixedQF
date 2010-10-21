ParseFmla <- function(fmla,data)
{
  fixe <- NULL;
  alea <- NULL;
  variables <- attr(terms.formula(fmla),"term.labels");
  
  for(i in 1:length(variables))
  {
    if(length(grep("\\|",variables[i]))>0)
    {
      ## C'est un effet aleatoire
      alea <- c(alea, variables[i]);
    }
    else
    {
      fixe <- c(fixe, variables[i]);
    }
  }

  ## On enleve les I()
  newvar <- vector();
  for(i in 1:length(variables))
  {
    if(length(grep('^I\\((.*)\\)$',variables[i]))>0)
    {
      strSansI <- sub('^I\\((.*)\\)$','\\1',variables[i]);
      varsInc <- attr(terms.formula(formula(paste('i~',strSansI))),"term.labels");
      newvar <- c(newvar, varsInc);
    }
    else
    {
      newvar <- c(newvar, variables[i]);
    }
  }
  
  varlist=vector();
  for(i in 1:length(newvar))
  {
    varlist <- c(varlist, sub('^.*\\| ','',newvar[i]));
  }
  
  ## Maitenant on crée les données croisés

  for(i in 1:length(varlist))
  {
    if(length(grep(':',varlist[i]))>0)
    {
      var <- varlist[i];
      v <- data[[sub(':.*$','',var)]];
      var <- sub('^[^:]*:','',var);
      while(nchar(var)>0)
      {
        v1 <- data[[sub(':.*$','',var)]];
        var <- sub('^[^:]*(:|$)','',var);
        v <- factor(v:v1);
      }
      data[[varlist[i]]] <- v;
    }
  }

  ## On crée la matrice d'influence des effets fixes
  
  X <- matrix(0,length(data[[1]]),0);
  nomColX = vector();

  if(attr(terms.formula(fmla),"intercept"))
  {
    X <- cbind(X, matrix(1,length(data[[1]]),1));
    nomColX <- c(nomColX,'intercept')
  }
  
  for(i in 1:length(fixe))
  {
    if(is.factor(data[[fixe[i]]]))
    {
      X <- cbind(X, nnet::class.ind(data[[fixe[i]]]));
      for(niveau in levels(data[[fixe[i]]]))
      {
        nomColX <- c(nomColX, paste(fixe[i],'=',niveau,sep=''));
      }
    }
    else
    {
      nomColX <- c(nomColX,fixe[i]);
      if(length(grep('^I\\((.*)\\)$',fixe[i]))>0)
      {
        strSansI <- sub('^I\\((.*)\\)$','\\1',fixe[i]);
        varsInc <- attr(terms.formula(formula(paste('i~',strSansI))),"term.labels");
        Colonne <- matrix(0,length(data[[1]]),1);
        for(k in 1:length(varsInc))
        {
          Colonne <- Colonne + data[[varsInc[k]]];
        }
        X <- cbind(X, Colonne);
      }
      else
      {
        X <- cbind(X, data[[fixe[i]]]);
      }
    }
  }


  ## On crée les matrices des effets aleatoires
  
  aleaMats <- list();
  nomAlea <- vector();
  
  for(i in 1:length(alea))
  {
    if(length(grep('^I\\((.*)\\)$',alea[i]))>0)
    {
      strSansI <- sub('^I\\((.*)\\)$','\\1',alea[i]);
      varsInc <- attr(terms.formula(formula(paste('i~',strSansI))),"term.labels");
      MatriceDeCovar <- matrix(0,length(data[[1]]),length(data[[1]]));
      for(k in 1:length(varsInc))
      {
        aleaname <- sub('^.*\\| ','',varsInc[k]);
        MatriceDApplication <- nnet::class.ind(data[[aleaname]]);
        MatriceDeCovar <- MatriceDeCovar + (MatriceDApplication %*% t(MatriceDApplication));
      }
    }
    else
    {
      aleaname <- sub('^.*\\| ','',alea[i]);
      MatriceDApplication <- nnet::class.ind(data[[aleaname]]);
      MatriceDeCovar <- MatriceDApplication %*% t(MatriceDApplication);
    }

    aleaMats <- c(aleaMats,list(MatriceDeCovar));
    nomAlea <- c(nomAlea,alea[i]);
  }

  ## On rajoute le bruit residuel
  aleaMats <- c(aleaMats,list(diag(1,length(data[[1]]))));
  nomAlea <- c(nomAlea, 'residual');


  list(X = X, namesX = nomColX, A = aleaMats, namesA = nomAlea);
}

