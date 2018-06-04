deCor <-
function(X)
{
  n<-nrow(X)
  q<-ncol(X)
  
  if(q>1)
  {
    Xnew<-matrix(NA,n,q)
    Xnew[,1]<-X[,1]
    for(j in 2:q)
    {
      xtmp<-X[,1:(j-1)]
      fit<-lm(X[,j]~xtmp)
      Xnew[,j]<-fit$residuals+fit$coefficients[1]
    }
  }else
  {
    Xnew<-X
  }
  return(Xnew)
}
