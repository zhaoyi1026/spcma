deCorM.X <-
function(M.tilde,X)
{
  n<-nrow(M.tilde)
  q<-ncol(M.tilde)
  
  if(q>1)
  {
    Mnew<-matrix(NA,n,q)
    Mnew[,1]<-M.tilde[,1]
    for(j in 2:q)
    {
      xtmp<-M.tilde[,1:(j-1)]
      fit<-lm(M.tilde[,j]~X+xtmp)
      Mnew[,j]<-fit$residuals+cbind(rep(1,n),X)%*%fit$coefficients[c(1,2)]
    }
  }else
  {
    Mnew<-M.tilde
  }
  return(Mnew)
}
