SPCA <-
function(X,M,adaptive=FALSE,var.per=0.8,n.pc=NULL,D=NULL,gamma=0,eps=1e-4,trace=TRUE,maxsteps=2000,lambda.tune=c("R2"),per.jump=0.7)
{
  n<-nrow(M)
  p<-ncol(M)
  
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  #==================================================
  # PCA on M~X residuals
  fit.m<-lm(M~X)
  E.m<-fit.m$residuals
  Sigma.m<-cov(E.m)
  svd.m<-svd(Sigma.m)
  if(adaptive)
  {
    n.pc<-sum(cumsum(svd.m$d)/sum(svd.m$d)<var.per)+1
  }else
  {
    if(is.null(n.pc))
    {
      n.pc<-p
    }
  }
  U<-svd.m$u[,1:n.pc]
  colnames(U)<-paste0("PC",1:n.pc)
  rownames(U)<-colnames(M)
  M.pc<-M%*%U
  colnames(M.pc)<-paste0("PC",1:n.pc)
  
  var.pc<-cumsum(svd.m$d[1:n.pc])/sum(svd.m$d)
  #==================================================
  
  
  #==================================================
  # Sparse PCA
  if(lambda.tune[1]=="R2")
  {
    # lambda chosen by adjusted total variance
    re.SPCA<-R2.flasso(E.m,U,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps,per.jump=per.jump)
    V<-re.SPCA$V
    SPCA.var.per.cum<-re.SPCA$var.per
  }
  W<-apply(V,2,function(x){return(x/sqrt(sum(x^2)))})
  colnames(V)=colnames(W)<-paste0("PC",1:n.pc)
  rownames(V)=rownames(W)<-colnames(M)
  #==================================================
  
  # U: original loading matrix
  # V: sparsify loading matrix
  # W: sparsify loading matrix with l2-norm 1
  re<-list(U=U,V=V,W=W,var.pc=var.pc,var.spc=SPCA.var.per.cum)
  return(re)
}
