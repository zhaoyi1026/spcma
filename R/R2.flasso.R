R2.flasso <-
function(E,U,D=NULL,gamma=0,eps=1e-4,maxsteps=2000,per.jump=0.7)
{
  p<-ncol(E)
  n.pc<-ncol(U)
  
  # PC
  E.pc<-E%*%U
  
  # % of variance
  var.total<-sum(diag(cov(E)))
  var.per=var.per.ind<-rep(NA,n.pc)
  for(j in 1:n.pc)
  {
    var.per.ind[j]<-var(E.pc[,j])/var.total
    var.per[j]<-sum(diag(cov(matrix(E.pc[,1:j],ncol=j))))/var.total
  }
  
  lambda.est<-rep(NA,n.pc)
  V<-matrix(NA,p,n.pc)
  var.per.new<-rep(NA,n.pc)
  # first PC lambda choice
  if(is.null(D))
  {
    out.tmp<-fusedlasso1d(y=E.pc[,1],X=E,gamma=gamma,eps=eps,maxsteps=maxsteps)  
  }else
  {
    out.tmp<-fusedlasso(y=E.pc[,1],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
  }
  var.per.tmp<-rep(NA,length(out.tmp$lambda))
  for(k in 1:length(out.tmp$lambda))
  {
    var.per.tmp[k]<-var(out.tmp$fit[,k])/var.total
  }
  var.per.diff.tmp<-var.per.tmp[2:length(var.per.tmp)]-var.per.tmp[1:(length(var.per.tmp)-1)]
  lambda.idx.tmp<-max(which(var.per.diff.tmp>quantile(var.per.diff.tmp,probs=per.jump)))+1
  # lambda.idx.tmp<-min(which(abs(var.per.tmp-var.per.ind[1])<var.diff))
  lambda.est[1]<-out.tmp$lambda[lambda.idx.tmp]
  V[,1]<-out.tmp$beta[,lambda.idx.tmp]
  var.per.new[1]<-var.per.tmp[lambda.idx.tmp]
  
  for(j in 2:n.pc)
  {
    Etmp<-deCor(cbind(E%*%V[,1:(j-1)],E.pc[,j]))
    
    if(is.null(D))
    {
      out.tmp<-fusedlasso1d(y=Etmp[,j],X=E,gamma=gamma,eps=eps,maxsteps=maxsteps)  
    }else
    {
      out.tmp<-fusedlasso(y=Etmp[,j],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
    }
    
    var.per.tmp=var.per.tol.tmp<-rep(NA,length(out.tmp$lambda))
    for(k in 1:length(out.tmp$lambda))
    {
      var.per.tmp[k]<-var(out.tmp$fit[,k])/var.total
      dtmp<-deCor(cbind(E%*%V[,1:(j-1)],out.tmp$fit[,k]))
      var.per.tol.tmp[k]<-sum(diag(cov(dtmp)))/var.total
    }
    var.per.diff.tmp<-var.per.tol.tmp[2:length(var.per.tol.tmp)]-var.per.tol.tmp[1:(length(var.per.tol.tmp)-1)]
    lambda.idx.tmp<-max(which(var.per.diff.tmp>quantile(var.per.diff.tmp,probs=per.jump)))+1
    lambda.est[j]<-out.tmp$lambda[lambda.idx.tmp]
    V[,j]<-out.tmp$beta[,lambda.idx.tmp]
    var.per.new[j]<-var.per.tol.tmp[lambda.idx.tmp]
  }
  
  re<-list(lambda=lambda.est,V=V,var.per=var.per.new)
  return(re)
}
