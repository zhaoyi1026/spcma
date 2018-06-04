mcma_PCA <-
function(X,M,Y,adaptive=FALSE,var.per=0.8,n.pc=NULL,boot=TRUE,sims=1000,boot.ci.type=c("bca","perc"),conf.level=0.95,p.adj.method=c("BH","bonferroni","BY"))
{
  n<-nrow(M)
  p<-ncol(M)
  
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  # PCA on M~X residuals
  fit.m<-lm(M~X)
  Sigma.m<-cov(fit.m$residuals)
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
  
  # run marginal mediation on PCs
  re.pc<-mcma_BK(X,M.pc,Y,sims=sims,boot=boot,boot.ci.type=boot.ci.type,conf.level=conf.level,p.adj.method=p.adj.method)
  re.pc$U<-U
  re.pc$var.per<-cumsum(svd.m$d[1:n.pc])/sum(svd.m$d)
  
  return(re.pc)
}
