spcma <-
function(X,M,Y,adaptive=FALSE,var.per=0.8,n.pc=NULL,D=NULL,gamma=0,eps=1e-4,maxsteps=2000,per.jump=0.7,
                     boot=TRUE,sims=1000,boot.ci.type=c("bca","perc"),conf.level=0.95,
                     p.adj.method=c("BH","bonferroni","BY"),PC.run=TRUE)
{
  n<-nrow(M)
  p<-ncol(M)
  
  re.SPCA<-SPCA(X,M,adaptive=adaptive,var.per=var.per,n.pc=n.pc,D=D,gamma=gamma,eps=eps,trace=FALSE,maxsteps=maxsteps,lambda.tune="R2",per.jump=per.jump)
  n.pc<-ncol(re.SPCA$U)
  
  #==================================================
  # run marginal mediaiton on sparse PCs
  M.spc<-M%*%re.SPCA$W
  M.spc.dCor<-deCorM.X(M.spc,X)
  re.spc<-mcma_BK(X,M.spc.dCor,Y,sims=sims,boot=boot,boot.ci.type=boot.ci.type,conf.level=conf.level,p.adj.method=p.adj.method)
  re.spc$W<-re.SPCA$W
  re.spc$var.per<-re.SPCA$var.spc
  #==================================================
  
  if(PC.run)
  {
    #==================================================
    # run marginal mediation on PCs
    M.pc<-M%*%re.SPCA$U
    re.pc<-mcma_BK(X,M.pc,Y,sims=sims,boot=boot,boot.ci.type=boot.ci.type,conf.level=conf.level,p.adj.method=p.adj.method)
    re.pc$U<-re.SPCA$U
    re.pc$var.per<-re.SPCA$var.pc
    #==================================================
    
    re<-list(PCA=re.pc,SPCA=re.spc)
  }else
  {
    re<-list(SPCA=re.spc)
  }
  
  return(re)
}
