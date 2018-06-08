mcma_BK <-
function(X,M,Y,sims=1000,boot=TRUE,boot.ci.type=c("bca","perc"),conf.level=0.95,p.adj.method=c("BH","bonferroni","BY"))
{
  n<-nrow(M)
  p<-ncol(M)
  
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  alpha<-matrix(NA,p,5)
  colnames(alpha)<-c("Estimate","pvalue","LB","UB","adjpv")
  rownames(alpha)<-colnames(M)
  beta=gamma=IE<-alpha
  
  # boot==TRUE
  IE.sims<-matrix(NA,sims,p)
  TE.sims<-matrix(NA,sims,p)
  DE.sims<-rep(NA,sims)
  # boot==FALSE
  IE.se<-rep(NA,p)
  for(j in 1:p)
  {
    dat.tmp<-data.frame(X=X,M=M[,j],Y=Y)
    
    fit.m<-lm(M~X,data=dat.tmp)
    fit.y<-lm(Y~X+M,data=dat.tmp)
    
    alpha[j,1:4]<-c(coef(fit.m)[2],summary(fit.m)$coefficients[2,4],confint(fit.m,level=conf.level)[2,])
    gamma[j,1:4]<-c(coef(fit.y)[2],summary(fit.y)$coefficients[2,4],confint(fit.y,level=conf.level)[2,])
    beta[j,1:4]<-c(coef(fit.y)[3],summary(fit.y)$coefficients[3,4],confint(fit.y,level=conf.level)[3,])
    
    if(boot)
    {
      re.med<-mediate(fit.m,fit.y,treat="X",mediator="M",sims=sims,boot=boot,boot.ci.type=boot.ci.type[1],conf.level=conf.level)
      
      IE[j,1:4]<-c(re.med$d1,re.med$d1.p,re.med$d1.ci)
      gamma[j,1:4]<-c(re.med$z1,re.med$z1.p,re.med$z1.ci) 
      
      IE.sims[,j]<-re.med$d1.sims
      TE.sims[,j]<-re.med$d1.sims+re.med$z1.sims
      IE.se[j]<-sd(re.med$d1.sims)
    }else
    {
      IE[j,1]<-alpha[j,1]*beta[j,1]
      IE.se[j]<-sqrt((alpha[j,1]*summary(fit.y)$coefficients[3,2])^2+(summary(fit.m)$coefficients[2,2]*beta[j,1])^2)
      IE[j,2]<-2*pnorm(abs(IE[j,1]/IE.se[j]),lower.tail=FALSE)
      IE[j,c(3,4)]<-c(IE[j,1]-IE.se[j]*qnorm((1-conf.level)/2,lower.tail=FALSE),IE[j,1]+IE.se[j]*qnorm((1-conf.level)/2,lower.tail=FALSE))
    }
  }
  alpha[,5]<-p.adjust(alpha[,2],method=p.adj.method[1])
  beta[,5]<-p.adjust(beta[,2],method=p.adj.method[1])
  gamma[,5]<-p.adjust(gamma[,2],method=p.adj.method[1])
  IE[,5]<-p.adjust(IE[,2],method=p.adj.method[1])
  
  IE.total<-matrix(NA,1,4)
  colnames(IE.total)<-c("Estiamte","pvalue","LB","UB")
  rownames(IE.total)<-"Total"
  DE<-matrix(NA,1,4)
  colnames(DE)<-c("Estimate","pvalue","LB","UB")
  rownames(DE)<-"DE"
  IE.total[1,1]<-sum(IE[,1])
  IE.total.se<-sqrt(sum(IE.se^2))
  IE.total[1,2]<-2*pnorm(abs(IE.total[1,1]/IE.total.se),lower.tail=FALSE)
  if(boot)
  {
    DE.sims<-apply(TE.sims,1,mean,na.rm=TRUE)-apply(IE.sims,1,sum,na.rm=TRUE)
    DE[1,1]<-mean(DE.sims,na.rm=TRUE)
    DE[1,2]<-2*pnorm(abs(DE[1,1]/sd(DE.sims,na.rm=TRUE)),lower.tail=FALSE)
    
    if(boot.ci.type[1]=="bca")
    {
      IE.total[1,c(3,4)]<-BC.CI(apply(IE.sims,1,sum),sims=sims,conf.level=conf.level)
      DE[1,c(3,4)]<-BC.CI(DE.sims,sims=sims,conf.level=conf.level)
    }else
    {
      IE.total[1,c(3,4)]<-quantile(apply(IE.sims,1,sum),probs=c((1-conf.level)/2,1-(1-conf.level)/2))
      DE[1,c(3,4)]<-quantile(DE.sims,probs=c((1-conf.level)/2,1-(1-conf.level)/2))
    }
  }else
  {
    IE.total[1,c(3,4)]<-c(IE.total[1,1]-IE.total.se*qnorm((1-conf.level)/2,lower.tail=FALSE),IE.total[1,1]+IE.total.se*qnorm((1-conf.level)/2,lower.tail=FALSE))
    
    # total effect model
    fit<-lm(Y~X)
    TE.se<-sqrt(summary(fit)$cov.unscaled[2,2]*(summary(fit)$sigma)^2)
    # direct effect
    DE[1,1]<-fit$coefficients[2]-IE.total[1,1]
    DE.se<-sqrt(IE.total.se^2+TE.se^2)
    DE[1,2]<-2*pnorm(abs(DE[1,1]/DE.se),lower.tail=FALSE)
    DE[1,c(3,4)]<-c(DE[1,1]-DE.se*qnorm((1-conf.level)/2,lower.tail=FALSE),DE[1,1]+DE.se*qnorm((1-conf.level)/2,lower.tail=FALSE))
  }
  
  re<-list(IE=IE,DE=DE,alpha=alpha,beta=beta,gamma=gamma,IE.total=IE.total)
  return(re)
}
