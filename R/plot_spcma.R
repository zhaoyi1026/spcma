plot_spcma <-
function(object,plot.coef=c("alpha","beta","IE"),cex.lab=1,cex.axis=1,pt.cex=1,...)
{
  plot.idx<-which(names(object)==plot.coef)
  
  out<-as.matrix(object[[plot.idx]][,c(1,3,4)])
  colnames(out)<-c("Estimate","LB","UB")
  
  K<-nrow(out)
  
  sig<-as.numeric(out[,2]*out[,3]>0)
  neg<-as.numeric(out[,1]<0)*sig*3
  pos<-as.numeric(out[,1]>0)*sig*2
  pt.type<-(1-sig)+pos+neg
  col.tmp<-c(1,2,4)[pt.type]
  pt.tmp<-c(19,17,15)[pt.type]
  
  par(mar=c(5,5,3,3))
  plot(range(1-0.5/K,K+0.5/K),range(out[,c(2,3)],na.rm=TRUE),type="n",xaxt="n",xlab="",ylab=plot.coef,cex.lab=cex.lab,cex.axis=cex.axis)
  axis(side=1,at=1:K,labels=rownames(out),cex.axis=cex.axis)
  abline(h=0,lty=2,lwd=2,col=8)
  points(1:K,out[,1],pch=pt.tmp,col=col.tmp,cex=pt.cex)
  for(j in 1:K)
  {
    lines(rep(j,2),out[j,c(2,3)],lty=2,lwd=1,col=col.tmp[j])
    lines(c(j-0.3/K,j+0.3/K),rep(out[j,2],2),lty=1,lwd=1,col=col.tmp[j])
    lines(c(j-0.3/K,j+0.3/K),rep(out[j,3],2),lty=1,lwd=1,col=col.tmp[j])
  }
}
