# Function to infer clusters separation using the sillouthe index (from Taskesen et al ).
dbscan_SH <- function(data,eps=NULL,showplot=F,prop_outliers=.1,eps_res=500,eps_range=NULL){
  data<-data.frame(data)
  if(is.null(eps)){
    cat("Optimising eps: ")
    if(is.null(eps_range)){
      eps_scale<-mean(apply(data,2,stats::sd)) # makes the search scale independent
      epsvec<-seq(0,4,length.out=eps_res)*eps_scale # space to search for eps parameter
    } else epsvec<-seq(eps_range[1],eps_range[2],length.out=eps_res)
    silvec<-numeric(length(epsvec))
    for(i in 1:length(epsvec)){
      eps<-epsvec[i]
      DBcl<-dbscan::dbscan(data,eps)
      cl<-DBcl$cluster
      cat(".")
      if(all(cl==1)) break else if(max(cl)==1) silvec[i]<-0 else
        if(all(cl==0)) silvec[i]<-0 else
          if(mean(cl==0)>prop_outliers) silvec[i]<-0 else{
            S<-cluster::silhouette(x=cl[cl!=0],dist=stats::dist(data[cl!=0,])) # exclude the 0's
            silvec[i]<-summary(S)$avg.width
          }
    }
    cat("\n")
    if(showplot){
      if(ncol(data)==2) graphics::par(mfrow=c(1,2))
      end<-length(silvec)-which(cumsum(rev(silvec))>0)[1]+10
      graphics::plot(epsvec[1:end],silvec[1:end],xlab="eps value",ylab="silhouette score")
    }
    eps<-epsvec[which.max(silvec)]
    if(showplot){
      graphics::abline(h=max(silvec),lty=2)
      graphics::abline(v=eps,lty=2)
    }
  }
  DBcl<-dbscan::dbscan(data,eps)$cluster
  if(showplot){
    graphics::plot(data,col=c("lightgrey",sample(grDevices::rainbow(max(DBcl),v=.8)))[DBcl+1],pch=16)
    graphics::par(mfrow=c(1,1))
  }
  cat("Used eps: ",eps,"\n")
  return(list(cluster=DBcl,eps=eps,SIL=max(silvec)))
}
