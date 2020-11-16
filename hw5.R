##Method 3
rm(list=ls())
library("TTR")
library("FITSio")
#ori=readFrameFromFITS("cB58_Lyman_break.fit")


target<-readFrameFromFITS("cB58_Lyman_break.fit")
scaler=sd(10^(target$LOGLAM))
pre_processs<-function(spec){
  spec_1=spec[spec$or_mask==0,]
  spec_1=spec
  fluxl=spec_1$FLUX
  wvl=10^spec_1$LOGLAM
  scaled_flux=(fluxl-mean(fluxl))/sd(fluxl)
  smoothed_flux<-SMA(scaled_flux,n=6)
  spec_set<-data.frame(flux=smoothed_flux,lam=wvl)
  spec_set<-na.omit(spec_set)
  spec_set
}

find_der0<-function(series,x_series,scaler){
  diffe<-diff(series)
  diffe<-na.omit(diffe)
  n=1
  end=length(diffe)-2
  points<-rep(NA,end+2)
  point<-rep(NA,end+2)
  for (i in c(2:end)){
    one_step=diffe[i]*diffe[i+1]
    #two_step=(diffe[i-1])*(diffe[i+2])
    if (one_step<0){
      points[n]<-x_series[i+1]
      point[n]<-series[i+1]
    }
    n=n+1
  }
  pointt=data.frame(x=points,y=scaler*point)
  pointt<-na.omit(pointt)#}
  pointt
}

spec_set_d<-pre_processs(spec=target)
x1_d<-find_der0(series = spec_set_d$flux,x_series=spec_set_d$lam,scaler = scaler)
ranked_d<-x1_d
#ranked_d<-get_simplified(x1_d)
#target_vec<-get_vectors(ranked_d)
ori2=ranked_d[,2]
n=length(ori2)


#files=list.files("data/")
#data=c()
#for (i in 1:5){
#  data[[i]]<-readFrameFromFITS(paste0("data/",files[i]))
#}

pre_process<-function(spect){
  spec<-na.omit(spect)
  spec[spec[,4]!=0,1]=median(spec[spec[,4]==0,1])
  fluxl<-spec$flux
  wvl<-10^spec$loglam
  scaled_flux<-(fluxl-mean(fluxl,na.rm = T))/sd(fluxl,na.rm = T)
  scaled_flux<-na.omit(scaled_flux)
  smoothed_flux<-SMA(scaled_flux,n=6)
  spec_set<-data.frame(flux=smoothed_flux,lam=wvl)
  spec_set<-na.omit(spec_set)
  spec_set
}
#find 3 closest
r=data.frame()
for(i in 1:3) r[[i]]=as.numeric()

name=list.files("data/")
#spectrumID<-rep(NA,length(names))
#distance<-rep(NA,length(names))
#i<-rep(NA,length(names))
#distance_rank=data.frame(spectrumID=spectrumID,distance=distance,i=i)
for(j in name[1:50]){
  dat=readFrameFromFITS(paste("data/",j,sep=""))
  m=dim(dat)[1]
  #dat[dat[,4]!=0,1]=NA
  #dat<-na.omit(dat)
  if (m-n>2){
    x.inv <- try(pre_process(spect=dat), silent=TRUE)
    if ('try-error' %in% class(x.inv)) next
    spect<-pre_process(spect=dat)
    derived<-find_der0(series = spect$flux,x_series=spect$lam,scaler = scaler)
    dat<-derived
  corr=rep(0,m-n+1)
  for(i in 1:(m-n+1)){
    dat2=dat[i:(n+i-1),2]
    #dat2[dat2[,4]!=0,1]=mean(dat2[dat2[,4]==0,1]) #replace and_mark!=0 with the mean of and_mark==0
    
    #dats=(dat2[,1]-mean(dat2[,1],na.rm = T))/sd(dat2[,1],na.rm = T)
    minus=dat2-ori2
    minus<-na.omit(minus)
    sum=0
    l=length(ori2)
    for (k in 1:l){
      sum=sum+(minus[k])^2
    }
    corr[i]=sum
  }
  r[j,1]=min(abs(corr[!is.na(corr)]))
  r[j,2]=j
  r[j,3]=which(abs(corr)==min(abs(corr[!is.na(corr)])))[1]-1 }
  else{
    r[j,1]=NA
    r[j,2]=NA
    r[j,3]=NA
  }
}

#export hw2.csv
r0=r[order(r[,1],decreasing=FALSE),]
rownames(r0)=1:dim(r0)[1]
colnames(r0)=c("distance","spectrumID","i")
#write.csv(r0,"hw2_2.csv",row.names=FALSE)