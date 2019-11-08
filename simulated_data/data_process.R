## Initializing matrices
  errs_com_mt1<- matrix(NA, 5, 6)
  colnames(errs_com_mt1)= c("t1_t1","t1_t2","t2_t1","t2_t2","both_t1","both_t2")
  errs_com_mt2<- matrix(NA, 5, 6)
  colnames(errs_com_mt2)= c("t1_t1","t1_t2","t2_t1","t2_t2","both_t1","both_t2")
  errs_com_mtb<- matrix(NA,5, 6)
  colnames(errs_com_mtb)= c("size","size_thr","par_size","patt","patt_thr","par_patt")
  errs_com_msim<- matrix(NA,5, 6)
  colnames(errs_com_msim)= c("t1_t1","t1_t2","t2_t1","t2_t2","both_t1","both_t2")
  errs_com_ms2<- matrix(NA, 5, 6)
  colnames(errs_com_ms2)= c("t1_t1","t1_t2","t2_t1","t2_t2","both_t1","both_t2")
  s7= 0
## Writing the matrices
  write.csv(errs_com_msim, file="errs_trait_summary.csv")
  write.csv(errs_com_mt1, file="errs_trait_1.csv")
  write.csv(errs_com_mt2, file="errs_trait_2.csv")
  write.csv(errs_com_mtb, file="errs_trait_both.csv")
  if(s7==1) write.csv(errs_com_ms2, file="errs_trait_acrj.csv")

## For trait 1 only
raw<-read.csv("summary.txt",sep="\t")
for(i in 1:5){
  start= (1+100*(i-1))
  end= 100*i
  errs_com_mt1[i,1]= format(round(mean(raw$avg_size[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt1[i,2]= format(round(mean(raw$avg_size_thr[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt1[i,3]= format(round(mean(raw$avg_par_size[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt1[i,4]= format(round(mean(raw$avg_patt[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt1[i,5]= format(round(mean(raw$avg_patt_thr[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt1[i,6]= format(round(mean(raw$avg_par_patt[start:end], na.rm= T),4), nsmall=4)
}
for(i in 1:5){
  start= (1+100*(i-1))
  end= 100*i
  errs_com_msim[i,1]= format(round(mean(raw$typ1_err[start:end], na.rm= T),5), nsmall=5)
  errs_com_msim[i,2]= format(round(mean(raw$typ2_err[start:end], na.rm= T),5), nsmall=5)
}
if(s7==1){
  for(i in 1:5){
    start= (1+100*(i-1))
    end= 100*i
    errs_com_ms2[i,1]= format(round(mean(raw$rej_rate[start:end], na.rm= T),4), nsmall=4)
    errs_com_ms2[i,2]= format(round(mean(raw$acp_rate[start:end], na.rm= T),4), nsmall=4)
  }
}

## For trait 2 only
raw<-read.csv("summary.txt",sep="\t")
for(i in 1:5){
  start= (1+100*(i-1))
  end= 100*i
  errs_com_mt2[i,1]= format(round(mean(raw$avg_size[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt2[i,2]= format(round(mean(raw$avg_size_thr[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt2[i,3]= format(round(mean(raw$avg_par_size[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt2[i,4]= format(round(mean(raw$avg_patt[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt2[i,5]= format(round(mean(raw$avg_patt_thr[start:end], na.rm= T),4), nsmall=4)
  errs_com_mt2[i,6]= format(round(mean(raw$avg_par_patt[start:end], na.rm= T),4), nsmall=4)
}
for(i in 1:5){
  start= (1+100*(i-1))
  end= 100*i
  errs_com_msim[i,3]= format(round(mean(raw$typ1_err[start:end], na.rm= T),5), nsmall=5)
  errs_com_msim[i,4]= format(round(mean(raw$typ2_err[start:end], na.rm= T),5), nsmall=5)
}
if(s7==1){
  for(i in 1:5){
    start= (1+100*(i-1))
    end= 100*i
    errs_com_ms2[i,3]= format(round(mean(raw$rej_rate[start:end], na.rm= T),4), nsmall=4)
    errs_com_ms2[i,4]= format(round(mean(raw$acp_rate[start:end], na.rm= T),4), nsmall=4)
  }
}

## For both traits 
raw<-read.csv("summary.txt",sep="\t")
for(i in 1:5){
  start= (1+100*(i-1))
  end= 100*i
  errs_com_mtb[i,1]= format(round(mean(raw$avg_size[start:end], na.rm= T),4), nsmall=4)
  errs_com_mtb[i,2]= format(round(mean(raw$avg_size_thr[start:end], na.rm= T),4), nsmall=4)
  errs_com_mtb[i,3]= format(round(mean(raw$avg_par_size[start:end], na.rm= T),4), nsmall=4)
  errs_com_mtb[i,4]= format(round(mean(raw$avg_patt[start:end], na.rm= T),4), nsmall=4)
  errs_com_mtb[i,5]= format(round(mean(raw$avg_patt_thr[start:end], na.rm= T),4), nsmall=4)
  errs_com_mtb[i,6]= format(round(mean(raw$avg_par_patt[start:end], na.rm= T),4), nsmall=4)
}
for(i in 1:5){
  start= (1+100*(i-1))
  end= 100*i
  errs_com_msim[i,5]= format(round(mean(raw$typ1_err[start:end], na.rm= T),5), nsmall=5)
  errs_com_msim[i,6]= format(round(mean(raw$typ2_err[start:end], na.rm= T),5), nsmall=5)
}
if(s7==1){
  for(i in 1:5){
    start= (1+100*(i-1))
    end= 100*i
    errs_com_ms2[i,5]= format(round(mean(raw$rej_rate[start:end], na.rm= T),4), nsmall=4)
    errs_com_ms2[i,6]= format(round(mean(raw$acp_rate[start:end], na.rm= T),4), nsmall=4)
  }
}

###############################################################################
## This seciton plots the phenotypic distributions of the desirables and undesirables, for Fig. 1
xx<- seq(5, 15, length.out = 1000)

dev.new()
host<- 10
diff<- 0.237
thr1<- 0.36
thr2<- 1.301
plot(xx, dnorm(xx, host, 0.5/1.2), type="l", ylim=c(0, 1.5), xlim=c(7.5, 14.5), ylab="relative probability")
points(xx, dnorm(xx, host+diff, 0.5*1.2), type="l", col="lightgray")
abline(v=(host+thr1), lty= 2)
abline(v=(host+thr2), lty= 1, col="orange")

xx<- seq(5, 15, length.out = 1000)

dev.new()
host<-  10
par<-   11.5
thr1<-  10.71
thr2<-  10.74
plot(xx, dnorm(xx, host, 0.5/1.2), type="l", ylim=c(0, 1.5), xlim=c(7.5, 14.5), ylab="relative probability")
points(xx, dnorm(xx, par, 0.5*1.2), type="l", lty=2)
abline(v=thr1, col="skyblue")
abline(v=thr2, col="orange")

dev.new()
plot(xx, dnorm(xx, 10, 0.5/1.2), type="l", ylim=c(0, 1.5), xlim=c(7.5, 15))
points(xx, dnorm(xx, 10+1, 0.5*1.2), type="l", col="lightgray")
abline(v=(10+0.48), lty=2)

dev.new()
plot(xx, dnorm(xx, 30, 1.5/1.2), type="l", ylim=c(0, 0.5), xlim=c(22.5, 45))
points(xx, dnorm(xx, 30+6, 1.5*1.2), type="l", col="lightgray")
abline(v=(30+2.513), lty=2)

dev.new()
plot(xx, dnorm(xx, 10, 0.5/1.2), type="l", ylim=c(0, 1.5), xlim=c(7.5, 15))
points(xx, dnorm(xx, 10, 0.5*1.2), type="l", col="lightgray")
points(xx, dnorm(xx, 10.5, 0.5*1.2), type="l", col="lightgray")
points(xx, dnorm(xx, 11, 0.5*1.2), type="l", col="lightgray")
points(xx, dnorm(xx, 11.5, 0.5*1.2), type="l", col="lightgray")
points(xx, dnorm(xx, 12, 0.5*1.2), type="l", col="lightgray")

