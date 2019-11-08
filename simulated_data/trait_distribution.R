## for t_1
raw<-read.csv("summary.txt",sep="\t")

h<- hist(raw$avg_par_size, breaks= 10, plot= F)
h$density= h$counts/sum(h$counts)* 2
dev.new()
plot(h, freq= F, xlim=c(7, 14), ylim=c(0,1), ylab="relative density", xlab="trait value")
points(xx, dnorm(xx, 10, 0.5/1.2), type="l", lty=1)
abline(v= mean(raw$avg_size_thr), col="skyblue")

## for t_2
raw<-read.csv("summary.txt",sep="\t")

h<- hist(raw$avg_par_patt, breaks= 10, plot= F)
h$density= h$counts/sum(h$counts)* 2
dev.new()
plot(h, freq= F, xlim=c(7, 14), ylim=c(0,1), ylab="relative density", xlab="trait value")
points(xx, dnorm(xx, 10, 0.5/1.2), type="l", lty=1)
abline(v= mean(raw$avg_patt_thr), col="skyblue")

## for t_both
raw<-read.csv("summary.txt",sep="\t")
length<- 500
sort<- matrix(NA, length, 6)
sort_remain<- matrix(NA, length, 6)
colnames(sort)<- c("thr1", "par1", "thr2", "par2", "typ1_err", "typ2_err")
colnames(sort_remain)<- c("thr1", "par1", "thr2", "par2", "typ1_err", "typ2_err")
## sorting the traits, making the one with smaller acceptnace threshold as trait 1
for (i in 1:length){
  tmp1<- raw$avg_par_size[i]- raw$avg_size_thr[i]
  tmp2<- raw$avg_par_patt[i]- raw$avg_patt_thr[i]
  if(raw$avg_size_thr[i]<= raw$avg_patt_thr[i]){
    # if(raw$avg_size_thr[i]<= 11){
    # if(tmp1> tmp2){
    sort[i,1]<- raw$avg_size_thr[i]
    sort[i,2]<- raw$avg_par_size[i]
    sort[i,3]<- raw$avg_patt_thr[i]
    sort[i,4]<- raw$avg_par_patt[i]
    sort[i,5]<- raw$typ1_err[i]
    sort[i,6]<- raw$typ2_err[i]
  }
  else{
    #if(raw$avg_patt_thr[i]<= 11){
    sort[i,1]<- raw$avg_patt_thr[i]
    sort[i,2]<- raw$avg_par_patt[i]
    sort[i,3]<- raw$avg_size_thr[i]
    sort[i,4]<- raw$avg_par_size[i]
    sort[i,5]<- raw$typ1_err[i]
    sort[i,6]<- raw$typ2_err[i]
    #}
  }
}
temp<- c(mean(sort[,1], na.rm=T), mean(sort[,2], na.rm=T), mean(sort[,3], na.rm=T), mean(sort[,4], na.rm=T),mean(sort[,5], na.rm=T),mean(sort[,6], na.rm=T))
temp

###
xx<- seq(0, 15, length.out = 1000)

h<- hist(as.data.frame(sort)$par1, breaks= 40, plot= F)
h$density= h$counts/sum(h$counts)* 2
dev.new()
plot(h, freq= F, xlim=c(7, 14), ylim=c(0,1), ylab="relative density", xlab="trait value")
points(xx, dnorm(xx, 10, 0.5/1.2), type="l", lty=1)
abline(v= temp[1], col="orange")

h<- hist(as.data.frame(sort)$par2, breaks= 40, plot= F)
h$density= h$counts/sum(h$counts)* 2
dev.new()
plot(h, freq= F, xlim=c(7, 14), ylim=c(0,1), ylab="relative density", xlab="trait value")
points(xx, dnorm(xx, 10, 0.5/1.2), type="l", lty=1)
abline(v= temp[3], col="orange")

