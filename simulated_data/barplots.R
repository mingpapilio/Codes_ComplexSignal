###### This plots Fig. 3 and Fig. S5a-d
sd_matrix<- matrix(NA, 5, 6)
rownames(sd_matrix)<- c("(0, 2)","(0.5, 1.5)","(1, 1)","(1.5, 0.5)","(2, 0)")
raw<-read.csv("summary.txt",sep="\t")
## ***Please repeat this section 3 times***
## ***making[i, 1-2] for trait 1 only (t_1 folder), [i, 3-4] for trait 2 only (t_2 folder), [i, 5-6] for both traits (t_both folder)***
for(i in 1:5){
  sd_matrix[i,3]<- sd(raw$typ1_err[((i-1)*100+1):(i*100)])
  sd_matrix[i,4]<- sd(raw$typ2_err[((i-1)*100+1):(i*100)])
}
## ******
## reverse the sequence of trait 2 only simulations
sd_matrix[,3]<- rev(sd_matrix[,3])
sd_matrix[,4]<- rev(sd_matrix[,4])
##
typ1_sd<- sd_matrix[,c(1,3,5)]
typ1_sd<- t(typ1_sd)
typ2_sd<- sd_matrix[,c(2,4,6)]
typ2_sd<- t(typ2_sd)
##
dev.new()
bar_centers<- barplot(typ1, border="white", beside=T, legend=c("trait 1 only", "trait 2 only", "both traits"), xlab="minimal difference (trait 1, trait 2)", ylim=c(0, 0.3))
segments(bar_centers, typ1- typ1_sd, bar_centers, typ1+typ1_sd)
arrows(bar_centers, typ1- typ1_sd, bar_centers, typ1+typ1_sd, lwd = 1.0, angle = 90,
       code = 3, length = 0.05)
##
dev.new()
bar_centers<- barplot(typ2, border="white", beside=T, legend=c("trait 1 only", "trait 2 only", "both traits"), xlab="minimal difference (trait 1, trait 2)", ylim=c(0, 1.0))
segments(bar_centers, typ2- typ2_sd, bar_centers, typ2+typ2_sd)
arrows(bar_centers, typ2- typ2_sd, bar_centers, typ2+typ2_sd, lwd = 1.0, angle = 90,
       code = 3, length = 0.05)


###### This plot Fig. S5e-h
surv_sd<- matrix(NA, 5, 6)
rownames(surv_sd)<- c("(0, 2)","(0.5, 1.5)","(1, 1)","(1.5, 0.5)","(2, 0)")
raw<-read.csv("summary.txt",sep="\t")
## ***Please repeat this section 3 times***
## ***making[i, 1-2] for trait 1 only (t_1 folder), [i, 3-4] for trait 2 only (t_2 folder), [i, 5-6] for both traits (t_both folder)***
for(i in 1:5){
  surv_sd[i,3]<- sd(raw$rej_rate[((i-1)*100+1):(i*100)])
  surv_sd[i,4]<- sd(raw$acp_rate[((i-1)*100+1):(i*100)])
}
## ******
## reverse the sequence of trait 2 only simulations
surv_sd[,3]<- rev(surv_sd[,3])
surv_sd[,4]<- rev(surv_sd[,4])
typ1_sd<- t(surv_sd[,c(1,3,5)])
typ2_sd<- t(surv_sd[,c(2,4,6)])
##
dev.new()
bar_centers<- barplot(1-typ1, border="white", beside=T, legend=c("trait 1 only", "trait 2 only", "both traits"), xlab="minimal difference (trait 1, trait 2)", ylim=c(0, 1.0))
segments(bar_centers, 1-typ1- typ1_sd, bar_centers, 1-typ1+typ1_sd)
arrows(bar_centers, 1-typ1- typ1_sd, bar_centers, 1-typ1+typ1_sd, lwd = 1.0, angle = 90,
       code = 3, length = 0.05)
##
dev.new()
bar_centers<- barplot(typ2, border="white", beside=T, legend=c("trait 1 only", "trait 2 only", "both traits"), xlab="minimal difference (trait 1, trait 2)", ylim=c(0, 1.0))
segments(bar_centers, typ2- typ2_sd, bar_centers, typ2+typ2_sd)
arrows(bar_centers, typ2- typ2_sd, bar_centers, typ2+typ2_sd, lwd = 1.0, angle = 90,
       code = 3, length = 0.05)


###### This plots Fig. 4 and S7
surv_sd<- matrix(NA, 1, 6)
acrj_sd<- matrix(NA, 1, 6)
## rownames(surv_sd)<- c("(0, 2)","(0.5, 1.5)","(1, 1)","(1.5, 0.5)","(2, 0)")
raw<-read.csv("summary.txt",sep="\t")
## ***Please repeat this section 3 times***
## ***making[i, 1-2] for trait 1 only (t_1 folder), [i, 3-4] for trait 2 only (t_2 folder), [i, 5-6] for both traits (t_both folder)***
surv_sd[1, 5]<- sd(raw$typ1_err)
surv_sd[1, 6]<- sd(raw$typ2_err)
acrj_sd[1, 5]<- sd(raw$rej_rate)
acrj_sd[1, 6]<- sd(raw$acp_rate)
## ******
## 
typ1_sd<- t(surv_sd[,c(1,3,5)])
typ2_sd<- t(surv_sd[,c(2,4,6)])
rej_sd<- t(acrj_sd[,c(1,3,5)])
acp_sd<- t(acrj_sd[,c(2,4,6)])
## 
raw<- read.csv("errs_trait_summary.csv")
raw<- read.csv("errs_trait_acrj.csv")
raw<- raw[,2:7]
raw<- t(raw)
raw<- as.matrix(raw)
typ1<- matrix(NA, 1, 3)
typ2<- matrix(NA, 1, 3)
for(i in 1:3){
  typ1[1, i]= mean(raw[((i-1)*2+1)])
  typ2[1, i]= mean(raw[((i-1)*2+2)])
}
##
dev.new()
bar_centers<- barplot(1-typ1, border="white", beside=T, legend=c("trait 1 only", "trait 2 only", "both traits"), xlab="minimal difference (trait 1, trait 2)", ylim=c(0, 1.0))
segments(bar_centers, 1-typ1- rej_sd, bar_centers, 1-typ1+rej_sd)
arrows(bar_centers, 1-typ1- rej_sd, bar_centers, 1-typ1+rej_sd, lwd = 1.0, angle = 90,
       code = 3, length = 0.05)
##
dev.new()
bar_centers<- barplot(typ2, border="white", beside=T, legend=c("trait 1 only", "trait 2 only", "both traits"), xlab="minimal difference (trait 1, trait 2)", ylim=c(0, 1.0))
segments(bar_centers, typ2- acp_sd, bar_centers, typ2+acp_sd)
arrows(bar_centers, typ2- acp_sd, bar_centers, typ2+acp_sd, lwd = 1.0, angle = 90,
       code = 3, length = 0.05)

