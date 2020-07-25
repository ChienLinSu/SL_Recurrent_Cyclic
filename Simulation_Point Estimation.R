############### Load R package ##########################

library(survival)
library(nleqslv)
library(parallel)

########## function to find the interval for a given time T #########

Find_interval=function(x){ test=c()
for(i in 1:dim(renewal_times)[1]){
test[i]=ifelse(findInterval(x, renewal_times[i,])==1,1,0)
} ## for i
loc=which(test==1)
return(renewal_times[loc,][1])
}## for function



########## Parameter setting ############################


Re.type.data=c()        ## save data

n=50                    ## number of individuals
c=22                    ## censoring variable ~Uniform (0,c)
sigma=1                 ## variance for random effect for Q
true_beta1=-0.3         ## coefficient for beta 1
true_beta2=-0.2         ## coefficient for beta 2
baseline_type1=0.15     ## baseline rate for type 1
baseline_type2=0.2      ## baseline rate for type 2


t2=2
time1_true=t2*baseline_type1
time2_true=t2*baseline_type2

renewal_times=matrix(c(0,2,2,4,4,6,6,8,8,10,10,12,12,14,14,16,16,18,18,20,20,22,22,24),ncol=2,byrow=TRUE)


##################### Generate recurrent event data with two event types  ######################################################################
##################### As mentioned in the simulation section of the paper, here we only consider one Binary covariate X1 to  ###################
##################### generate recurrent event data. However, to show how to incorporate other covariate                     ################### 
##################### into the data frame we generated, we thus simulate extra covariate X2 as an example. #####################################

Data=c()
p.id<-1

while (p.id<=n){

rr.rec1=c()
rr.rec2=c()

X1=rbinom(1,1,0.5)  ## covariate 1

X2=runif(1,0,1)     ## covariate 2: we don't incorporate this covariate into the event time generation
censor.time=runif(1,0,c)
Q=1
##Q=rgamma(1,shape=1/sigma,scale=sigma) 
T1<-T1.new<-0

while(T1.new<censor.time){
a=runif(1)
gene.function1=function(x){
exp(-1*Q*exp(sum(c(true_beta1)*c(X1)))*baseline_type1*((x-Find_interval(x))-(T1.new-Find_interval(T1.new))))-a
}
T1.new<-nleqslv(T1.new+0.1,gene.function1)$x


T1<-c(T1,T1.new)
}

rec.type1<-c(T1[T1<censor.time],censor.time) # observed type 1 recurrent events 

 for (j in 1:(length(rec.type1)-1)) {
      if (j==(length(rec.type1)-1)) {
        rr.rec1<-rbind(rr.rec1, c(p.id, rec.type1[j:(j+1)],0,X1,0,X2,0,1,censor.time)) 
      } 
      else {
        rr.rec1<-rbind(rr.rec1, c(p.id, rec.type1[j:(j+1)], 1, X1,0,X2,0, 1, censor.time))
      }
    }

rr.rec1.sort=rr.rec1[order(rr.rec1[,3]),] 

T2<-T2.new<-0
while(T2.new<censor.time){
b=runif(1)
gene.function2=function(x){
exp(-1*Q*exp(sum(c(true_beta2)*c(X1)))*baseline_type2*((x-Find_interval(x))-(T2.new-Find_interval(T2.new))))-b
}

T2.new<-nleqslv(T2.new+0.1,gene.function2)$x
T2<-c(T2,T2.new)
}

rec.type2<-c(T2[T2<censor.time],censor.time) # observed type 2 recurrent events 

 for (j in 1:(length(rec.type2)-1)) {
      if (j==(length(rec.type2)-1)) {
        rr.rec2<-rbind(rr.rec2, c(p.id, rec.type2[j:(j+1)],0,0,X1,0,X2,2,censor.time)) 
      } 
      else {
        rr.rec2<-rbind(rr.rec2, c(p.id, rec.type2[j:(j+1)], 1, 0,X1,0,X2, 2, censor.time))
      }
    }

total.rr.rec<-data.frame(rbind(rr.rec1,rr.rec2))
names(total.rr.rec)=c("p.id", "start.time", "end.time", "Censor", "X1_1", "X1_2","X2_1","X2_2", "type", "censoring")
    

Data=rbind(Data,total.rr.rec)
p.id <- p.id + 1

}  ##for while
Original_Data=Data


length(which(Data$type==1))/n ## average number of recurrent event of type 1
length(which(Data$type==2))/n ## average number of recurrent event of type 2


## use coxph to obtain marginal estimates
result<-coxph(Surv(start.time, end.time, Censor)~X1_1+X1_2+cluster(p.id)+strata(type), data=Data, robust=T)
as.numeric(summary(result)$coefficients[,1]) ## marginal estimates
as.numeric(summary(result)$coefficients[,4]) ## standard error 




Centime=c()
Cov1=c()
Cov2=c()
for(i in 1:n){
Centime[i]=Data[,10][which(Data[,1]==i)[1]]
Cov1[i]=Data[,5][which(Data[,9]==1 & Data[,1]==i)[1]]
Cov2[i]=Data[,6][which(Data[,9]==2 & Data[,1]==i)[1]]
}


######################### baseline mean function based on Cai and Schaubel(2004) ###############################

Cai_Type_baseline=function(t,type,p,covariate){ ##t:time; type: event type; p: dimension of covariate
if(t==0){
total=0
}else{

V1=unique(sort(c(0,Data[which(Data[,"type"]==type & Data[,"Censor"]==1 ),"end.time"])))
V1=V1[which(V1<t)]

XX=matrix(c(1:length(V1)),1,length(V1))
Ratio=mclapply(XX,function(h){ 
df=data.frame(x=ifelse(Data[which(Data[,"type"]==type),"end.time"]==V1[h]  & Data[which(Data[,"type"]==type),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(V1[h]<=Centime)*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))####
ratio=denominator/numerator
}
)
total=sum(unlist(Ratio))
}## for else
return(max(0,total))
}## for function


Cai_Type_baseline(t2,type=1,p=1,covariate=Cov1)
Cai_Type_baseline(t2,type=2,p=2,covariate=Cov1)




#########################################################################
##              Variance estimation for mean function                   #
#########################################################################

############# for A  ######################

Type1_covZ=list()
Type2_covZ=list()

for(i in 1:n){
Type1_covZ[[i]]=matrix(c(Cov1[i],0),2,1)
Type2_covZ[[i]]=matrix(c(0,Cov1[i]),2,1)
}


Type1.uncensored.time=round(Data[,"end.time"][which(Data[,"type"]==1 & Data[,"Censor"]==1)],6)
Type2.uncensored.time=round(Data[,"end.time"][which(Data[,"type"]==2 & Data[,"Censor"]==1)],6)

XX_type1=matrix(c(1:length(Type1.uncensored.time)),1,length(Type1.uncensored.time))
S0_type1=mclapply(XX_type1,function(h){
ans=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type1=unlist(S0_type1)

XX_type2=matrix(c(1:length(Type2.uncensored.time)),1,length(Type2.uncensored.time))
S0_type2=mclapply(XX_type2,function(h){
ans=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type2=unlist(S0_type2)


S1_type1=mclapply(XX_type1,function(h){
a1=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)


S1_type2=mclapply(XX_type2,function(h){
f1=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
f2=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(f1,f2),2,1)
}
)

S2_type1=mclapply(XX_type1,function(h){
a11=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*(unlist(lapply(Type1_covZ,'[[',1))^2)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a12=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*(unlist(lapply(Type1_covZ,'[[',1))*unlist(lapply(Type1_covZ,'[[',2)))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a22=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*(unlist(lapply(Type1_covZ,'[[',2))^2)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a11,a12,a12,a22),2,2,byrow=TRUE)
}
)

S2_type2=mclapply(XX_type2,function(h){
g11=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*(unlist(lapply(Type2_covZ,'[[',1))^2)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
g12=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*(unlist(lapply(Type2_covZ,'[[',1))*unlist(lapply(Type2_covZ,'[[',2)))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
g22=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*(unlist(lapply(Type2_covZ,'[[',2))^2)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(g11,g12,g12,g22),2,2,byrow=TRUE)
}
)


aq1=mclapply(XX_type1,function(k){
(((S2_type1[[k]]/S0_type1[k])-((S1_type1[[k]]/S0_type1[k])%*%t((S1_type1[[k]]/S0_type1[k])))))
}
)

aq2=mclapply(XX_type2,function(k){
(((S2_type2[[k]]/S0_type2[k])-((S1_type2[[k]]/S0_type2[k])%*%t((S1_type2[[k]]/S0_type2[k])))))
}
)

Hat_A_type1=Reduce('+', aq1)
Hat_A_type2=Reduce('+', aq2)
Hat_A_type=Hat_A_type1+Hat_A_type2



############### Calculate hat H in Cai and Schaubel (2004) #####################
H_type=function(t,type){

if(type==1){
Type1.uncensored.time=round(Data[,"end.time"][which(Data[,"type"]==1 & Data[,"Censor"]==1)],6)
Type1.less.t=Type1.uncensored.time[which(Type1.uncensored.time<=t)]
XX_type1=matrix(c(1:length(Type1.less.t)),1,length(Type1.less.t))

S0_type1=mclapply(XX_type1,function(h){
ans=sum(as.numeric(Type1.less.t[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type1=unlist(S0_type1)


S1_type1=mclapply(XX_type1,function(h){
a1=sum(as.numeric(Type1.less.t[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(Type1.less.t[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)

Ratio_type1=mclapply(XX_type1,function(h){ 
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==1),"end.time"],6)==round(Type1.less.t[h],6)  & Data[which(Data[,"type"]==1),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==1),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(Type1.less.t[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))####
ratio=denominator/numerator
}
)

Ans1=-1*sum((1/S0_type1)*unlist(lapply(S1_type1,'[[',1))*unlist(Ratio_type1))
Ans2=-1*sum((1/S0_type1)*unlist(lapply(S1_type1,'[[',2))*unlist(Ratio_type1))
Ans=matrix(c(Ans1,Ans2),2,1)

}else{

Type2.uncensored.time=round(Data[,"end.time"][which(Data[,"type"]==2 & Data[,"Censor"]==1)],6)
Type2.less.t=(Type2.uncensored.time[which(Type2.uncensored.time<=t)])

XX_type2=matrix(c(1:length(Type2.less.t)),1,length(Type2.less.t))
S0_type2=mclapply(XX_type2,function(h){
ans=sum(as.numeric(Type2.less.t[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type2=unlist(S0_type2)



S1_type2=mclapply(XX_type2,function(h){
f1=sum(as.numeric(Type2.less.t[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
f2=sum(as.numeric(Type2.less.t[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(f1,f2),2,1)
}
)

Ratio_type2=mclapply(XX_type2,function(h){ 
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==2),"end.time"],6)==round(Type2.less.t[h],6)  & Data[which(Data[,"type"]==2),"Censor"]==1,1,0),id=rep(1:n,as.vector(table(Data[which(Data[,"type"]==2),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(Type2.less.t[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))####
ratio=denominator/numerator
}
)

Ans1=-1*sum((1/S0_type2)*unlist(lapply(S1_type2,'[[',1))*unlist(Ratio_type2))
Ans2=-1*sum((1/S0_type2)*unlist(lapply(S1_type2,'[[',2))*unlist(Ratio_type2))
Ans=matrix(c(Ans1,Ans2),2,1)
}## for else

return(Ans)
}## for function

##H_type(t=2,type=1)




################################ Calculate sum U in Eq.(9) of Cai and Schaubel (2004)#######################################
########################################################################################

Sum_est_U=mclapply(matrix(c(1:n),1,n),function(i){

data_i=Data[which(Data[,"p.id"]==i),]
obs.time=data_i[which(data_i[,"type"]==1),"end.time"]
delta.i=data_i[which(data_i[,"type"]==1),"Censor"]
est_beta=as.numeric(summary(result)$coefficients[,1])
XX_type1_i=matrix(c(1:length(obs.time)),1,length(obs.time))

S0_type1_i=mclapply(XX_type1_i,function(h){
ans=sum(as.numeric(obs.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
}
)
S0_type1_i=unlist(S0_type1_i)


S1_type1_i=mclapply(XX_type1_i,function(h){
a1=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)


U_i_type1_part1=mclapply(matrix(c(1:length(delta.i)),1,length(delta.i)),function(s){
comp1=(Type1_covZ[[i]]-(S1_type1_i[[s]]/S0_type1_i[s]))*delta.i[s]
}## function s
)
U_i_type1_part1=Reduce('+', U_i_type1_part1)


XX_type1=matrix(c(1:length(Type1.uncensored.time)),1,length(Type1.uncensored.time))
S0_type1=mclapply(XX_type1,function(h){
ans=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type1=unlist(S0_type1)

S1_type1=mclapply(XX_type1,function(h){
a1=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)

Ratio_type1=mclapply(XX_type1,function(h){ 
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==1),"end.time"],6)==round(Type1.uncensored.time[h],6)  & Data[which(Data[,"type"]==1),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==1),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))####
ratio=denominator/numerator
}
)


U_i_type1_part2=mclapply(matrix(c(1:length(XX_type1)),1,length(XX_type1)),function(s){
comp1=(Type1_covZ[[i]]-(S1_type1[[s]]/S0_type1[s]))*as.numeric(Centime[i]>=Type1.uncensored.time[s])*as.numeric(exp(t(Type1_covZ[[i]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))*unlist(Ratio_type1)[s]
}## function s
)
U_i_type1_part2=Reduce('+', U_i_type1_part2)


##U_i_type1_part1-U_i_type1_part2

data_i=Data[which(Data[,"p.id"]==i),]
obs.time=data_i[which(data_i[,"type"]==2),"end.time"]
delta.i=data_i[which(data_i[,"type"]==2),"Censor"]
est_beta=as.numeric(summary(result)$coefficients[,1])
XX_type2_i=matrix(c(1:length(obs.time)),1,length(obs.time))

S0_type2_i=mclapply(XX_type2_i,function(h){
ans=sum(as.numeric(obs.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
}
)
S0_type2_i=unlist(S0_type2_i)


S1_type2_i=mclapply(XX_type2_i,function(h){
a1=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)

U_i_type2_part1=mclapply(matrix(c(1:length(delta.i)),1,length(delta.i)),function(s){

comp1=(Type2_covZ[[i]]-(S1_type2_i[[s]]/S0_type2_i[s]))*delta.i[s]
}## function s
)
U_i_type2_part1=Reduce('+', U_i_type2_part1)

XX_type2=matrix(c(1:length(Type2.uncensored.time)),1,length(Type2.uncensored.time))

S0_type2=mclapply(XX_type2,function(h){
ans=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type2=unlist(S0_type2)

S1_type2=mclapply(XX_type2,function(h){
a1=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)

Ratio_type2=mclapply(XX_type2,function(h){ 
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==2),"end.time"],6)==round(Type2.uncensored.time[h],6)  & Data[which(Data[,"type"]==2),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==2),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))####
ratio=denominator/numerator
}
)

U_i_type2_part2=mclapply(matrix(c(1:length(XX_type2)),1,length(XX_type2)),function(s){
comp1=(Type2_covZ[[i]]-(S1_type2[[s]]/S0_type2[s]))*as.numeric(Centime[i]>=Type2.uncensored.time[s])*as.numeric(exp(t(Type2_covZ[[i]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))*unlist(Ratio_type2)[s]
}## function s
)
U_i_type2_part2=Reduce('+', U_i_type2_part2)

U_sum=(U_i_type1_part1-U_i_type1_part2)+(U_i_type2_part1-U_i_type2_part2)

}## for function i
)##






#####################  Variance estimation for the mean function proposed by Cai and Schaubel (2004) ##########################
Cai_var=function(t,type,CovZ){### t:time; type: event type;  CovZ: covariate  
if(type==1){
type.obstime=Type1.uncensored.time
}else{
type.obstime=Type2.uncensored.time
}


Total=mclapply(matrix(c(1:n),1,n),function(i){
data_i=Data[which(Data[,"p.id"]==i),]
loc=which(data_i[which(data_i[,"type"]==type),"end.time"]<=t)
if(length(loc)==0){
C_part11=0
}else{
obs.time=data_i[which(data_i[,"type"]==type),"end.time"][loc]
delta.i=data_i[which(data_i[,"type"]==type),"Censor"][loc]
est_beta=as.numeric(summary(result)$coefficients[,1])
XX_type_i=matrix(c(1:length(obs.time)),1,length(obs.time))


S0_type_i=mclapply(XX_type_i,function(h){
ans=mean(as.numeric(obs.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(CovZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
}
)
S0_type_i=unlist(S0_type_i)

C_part11=mclapply(matrix(c(1:length(delta.i)),1,length(delta.i)),function(s){
(1/S0_type_i[s])*delta.i[s]
}## function s
)
}


type.obstime.lesst=type.obstime[which(type.obstime<=t)]

XX_type=matrix(c(1:length(type.obstime.lesst)),1,length(type.obstime.lesst))
S0_type=mclapply(XX_type,function(h){
ans=mean(as.numeric(type.obstime.lesst[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(CovZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type=unlist(S0_type)


Ratio_type.t=mclapply(XX_type,function(h){ 
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==type),"end.time"],6)==round(type.obstime.lesst[h],6)  & Data[which(Data[,"type"]==type),"Censor"]==1,1,0),id=rep(1:n,as.vector(table(Data[which(Data[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(type.obstime.lesst[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(CovZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
ratio=denominator/numerator
}
)


C_part12=mclapply(matrix(c(1:length(XX_type)),1,length(XX_type)),function(s){
comp1=(1/S0_type[s])*as.numeric(Centime[i]>=type.obstime.lesst[s])*as.numeric(exp(t(CovZ[[i]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))*unlist(Ratio_type.t)[s]
}## function s
)
C_part12=Reduce('+', C_part12)

Part1_C=sum(unlist(C_part11))-sum(unlist(C_part12))
total_sum=(Part1_C+t(H_type(t,type))%*%solve((1/n)*Hat_A_type)%*%Sum_est_U[[i]])^2
}## for function Total
)
return(mean(unlist(Total)))
}


#################### Cyclic Mean function and variance estimation ##################
####################################################################################

Type1.uncensored.time=round(Data[,"end.time"][which(Data[,"type"]==1 & Data[,"Censor"]==1)],6)
Type2.uncensored.time=round(Data[,"end.time"][which(Data[,"type"]==2 & Data[,"Censor"]==1)],6)

Type1.delta=Data[,4][which(Data[,9]==1)]
Type2.delta=Data[,4][which(Data[,9]==2)]

Type1.left.time=c()
Type2.left.time=c()

Type1.Difference.v=c()
Type2.Difference.v=c()





#### locate the failure times into renewal intervals for Type 1 and Type 2 

for(i in 1:length(Type1.uncensored.time)){
Type1.left.time[i]=Find_interval(Type1.uncensored.time[i])
Type1.Difference.v[i]=Type1.uncensored.time[i]-Type1.left.time[i]
}

for(i in 1:length(Type2.uncensored.time)){
Type2.left.time[i]=Find_interval(Type2.uncensored.time[i])
Type2.Difference.v[i]=Type2.uncensored.time[i]-Type2.left.time[i]
}



Cyclic.Type.mean.function=function(u,diffv,type,p,covariate){
if(u==0){
total=0
}else{
V=round(sort(unique(diffv[which(diffv<=u)])),6)


YY=matrix(c(1:length(V)),1,length(V))
final=mclapply(YY,function(g){

uncensored.time=renewal_times[,1][-length(renewal_times[,1])]+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==type),"end.time"],6)==round(uncensored.time[h],6)  & Data[which(Data[,"type"]==type),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime)*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))####
ratio=c(denominator,numerator)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,2])
rr=K1/K2
}## function fot g
)##mclapply
total=sum(unlist(final))
}##for else
return(total)
} ## for function


Cyclic.Type.mean.function(u=t2,diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1)
Cyclic.Type.mean.function(u=t2,diffv=Type2.Difference.v,type=2,p=2,covariate=Cov2)




###################Calculate C_j(u;beta)#######
###Cyclic_var(u=2,diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1)

C_type=function(u,diffv,type,p,covariate){ ##0<u<2
if(u==0){
total=0
}else{
V=round(sort(unique(diffv[which(diffv<=u)])),6)
YY=matrix(c(1:length(V)),1,length(V))
final=mclapply(YY,function(g){
uncensored.time=renewal_times[,1][-length(renewal_times[,1])]+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==type),"end.time"],6)==round(uncensored.time[h],6)  & Data[which(Data[,"type"]==type),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime)*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))
denominator2=sum(as.numeric(uncensored.time[h]<=Centime)*covariate*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))
ratio=c(denominator,numerator,denominator2)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),11,3,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),11,3,byrow=TRUE)[,2])
K3=sum(matrix(unlist(Ratio),11,3,byrow=TRUE)[,3])

rr=(K1*K3)/(K2)^2
} ## function fot g
) ## for mclapply
total=sum(unlist(final))
} ## for else
return(matrix(c(total,0),2,1))
} ## for function

#C_type(u=2,diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1)
#C_type(u=2,diffv=Type2.Difference.v,type=2,p=2,covariate=Cov2)

#########################
##Cyclic_var(u=2,diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1)



############## Calculate sum U ####################

Cyclic_Sum_est_U=mclapply(matrix(c(1:n),1,n),function(i){

data_i=Data[which(Data[,"p.id"]==i),]
obs.time=data_i[which(data_i[,"type"]==1),"end.time"]
delta.i=data_i[which(data_i[,"type"]==1),"Censor"]
est_beta=as.numeric(summary(result)$coefficients[,1])
XX_type1_i=matrix(c(1:length(obs.time)),1,length(obs.time))

S0_type1_i=mclapply(XX_type1_i,function(h){
ans=sum(as.numeric(obs.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
}
)
S0_type1_i=unlist(S0_type1_i)


S1_type1_i=mclapply(XX_type1_i,function(h){
a1=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)


U_i_type1_part1=mclapply(matrix(c(1:length(delta.i)),1,length(delta.i)),function(s){
comp1=(Type1_covZ[[i]]-(S1_type1_i[[s]]/S0_type1_i[s]))*delta.i[s]
}## function s
)
U_i_type1_part1=Reduce('+', U_i_type1_part1)

XX_type1=matrix(c(1:length(Type1.uncensored.time)),1,length(Type1.uncensored.time))
S0_type1=mclapply(XX_type1,function(h){
ans=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type1=unlist(S0_type1)
S1_type1=mclapply(XX_type1,function(h){
a1=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(Type1.uncensored.time[h]<=Centime)*unlist(lapply(Type1_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type1_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)

V.type1=round(Type1.Difference.v,6)
YY=matrix(c(1:length(V.type1)),1,length(V.type1))
final.type1=mclapply(YY,function(g){
uncensored.time=renewal_times[,1][-length(renewal_times[,1])]+V.type1[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==1),"end.time"],6)==round(uncensored.time[h],6)  & Data[which(Data[,"type"]==1),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==1),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime)*exp((as.numeric(summary(result)$coefficients[,1][1])*Cov1)))####
ratio=c(denominator,numerator)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,2])
rr=K1/K2
}## function fot g
)##mclapply

U_i_type1_part2=mclapply(matrix(c(1:length(XX_type1)),1,length(XX_type1)),function(s){
comp1=(Type1_covZ[[i]]-(S1_type1[[s]]/S0_type1[s]))*as.numeric(Centime[i]>=Type1.uncensored.time[s])*as.numeric(exp(t(Type1_covZ[[i]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))*unlist(final.type1)[s]
}## function s
)
U_i_type1_part2=Reduce('+', U_i_type1_part2)
data_i=Data[which(Data[,"p.id"]==i),]
obs.time=data_i[which(data_i[,"type"]==2),"end.time"]
delta.i=data_i[which(data_i[,"type"]==2),"Censor"]
est_beta=as.numeric(summary(result)$coefficients[,1])
XX_type2_i=matrix(c(1:length(obs.time)),1,length(obs.time))
S0_type2_i=mclapply(XX_type2_i,function(h){
ans=sum(as.numeric(obs.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
}
)
S0_type2_i=unlist(S0_type2_i)
S1_type2_i=mclapply(XX_type2_i,function(h){
a1=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(obs.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)

U_i_type2_part1=mclapply(matrix(c(1:length(delta.i)),1,length(delta.i)),function(s){
comp1=(Type2_covZ[[i]]-(S1_type2_i[[s]]/S0_type2_i[s]))*delta.i[s]
}## function s
)
U_i_type2_part1=Reduce('+', U_i_type2_part1)

XX_type2=matrix(c(1:length(Type2.uncensored.time)),1,length(Type2.uncensored.time))
S0_type2=mclapply(XX_type2,function(h){
ans=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1))
)))
}
)
S0_type2=unlist(S0_type2)
S1_type2=mclapply(XX_type2,function(h){
a1=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',1))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
a2=sum(as.numeric(Type2.uncensored.time[h]<=Centime)*unlist(lapply(Type2_covZ,'[[',2))*exp((apply(matrix(c(1:n),1,n),2,function(w)t(Type2_covZ[[w]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))))
ans=matrix(c(a1,a2),2,1)
}
)


V.type2=round(((Type2.Difference.v)),6)
YY=matrix(c(1:length(V.type2)),1,length(V.type2))
final.type2=mclapply(YY,function(g){
uncensored.time=renewal_times[,1][-length(renewal_times[,1])]+V.type2[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==2),"end.time"],6)==round(uncensored.time[h],6)  & Data[which(Data[,"type"]==2),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==2),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime)*exp((as.numeric(summary(result)$coefficients[,1][2])*Cov1)))####
ratio=c(denominator,numerator)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,2])
rr=K1/K2
}## function fot g
)##mclapply
U_i_type2_part2=mclapply(matrix(c(1:length(XX_type2)),1,length(XX_type2)),function(s){
comp1=(Type2_covZ[[i]]-(S1_type2[[s]]/S0_type2[s]))*as.numeric(Centime[i]>=Type2.uncensored.time[s])*as.numeric(exp(t(Type2_covZ[[i]])%*%matrix(c(as.numeric(summary(result)$coefficients[,1])),2,1)))*unlist(final.type2)[s]
}## function s
)
U_i_type2_part2=Reduce('+', U_i_type2_part2)

U_sum=(U_i_type1_part1-U_i_type1_part2)+(U_i_type2_part1-U_i_type2_part2)

}## for function i
)##


Cyclic_var=function(u,diffv,type,p,covariate){
Total_cyclic=mclapply(matrix(c(1:n),1,n),function(i){
data_i=Data[which(Data[,"p.id"]==i),]
obs.time=data_i[which(data_i[,"type"]==type),"end.time"]
delta.i=data_i[which(data_i[,"type"]==type),"Censor"]
Type.Difference.v=c()
for(i in 1:length(obs.time)){
Type.Difference.v[i]=obs.time[i]-Find_interval(obs.time[i])
}
est_beta=as.numeric(summary(result)$coefficients[,1])
XX_type_i=matrix(c(1:length(Type.Difference.v)),1,length(Type.Difference.v))
V=round(sort(unique(Type.Difference.v[which(Type.Difference.v<=u)])),6)
YY=matrix(c(1:length(V)),1,length(V))
final_part1=mclapply(YY,function(g){
uncensored.time=renewal_times[,1][-length(renewal_times[,1])]+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(data_i[which(data_i[,"type"]==type),"end.time"],6)==round(uncensored.time[h],6)  & data_i[which(data_i[,"type"]==type),"Censor"]==1,1,0),id=rep(1:1,as.vector(table(data_i[which(data_i[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime)*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))####
ratio=c(denominator,numerator)
}## for function
)## for mclappy
K1=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,1]) ##
K2=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,2])
rr=(K1)/((1/n)*K2)
}## function fot g
)

#### part2

V_total=round(sort(unique(diffv[which(diffv<=u)])),6)

QQ=matrix(c(1:length(V_total)),1,length(V_total))
final_part2=mclapply(QQ,function(g){
uncensored.time=renewal_times[,1][-length(renewal_times[,1])]+V_total[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Data[which(Data[,"type"]==type),"end.time"],6)==round(uncensored.time[h],6)  & Data[which(Data[,"type"]==type),"Censor"]==1,1,0),id=rep(1:length(Centime),as.vector(table(Data[which(Data[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime)*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))####
denominator2=as.numeric(uncensored.time[h]<=Centime[i])*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate[i]))
ratio=c(denominator,numerator,denominator2)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),11,3,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),11,3,byrow=TRUE)[,2])
K3=sum(matrix(unlist(Ratio),11,3,byrow=TRUE)[,3])

rr=(n*K1*K3)/(K2)^2
}## function fot g
)##mclapply

Part_Cyclic_1=sum(unlist(final_part1))-sum(unlist(final_part2))
total_sum=(Part_Cyclic_1+t(C_type(u,diffv,type,p,covariate))%*%solve((1/n)*Hat_A_type)%*%Cyclic_Sum_est_U[[i]])^2
}## for function Total
)

return(mean(unlist(Total_cyclic)))


} ## for Part1



## Standard error of Cai and Schauble's estimator

sqrt((1/n)*Cai_var(t=2,type=1,CovZ=Type1_covZ))  ## for type 1
sqrt((1/n)*Cai_var(t=2,type=2,CovZ=Type2_covZ))  ## for type 2

###Standard error of our proposed cyclic mean function estimator

sqrt((1/n)*(Cyclic_var(u=2,diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1))) ## for type 1
sqrt((1/n)*(Cyclic_var(u=2,diffv=Type2.Difference.v,type=2,p=2,covariate=Cov1))) ## for type 2
















