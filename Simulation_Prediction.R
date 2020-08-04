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

Re.type.data=c()  ##¡@save data

n=100                ## number of individuals
c=22                 ## censor variable uniform (0,c)
true_beta1=-0.3      ## coefficient for beta 1
true_beta2=-0.2      ## coefficient for beta 2
baseline_type1=0.15  ## baseline rate for type 1
baseline_type2=0.2   ## baseline rate for type 2

t2=2
time1_true=t2*baseline_type1
time2_true=t2*baseline_type2

renewal_times=matrix(c(0,2,2,4,4,6,6,8,8,10,10,12,12,14,14,16,16,18,18,20,20,22,22,24),ncol=2,byrow=TRUE)


##################### Generate recurrent event data with two event types ############################################


Data=c()
p.id<-1

while (p.id<=n){

rr.rec1=c()
rr.rec2=c()

X1=rbinom(1,1,0.5)  ## covariate 1

X2=runif(1,0,1)     ## covariate 2
censor.time=runif(1,0,c)
Q=1
##sigma=1         ## variance for random effect for Q
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






############## Create training and validation datasets

Data_training=subset(Data,p.id<=((n/2)))
Data_validation=subset(Data,p.id>=((n/2)+1))

## use coxph to obtain marginal estimates
result<-coxph(Surv(start.time, end.time, Censor)~X1_1+X1_2+cluster(p.id)+strata(type), data=Data_training, robust=T)
as.numeric(summary(result)$coefficients[,1]) ## marginal estimates
as.numeric(summary(result)$coefficients[,4]) ## standard error 





Centime=c()
Cov1=c()

for(i in 1:n){
Centime[i]=Data[,10][which(Data[,1]==i)[1]]
Cov1[i]=Data[,5][which(Data[,9]==1 & Data[,1]==i)[1]]
}

Centime_training=Centime[1:(n/2)]
Cov1_training=Cov1[1:(n/2)]
Centime_validation=Centime[((n/2)+1):n]
Cov1_validation=Cov1[((n/2)+1):n]


### True recurrent number of each subject for Type 1 and Type 2 events
Type1_number=c()
Type2_number=c()
for(i in 1:n){
Type1_number[i]=length(which(Data[,"p.id"]==i & Data[,"Censor"]==1 & Data[,"type"]==1))
Type2_number[i]=length(which(Data[,"p.id"]==i & Data[,"Censor"]==1 & Data[,"type"]==2))
}
################baseline mean function based on Cai and Shaubel ###############################

Cai_Type_baseline=function(t,type,p,covariate){ ##p=1:type 1

if(t==0){
total=0
}else{

V1=unique(sort(c(0,Data_training[which(Data_training[,"type"]==type & Data_training[,"Censor"]==1 ),"end.time"])))
V1=V1[which(V1<t)]

XX=matrix(c(1:length(V1)),1,length(V1))
Ratio=mclapply(XX,function(h){ 
df=data.frame(x=ifelse(Data_training[which(Data_training[,"type"]==type),"end.time"]==V1[h]  & Data_training[which(Data_training[,"type"]==type),"Censor"]==1,1,0),id=rep(1:length(Centime_training),as.vector(table(Data_training[which(Data_training[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(V1[h]<=Centime_training)*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))####
ratio=denominator/numerator
}
)
total=sum(unlist(Ratio))
}## for else
return(max(0,total))
}## for function


## Prediction based on Cai and Schaubel's estimator
Pred_Cai_type1=exp(Cov1_validation*(as.numeric(summary(result)$coefficients[,1][1])))*as.numeric(mclapply(matrix(1:(n/2),1,n/2),function(k){Cai_Type_baseline(Centime_validation[k],type=1,p=1,covariate=Cov1_training)
}))
Pred_Cai_type2=exp(Cov1_validation*(as.numeric(summary(result)$coefficients[,1][2])))*as.numeric(mclapply(matrix(1:(n/2),1,n/2),function(k){Cai_Type_baseline(Centime_validation[k],type=2,p=2,covariate=Cov1_training)
}))

## Average square error for Type 1 and Type 2
mean((Pred_Cai_type1-Type1_number[((n/2)+1):n])^2)
mean((Pred_Cai_type2-Type2_number[((n/2)+1):n])^2)


############## Prediction based on our proposed method ###########################
###################################################################################

Type1.uncensored.time=round(Data_training[,"end.time"][which(Data_training[,"type"]==1 & Data_training[,"Censor"]==1)],6)
Type2.uncensored.time=round(Data_training[,"end.time"][which(Data_training[,"type"]==2 & Data_training[,"Censor"]==1)],6)

Type1.delta=Data_training[,4][which(Data_training[,9]==1)]
Type2.delta=Data_training[,4][which(Data_training[,9]==2)]

Type1.left.time=c()
Type2.left.time=c()

Type1.Difference.v=c()
Type2.Difference.v=c()


for(i in 1:length(Type1.uncensored.time)){
Type1.left.time[i]=Find_interval(Type1.uncensored.time[i])
Type1.Difference.v[i]=Type1.uncensored.time[i]-Type1.left.time[i]
}

for(i in 1:length(Type2.uncensored.time)){
Type2.left.time[i]=Find_interval(Type2.uncensored.time[i])
Type2.Difference.v[i]=Type2.uncensored.time[i]-Type2.left.time[i]
}



Cyclic.Type.mean.function=function(u,diffv,type,p,covariate){ ##0<u<2
if(u==0){
total=0

}else{
V=round(sort(unique(diffv[which(diffv<=u)])),6)
YY=matrix(c(1:length(V)),1,length(V))
final=mclapply(YY,function(g){
uncensored.time=renewal_times[,1][-length(renewal_times[,1])]+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Data_training[which(Data_training[,"type"]==type),"end.time"],6)==round(uncensored.time[h],6)  & Data_training[which(Data_training[,"type"]==type),"Censor"]==1,1,0),id=rep(1:length(Centime_training),as.vector(table(Data_training[which(Data_training[,"type"]==type),"p.id"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime_training)*exp((as.numeric(summary(result)$coefficients[,1][p])*covariate)))
ratio=c(denominator,numerator)
}## for function
)## for mclappy
K1=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),11,2,byrow=TRUE)[,2])
rr=K1/K2
}
)
total=sum(unlist(final))
}##for else
return(total)
} ## for function

Prediction_result=mclapply(matrix(c(1:(n/2)),1,(n/2)),function(i){
if(Find_interval(Centime_validation[i])==0){
Pred_Cyclic_type1=exp((as.numeric(summary(result)$coefficients[,1][1])*Cov1_validation[i]))*Cyclic.Type.mean.function(u=Centime_validation[i],diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1_training)
Pred_Cyclic_type2=exp((as.numeric(summary(result)$coefficients[,1][2])*Cov1_validation[i]))*Cyclic.Type.mean.function(u=Centime_validation[i],diffv=Type2.Difference.v,type=2,p=2,covariate=Cov1_training)
aa=as.vector(c(Pred_Cyclic_type1,Pred_Cyclic_type2))
}else{
Diff.v=Centime_validation[i]-Find_interval(Centime_validation[i])
muplier=Find_interval(Centime_validation[i])/2
Pred_Cyclic_type1=exp((as.numeric(summary(result)$coefficients[,1][1])*Cov1_validation[i]))*((muplier*Cyclic.Type.mean.function(u=2,diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1_training))+Cyclic.Type.mean.function(u=Diff.v,diffv=Type1.Difference.v,type=1,p=1,covariate=Cov1_training))
Pred_Cyclic_type2=exp((as.numeric(summary(result)$coefficients[,1][2])*Cov1_validation[i]))*((muplier*Cyclic.Type.mean.function(u=2,diffv=Type2.Difference.v,type=2,p=2,covariate=Cov1_training))+Cyclic.Type.mean.function(u=Diff.v,diffv=Type2.Difference.v,type=2,p=2,covariate=Cov1_training))
aa=as.vector(c(Pred_Cyclic_type1,Pred_Cyclic_type2))
}})


## Prediction based on our proposed estimator
Pred_Cyclic_type1=unlist(lapply(Prediction_result,'[[',1))
Pred_Cyclic_type2=unlist(lapply(Prediction_result,'[[',2))

## Average square error for Type 1 and Type 2
mean((Pred_Cyclic_type1-Type1_number[((n/2)+1):n])^2)
mean((Pred_Cyclic_type2-Type2_number[((n/2)+1):n])^2)


 