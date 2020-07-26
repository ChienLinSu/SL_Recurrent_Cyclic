############################# Load R package ##########################################

library(survival)
library(ggplot2)
library(parallel)


############## The real data set can be downloaded from                      ############################
############## https://projecteuclid.org/euclid.aoas/1346418578\#supplemental############################
############## (Wang et al. (2012, Annals of Applied Statistics) )           ############################

## please read the dataset completedata.txt into R in your computer
data <- read.table("F:/Postdoc at McGill/Lin Feng Chang/completedata.txt", header=T)

DD=as.data.frame(table(cut(data[,2], breaks=quantile(data[,2], probs =c(0,0.2,0.4,0.6,0.8,1))), cut(data[,3], breaks=quantile(data[,3], probs =c(0,0.3,0.5,0.7,0.9,1)))))

x_breaks=quantile(data[,2], probs =c(0,0.2,0.4,0.6,0.8,1))
y_breaks=quantile(data[,3], probs =c(0,0.15,0.3,0.5,0.6,0.8,1))

points <- cbind(data[,2], data[,3])
Cordo=data.frame(t1=findInterval(points[,1], x_breaks),t2=findInterval(points[,2], y_breaks))

k=0
New_data=c()
for(i in 1:6){ ## 6 intervals for y
for(j in 1:5){ ## 5 intervals for x
k=k+1
location=which(Cordo[,1]==j & Cordo[,2]==i) 
New_data=rbind(New_data,cbind(rep(k,length(location)),data[location,]))

}## for j
}## for i

colnames(New_data)[1] <- "Area"
delete1=which(ifelse(is.na(New_data[,8]),1,0)==1)
New_data=New_data[-delete1,]
delete2=which(ifelse(is.na(New_data[,9]),1,0)==1)
New_data=New_data[-delete2,]
table(New_data[,1])



################################## Define variables ######################################
##########################################################################################


################### Fire type ##################


Fire_type=c()
for(k in 1:length(New_data[,7])){
if(New_data[,7][k]==1){
Fire_type[k]=1  ##1:structure fire
}else if(New_data[,7][k]==3){
Fire_type[k]=2 ## 2:vegetation fire
}else{
Fire_type[k]=3
}
}

##################### Alarm ####################


Alarm=c()

for(k in 1:length(New_data[,6])){
if(New_data[,6][k]==1){
Alarm[k]=1
}else{
Alarm[k]=0
}
}

################### Heatsource #################

Heatsource=c()
for(k in 1:length(New_data[,8])){
if(New_data[,8][k]==5){
Heatsource[k]=1
}else if(New_data[,8][k]==7){
Heatsource[k]=2
}else{
Heatsource[k]=3
}
}



############  Objignited ####################

Objignited=c()
for(k in 1:length(New_data[,9])){
if(New_data[,9][k]==7){  ##7:outdoor item
Objignited[k]=1
}else {
Objignited[k]=0
}
}

############ Days ##########################

Dayweek=New_data[,12]
Urban=New_data[,5]

Label=New_data[,13]

Heatsource_Hot=ifelse(Heatsource==1,1,0)
Heatsource_Cigarette=ifelse(Heatsource==2,1,0)
Monday=ifelse(Dayweek==1,1,0)
Tuesday=ifelse(Dayweek==2,1,0)
Wednesday=ifelse(Dayweek==3,1,0)
Thursday=ifelse(Dayweek==4,1,0)
Friday=ifelse(Dayweek==5,1,0)
Saturday=ifelse(Dayweek==6,1,0)
Sunday=ifelse(Dayweek==7,1,0)

Arrival_time=c()

for(k in 1:length(New_data[,2])){
if(New_data[k,2]>=2004 & New_data[k,2]<2005){
Arrival_time[k]=New_data[k,10]
}else if(New_data[k,2]>=2005 & New_data[k,2]<2006){
Arrival_time[k]=New_data[k,10]+366
}else if(New_data[k,2]>=2006 & New_data[k,2]<2007){
Arrival_time[k]=New_data[k,10]+(365+366)
}else{
Arrival_time[k]=New_data[k,10]+(366+365+365)
}


}## for k

Final_Data=data.frame(Area=New_data[,1],year=New_data[,2],Arrival_time,Fire_type,Urban,Alarm,Heatsource_Hot,Heatsource_Cigarette,Objignited,Tuesday,Wednesday,Thursday,Friday,Saturday,Sunday,Label,X_axis=New_data[,3],Y_axis=New_data[,4])

location=which(Final_Data[,4]==3) ### delete type 3

Final_Data=Final_Data[-location,]

final_data=c()
for(i in 1:30){

rr.rec1=c()
T1=0
censor.time=1461

location1=which(Final_Data[,1]==i & Final_Data[,4]==1) ## 1:event type 1

rec.type1=c(T1,Final_Data[location1,3],censor.time)
Type1_data=Final_Data[location1,]


for(j in 1:(length(rec.type1)-1)){
if (j==(length(rec.type1)-1)) {
rr.rec1=rbind(rr.rec1,c(i,rec.type1[j:(j+1)],0,1,Type1_data[j-1,5],0,Type1_data[j-1,6],0,Type1_data[j-1,7],0,Type1_data[j-1,8],0,Type1_data[j-1,9],0,Type1_data[j-1,10],0,Type1_data[j-1,11],0,Type1_data[j-1,12],0,Type1_data[j-1,13],0,Type1_data[j-1,14],0,Type1_data[j-1,15],0,Type1_data[j-1,16],0))

}else{

rr.rec1=rbind(rr.rec1,c(i,rec.type1[j:(j+1)],1,1,Type1_data[j,5],0,Type1_data[j,6],0,Type1_data[j,7],0,Type1_data[j,8],0,Type1_data[j,9],0,Type1_data[j,10],0,Type1_data[j,11],0,Type1_data[j,12],0,Type1_data[j,13],0,Type1_data[j,14],0,Type1_data[j,15],0,Type1_data[j,16],0))
}

}## for j

if(length(rr.rec1)==17){
rr.rec1=cbind(rr.rec1,matrix(rep(0,(29-length(rr.rec1))),1,(29-length(rr.rec1))))
}

rr.rec2=c()
T2=0
censor.time=1461


location2=which(Final_Data[,1]==i & Final_Data[,4]==2)

rec.type2=c(T2,Final_Data[location2,3],censor.time)
Type2_data=Final_Data[location2,]


for(j in 1:(length(rec.type2)-1)){
if (j==(length(rec.type2)-1)) {
rr.rec2=rbind(rr.rec2,c(i,rec.type2[j:(j+1)],0,2,0,Type2_data[j-1,5],0,Type2_data[j-1,6],0,Type2_data[j-1,7],0,Type2_data[j-1,8],0,Type2_data[j-1,9],0,Type2_data[j-1,10],0,Type2_data[j-1,11],0,Type2_data[j-1,12],0,Type2_data[j-1,13],0,Type2_data[j-1,14],0,Type2_data[j-1,15],0,Type2_data[j-1,16]))

}else{

rr.rec2=rbind(rr.rec2,c(i,rec.type2[j:(j+1)],1,2,0,Type2_data[j,5],0,Type2_data[j,6],0,Type2_data[j,7],0,Type2_data[j,8],0,Type2_data[j,9],0,Type2_data[j,10],0,Type2_data[j,11],0,Type2_data[j,12],0,Type2_data[j,13],0,Type2_data[j,14],0,Type2_data[j,15],0,Type2_data[j,16]))
}

}## for j


if(length(rr.rec2)==17){
rr.rec2=cbind(rr.rec2,matrix(rep(0,(29-length(rr.rec2))),1,(29-length(rr.rec2))))

}

total.rr.rec<-data.frame(rbind(rr.rec1,rr.rec2))
names(total.rr.rec)=c("Area", "start.time", "end.time", "Censor","Type","Urban_1","Urban_2","Alarm_1","Alarm_2","Heat_Hot_1","Heat_Hot_2","Heat_Cig_1","Heat_Cig_2","Obji_1","Obji_2","Tues_1","Tues_2","Wed_1","Wed_2","Thur_1","Thur_2","Fri_1","Fri_2","Sat_1","Sat_2","Sun_1","Sun_2","Label_1","Label_2")
 
final_data=rbind(final_data,total.rr.rec)
} ## for i 

final_data=final_data[-which((final_data[,3]- final_data[,2])==0),]
final_data$weekend_1=ifelse(final_data$Sat_1==1 |final_data$Sun_1==1,1,0) 
final_data$weekend_2=ifelse(final_data$Sat_2==1 |final_data$Sun_2==1,1,0) 
final_data_sub <- subset(final_data, Censor>0)



###################################### Results for Table 4 in the paper ######################################################################
###################################################################################################################################3


## summary statistics for urban

##count
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Urban_1"]==1))
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Urban_1"]==0))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Urban_2"]==1))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Urban_2"]==0))

## occurrence rate per year
final_data_sub_ty1=subset(final_data_sub,Type==1)
365*111/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Urban_1==1])
365*19/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Urban_1==0])

final_data_sub_ty2=subset(final_data_sub,Type==2)
365*83/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Urban_2==1])
365*76/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Urban_2==0])

## Poisson regression
df1=data.frame(count=c(111,19),type=c(1,0),Length=c(24454,6181))
df2=data.frame(count=c(83,76),type=c(1,0),Length=c(20675,12960))

summary(glm(count ~ type, offset = log(Length), family=poisson, data=df1))
summary(glm(count ~ type, offset = log(Length), family=poisson, data=df2))



######summary statistics for alarm

## count
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Alarm_1"]==1))
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Alarm_1"]==0))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Alarm_2"]==1))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Alarm_2"]==0))

## occurrence rate per year
365*114/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Alarm_1==1])
365*16/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Alarm_1==0])
365*146/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Alarm_2==1])
365*13/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Alarm_2==0])

## Poisson regression
df1=data.frame(count=c(114,16),type=c(1,0),Length=c(27593,3042))
df2=data.frame(count=c(146,13),type=c(1,0),Length=c(31438,2197))
summary(glm(count ~ type, offset = log(Length), family=poisson, data=df1))
summary(glm(count ~ type, offset = log(Length), family=poisson, data=df2))



###### summary statistics for heatsources

## count
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Heat_Hot_1"]==1))
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Heat_Cig_1"]==1))
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Heat_Hot_1"]==0 & final_data_sub[,"Heat_Cig_1"]==0))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Heat_Hot_2"]==1))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Heat_Cig_2"]==1))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Heat_Hot_2"]==0 & final_data_sub[,"Heat_Cig_2"]==0))




## occurrence rate per year
365*53/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Heat_Hot_1==1])
365*30/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Heat_Cig_1==1])
365*47/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Heat_Hot_1==0 &¡@final_data_sub_ty1$Heat_Cig_1==0])
365*16/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Heat_Hot_2==1])
365*88/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Heat_Cig_2==1])
365*55/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Heat_Hot_2==0 & final_data_sub_ty2$Heat_Cig_2==0])

## Poisson regression
df1=data.frame(count=c(53,30,47),type=c(1,2,0),Length=c(13911,5343,11381))
df2=data.frame(count=c(16,88,55),type=c(1,2,0),Length=c(1931,22262,9442))

summary(glm(count ~ as.factor(type), offset = log(Length), family=poisson, data=df1))
summary(glm(count ~ as.factor(type), offset = log(Length), family=poisson, data=df2))





### summary statistics for Object

## count
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Obji_1"]==1))
length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"Obji_1"]==0))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Obji_2"]==1))
length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"Obji_2"]==0))

## occurrence rate per year
365*41/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Obji_1==1])
365*89/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$Obji_1==0])
365*151/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Obji_2==1])
365*8/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$Obji_2==0])

## Poisson regression
df1=data.frame(count=c(41,89),type=c(1,0),Length=c(7904,22731))
df2=data.frame(count=c(151,8),type=c(1,0),Length=c(31624,2011))
summary(glm(count ~ type, offset = log(Length), family=poisson, data=df1))
summary(glm(count ~ type, offset = log(Length), family=poisson, data=df2))



### summary statistics for Day/Weeks

## count
#length(which( final_data_sub[,"Type"]==1 & final_data_sub[,"weekend_1"]==1))
#length(which(final_data_sub[,"Type"]==1 & final_data_sub[,"weekend_1"]==0))
#length(which(final_data_sub[,"Type"]==2 & final_data_sub[,"weekend_2"]==1))
#length(which( final_data_sub[,"Type"]==2 & final_data_sub[,"weekend_2"]==0))

## occurrence rate per year
#365*44/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$weekend_1==1])
#365*86/sum((final_data_sub_ty1$end.time-final_data_sub_ty1$start.time)[final_data_sub_ty1$weekend_1==0])
#365*61/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$weekend_2==1])
#365*98/sum((final_data_sub_ty2$end.time-final_data_sub_ty2$start.time)[final_data_sub_ty2$weekend_2==0])

## Poisson regression model
#df1=data.frame(count=c(44,86),type=c(1,0),Length=c(7020,23587))
#df2=data.frame(count=c(61,98),type=c(1,0),Length=c(13063,20544))
#summary(glm(count ~ type, offset = log(Length), family=poisson, data=df1))
#summary(glm(count ~ type, offset = log(Length), family=poisson, data=df2))


######################################### Fitted marginal results for Table 5 in the paper ###############################
##########################################################################################################################

fit_result<-coxph(Surv(start.time, end.time, Censor)~Urban_1+Urban_2+Alarm_1+Alarm_2+Heat_Hot_1+Heat_Hot_2+Heat_Cig_1+Heat_Cig_2+Obji_1+Obji_2+cluster(Area)+strata(Type), data=final_data_sub, robust=T)
fit_result





############################################## Figure 2 in the paper ##################################################
#######################################################################################################################

### plot distribution of structure fire

event.time.structure=final_data_sub[,"end.time"][which(final_data_sub[,"Type"]==1 & final_data_sub[,"Censor"]==1) ]


Year=c()
New.time=c()
for(k in 1:length(event.time.structure)){
if(event.time.structure[k]<=366){
Year[k]=2004
New.time[k]=event.time.structure[k]
}else if(366<event.time.structure[k] & event.time.structure[k]<=731){
Year[k]=2005
New.time[k]=event.time.structure[k]-366
}else if(731<event.time.structure[k] & event.time.structure[k]<=1096){
Year[k]=2006
New.time[k]=event.time.structure[k]-731
}else{
Year[k]=2007
New.time[k]=event.time.structure[k]-1096
}
}

df_stru=data.frame(Year,event.time.structure,New.time)

new_p1=ggplot(df_stru, aes(x=New.time)) + 
 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=50,
                   colour="black", fill="white") +
  geom_density(alpha=0.5, fill="#FF6666")+ 
 facet_grid(. ~ Year)+
  labs(x="Month",y="Density") +
ggtitle("Distribution of Structure Fire")+
 theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = c(31, 90, 151,212,273,334), labels = c('1', '3', '5', '7', '9', '11'))+
ylim(0,0.009)





###### plot distribution of vegetation fire


event.time.vegetatiion=final_data_sub[,"end.time"][which(final_data_sub[,"Type"]==2 & final_data_sub[,"Censor"]==1) ]

Year=c()
New.time=c()
for(k in 1:length(event.time.vegetatiion)){
if(event.time.vegetatiion[k]<=366){
Year[k]=2004
New.time[k]=event.time.vegetatiion[k]
}else if(366<event.time.vegetatiion[k] & event.time.vegetatiion[k]<=731){
Year[k]=2005
New.time[k]=event.time.vegetatiion[k]-366
}else if(731<event.time.vegetatiion[k] & event.time.vegetatiion[k]<=1096){
Year[k]=2006
New.time[k]=event.time.vegetatiion[k]-731
}else{
Year[k]=2007
New.time[k]=event.time.vegetatiion[k]-1096
}
}

df_vege=data.frame(Year,event.time.vegetatiion,New.time)

new_p2=ggplot(df_vege, aes(x=New.time)) + 
 geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=50,
                   colour="black", fill="white") +
  geom_density(alpha=0.5, fill="#FF6666")+ 
 facet_grid(. ~ Year)+
  labs(x="Month",y="Density") +
ggtitle("Distribution of Vegetation Fire")+
 theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = c(31, 90, 151,212,273,334), labels = c('1', '3', '5', '7', '9', '11'))+
ylim(0,0.009)



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




multiplot(new_p1, new_p2, cols=2)



############################################# Plot Figure 3 in the paper ###############################################
########################################################################################################################

Structure_data=subset(final_data_sub,Type==1)
Vegetation_data=subset(final_data_sub,Type==2)

Structure_Centime=c()
Vegetation_Centime=c()

for(i in 1:length(table(Structure_data$Area))){
location=which((Structure_data[,1])==as.numeric(names(table(Structure_data[,1])))[i])
Structure_Centime=c(Structure_Centime,max(as.numeric(as.vector(Structure_data[location,3]))))
}


for(i in 1:length(table(Vegetation_data$Area))){
location=which((Vegetation_data[,1])==as.numeric(names(table(Vegetation_data[,1])))[i])
Vegetation_Centime=c(Vegetation_Centime,max(as.numeric(as.vector(Vegetation_data[location,3]))))
}



Urb_structure=c()
Alarm_structure=c()
Heat_Hot_structure=c()
Heat_Cig_structure=c()
Obji_structure=c()
Weekend_structure=c()


for(i in 1:length(table(Structure_data$Area))){
location=which(Structure_data[,1]==as.numeric(names(table(Structure_data$Area)))[i])
Urb_structure=c(Urb_structure,Structure_data[location,"Urban_1"][1])
Alarm_structure=c(Alarm_structure,Structure_data[location,"Alarm_1"][1])
Heat_Hot_structure=c(Heat_Hot_structure,Structure_data[location,"Heat_Hot_1"][1])
Heat_Cig_structure=c(Heat_Cig_structure,Structure_data[location,"Heat_Cig_1"][1])
Obji_structure=c(Obji_structure,Structure_data[location,"Obji_1"][1])
Weekend_structure=c(Weekend_structure,Structure_data[location,"weekend_1"][1])
}


Urb_vegetation=c()
Alarm_vegetation=c()
Heat_Hot_vegetation=c()
Heat_Cig_vegetation=c()
Obji_vegetation=c()
Weekend_vegetation=c()


for(i in 1:length(table(Vegetation_data$Area))){
location=which(Vegetation_data[,1]==as.numeric(names(table(Vegetation_data$Area)))[i])
Urb_vegetation=c(Urb_vegetation,Vegetation_data[location,"Urban_2"][1])
Alarm_vegetation=c(Alarm_vegetation,Vegetation_data[location,"Alarm_2"][1])
Heat_Hot_vegetation=c(Heat_Hot_vegetation,Vegetation_data[location,"Heat_Hot_2"][1])
Heat_Cig_vegetation=c(Heat_Cig_vegetation,Vegetation_data[location,"Heat_Cig_2"][1])
Obji_vegetation=c(Obji_vegetation,Vegetation_data[location,"Obji_2"][1])
Weekend_vegetation=c(Weekend_vegetation,Vegetation_data[location,"weekend_2"][1])
}

Str_restart.time=c()
Str_reend.time=c()

for(i in 1:length(Structure_data[,2])){
if(Structure_data[i,2]<=366){
Str_restart.time[i]=Structure_data[i,2]
}else if(366<Structure_data[i,2] & Structure_data[i,2]<=731){
Str_restart.time[i]=Structure_data[i,2]-366
}else if(731<Structure_data[i,2] & Structure_data[i,2]<=1096){
Str_restart.time[i]=Structure_data[i,2]-731
}else if(1096<Structure_data[i,2] & Structure_data[i,2]<=1461){
Str_restart.time[i]=Structure_data[i,2]-1096
}
}# for loop


for(i in 1:length(Structure_data[,3])){
if(Structure_data[i,3]<=366){
Str_reend.time[i]=Structure_data[i,3]
}else if(366<Structure_data[i,3] & Structure_data[i,3]<=731){
Str_reend.time[i]=Structure_data[i,3]-366
}else if(731<Structure_data[i,3] & Structure_data[i,3]<=1096){
Str_reend.time[i]=Structure_data[i,3]-731
}else if(1096<Structure_data[i,3] & Structure_data[i,3]<=1461){
Str_reend.time[i]=Structure_data[i,3]-1096
}
}# for loop


Veg_restart.time=c()
Veg_reend.time=c()

for(i in 1:length(Vegetation_data[,2])){
if(Vegetation_data[i,2]<=366){
Veg_restart.time[i]=Vegetation_data[i,2]
}else if(366<Vegetation_data[i,2] & Vegetation_data[i,2]<=731){
Veg_restart.time[i]=Vegetation_data[i,2]-366
}else if(731<Vegetation_data[i,2] & Vegetation_data[i,2]<=1096){
Veg_restart.time[i]=Vegetation_data[i,2]-731
}else if(1096<Vegetation_data[i,2] & Vegetation_data[i,2]<=1461){
Veg_restart.time[i]=Vegetation_data[i,2]-1096
}
}# for loop


for(i in 1:length(Vegetation_data[,3])){
if(Vegetation_data[i,3]<=366){
Veg_reend.time[i]=Vegetation_data[i,3]
}else if(366<Vegetation_data[i,3] & Vegetation_data[i,3]<=731){
Veg_reend.time[i]=Vegetation_data[i,3]-366
}else if(731<Vegetation_data[i,3] & Vegetation_data[i,3]<=1096){
Veg_reend.time[i]=Vegetation_data[i,3]-731
}else if(1096<Vegetation_data[i,3] & Vegetation_data[i,3]<=1461){
Veg_reend.time[i]=Vegetation_data[i,3]-1096
}
}# for loop



renewal_times=c(0,366,731,1096)
Centime_structure=Structure_Centime


## construct the cyclic mean function for structure fire

Cyclic.Structure.mean.function=function(u){ 
if(u==0){
total=0

}else{

V=round(sort(unique(Str_reend.time[which(Str_reend.time<=u)])),6)
YY=matrix(c(1:length(V)),1,length(V))
final=mclapply(YY,function(g){

uncensored.time=renewal_times+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Structure_data[which(Structure_data[,"Type"]==1),"end.time"],6)==round(uncensored.time[h],6)  & Structure_data[which(Structure_data[,"Type"]==1),"Censor"]==1,1,0),id=rep(1:length(table(Structure_data[,"Area"])),as.vector(table(Structure_data[which(Structure_data[,"Type"]==1),"Area"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime_structure)*exp((fit_result$coef[1]*Urb_structure)+(fit_result$coef[3]*Alarm_structure)+(fit_result$coef[5]*Heat_Hot_structure)+(fit_result$coef[7]*Heat_Cig_structure)+(fit_result$coef[9]*Obji_structure)
))
ratio=c(denominator,numerator)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,2])
rr=K1/K2
}## function fot g
)##mclapply
total=sum(unlist(final))
}##for else
return(total)
} ## for function



## construct the cyclic mean function for vegetation fire

Centime_vegetation=Vegetation_Centime
Cyclic.Vegetation.mean.function=function(u){ 
if(u==0){
total=0

}else{

V=round(sort(unique(Veg_reend.time[which(Veg_reend.time<=u)])),6)
YY=matrix(c(1:length(V)),1,length(V))
final=mclapply(YY,function(g){

uncensored.time=renewal_times+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Vegetation_data[which(Vegetation_data[,"Type"]==2),"end.time"],6)==round(uncensored.time[h],6)  & Vegetation_data[which(Vegetation_data[,"Type"]==2),"Censor"]==1,1,0),id=rep(1:length(table(Vegetation_data[,"Area"])),as.vector(table(Vegetation_data[which(Vegetation_data[,"Type"]==2),"Area"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime_vegetation)*exp((fit_result$coef[2]*Urb_vegetation)+(fit_result$coef[4]*Alarm_vegetation)+(fit_result$coef[6]*Heat_Hot_vegetation)+(fit_result$coef[8]*Heat_Cig_vegetation)+(fit_result$coef[10]*Obji_vegetation)
))
ratio=c(denominator,numerator)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,2])
rr=K1/K2
}##  
)## mclapply
total=sum(unlist(final))
}##for else
return(total)
} ## for function



##¡@calculate the values of cyclic mean functions for structure and vegetation fires
xx=seq(0,366,by=1)
type1_stru_cyclic=c()
type2_vege_cyclic=c()

for(k in 1:length(xx)){
type1_stru_cyclic[k]=Cyclic.Structure.mean.function(xx[k])
type2_vege_cyclic[k]=Cyclic.Vegetation.mean.function(xx[k])
cat("number=",k,"\n")
}



## Construct the mean function of structure fire based on the estimator proposed by Cai and Schaubel(2004)
Structure_mean_function=function(t){

if(t==0){
total=0

}else{

V1=unique(sort(c(0,Structure_data[,3][which(Structure_data[,3]<=t )])))
ttr=c()
total=0
for(h in 1:length(V1)){
df=data.frame(x=ifelse(Structure_data[,3]==V1[h] & Structure_data[,4]==1,1,0),id=rep(1:length(Structure_Centime),as.vector(table(Structure_data[,1]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(V1[h]<=Structure_Centime)*exp((fit_result$coef[1]*Urb_structure)+(fit_result$coef[3]*Alarm_structure)+(fit_result$coef[5]*Heat_Hot_structure)+(fit_result$coef[7]*Heat_Cig_structure)+(fit_result$coef[9]*Obji_structure)
)) 
total=total+(denominator/numerator)

} ## for h
}

return(max(0,total))
}## for function


## Construct the mean function of vegetation fire based on the estimator proposed by Cai and Schaubel(2004)
Vegetation_mean_function=function(t){

if(t==0){
total=0

}else{

V1=unique(sort(c(0,Vegetation_data[,3][which(Vegetation_data[,3]<=t )])))
ttr=c()
total=0
for(h in 1:length(V1)){
df=data.frame(x=ifelse(Vegetation_data[,3]==V1[h] & Vegetation_data[,4]==1,1,0),id=rep(1:length(Vegetation_Centime),as.vector(table(Vegetation_data[,1]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(V1[h]<=Vegetation_Centime)*exp((fit_result$coef[2]*Urb_vegetation)+(fit_result$coef[4]*Alarm_vegetation)+(fit_result$coef[6]*Heat_Hot_vegetation)+(fit_result$coef[8]*Heat_Cig_vegetation)+(fit_result$coef[10]*Obji_vegetation)
)) 
total=total+(denominator/numerator)
} ## for h
}
return(max(0,total))
}## for function


xx=seq(0,1461,by=1)
type1_stru=c()
type2_vege=c()

for(k in 1:length(xx)){
type1_stru[k]=Structure_mean_function(xx[k])
type2_vege[k]=Vegetation_mean_function(xx[k])
cat("number=",k,"\n")
}



####Plot mean functions 

xx_1=seq(0,366,by=1)
xx_2=seq(366,731,by=1)
xx_3=seq(731,1096,by=1)
xx_4=seq(1096,1461,by=1)

Cai_outcome_structure_1=data.frame(time=xx_1,type=type1_stru[1:367])
Cai_outcome_structure_2=data.frame(time=xx_2,type=type1_stru[367:732]-type1_stru[367])
Cai_outcome_structure_3=data.frame(time=xx_3,type=type1_stru[732:1097]-type1_stru[732])
Cai_outcome_structure_4=data.frame(time=xx_4,type=type1_stru[1097:1462]-type1_stru[1097])

Cai_outcome_vegetation_1=data.frame(time=xx_1,type=type2_vege[1:367])
Cai_outcome_vegetation_2=data.frame(time=xx_2,type=type2_vege[367:732]-type2_vege[367])
Cai_outcome_vegetation_3=data.frame(time=xx_3,type=type2_vege[732:1097]-type2_vege[732])
Cai_outcome_vegetation_4=data.frame(time=xx_4,type=type2_vege[1097:1462]-type2_vege[1097])


outcome_structure_1=data.frame(time=xx_1,type=type1_stru_cyclic)
outcome_vegetation_1=data.frame(time=xx_1,type=type2_vege_cyclic)

outcome_structure_2=data.frame(time=xx_2,type=type1_stru_cyclic[1:366])
outcome_vegetation_2=data.frame(time=xx_2,type=type2_vege_cyclic[1:366])

outcome_structure_3=data.frame(time=xx_3,type=type1_stru_cyclic[1:366])
outcome_vegetation_3=data.frame(time=xx_3,type=type2_vege_cyclic[1:366])

outcome_structure_4=data.frame(time=xx_4,type=type1_stru_cyclic[1:366])
outcome_vegetation_4=data.frame(time=xx_4,type=type2_vege_cyclic[1:366])

visuals=rbind(Cai_outcome_structure_1,Cai_outcome_vegetation_1,Cai_outcome_structure_2,Cai_outcome_vegetation_2,Cai_outcome_structure_3,
Cai_outcome_vegetation_3,Cai_outcome_structure_4,Cai_outcome_vegetation_4,outcome_structure_1,outcome_vegetation_1,outcome_structure_2,
outcome_vegetation_2,outcome_structure_3,outcome_vegetation_3,outcome_structure_4,outcome_vegetation_4)

visuals$Fire_type=c(rep("Structure",length(xx_1)),rep("Vegetation",length(xx_1)),rep("Structure",length(xx_2)),rep("Vegetation",length(xx_2))
,rep("Structure",length(xx_3)),rep("Vegetation",length(xx_3)),rep("Structure",length(xx_4)),rep("Vegetation",length(xx_4)),rep("Cyclic_Structure",length(xx_1)),
rep("Cyclic_Vegetation",length(xx_1)),rep("Cyclic_Structure",length(xx_2)),rep("Cyclic_Vegetation",length(xx_2)),rep("Cyclic_Structure",length(xx_3)),
rep("Cyclic_Vegetation",length(xx_3)),rep("Cyclic_Structure",length(xx_4)),rep("Cyclic_Vegetation",length(xx_4))) 

visuals$Fire_type=c(rep("Non-cyclic structure",length(xx_1)),rep("Non-cyclic vegetation",length(xx_1)),rep("Non-cyclic structure",length(xx_2)),rep("Non-cyclic vegetation",length(xx_2))
,rep("Non-cyclic structure",length(xx_3)),rep("Non-cyclic vegetation",length(xx_3)),rep("Non-cyclic structure",length(xx_4)),rep("Non-cyclic vegetation",length(xx_4)),rep("Cyclic structure",length(xx_1)),
rep("Cyclic vegetation",length(xx_1)),rep("Cyclic structure",length(xx_2)),rep("Cyclic vegetation",length(xx_2)),rep("Cyclic structure",length(xx_3)),
rep("Cyclic vegetation",length(xx_3)),rep("Cyclic structure",length(xx_4)),rep("Cyclic vegetation",length(xx_4))) 

LegendTitle = "Fire type"
 ggplot(visuals, aes(time,type,group=Fire_type,col=Fire_type)) + 
 geom_line(aes(linetype=Fire_type))+

scale_linetype_manual(name = LegendTitle,values=c(1,1,4,4)) +
scale_color_manual(name = LegendTitle,values=c("blue","red","blue","red" ))+
theme_classic()+
geom_vline(aes(xintercept=0), colour="gray69", linetype="dashed")+
geom_vline(aes(xintercept=366), colour="gray69", linetype="dashed")+
geom_vline(aes(xintercept=731), colour="gray69", linetype="dashed")+
geom_vline(aes(xintercept=1096), colour="gray69", linetype="dashed")+
geom_vline(aes(xintercept=1462), colour="gray69", linetype="dashed")+
scale_x_continuous(name="Time (year)",breaks=c(0,366,731,1096,1462),
        labels=c(2004, 2005, 2006,2007,2008))+
scale_y_continuous(name="Mean function")




#####################  Figure 4: Prediction figure #############################################################
################################################################################################################

fit_result_pred<-coxph(Surv(start.time, end.time, Censor)~Urban_1+Urban_2+cluster(Area)+strata(Type), data=final_data_sub, robust=T)

Str_restart.time=c()
Str_reend.time=c()

for(i in 1:length(Structure_data[,2])){
if(Structure_data[i,2]<=366){
Str_restart.time[i]=Structure_data[i,2]

}else if(366<Structure_data[i,2] & Structure_data[i,2]<=731){
Str_restart.time[i]=Structure_data[i,2]-366

}else if(731<Structure_data[i,2] & Structure_data[i,2]<=1096){
Str_restart.time[i]=Structure_data[i,2]-731

}else if(1096<Structure_data[i,2] & Structure_data[i,2]<=1461){
Str_restart.time[i]=Structure_data[i,2]-1096

}

}# for loop


for(i in 1:length(Structure_data[,3])){
if(Structure_data[i,3]<=366){

Str_reend.time[i]=Structure_data[i,3]


}else if(366<Structure_data[i,3] & Structure_data[i,3]<=731){

Str_reend.time[i]=Structure_data[i,3]-366

}else if(731<Structure_data[i,3] & Structure_data[i,3]<=1096){

Str_reend.time[i]=Structure_data[i,3]-731

}else if(1096<Structure_data[i,3] & Structure_data[i,3]<=1461){

Str_reend.time[i]=Structure_data[i,3]-1096
}

}# for loop


Veg_restart.time=c()
Veg_reend.time=c()


for(i in 1:length(Vegetation_data[,2])){
if(Vegetation_data[i,2]<=366){
Veg_restart.time[i]=Vegetation_data[i,2]

}else if(366<Vegetation_data[i,2] & Vegetation_data[i,2]<=731){
Veg_restart.time[i]=Vegetation_data[i,2]-366

}else if(731<Vegetation_data[i,2] & Vegetation_data[i,2]<=1096){
Veg_restart.time[i]=Vegetation_data[i,2]-731

}else if(1096<Vegetation_data[i,2] & Vegetation_data[i,2]<=1461){
Veg_restart.time[i]=Vegetation_data[i,2]-1096

}

}# for loop


for(i in 1:length(Vegetation_data[,3])){
if(Vegetation_data[i,3]<=366){

Veg_reend.time[i]=Vegetation_data[i,3]

}else if(366<Vegetation_data[i,3] & Vegetation_data[i,3]<=731){

Veg_reend.time[i]=Vegetation_data[i,3]-366

}else if(731<Vegetation_data[i,3] & Vegetation_data[i,3]<=1096){

Veg_reend.time[i]=Vegetation_data[i,3]-731

}else if(1096<Vegetation_data[i,3] & Vegetation_data[i,3]<=1461){

Veg_reend.time[i]=Vegetation_data[i,3]-1096
}

}# for loop



Cyclic.Structure.mean.function=function(u){ 
if(u==0){
total=0

}else{

V=round(sort(unique(Str_reend.time[which(Str_reend.time<=u)])),6)
YY=matrix(c(1:length(V)),1,length(V))
final=mclapply(YY,function(g){

uncensored.time=renewal_times+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Structure_data[which(Structure_data[,"Type"]==1),"end.time"],6)==round(uncensored.time[h],6)  & Structure_data[which(Structure_data[,"Type"]==1),"Censor"]==1,1,0),id=rep(1:length(table(Structure_data[,"Area"])),as.vector(table(Structure_data[which(Structure_data[,"Type"]==1),"Area"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime_structure)*exp((fit_result_pred$coef[1]*Urb_structure)))####
ratio=c(denominator,numerator)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,1]) 
K2=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,2])
rr=K1/K2
}## function fot g
)##mclapply
total=sum(unlist(final))
}##for else
return(total)
} ## for function



Cyclic.Vegetation.mean.function=function(u){ 
if(u==0){
total=0

}else{

V=round(sort(unique(Veg_reend.time[which(Veg_reend.time<=u)])),6)
YY=matrix(c(1:length(V)),1,length(V))
final=mclapply(YY,function(g){

uncensored.time=renewal_times+V[g] 
XX=matrix(c(1:length(uncensored.time)),1,length(uncensored.time))
Ratio=mclapply(XX,function(h){
df=data.frame(x=ifelse(round(Vegetation_data[which(Vegetation_data[,"Type"]==2),"end.time"],6)==round(uncensored.time[h],6)  & Vegetation_data[which(Vegetation_data[,"Type"]==2),"Censor"]==1,1,0),id=rep(1:length(table(Vegetation_data[,"Area"])),as.vector(table(Vegetation_data[which(Vegetation_data[,"Type"]==2),"Area"]))))
denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
numerator=sum(as.numeric(uncensored.time[h]<=Centime_vegetation)*exp((fit_result_pred$coef[2]*Urb_vegetation)))
ratio=c(denominator,numerator)
}## for function
)## for mclappy

K1=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,1]) ##
K2=sum(matrix(unlist(Ratio),4,2,byrow=TRUE)[,2])
rr=K1/K2
}## function fot g
)##mclapply
total=sum(unlist(final))
}##for else
return(total)
} ## for function

xx=seq(0,366,by=1)

type1_stru_cyclic_urban=c()
type1_stru_cyclic_rural=c()

type2_vege_cyclic_urban=c()
type2_vege_cyclic_rural=c()

for(k in 1:length(xx)){

type1_stru_cyclic_urban[k]=Cyclic.Structure.mean.function(xx[k])*exp(fit_result_pred$coef[1]*1)
type1_stru_cyclic_rural[k]=Cyclic.Structure.mean.function(xx[k])*exp(fit_result_pred$coef[1]*0) 
type2_vege_cyclic_urban[k]=Cyclic.Vegetation.mean.function(xx[k])*exp(fit_result_pred$coef[2]*1)
type2_vege_cyclic_rural[k]=Cyclic.Vegetation.mean.function(xx[k])*exp(fit_result_pred$coef[2]*0)
cat("number=",k,"\n")
}

months=c( "January","March","May","July","September","November")
plot(xx,type1_stru_cyclic_urban,type="l",lty=3,col="blue",xlim=c(0,366),xlab="Months in 2008",ylab="Prediction of Cumulative Number of Fires",xaxt='n')
axis(1, at=c(31,91,152,213,274,335), labels =months)
lines(xx,type1_stru_cyclic_rural,col="blue",lty=2)
lines(xx,type2_vege_cyclic_urban,col="2",lty=1)
lines(xx,type2_vege_cyclic_rural,col="2",lty=5)
legend("topleft", legend=c("Structure fire (Urban)", "Structure fire (Rural)","Vegetation fire (Urban)","Vegetation fire (Rural)"),
       col=c("blue","blue","2", "2"), lty=c(3,2,1,5))











