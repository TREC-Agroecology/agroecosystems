setwd ("/Users/sreader/Documents/GitHub/agroecosystems/data")
BD<-read.csv('soil-BD-summer-2019-ECHO.csv', header=T, as.is = T)
head(BD)
##we lost one sample (181 post) So here I take out 181 post then pre
vec=is.na(BD$SoilBD_gcm3)
BDnew=BD[!vec,]
BDnew2=BDnew[-37,]

#T.test
Pre<-subset(BDnew2, Time=='pre')
Post<-subset(BDnew2, Time=='post')
t.test(x=Pre$SoilBD_gcm3, y=Post$SoilBD_gcm3, alternative="two.sided", paired=T)

##ANOVAs
BD<-BDnew2$SoilBD_gcm3
Time<-BDnew2$Time
CropTrt<-BDnew2$CropTrt
BDaov<-aov(BD~CropTrt, data = BDnew2)
summary(BDaov)
BDaov2<-aov(BD~CropTrt+Time, data=BDnew2)
summary(BDaov2)
BDaov3<-aov(BD~CropTrt*Time, data=BDnew2)
summary(BDaov3)
BDaov4<-aov(BD~CropTrt*SiteName, data=BDnew2)
summary(BDaov4)
BDaov5<-aov(BD~CropTrt+SiteName, data=BDnew2)
summary(BDaov5)
BDaov6<-aov(BD~CropTrt+SiteName+Time, data=BDnew2)
summary(BDaov6)
BDaov7<-aov(BD~CropTrt+SiteName*Time, data=BDnew2)
summary(BDaov7)

#Some plots
par (mfrow=c(2,2))
boxplot(BD~CropTrt, data = BDnew2)
boxplot(BD~Time, data=BDnew2)
boxplot (BD~Time*CropTrt)
boxplot (BD~SiteName, data=BDnew2)
