################################################################################################
################################################################################################

setwd("C:/Jonathan/Dropbox/BLUE/Zurich/LPI/LivingPlanetIndexDatabase_2025-02-01_19.19.51")
data<-read.csv("LPD_2024_public.csv", header=T)

Af1<-subset(data, Region == "Africa")
Af2<-subset(Af1, Included.in.LPR2024 == 1)
Af3<-subset(Af2, Class == "Mammalia")

#ignoring units for now
#unique(Af$Units)

# subset 5+ counts
obs.cnt<-NULL
for (n in 1:dim(Af3)[1]){
obs.cnt[n]<-sum(Af3[n, 31:101]!="NULL")
}

Af4<-Af3[obs.cnt>4,]

# subset species 10+ populations
sp<-unique(Af4$Binomial)#

sp.cnt<-NULL
for(x in 1:length(sp)){
Af.tmp<-subset(Af4, Binomial == sp[x])
sp.cnt[x]<-dim(Af.tmp)[1]
}

keep.sp<-sp[sp.cnt>9]
Af5<-Af4[Af4$Binomial %in% keep.sp,]
#14 sp remaining
