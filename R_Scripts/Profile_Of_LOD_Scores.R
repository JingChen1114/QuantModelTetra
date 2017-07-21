# R script written by Dr Jing Chen and Dr Lindsey Leach for plotting the LOD score profiles under
# 12 different parental genotype configurations at a putative QTL for a quantitative trait (See Figure S1). 
LOD=read.table("LOD_scores.txt",head=F)
x=LOD[,1]
x1=seq(1:50)
x2=seq(1:50)
x3=seq(1:50)
x4=seq(1:50)
x5=seq(1:50)
x6=seq(1:50)
x7=seq(1:50)
x8=seq(1:50)
x9=seq(1:50)
x10=seq(1:50)
x11=seq(1:50)
x12=seq(1:50)
for(i in seq(1:50)){
  x1[i]=x[i]
}
for(i in seq(1:50)){
  j=50+i
  x2[i]=x[j]
}
for(i in seq(1:50)){
  j=100+i
  x3[i]=x[j]
}
for(i in seq(1:50)){
  j=150+i
  x4[i]=x[j]
}
for(i in seq(1:50)){
  j=200+i
  x5[i]=x[j]
}
for(i in seq(1:50)){
  j=250+i
  x6[i]=x[j]
}
for(i in seq(1:50)){
  j=300+i
  x7[i]=x[j]
}
for(i in seq(1:50)){
  j=350+i
  x8[i]=x[j]
}
for(i in seq(1:50)){
  j=400+i
  x9[i]=x[j]
}
for(i in seq(1:50)){
  j=450+i
  x10[i]=x[j]
}
for(i in seq(1:50)){
  j=500+i
  x11[i]=x[j]
}
for(i in seq(1:50)){
  j=550+i
  x12[i]=x[j]
}

y=seq(1:50)
for(i in seq(1:50)){
  y[i]=0.005*(i-1)
}

par(mfrow=c(3,4))
plot(y,x1,type="l",main="Model[1]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=1,ylim=c(0,7))
plot(y,x2,type="l",main="Model[2]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x3,type="l",main="Model[3]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x4,type="l",main="Model[4]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x5,type="l",main="Model[5]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x6,type="l",main="Model[6]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x7,type="l",main="Model[7]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x8,type="l",main="Model[8]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x9,type="l",main="Model[9]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x10,type="l",main="Model[10]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x11,type="l",main="Model[11]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))
plot(y,x12,type="l",main="Model[12]",ylab="LOD scores",xlab="Coefficient of Double Reduction",col="red",lwd=2,ylim=c(0,7))