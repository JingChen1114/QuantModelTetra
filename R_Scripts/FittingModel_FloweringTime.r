#R script written by Dr Lindsey Leach and Dr Jing Chen to produce figure 1e.
ph=read.table("data_in_FloweringTime.txt",head=T)
x=ph[,1]
hist(x,xlim=range(x),main="Flowering Time",xlab="",nclass=30,freq=F)
max1=max(x)
min1=min(x)
x1=seq(min1,max1,0.01)
#offspring probability distribution for the 5 genotypes (from data_out.txt) 
f0=0.0342
f1=0.2331
f2=0.4654
f3=0.2331
f4=0.0342
#offspring genotypic values for the five genotypes estimated using the orthogonal model
g0=25.4637
g1=28.4914
g2=32.3313
g3=40.1198
g4=49.9044
e=2.9647
y1=(1/(((2*pi)^(1/2))*e))*(f0*exp(-((x1-g0)^2)/(2*(e^2)))+
    f1*exp(-((x1-g1)^2)/(2*(e^2)))+f2*exp(-((x1-g2)^2)/(2*
    (e^2)))+f3*exp(-((x1-g3)^2)/(2*(e^2)))+f4*exp(-((x1-g4)^2)
    /(2*(e^2))))

lines(x1,y1,type="l",col="red",lwd=2)

y10=(1/(((2*pi)^(1/2))*e))*(f0*exp(-((x1-g0)^2)/(2*(e^2))))
lines(x1,y10,type="p",col="blue",cex=0.1)

y11=(1/(((2*pi)^(1/2))*e))*(f1*exp(-((x1-g1)^2)/(2*(e^2))))
lines(x1,y11,type="p",col="blue",cex=0.1)

y12=(1/(((2*pi)^(1/2))*e))*(f2*exp(-((x1-g2)^2)/(2*(e^2))))
lines(x1,y12,type="p",col="blue",cex=0.1)

y13=(1/(((2*pi)^(1/2))*e))*(f3*exp(-((x1-g3)^2)/(2*(e^2))))
lines(x1,y13,type="p",col="blue",cex=0.1)

y14=(1/(((2*pi)^(1/2))*e))*(f4*exp(-((x1-g4)^2)/(2*(e^2))))
lines(x1,y14,type="p",col="blue",cex=0.1)