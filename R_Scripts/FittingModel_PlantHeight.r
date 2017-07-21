#R script written by Dr Lindsey Leach and Dr Jing Chen to produce figure 1f.
ph=read.table("data_in_PlantHeight.txt",head=T)
x=ph[,1]
hist(x,xlim=range(x),main="Plant Height",xlab="",nclass=30,freq=F)
max1=max(x)
min1=min(x)
x1=seq(min1,max1,0.01)
#offspring probability distribution for the 5 genotypes (from data_out.txt) 
f2=0.0413
f3=0.4175
f4=0.5413
#offspring genotypic values for the five genotypes estimated using the orthogonal model
g2=22.6260
g3=40.6551
g4=51.8400
e=6.4031
y1=(1/(((2*pi)^(1/2))*e))*(f2*exp(-((x1-g2)^2)/(2*
    (e^2)))+f3*exp(-((x1-g3)^2)/(2*(e^2)))+f4*exp(-((x1-g4)^2)
    /(2*(e^2))))

lines(x1,y1,type="l",col="red",lwd=2)


y12=(1/(((2*pi)^(1/2))*e))*(f2*exp(-((x1-g2)^2)/(2*(e^2))))
lines(x1,y12,type="p",col="blue",cex=0.1)

y13=(1/(((2*pi)^(1/2))*e))*(f3*exp(-((x1-g3)^2)/(2*(e^2))))
lines(x1,y13,type="p",col="blue",cex=0.1)

y14=(1/(((2*pi)^(1/2))*e))*(f4*exp(-((x1-g4)^2)/(2*(e^2))))
lines(x1,y14,type="p",col="blue",cex=0.1)