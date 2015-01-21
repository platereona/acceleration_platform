library(mvtnorm)
library(scatterplot3d)
data<-matrix(scan("~/Dropbox/kiryu/mixem/algo/output.txt"),6,119)
x1<-seq(from=0,to=1,length=140)
x2<-x1
for(i in 1:17){
m1<-c(data[1,(i-1)*7+1],data[2,(i-1)*7+1])
m2<-c(data[1,(i-1)*7+2],data[2,(i-1)*7+2])
m3<-c(data[1,(i-1)*7+3],data[2,(i-1)*7+3])
m4<-c(data[1,(i-1)*7+4],data[2,(i-1)*7+4])
m5<-c(data[1,(i-1)*7+5],data[2,(i-1)*7+5])
m6<-c(data[1,(i-1)*7+6],data[2,(i-1)*7+6])
m7<-c(data[1,(i-1)*7+7],data[2,(i-1)*7+7])
s1<-matrix(c(data[3,(i-1)*7+1],data[5,(i-1)*7+1],data[5,(i-1)*7+1],data[4,(i-1)*7+1]),2,2)
s2<-matrix(c(data[3,(i-1)*7+2],data[5,(i-1)*7+2],data[5,(i-1)*7+2],data[4,(i-1)*7+2]),2,2)
s3<-matrix(c(data[3,(i-1)*7+3],data[5,(i-1)*7+3],data[5,(i-1)*7+3],data[4,(i-1)*7+3]),2,2)
s4<-matrix(c(data[3,(i-1)*7+4],data[5,(i-1)*7+4],data[5,(i-1)*7+4],data[4,(i-1)*7+4]),2,2)
s5<-matrix(c(data[3,(i-1)*7+5],data[5,(i-1)*7+5],data[5,(i-1)*7+5],data[4,(i-1)*7+5]),2,2)
s6<-matrix(c(data[3,(i-1)*7+6],data[5,(i-1)*7+6],data[5,(i-1)*7+6],data[4,(i-1)*7+6]),2,2)
s7<-matrix(c(data[3,(i-1)*7+7],data[5,(i-1)*7+7],data[5,(i-1)*7+7],data[4,(i-1)*7+7]),2,2)
pi1<-data[6,(i-1)*7+1]
pi2<-data[6,(i-1)*7+2]
pi3<-data[6,(i-1)*7+3]
pi4<-data[6,(i-1)*7+4]
pi5<-data[6,(i-1)*7+5]
pi6<-data[6,(i-1)*7+6]
pi7<-data[6,(i-1)*7+7]
f1<-function(x1, x2){dmvnorm(matrix(c(x1, x2), ncol=2), mean=m1, sigma=s1)}
f2<-function(x1, x2){dmvnorm(matrix(c(x1, x2), ncol=2), mean=m2, sigma=s2)}
f3<-function(x1, x2){dmvnorm(matrix(c(x1, x2), ncol=2), mean=m3, sigma=s3)}
f4<-function(x1, x2){dmvnorm(matrix(c(x1, x2), ncol=2), mean=m4, sigma=s4)}
f5<-function(x1, x2){dmvnorm(matrix(c(x1, x2), ncol=2), mean=m5, sigma=s5)}
f6<-function(x1, x2){dmvnorm(matrix(c(x1, x2), ncol=2), mean=m6, sigma=s6)}
f7<-function(x1, x2){dmvnorm(matrix(c(x1, x2), ncol=2), mean=m7, sigma=s7)}
z1<-outer(x1, x2, f1)
z2<-outer(x1, x2, f2)
z3<-outer(x1, x2, f3)
z4<-outer(x1, x2, f4)
z5<-outer(x1, x2, f5)
z6<-outer(x1, x2, f6)
z7<-outer(x1, x2, f7)
x<-seq(from=0,to=1,length=140)
y<-seq(from=0,to=1,length=140)
z<-pi1*z1+pi2*z2+pi3*z3+pi4*z4+pi5*z5+pi6*z6+pi7*z7
png(sprintf("~/Dropbox/kiryu/mixem/algo/persp%d.png",i), width = 500, height = 500)
persp(x, y, z, theta = 330, phi = 45, expand = 0.5, col = "pink",ticktype="simple",xlab="X",ylab="Y",zlab="",main=sprintf("Repeat Count = %d",i))
dev.off()
}
