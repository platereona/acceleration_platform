data<-matrix(scan("~/Desktop/moris/kmeans/output.txt"),3,10000)
png("~/Desktop/moris/kmeans/output.png", width = 500, height = 500)
plot(data[1,],data[2,],col=rainbow(10)[data[3,]],xlab="X1",ylab="X2",main="Lloyd k-means(n=10000,d=2,k=10)")