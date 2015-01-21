data<-matrix(scan("~/Dropbox/kiryu/mixem/algo/result.txt"),3,200)
png("~/Dropbox/kiryu/mixem/algo/classter.png", width = 500, height = 500)
plot(data[1,],data[2,],col=rainbow(7)[data[3,]],xlab="X",ylab="Y")