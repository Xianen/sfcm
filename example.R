#cmp_fcm_sfcm.R

library(e1071);
source("SparseFuzzyClustering.R");
source("CMeansSparseCluster.R");
source("calCER.R");
source("cvt_labels.R");

###generate simulation data###
n<- 300; #number of samples 
K<- 3; #number of classes
q<- 50; #actual number of dimensions
p<- 1000; #total number of dimensions 
m<- 1.1; #degree of fuzzification in cmeans 
u<- 3; #offset of cluster centers
sfcm_w1bnd<- 7; #w1 bound for sfcm

x <- matrix(rnorm(n*p), ncol=p);
c1start<- 1;
c1end<- n/3;
c2start<- c1end+1;
c2end<- 2*n/3;

x[c1start: c1end, 1: q]<- x[c1start: c1end, 1: q]+u;
x[c2start: c2end, 1: q]<- x[c2start: c2end, 1: q]-u;		
#the rest are 3rd class

x <- scale(x, TRUE, TRUE);
truelabels=c(rep(1,n/3), rep(2,n/3), rep(3,n/3)); 

###parameter tuning###
#this is not a necessary step, remove the comment below if you want to tune parameters, it will take a little bit long time;
#start of tuning
#source("CMeansSparseCluster.permute.R");
#TUNE<- CMeansSparseCluster.permute(x,K=K,nvals=50); 
#print(TUNE$bestw);
#pdf("Rplot_gap_stat.pdf");
#plot(TUNE$wbounds, type="o", TUNE$gaps, xlab="w1bound",ylab="gap statistics", main="gap statistics");
#dev.off();
#end of tuning


###compare cmeans and sfcm###
#1. fuzzy c means
fcm.out<- cmeans(x, K, m=m);
fcm.Cer<- calCER(fcm.out$cluster, truelabels);		
	
#2. sparse fuzzy c means
sfcm.out<-CMeansSparseCluster(x, K, sfcm_w1bnd); 
sfcm.WNonZerosNum<- sum(sfcm.out$ws!=0);
sfcm.Cer<- calCER(sfcm.out$Cs, truelabels);


###print and plot results###
#print our results
cat("fcm.Cer: ",fcm.Cer,"\n");
cat("sfcm.Cer: ",sfcm.Cer,"\nsfcm.WNonZerosNum", sfcm.WNonZerosNum,"\n");		

#plot results, the first two dimension
x1<- x[,1];
x2<- x[,2];

#align class no.
fcm_labels<- cvt_labels(fcm.out$cluster);
sfcm_labels<- cvt_labels(sfcm.out$Cs);

pdf(file="Rplot.pdf");
layout(matrix(1:3, 1,3));
plot(x1,x2, main="truelabels", sub="(a)",pch=truelabels,col=truelabels);
plot(x1,x2, main="fcm",sub="(b)", pch=fcm_labels, col=fcm_labels);
plot(x1,x2, main="sfcm", sub="(c)", pch=sfcm_labels, col=sfcm_labels);

layout(1);
plot(sfcm.out$ws, main="feature weights", xlab="feature index", ylab="weight value");
dev.off();
