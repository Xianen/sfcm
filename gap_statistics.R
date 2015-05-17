source("CMeansSparseCluster.permute.R");
K<- 3;
TUNE<- CMeansSparseCluster.permute(x,K=K,nvals=50); print(TUNE$bestw);
pdf("Rplot_gap_stat.pdf");
#tiff(filename="fig_1.tif",width=800, height=800,res=150);
plot(TUNE$wbounds, type="o", TUNE$gaps, xlab="w1bound",ylab="gap statistics", main="gap statistics");
#dev.off();

#tiff(filename="fig_1_b.tif",width=800, height=800,res=150);
 plot(sfcm.out$ws, main="feature weights", xlab="feature index", ylab="weight value");
 dev.off();
 
