#The code is revised from Witten¡¯s sparse clustering framework, 
#which is available at http://cran.r-project.org/web/packages/sparcl/index.html.

#SparseFuzzyClustering.R

GetFWCSS<- function(x, U, ws=NULL){
	wcss.perfeature <- numeric(ncol(x))
  	K=ncol(U);  
  	ms_sum<- apply(U^m,2,sum); #ms:membership
  
  	for(k in 1:K){
	  	ms_k<- diag(U[,k])^m;
	  	center_k<- apply(ms_k %*% x,2,sum)/ms_sum[k];
	  	wcss.perfeature<- wcss.perfeature+ apply(ms_k %*% (x-matrix(rep(1,n)) %*% center_k)^2, 2, sum);  	
  	}
  
	#print(wcss.perfeature);
 	bcss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.perfeature
  	if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature))
  	if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}

UpdateUs <- function(x, K, ws, U){
  	x <- x[,ws!=0]
  	z <- sweep(x, 2, sqrt(ws[ws!=0]), "*")
  	nrowz <- nrow(z)
  	
  	ms_sum<- apply(U,2,sum);  #ms: membership
  	mus <- NULL;
  	for(k in 1:K){
  		ms_k<- diag(U[,k])^m;
  		center_k<- apply(ms_k %*% z, 2, sum)/ms_sum[k];
  		mus<-rbind(mus,center_k);
  	}  	
  	#print(mus);

  	if(is.null(U)){
  		cm<- cmeans(x, centers=K, m=m);		
  	} else {
  		cm <- cmeans(z, centers=mus,m=m);    
  	}
  	return(cm$membership);  	 
}


UpdateWs <- function(x, U, l1bound){
  	wcss.perfeature <- GetFWCSS(x, U)$wcss.perfeature
  	
  	tss.perfeature <-    apply(scale(x,center=TRUE, scale=FALSE)^2, 2, sum)            
  	#GetFWCSS(x, rep(1, nrow(x)))$wcss.perfeature
  	#print(wcss.perfeature);
  	#print(tss.perfeature);
  	
  	#print("-wcss.perfeature+tss.perfeature:\n");
  	#print(tss.perfeature);
  	
  	lam <- BinarySearch(-wcss.perfeature+tss.perfeature, l1bound)
  	ws.unscaled <- soft(-wcss.perfeature+tss.perfeature,lam)
  	return(ws.unscaled/l2n(ws.unscaled))
}

################################################################
BinarySearch <- function(argu,sumabs){
  	if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  	lam1 <- 0
  	lam2 <- max(abs(argu))-1e-5
  	iter <- 1
  	while(iter<=15 && (lam2-lam1)>(1e-4)){
    		su <- soft(argu,(lam1+lam2)/2)
    		if(sum(abs(su/l2n(su)))<sumabs){
     		lam2 <- (lam1+lam2)/2
    		} else {
      		lam1 <- (lam1+lam2)/2
    		}
    		iter <- iter+1
  		}
  	return((lam1+lam2)/2)
}

soft <- function(x,d){
  	return(sign(x)*pmax(0, abs(x)-d))
}

l2n <- function(vec){
  	return(sqrt(sum(vec^2)))
}
