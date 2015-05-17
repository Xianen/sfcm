#The code is revised from Witten¡¯s sparse clustering framework, 
#which is available at http://cran.r-project.org/web/packages/sparcl/index.html.
`CMeansSparseCluster` <-
function(x, K=NULL, wbounds=NULL, nstart=20, silent=FALSE, maxiter=6, centers=NULL){
  	# The criterion is : minimize_{w, C} sum_j w_j (FWCSS_j - TSS_j) s.t. ||w||_2=1, ||w||_1<=s, w_j>=0
  	# x is the data, nxp
  	# K is the number of clusters desired
  	# wbounds is a vector of L1 constraints on w, of the form  sum(abs(w))<=wbounds[i]
 	if(is.null(K) && is.null(centers)) stop("Must provide either K or centers.")
 	if(!is.null(K) && !is.null(centers)){
    		if(nrow(centers)!=K) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    		if(nrow(centers)==K) K <- NULL
  	}
  	if(!is.null(centers) && ncol(centers)!=ncol(x)) stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  	if(is.null(wbounds)) wbounds <- seq(1.1, sqrt(ncol(x)), len=20)
  	if(min(wbounds)<=1) stop("wbounds should be greater than 1")
 	wbounds <- c(wbounds) # In case wbounds is a single number, turn it into a vector
  	
  	out <- list()
  	if(!is.null(K)) U <-cmeans(x, centers=K,m=m)$membership;
  	if(is.null(K))  U <-cmeans(x, centers=centers, m=m)$membership; 
  	#print(U);
  	for(i in 1:length(wbounds)){
    		if(length(wbounds)>1 && !silent) cat(i,fill=FALSE)
    		ws <- rep(1/sqrt(ncol(x)), ncol(x)) # Start with equal weights on each feature
    		ws.old <- rnorm(ncol(x))
    		store.bcss.ws <- NULL
    		niter <- 0
    		while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
      		if(!silent) cat(niter, fill=FALSE)
      		niter <- niter+1
      		ws.old <- ws
      		if(!is.null(K)){
        			if(niter>1)U<- UpdateUs(x, K, ws, U);   # if niter=1, no need to update!!
      		} else {
        		if(niter>1) U<- UpdateUs(x, nrow(centers), ws, U);  # if niter=1, no need to update!!
      		}
      		ws<- UpdateWs(x, U, wbounds[i]);  
      		store.bcss.ws <- c(store.bcss.ws, sum(GetFWCSS(x, U)$bcss.perfeature*ws))
    		}

    		Cs<- apply(U,1, which.max);
    		out[[i]] <- list(ws=ws, U=U, Cs=Cs, wcss=GetFWCSS(x, U, ws), crit=store.bcss.ws, wbound=wbounds[i]) 
  	}
  	if(!silent) cat(fill=TRUE)
  	if(length(wbounds)==1){
    		out <- out[[1]]
    		class(out) <- "cmeanssparse"
    	return(out)
  	}
  	class(out) <- "multicmeanssparse"
  	return(out)
}

plot.multicmeanssparse <- function(x,...){
  N <- length(x)
 # par(mfrow=c(ceiling(N/2),2))
 	#x11(); par(mfrow=c(2,2));
  for(i in 1:N){
    if(round((i-1)/4)*4==i-1)  {x11(); par(mfrow=c(2,2));}
    plot(x[[i]]$ws, main=paste("Wbound is ", sep="", round(x[[i]]$wbound,3)), xlab="Feature Index", ylab="Wj")
  }
}

plot.cmeanssparse <- function(x,...){
  plot(x$ws, main=paste("Wbound is ", sep="", round(x$wbound,3)), xlab="Feature Index", ylab="Wj")
}

PrintIt <- function(x){
  cat("Number of non-zero weights: ", sum(x$ws!=0), fill=TRUE)
  cat("Sum of weights: ", sum(x$ws), fill=TRUE)
  cat("Clustering: ", x$Cs, fill=TRUE)
  cat("Membership: ", x$U, fill=TRUE)
  cat(fill=TRUE)
}

print.cmeanssparse <- function(x,...){
  cat("Wbound is ", x$wbound, ":", fill=TRUE)
  PrintIt(x)
}

print.multicmeanssparse <- function(x,...){
  for(i in 1:length(x)){
  	cat("Wbound is ", x[[i]]$wbound, ":", fill=TRUE)
   	PrintIt(x[[i]])
  }  
}

statresults<- function(w1bounds, sfcmsout, truelabels ){
	cat("w1bounds    rightrates    Number of non-zero weights   Sum of weights \n");
	cat("_________________________________________________________________\n");
	results<- NULL;
	for(i in 1:length(sfcmsout)){
	#print(i);
		Cs<-sfcmsout[[i]]$Cs;
		errornums<- 0;
		for(j in 1:length(Cs)){
			if(Cs[j]!=truelabels[j]) errornums<- errornums+1;
		}
		rightrate=1-errornums/length(Cs);
		results<- rbind(results, c(w1bounds[i], rightrate,  sum(sfcmsout[[i]]$ws!=0),  sum(sfcmsout[[i]]$ws) ));
		#print(w1bounds[i]);		print(rightrate);		print(sum(sfcmsout[[i]]$ws!=0));		print(sum(sfcmsout[[i]]$ws));
		#print(c(w1bounds[i], rightrate,  sum(sfcmsout[[i]]$ws!=0),  sum(sfcmsout[[i]]$ws) ));
	}
	print(results);
}
