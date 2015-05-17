calCER<- function(Cs1, Cs2){

	if(length(Cs1)!= length(Cs2))  stop("calCER: Cs1 and Cs2 are different length!");
	ErrNum<- 0;
	n<- length(Cs1);
	
	for(i in 1:(n-1)){
		for(j in( i+1):n){
			if(Cs1[i]==Cs1[j] && Cs2[i]!=Cs2[j]) ErrNum<- ErrNum+1;
			if(Cs1[i]!=Cs1[j] && Cs2[i]==Cs2[j]) ErrNum<- ErrNum+1;
			
		}
	}
	
	out<- ErrNum/(n*(n-1)/2);
	return(out);
}
