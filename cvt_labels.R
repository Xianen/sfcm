cvt_labels<- function(labels){
#re-assign labels to make them start from 1
maxlabel<- max(labels);
map<- seq(-1,-1,len=maxlabel);
len<- length(labels);
curr_label<- 1;
for(i in 1:len){
	if(-1==map[labels[i]]) {
		map[labels[i]]<- curr_label;
		curr_label<- curr_label+1;
	}
}
for(i in 1: len){
	labels[i]<- map[labels[i]];
}

return(labels);
}
