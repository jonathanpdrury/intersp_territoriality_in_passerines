add.outgroup <- function(phy, root.node=NULL){ #from Sergio Sanchez https://stat.ethz.ch/pipermail/r-sig-phylo/2015-February/003886.html
	require(ape)
	if (is.null(phy$edge.length)){
	tt <- rtree(2)
	tt$edge.length <- NULL
	tt$tip.label <- c("outgroup","drop")
	phy$root.edge <- 1
	if (!is.null(root.node)){
	rphy <- root(phy, node=root.node)
	ot <- bind.tree(rphy, tt, position=1)
	} else {
	ot <- bind.tree(phy, tt, position=1)
	}
	ot <- bind.tree(phy, tt, position=1)
	res <- drop.tip(ot, "drop")
	return(res)
	} else if (is.ultrametric(phy)){
	th <- as.numeric(sort(branching.times(phy), decreasing=TRUE))[1]
	re <- th/5
	phy$root.edge <- re
	tt <- rtree(2)
	tt$edge.length <- c(0,0)
	tt$tip.label <- c("outgroup","drop")
	tt$root.edge <- th + re
	ot <- bind.tree(phy, tt, position=re)
	res <- drop.tip(ot, "drop")
	return(res)
	} else {
	tl <- max(phy$edge.length)
	re <- tl/3
	tt <- rtree(2)
	tt$edge.length <- c(0,0)
	tt$tip.label <- c("outgroup","drop")
	tt$root.edge <- tl
	if (!is.null(root.node)){
	rphy <- root(phy, node=root.node)
	rphy$root.edge <- re
	ot <- bind.tree(rphy, tt, position=re)
	} else {
	phy$root.edge <- re
	ot <- bind.tree(phy, tt, position=re)
	}
	res <- drop.tip(ot, "drop")
	return(res)
	}
}

random.effect.sorting<-function(data,counter.max=1e6,seed=Sys.time()){
	if(class(data)!="matrix"){data<-as.matrix(data)}
	check<-unique(c(levels(as.factor(data[,1])),levels(as.factor(data[,2])))) ##arrange the orders of sp1 and sp2 so that they are as evenly distributed as possible
	ord=2
	counter=0
	half=dim(data)[1]/2
	half.samp=sample(1:dim(data)[1],half)
	data[half.samp,]<-data[half.samp,c(2,1)]
	set.seed(seed)
	while(ord[1]>1 && counter < counter.max){
		counter=counter+1
		t1<-table(data[,1])
		t2<-table(data[,2])
		t1<-t1[match(check,names(t1))]
		names(t1)<-check
		t2<-t2[match(check,names(t2))]
		names(t2)<-check
		#int<-data.frame(check,as.vector(t1),as.vector(t2))
		#colnames(int)<-c("check","t1","t2")	
		t1[is.na(t1)] <- 0
		t2[is.na(t2)] <- 0
		diff<-t1-t2
		ord<-sort(abs(diff),decreasing=TRUE)
		for(j in 1:4){
			if(ord[j]>1){
				rowx<-sample(seq(1:(dim(data[which(data[,1]==names(ord[j])| data[,2]==names(ord[j])),])[1])),(ord[j]/2))
				data[which(data[,1]==names(ord[j])| data[,2]==names(ord[j])),][rowx,]<-data[which(data[,1]==names(ord[j])| data[,2]==names(ord[j])),][rowx,c(2,1)]
			}
		}	
	}	
	if(counter==counter.max){print("hit counter max")}
	return(list(sorted.data=data,asymmetry=ord))
}	
