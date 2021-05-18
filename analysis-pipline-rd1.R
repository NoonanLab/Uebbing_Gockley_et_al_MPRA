### FUNCTIONS AND PACKAGES ###
library(tidyr)
paste5<-function(...,sep=" ",collapse=NULL,na.rm=F){
	if(na.rm==F){
		paste(...,sep=sep,collapse=collapse)
	}else{
		if(na.rm==T){
			paste.na<-function(x,sep){
				x<-gsub("^\\s+|\\s+$","",x)
				ret<-paste(na.omit(x),collapse=sep)
				is.na(ret)<-ret==""
				return(ret)
			}
			df<-data.frame(...,stringsAsFactors=F)
			ret<-apply(df,1,FUN=function(x)paste.na(x,sep))
			if(is.null(collapse)){
				ret
			}else{
				paste.na(ret,sep=collapse)
}}}}
# http://stackoverflow.com/questions/13673894/suppress-nas-in-paste
"%out%"<-function(x,table){match(x,table,nomatch=0)==0}

### INPUT DATA ###
MPRAdata1<-read.delim(gzfile("data/COUNTS_ANNOTATED.txt.gz"),stringsAsFactors=F) # primary input data
if(file.exists("data/MPRAsum1.tsv")){ # If file available, load; else create
	MPRAsum1<-read.delim("data/MPRAsum1.tsv",stringsAsFactors=F)
}else{
	
	### NORMALIZATION AND FILTERS ###
	lib.sizes<-colSums(MPRAdata1[,10:19])
	MPRAnorm1<-MPRAdata1
	for(i in 1:10){
		MPRAnorm1[,(9+i)]<-log2((MPRAnorm1[,(9+i)]+.1)*1e6/lib.sizes[i])
	}
	# Find Poisson error cutoff #
	par(mfrow=c(2,2),las=1)
	plot(density(MPRAnorm1$Rep_1_2_pDNA_Tag_Counts),main="Rep 12",xlab=expression(paste(log[2]," pDNA counts")),ylab="Relative amount",yaxt="n")
	abline(v=-5.25,col=2)
	plot(density(MPRAnorm1$Rep_1_3_pDNA_Tag_Counts),main="Rep 13",xlab=expression(paste(log[2]," pDNA counts")),ylab="Relative amount",yaxt="n")
	abline(v=-5.25,col=2)
	plot(density(MPRAnorm1$Rep_2_2_pDNA_Tag_Counts),main="Rep 22",xlab=expression(paste(log[2]," pDNA counts")),ylab="Relative amount",yaxt="n")
	abline(v=-5.25,col=2)
	plot(density(MPRAnorm1$Rep_2_3_pDNA_Tag_Counts),main="Rep 23",xlab=expression(paste(log[2]," pDNA counts")),ylab="Relative amount",yaxt="n")
	abline(v=-5.25,col=2)
	
	MPRAnorm1<-MPRAnorm1[rowSums(MPRAnorm1[,c(11,16:19)])/5 >= -5.25,] # cutoff for Poisson error (only pDNA!)
	MPRAsum1<-cbind(MPRAnorm1[,c(1,2,4:19)],Lin1_cDNA_Tag_Counts=rowSums(MPRAnorm1[,12:13])/2,Lin2_cDNA_Tag_Counts=rowSums(MPRAnorm1[,14:15])/2,
		Lin1_pDNA_Tag_Counts=rowSums(MPRAnorm1[,16:17])/2,Lin2_pDNA_Tag_Counts=rowSums(MPRAnorm1[,18:19])/2)
	rm(list=c("i","lib.sizes"))
	write.table(MPRAsum1,"data/MPRAsum1.tsv",row.names=F,sep="\t",quote=F)
}

if(file.exists("data/MPRAfrag1.tsv")){ # If file available load; else create
	MPRAfrag1<-read.delim("data/MPRAfrag1.tsv",stringsAsFactors=F)
	names2List<-load("data/fragnamelist.RData")
}else{
	### SORT BARCODES INTO FRAGMENTS ###
	Rep12_cDNAdata<-MPRAsum1[,c(2,11)] %>% nest(-Alignment)
	lens<-double(length(Rep12_cDNAdata$data))
	names(lens)<-Rep12_cDNAdata$Alignment
	for(i in 1:length(lens)){
		lens[i]<-dim(Rep12_cDNAdata$data[[i]])[1]
	}
	Rep12_cDNAsortedVals<-unlist(Rep12_cDNAdata$data)
	names(Rep12_cDNAsortedVals)<-rep(Rep12_cDNAdata$Alignment,times=lens)
	Rep12_cDNAmedian<-tapply(Rep12_cDNAsortedVals,names(Rep12_cDNAsortedVals),median)
	
	Rep13_cDNAdata<-MPRAsum1[,c(2,12)] %>% nest(-Alignment)
	Rep13_cDNAsortedVals<-unlist(Rep13_cDNAdata$data)
	names(Rep13_cDNAsortedVals)<-rep(Rep13_cDNAdata$Alignment,times=lens)
	Rep13_cDNAmedian<-tapply(Rep13_cDNAsortedVals,names(Rep13_cDNAsortedVals),median)
	
	Rep22_cDNAdata<-MPRAsum1[,c(2,13)] %>% nest(-Alignment)
	Rep22_cDNAsortedVals<-unlist(Rep22_cDNAdata$data)
	names(Rep22_cDNAsortedVals)<-rep(Rep22_cDNAdata$Alignment,times=lens)
	Rep22_cDNAmedian<-tapply(Rep22_cDNAsortedVals,names(Rep22_cDNAsortedVals),median)
	
	Rep23_cDNAdata<-MPRAsum1[,c(2,14)] %>% nest(-Alignment)
	Rep23_cDNAsortedVals<-unlist(Rep23_cDNAdata$data)
	names(Rep23_cDNAsortedVals)<-rep(Rep23_cDNAdata$Alignment,times=lens)
	Rep23_cDNAmedian<-tapply(Rep23_cDNAsortedVals,names(Rep23_cDNAsortedVals),median)
	
	Lin1_cDNAdata<-MPRAsum1[,c(2,19)] %>% nest(-Alignment)
	Lin1_cDNAsortedVals<-unlist(Lin1_cDNAdata$data)
	names(Lin1_cDNAsortedVals)<-rep(Lin1_cDNAdata$Alignment,times=lens)
	Lin1_cDNAmedian<-tapply(Lin1_cDNAsortedVals,names(Lin1_cDNAsortedVals),median)
	
	Lin2_cDNAdata<-MPRAsum1[,c(2,20)] %>% nest(-Alignment)
	Lin2_cDNAsortedVals<-unlist(Lin2_cDNAdata$data)
	names(Lin2_cDNAsortedVals)<-rep(Lin2_cDNAdata$Alignment,times=lens)
	Lin2_cDNAmedian<-tapply(Lin2_cDNAsortedVals,names(Lin2_cDNAsortedVals),median)
	
	Rep12_pDNAdata<-MPRAsum1[,c(2,15)] %>% nest(-Alignment)
	Rep12_pDNAsortedVals<-unlist(Rep12_pDNAdata$data)
	names(Rep12_pDNAsortedVals)<-rep(Rep12_pDNAdata$Alignment,times=lens)
	Rep12_pDNAmedian<-tapply(Rep12_pDNAsortedVals,names(Rep12_pDNAsortedVals),median)
	
	Rep13_pDNAdata<-MPRAsum1[,c(2,16)] %>% nest(-Alignment)
	Rep13_pDNAsortedVals<-unlist(Rep13_pDNAdata$data)
	names(Rep13_pDNAsortedVals)<-rep(Rep13_pDNAdata$Alignment,times=lens)
	Rep13_pDNAmedian<-tapply(Rep13_pDNAsortedVals,names(Rep13_pDNAsortedVals),median)
	
	Rep22_pDNAdata<-MPRAsum1[,c(2,17)] %>% nest(-Alignment)
	Rep22_pDNAsortedVals<-unlist(Rep22_pDNAdata$data)
	names(Rep22_pDNAsortedVals)<-rep(Rep22_pDNAdata$Alignment,times=lens)
	Rep22_pDNAmedian<-tapply(Rep22_pDNAsortedVals,names(Rep22_pDNAsortedVals),median)
	
	Rep23_pDNAdata<-MPRAsum1[,c(2,18)] %>% nest(-Alignment)
	Rep23_pDNAsortedVals<-unlist(Rep23_pDNAdata$data)
	names(Rep23_pDNAsortedVals)<-rep(Rep23_pDNAdata$Alignment,times=lens)
	Rep23_pDNAmedian<-tapply(Rep23_pDNAsortedVals,names(Rep23_pDNAsortedVals),median)
	
	Lin1_pDNAdata<-MPRAsum1[,c(2,21)] %>% nest(-Alignment)
	Lin1_pDNAsortedVals<-unlist(Lin1_pDNAdata$data)
	names(Lin1_pDNAsortedVals)<-rep(Lin1_pDNAdata$Alignment,times=lens)
	Lin1_pDNAmedian<-tapply(Lin1_pDNAsortedVals,names(Lin1_pDNAsortedVals),median)
	
	Lin2_pDNAdata<-MPRAsum1[,c(2,22)] %>% nest(-Alignment)
	Lin2_pDNAsortedVals<-unlist(Lin2_pDNAdata$data)
	names(Lin2_pDNAsortedVals)<-rep(Lin2_pDNAdata$Alignment,times=lens)
	Lin2_pDNAmedian<-tapply(Lin2_pDNAsortedVals,names(Lin2_pDNAsortedVals),median)
	
	### pDNA AND cDNA CORRELATIONS BETWEEN REPLICATES ###
	cor.test(Rep12_cDNAsortedVals,Rep13_cDNAsortedVals,method="spearman") # rho = 0.8303966, P = 0
	cor.test(Rep12_cDNAsortedVals,Rep22_cDNAsortedVals,method="spearman") # rho = 0.8522068, P = 0
	cor.test(Rep12_cDNAsortedVals,Rep23_cDNAsortedVals,method="spearman") # rho = 0.8110765, P = 0
	cor.test(Rep13_cDNAsortedVals,Rep22_cDNAsortedVals,method="spearman") # rho = 0.8358736, P = 0
	cor.test(Rep13_cDNAsortedVals,Rep23_cDNAsortedVals,method="spearman") # rho = 0.8091383, P = 0
	cor.test(Rep22_cDNAsortedVals,Rep23_cDNAsortedVals,method="spearman") # rho = 0.8272762, P = 0
	
	cor.test(Rep12_pDNAsortedVals,Rep13_pDNAsortedVals,method="spearman") # rho = 0.8839435, P = 0
	cor.test(Rep12_pDNAsortedVals,Rep22_pDNAsortedVals,method="spearman") # rho = 0.8947083, P = 0
	cor.test(Rep12_pDNAsortedVals,Rep23_pDNAsortedVals,method="spearman") # rho = 0.882196, P = 0
	cor.test(Rep13_pDNAsortedVals,Rep22_pDNAsortedVals,method="spearman") # rho = 0.874903, P = 0
	cor.test(Rep13_pDNAsortedVals,Rep23_pDNAsortedVals,method="spearman") # rho = 0.8773553, P = 0
	cor.test(Rep22_pDNAsortedVals,Rep23_pDNAsortedVals,method="spearman") # rho = 0.8727751, P = 0
	
	cor.test(Rep12_cDNAmedian,Rep13_cDNAmedian,method="spearman") # rho = 0.8542519, P = 0
	cor.test(Rep12_cDNAmedian,Rep22_cDNAmedian,method="spearman") # rho = 0.8678822, P = 0
	cor.test(Rep12_cDNAmedian,Rep23_cDNAmedian,method="spearman") # rho = 0.8404657, P = 0
	cor.test(Rep13_cDNAmedian,Rep22_cDNAmedian,method="spearman") # rho = 0.8593541, P = 0
	cor.test(Rep13_cDNAmedian,Rep23_cDNAmedian,method="spearman") # rho = 0.8415389, P = 0
	cor.test(Rep22_cDNAmedian,Rep23_cDNAmedian,method="spearman") # rho = 0.8509453, P = 0
	
	cor.test(Rep12_pDNAmedian,Rep13_pDNAmedian,method="spearman") # rho = 0.8803396, P = 0
	cor.test(Rep12_pDNAmedian,Rep22_pDNAmedian,method="spearman") # rho = 0.888206, P = 0
	cor.test(Rep12_pDNAmedian,Rep23_pDNAmedian,method="spearman") # rho = 0.8813426, P = 0
	cor.test(Rep13_pDNAmedian,Rep22_pDNAmedian,method="spearman") # rho = 0.8732541, P = 0
	cor.test(Rep13_pDNAmedian,Rep23_pDNAmedian,method="spearman") # rho = 0.8762132, P = 0
	cor.test(Rep22_pDNAmedian,Rep23_pDNAmedian,method="spearman") # rho = 0.8725649, P = 0
	
	### ACTIVITY CALCULATION WITH 12 BARCODE TAG CUTOFF ###
	Rep12_cDNAsortedVals12<-Rep12_cDNAsortedVals[names(Rep12_cDNAsortedVals) %in% names(lens[lens>11])]
	Rep12_pDNAsortedVals12<-Rep12_pDNAsortedVals[names(Rep12_pDNAsortedVals) %in% names(lens[lens>11])]
	Rep12_activitySortedVals12<-Rep12_cDNAsortedVals12-Rep12_pDNAsortedVals12
	sds<-double(length(unique(names(Rep12_activitySortedVals12))))
	names(sds)<-unique(names(Rep12_activitySortedVals12))
	sds<-tapply(Rep12_activitySortedVals12,names(Rep12_activitySortedVals12),sd)
	Rep12_activitySortedVals12sd<-Rep12_activitySortedVals12[names(Rep12_activitySortedVals12) %in% names(sds[sds!=0])]
	
	Rep13_cDNAsortedVals12<-Rep13_cDNAsortedVals[names(Rep13_cDNAsortedVals) %in% names(lens[lens>11])]
	Rep13_pDNAsortedVals12<-Rep13_pDNAsortedVals[names(Rep13_pDNAsortedVals) %in% names(lens[lens>11])]
	Rep13_activitySortedVals12<-Rep13_cDNAsortedVals12-Rep13_pDNAsortedVals12
	sds<-double(length(unique(names(Rep13_activitySortedVals12))))
	names(sds)<-unique(names(Rep13_activitySortedVals12))
	sds<-tapply(Rep13_activitySortedVals12,names(Rep13_activitySortedVals12),sd)
	Rep13_activitySortedVals12sd<-Rep13_activitySortedVals12[names(Rep13_activitySortedVals12) %in% names(sds[sds!=0])]
	
	Rep22_cDNAsortedVals12<-Rep22_cDNAsortedVals[names(Rep22_cDNAsortedVals) %in% names(lens[lens>11])]
	Rep22_pDNAsortedVals12<-Rep22_pDNAsortedVals[names(Rep22_pDNAsortedVals) %in% names(lens[lens>11])]
	Rep22_activitySortedVals12<-Rep22_cDNAsortedVals12-Rep22_pDNAsortedVals12
	sds<-double(length(unique(names(Rep22_activitySortedVals12))))
	names(sds)<-unique(names(Rep22_activitySortedVals12))
	sds<-tapply(Rep22_activitySortedVals12,names(Rep22_activitySortedVals12),sd)
	Rep22_activitySortedVals12sd<-Rep22_activitySortedVals12[names(Rep22_activitySortedVals12) %in% names(sds[sds!=0])]
	
	Rep23_cDNAsortedVals12<-Rep23_cDNAsortedVals[names(Rep23_cDNAsortedVals) %in% names(lens[lens>11])]
	Rep23_pDNAsortedVals12<-Rep23_pDNAsortedVals[names(Rep23_pDNAsortedVals) %in% names(lens[lens>11])]
	Rep23_activitySortedVals12<-Rep23_cDNAsortedVals12-Rep23_pDNAsortedVals12
	sds<-double(length(unique(names(Rep23_activitySortedVals12))))
	names(sds)<-unique(names(Rep23_activitySortedVals12))
	sds<-tapply(Rep23_activitySortedVals12,names(Rep23_activitySortedVals12),sd)
	Rep23_activitySortedVals12sd<-Rep23_activitySortedVals12[names(Rep23_activitySortedVals12) %in% names(sds[sds!=0])]
	
	Lin1_cDNAsortedVals12<-Lin1_cDNAsortedVals[names(Lin1_cDNAsortedVals) %in% names(lens[lens>11])]
	Lin1_pDNAsortedVals12<-Lin1_pDNAsortedVals[names(Lin1_pDNAsortedVals) %in% names(lens[lens>11])]
	Lin1_activitySortedVals12<-Lin1_cDNAsortedVals12-Lin1_pDNAsortedVals12
	sds<-double(length(unique(names(Lin1_activitySortedVals12))))
	names(sds)<-unique(names(Lin1_activitySortedVals12))
	sds<-tapply(Lin1_activitySortedVals12,names(Lin1_activitySortedVals12),sd)
	Lin1_activitySortedVals12sd<-Lin1_activitySortedVals12[names(Lin1_activitySortedVals12) %in% names(sds[sds!=0])]
	
	Lin2_cDNAsortedVals12<-Lin2_cDNAsortedVals[names(Lin2_cDNAsortedVals) %in% names(lens[lens>11])]
	Lin2_pDNAsortedVals12<-Lin2_pDNAsortedVals[names(Lin2_pDNAsortedVals) %in% names(lens[lens>11])]
	Lin2_activitySortedVals12<-Lin2_cDNAsortedVals12-Lin2_pDNAsortedVals12
	sds<-double(length(unique(names(Lin2_activitySortedVals12))))
	names(sds)<-unique(names(Lin2_activitySortedVals12))
	sds<-tapply(Lin2_activitySortedVals12,names(Lin2_activitySortedVals12),sd)
	Lin2_activitySortedVals12sd<-Lin2_activitySortedVals12[names(Lin2_activitySortedVals12) %in% names(sds[sds!=0])]
	
	# Create list of fragments present with >11 barcodes in both species (needed for differential activity testing) #
	chimpActSortedVals12<-Rep12_activitySortedVals12[grep("_Chimp",names(Rep12_activitySortedVals12))]
	names(chimpActSortedVals12)<-sub("_Chimp","",names(chimpActSortedVals12))
	humanActSortedVals12<-Rep12_activitySortedVals12[grep("_Chimp",names(Rep12_activitySortedVals12),invert=T)]
	nDummy<-c(unique(names(chimpActSortedVals12)),unique(names(humanActSortedVals12)))
	names2List<-nDummy[duplicated(nDummy)]
	save(names2List,file="data/fragnamelist.RData")
	rm(list=c("chimpActSortedVals12","humanActSortedVals12","nDummy"))

	### FRAGMENT MEDIAN ###
	Rep12_cDNAmedian12<-tapply(Rep12_cDNAsortedVals12,names(Rep12_cDNAsortedVals12),median)
	Rep13_cDNAmedian12<-tapply(Rep13_cDNAsortedVals12,names(Rep13_cDNAsortedVals12),median)
	Rep22_cDNAmedian12<-tapply(Rep22_cDNAsortedVals12,names(Rep22_cDNAsortedVals12),median)
	Rep23_cDNAmedian12<-tapply(Rep23_cDNAsortedVals12,names(Rep23_cDNAsortedVals12),median)
	Lin1_cDNAmedian12<-tapply(Lin1_cDNAsortedVals12,names(Lin1_cDNAsortedVals12),median)
	Lin2_cDNAmedian12<-tapply(Lin2_cDNAsortedVals12,names(Lin2_cDNAsortedVals12),median)
	
	Rep12_pDNAmedian12<-tapply(Rep12_pDNAsortedVals12,names(Rep12_pDNAsortedVals12),median)
	Rep13_pDNAmedian12<-tapply(Rep13_pDNAsortedVals12,names(Rep13_pDNAsortedVals12),median)
	Rep22_pDNAmedian12<-tapply(Rep22_pDNAsortedVals12,names(Rep22_pDNAsortedVals12),median)
	Rep23_pDNAmedian12<-tapply(Rep23_pDNAsortedVals12,names(Rep23_pDNAsortedVals12),median)
	Lin1_pDNAmedian12<-tapply(Lin1_pDNAsortedVals12,names(Lin1_pDNAsortedVals12),median)
	Lin2_pDNAmedian12<-tapply(Lin2_pDNAsortedVals12,names(Lin2_pDNAsortedVals12),median)
	
	Rep12_actMedian<-tapply(Rep12_cDNAsortedVals-Rep12_pDNAsortedVals,names(Rep12_pDNAsortedVals),median)
	Rep13_actMedian<-tapply(Rep13_cDNAsortedVals-Rep13_pDNAsortedVals,names(Rep13_pDNAsortedVals),median)
	Rep22_actMedian<-tapply(Rep22_cDNAsortedVals-Rep22_pDNAsortedVals,names(Rep22_pDNAsortedVals),median)
	Rep23_actMedian<-tapply(Rep23_cDNAsortedVals-Rep23_pDNAsortedVals,names(Rep23_pDNAsortedVals),median)
	Lin1_actMedian<-tapply(Lin1_cDNAsortedVals-Lin1_pDNAsortedVals,names(Lin1_pDNAsortedVals),median)
	Lin2_actMedian<-tapply(Lin2_cDNAsortedVals-Lin2_pDNAsortedVals,names(Lin2_pDNAsortedVals),median)
	
	### cDNA AND pDNA CORRELATIONS BETWEEN REPLICATES AFTER 12 TAG CUTOFF ###
	cor.test(Rep12_cDNAsortedVals12,Rep13_cDNAsortedVals12,method="spearman") # rho = 0.8302218, P = 0
	cor.test(Rep12_cDNAsortedVals12,Rep22_cDNAsortedVals12,method="spearman") # rho = 0.8520344, P = 0
	cor.test(Rep12_cDNAsortedVals12,Rep23_cDNAsortedVals12,method="spearman") # rho = 0.810885, P = 0
	cor.test(Rep13_cDNAsortedVals12,Rep22_cDNAsortedVals12,method="spearman") # rho = 0.8356911, P = 0
	cor.test(Rep13_cDNAsortedVals12,Rep23_cDNAsortedVals12,method="spearman") # rho = 0.8089392, P = 0
	cor.test(Rep22_cDNAsortedVals12,Rep23_cDNAsortedVals12,method="spearman") # rho = 0.8270934, P = 0
	
	cor.test(Rep12_pDNAsortedVals12,Rep13_pDNAsortedVals12,method="spearman") # rho = 0.883811, P = 0
	cor.test(Rep12_pDNAsortedVals12,Rep22_pDNAsortedVals12,method="spearman") # rho = 0.8945874, P = 0
	cor.test(Rep12_pDNAsortedVals12,Rep23_pDNAsortedVals12,method="spearman") # rho = 0.8820427, P = 0
	cor.test(Rep13_pDNAsortedVals12,Rep22_pDNAsortedVals12,method="spearman") # rho = 0.8747602, P = 0
	cor.test(Rep13_pDNAsortedVals12,Rep23_pDNAsortedVals12,method="spearman") # rho = 0.8772182, P = 0
	cor.test(Rep22_pDNAsortedVals12,Rep23_pDNAsortedVals12,method="spearman") # rho = 0.8726161, P = 0
	
	cor.test(Rep12_cDNAsortedVals12,Rep12_pDNAsortedVals12,method="spearman") # rho = 0.8264091, P = 0
	cor.test(Rep13_cDNAsortedVals12,Rep13_pDNAsortedVals12,method="spearman") # rho = 0.8105938, P = 0
	cor.test(Rep22_cDNAsortedVals12,Rep22_pDNAsortedVals12,method="spearman") # rho = 0.8426397, P = 0
	cor.test(Rep23_cDNAsortedVals12,Rep23_pDNAsortedVals12,method="spearman") # rho = 0.8134857, P = 0
	
	cor.test(Rep12_cDNAmedian12,Rep12_pDNAmedian12,method="spearman") # rho = 0.8031211, P = 0
	cor.test(Rep13_cDNAmedian12,Rep13_pDNAmedian12,method="spearman") # rho = 0.7875053, P = 0
	cor.test(Rep22_cDNAmedian12,Rep22_pDNAmedian12,method="spearman") # rho = 0.8067808, P = 0
	cor.test(Rep23_cDNAmedian12,Rep23_pDNAmedian12,method="spearman") # rho = 0.7964124, P = 0
	
	### ACTIVITY t TEST ###
	# Calculate My (activity distribution maximum) for t test #
	Rep12_ActMy<-density(Rep12_activitySortedVals12sd)$x[which.max(density(Rep12_activitySortedVals12)$y)]
	Rep13_ActMy<-density(Rep13_activitySortedVals12sd)$x[which.max(density(Rep13_activitySortedVals12)$y)]
	Rep22_ActMy<-density(Rep22_activitySortedVals12sd)$x[which.max(density(Rep22_activitySortedVals12)$y)]
	Rep23_ActMy<-density(Rep23_activitySortedVals12sd)$x[which.max(density(Rep23_activitySortedVals12)$y)]
	Lin1_ActMy<-density(Lin1_activitySortedVals12sd)$x[which.max(density(Lin1_activitySortedVals12)$y)]
	Lin2_ActMy<-density(Lin2_activitySortedVals12sd)$x[which.max(density(Lin2_activitySortedVals12)$y)]
	
	# Activity t test #
	Rep12_tcvpRes<-double(length(unique(names(Rep12_activitySortedVals12sd))))
	Rep12_tcvpRes<-tapply(Rep12_activitySortedVals12sd-Rep12_ActMy,names(Rep12_activitySortedVals12sd),t.test,alternative="greater")
	Rep12_tcvpP<-double(length(Rep12_tcvpRes))
	for(i in 1:length(Rep12_tcvpP)) Rep12_tcvpP[i]<-Rep12_tcvpRes[[i]]$p.value
	names(Rep12_tcvpP)<-names(Rep12_tcvpRes)
	
	Rep13_tcvpRes<-double(length(Rep13_activitySortedVals12sd))
	Rep13_tcvpRes<-tapply(Rep13_activitySortedVals12sd-Rep13_ActMy,names(Rep13_activitySortedVals12sd),t.test,alternative="greater")
	Rep13_tcvpP<-double(length(Rep13_tcvpRes))
	for(i in 1:length(Rep13_tcvpP)) Rep13_tcvpP[i]<-Rep13_tcvpRes[[i]]$p.value
	names(Rep13_tcvpP)<-names(Rep13_tcvpRes)
	
	Rep22_tcvpRes<-double(length(Rep22_activitySortedVals12sd))
	Rep22_tcvpRes<-tapply(Rep22_activitySortedVals12sd-Rep22_ActMy,names(Rep22_activitySortedVals12sd),t.test,alternative="greater")
	Rep22_tcvpP<-double(length(Rep22_tcvpRes))
	for(i in 1:length(Rep22_tcvpP)) Rep22_tcvpP[i]<-Rep22_tcvpRes[[i]]$p.value
	names(Rep22_tcvpP)<-names(Rep22_tcvpRes)
	
	Rep23_tcvpRes<-double(length(Rep23_activitySortedVals12sd))
	Rep23_tcvpRes<-tapply(Rep23_activitySortedVals12sd-Rep23_ActMy,names(Rep23_activitySortedVals12sd),t.test,alternative="greater")
	Rep23_tcvpP<-double(length(Rep23_tcvpRes))
	for(i in 1:length(Rep23_tcvpP)) Rep23_tcvpP[i]<-Rep23_tcvpRes[[i]]$p.value
	names(Rep23_tcvpP)<-names(Rep23_tcvpRes)
	
	Lin1_tcvpRes<-double(length(Lin1_activitySortedVals12sd))
	Lin1_tcvpRes<-tapply(Lin1_activitySortedVals12sd-Lin1_ActMy,names(Lin1_activitySortedVals12sd),t.test,alternative="greater")
	Lin1_tcvpP<-double(length(Lin1_tcvpRes))
	for(i in 1:length(Lin1_tcvpP)) Lin1_tcvpP[i]<-Lin1_tcvpRes[[i]]$p.value
	names(Lin1_tcvpP)<-names(Lin1_tcvpRes)
	
	Lin2_tcvpRes<-double(length(Lin2_activitySortedVals12sd))
	Lin2_tcvpRes<-tapply(Lin2_activitySortedVals12sd-Lin2_ActMy,names(Lin2_activitySortedVals12sd),t.test,alternative="greater")
	Lin2_tcvpP<-double(length(Lin2_tcvpRes))
	for(i in 1:length(Lin2_tcvpP)) Lin2_tcvpP[i]<-Lin2_tcvpRes[[i]]$p.value
	names(Lin2_tcvpP)<-names(Lin2_tcvpRes)
	
	### COLLECT FRAGMENT TABLE ###
	MPRAfrag1<-unique(MPRAsum1[,2:8])
	MPRAfrag1<-merge(MPRAfrag1,Rep12_cDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep12_cDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep13_cDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep13_cDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep22_cDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep22_cDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep23_cDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep23_cDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Lin1_cDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin1_cDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Lin2_cDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin2_cDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep12_pDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep12_pDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep13_pDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep13_pDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep22_pDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep22_pDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep23_pDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep23_pDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Lin1_pDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin1_pDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Lin2_pDNAmedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin2_pDNAmedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep12_actMedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep12_actMedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep13_actMedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep13_actMedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep22_actMedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep22_actMedian"
	MPRAfrag1<-merge(MPRAfrag1,Rep23_actMedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep23_actMedian"
	MPRAfrag1<-merge(MPRAfrag1,Lin1_actMedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin1_actMedian"
	MPRAfrag1<-merge(MPRAfrag1,Lin2_actMedian,by.x="Alignment",by.y=0)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin2_actMedian"
	
	MPRAfrag1<-merge(MPRAfrag1,Rep12_tcvpP,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep12_ttest_P"
	MPRAfrag1<-merge(MPRAfrag1,Rep13_tcvpP,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep13_ttest_P"
	MPRAfrag1<-merge(MPRAfrag1,Rep22_tcvpP,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep22_ttest_P"
	MPRAfrag1<-merge(MPRAfrag1,Rep23_tcvpP,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Rep23_ttest_P"
	MPRAfrag1<-merge(MPRAfrag1,Lin1_tcvpP,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin1_ttest_P"
	MPRAfrag1<-merge(MPRAfrag1,Lin2_tcvpP,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"Lin2_ttest_P"
	
	MPRAfrag1<-merge(MPRAfrag1,lens,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfrag1)[ncol(MPRAfrag1)]<-"noTags"

	rm(list=c("lens","Lin1_activitySortedVals12","Lin1_activitySortedVals12sd","Lin1_actMedian","Lin1_ActMy","Lin1_cDNAdata",
		"Lin1_cDNAmedian","Lin1_cDNAmedian12","Lin1_cDNAsortedVals","Lin1_cDNAsortedVals12","Lin1_pDNAdata","Lin1_pDNAmedian","Lin1_pDNAmedian12",
		"Lin1_pDNAsortedVals","Lin1_pDNAsortedVals12","Lin1_tcvpP","Lin1_tcvpRes","Lin2_activitySortedVals12","Lin2_activitySortedVals12sd",
		"Lin2_actMedian","Lin2_ActMy","Lin2_cDNAdata","Lin2_cDNAmedian","Lin2_cDNAmedian12","Lin2_cDNAsortedVals","Lin2_cDNAsortedVals12",
		"Lin2_pDNAdata","Lin2_pDNAmedian","Lin2_pDNAmedian12","Lin2_pDNAsortedVals","Lin2_pDNAsortedVals12","Lin2_tcvpP","Lin2_tcvpRes",
		"Rep12_activitySortedVals12","Rep12_activitySortedVals12sd","Rep12_actMedian","Rep12_ActMy","Rep12_cDNAdata","Rep12_cDNAmedian",
		"Rep12_cDNAmedian12","Rep12_cDNAsortedVals","Rep12_cDNAsortedVals12","Rep12_pDNAdata","Rep12_pDNAmedian","Rep12_pDNAmedian12",
		"Rep12_pDNAsortedVals","Rep12_pDNAsortedVals12","Rep12_tcvpP","Rep12_tcvpRes","Rep13_activitySortedVals12","Rep13_activitySortedVals12sd",
		"Rep13_actMedian","Rep13_ActMy","Rep13_cDNAdata","Rep13_cDNAmedian","Rep13_cDNAmedian12","Rep13_cDNAsortedVals","Rep13_cDNAsortedVals12",
		"Rep13_pDNAdata","Rep13_pDNAmedian","Rep13_pDNAmedian12","Rep13_pDNAsortedVals","Rep13_pDNAsortedVals12","Rep13_tcvpP","Rep13_tcvpRes",
		"Rep22_activitySortedVals12","Rep22_activitySortedVals12sd","Rep22_actMedian","Rep22_ActMy","Rep22_cDNAdata","Rep22_cDNAmedian",
		"Rep22_cDNAmedian12","Rep22_cDNAsortedVals","Rep22_cDNAsortedVals12","Rep22_pDNAdata","Rep22_pDNAmedian","Rep22_pDNAmedian12",
		"Rep22_pDNAsortedVals","Rep22_pDNAsortedVals12","Rep22_tcvpP","Rep22_tcvpRes","Rep23_activitySortedVals12","Rep23_activitySortedVals12sd",
		"Rep23_actMedian","Rep23_ActMy","Rep23_cDNAdata","Rep23_cDNAmedian","Rep23_cDNAmedian12","Rep23_cDNAsortedVals","Rep23_cDNAsortedVals12",
		"Rep23_pDNAdata","Rep23_pDNAmedian","Rep23_pDNAmedian12","Rep23_pDNAsortedVals","Rep23_pDNAsortedVals12","Rep23_tcvpP","Rep23_tcvpRes","sds"))
	write.table(MPRAfrag1,"data/MPRAfrag1.tsv",quote=F,row.names=F,sep="\t")
}

### FRAGMENT ACTIVITY CORRELATIONS BETWEEN REPLICATES ###
cor.test(MPRAfrag1$Rep12_actMedian[!is.na(MPRAfrag1$Rep12_ttest_P)],MPRAfrag1$Rep13_actMedian[!is.na(MPRAfrag1$Rep13_ttest_P)]) # rho = 0.7648721, P = 0
cor.test(MPRAfrag1$Rep12_actMedian[!is.na(MPRAfrag1$Rep12_ttest_P)],MPRAfrag1$Rep22_actMedian[!is.na(MPRAfrag1$Rep22_ttest_P)]) # rho = 0.7767904, P = 0
cor.test(MPRAfrag1$Rep12_actMedian[!is.na(MPRAfrag1$Rep12_ttest_P)],MPRAfrag1$Rep23_actMedian[!is.na(MPRAfrag1$Rep23_ttest_P)]) # rho = 0.7375057, P = 0
cor.test(MPRAfrag1$Rep13_actMedian[!is.na(MPRAfrag1$Rep13_ttest_P)],MPRAfrag1$Rep22_actMedian[!is.na(MPRAfrag1$Rep22_ttest_P)]) # rho = 0.7667899, P = 0
cor.test(MPRAfrag1$Rep13_actMedian[!is.na(MPRAfrag1$Rep13_ttest_P)],MPRAfrag1$Rep23_actMedian[!is.na(MPRAfrag1$Rep23_ttest_P)]) # rho = 0.7340809, P = 0
cor.test(MPRAfrag1$Rep22_actMedian[!is.na(MPRAfrag1$Rep22_ttest_P)],MPRAfrag1$Rep23_actMedian[!is.na(MPRAfrag1$Rep23_ttest_P)]) # rho = 0.7523065, P = 0

nrow(MPRAfrag1) # 97116 total recovered fragments
nrow(MPRAfrag1[MPRAfrag1$noTags>11,]) # 78487 fragments with >11 barcode tags

### PERMISSIVE COLLECTION OF ACTIVE FRAGMENTS FOR ROUND TWO ###
if(file.exists("data/MPRAactive1_1.tsv")){ # If file available, load; else create
	MPRAactive1_1<-read.delim("data/MPRAactive1_1.tsv",stringsAsFactors=F)
}else{
# Adjusted P value <0.1 and active in at least two replicates OR both lineage OR a Mann-Whitney U test between lineages
	LinActive<-MPRAfrag1[p.adjust(MPRAfrag1$Lin1_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Lin2_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Lin1_ttest_P),]
	LinActive<-cbind(LinActive,LinActive=rep("Lin1.2",nrow(LinActive)))
	LinActive$LinActive<-as.character(LinActive$LinActive)
	
	allActive<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep12<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep13<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")>.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep22<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep23<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")>.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep12.13<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")>.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep12.22<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep12.23<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")>.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep13.22<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")>.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep13.23<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")>.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")>.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	no_Rep22.23<-MPRAfrag1[p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.1 & p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.1 
		& p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")>.1 & p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")>.1 & !is.na(MPRAfrag1$Rep12_ttest_P),]
	
	# Mann-Whitney U test results #
	MWU_active<-read.delim("data/MWU_active_fragments.lst",header=F,stringsAsFactors=F)
	MWU_active<-as.vector(t(MWU_active))
	MWU_active<-MPRAfrag1[MPRAfrag1$Alignment %in% MWU_active,]
	MWU_active<-cbind(MWU_active,MWUactive=rep("MWU",nrow(MWU_active)))
	MWU_active$WilcActive<-as.character(MWU_active$MWUactive)
	MWU_active<-MWU_active[MWU_active$noTags>11,]
	
	### COLLECT ALL ACTIVE FRAGMENTS ###
	MPRAactive1_1<-rbind(cbind(allActive,active=rep("allReps",nrow(allActive))),cbind(no_Rep12,active=rep("Rep13.22.23",nrow(no_Rep12))),
		cbind(no_Rep13,active=rep("Rep12.22.23",nrow(no_Rep13))),cbind(no_Rep22,active=rep("Rep12.13.23",nrow(no_Rep22))),
		cbind(no_Rep23,active=rep("Rep12.13.22",nrow(no_Rep23))),cbind(no_Rep12.13,active=rep("Rep22.23",nrow(no_Rep12.13))),
		cbind(no_Rep12.22,active=rep("Rep13.23",nrow(no_Rep12.22))),cbind(no_Rep12.23,active=rep("Rep13.22",nrow(no_Rep12.23))),
		cbind(no_Rep13.22,active=rep("Rep12.23",nrow(no_Rep13.22))),cbind(no_Rep13.23,active=rep("Rep12.22",nrow(no_Rep13.23))),
		cbind(no_Rep22.23,active=rep("Rep12.13",nrow(no_Rep22.23))))
	MPRAactive1_1$active<-as.character(MPRAactive1_1$active)
	MPRAactive1_1<-merge(MPRAactive1_1,LinActive,all=T)
	MPRAactive1_1<-merge(MPRAactive1_1,MWU_active,all=T)
	actDummy<-double(nrow(MPRAactive1_1))
	for(i in 1:nrow(MPRAactive1_1)){
		actDummy[i]<-paste5(MPRAactive1_1[i,33],MPRAactive1_1[i,34],MPRAactive1_1[i,35],sep=".",na.rm=T)
	}
	MPRAactive1_1<-cbind(MPRAactive1_1[,1:32],active=actDummy)
	MPRAactive1_1$active<-as.character(MPRAactive1_1$active)
	write.table(MPRAactive1_1,"data/MPRAactive1_1.tsv",row.names=F,sep="\t",quote=F)
	MPRAdataDummy<-unique(MPRAdata1[,c(2,4:6)])
	bedAct<-cbind(MPRAactive1_1[,c(1:4,6,7)])
	for(i in 1:nrow(bedAct)){
		if(grepl("Chimp",bedAct$Alignment[i])==T){
			bedAct$Chr[i]<-MPRAdataDummy$Chr[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
			bedAct$Start[i]<-MPRAdataDummy$Start[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
			bedAct$Stop[i]<-MPRAdataDummy$Stop[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
	}}
	bedAct<-unique(bedAct[,-1])
	bedAct$Start<-bedAct$Start-1
	write.table(bedAct,"data/MPRAactive1_1.bed",row.names=F,col.names=F,quote=F,sep="\t")
	rm(list=c("actDummy","allActive","bedAct","i","LinActive","MPRAdataDummy","MWU_active","no_Rep12","no_Rep13","no_Rep22","no_Rep23","no_Rep12.13",
		"no_Rep12.22","no_Rep12.23","no_Rep13.22","no_Rep13.23","no_Rep22.23"))
}

### COMPARE ORTHOLOGOUS HUMAN AND CHIMPANZEE FRAGMENTS: FOLDED FRAGMENT TABLE ###
if(file.exists("data/MPRAfold1.tsv")){ # If available load; else write
	MPRAfold1<-read.delim("data/MPRAfold1.tsv",stringsAsFactors=F)
}else{
	MPRAfold1<-MPRAfrag1
	MPRAfold1$Alignment<-sub('_Chimp','',MPRAfrag1$Alignment)
	MPRAfold1<-merge(MPRAfold1[MPRAfold1$Species=='hg19',-5],by="Alignment",MPRAfold1[MPRAfold1$Species=='PanTro2',c(1,8:ncol(MPRAfold1))])
	colnames(MPRAfold1)<-sub('.x','_hs',colnames(MPRAfold1),fixed=T)
	colnames(MPRAfold1)<-sub('.y','_pt',colnames(MPRAfold1),fixed=T)
	
	MPRAfoldActive1<-MPRAfold1[MPRAfold1$Alignment %in% unique(sub('_Chimp','',MPRAactive1_1$Alignment)),]
	MPRAfoldActive1<-cbind(MPRAfoldActive1,active_hs=double(nrow(MPRAfoldActive1)),active_pt=double(nrow(MPRAfoldActive1)))
	for(i in 1:nrow(MPRAactive1_1)){
		alig<-sub('_Chimp','',MPRAactive1_1$Alignment[i])
		if(grepl('_Chimp',MPRAactive1_1[i,1])){
			MPRAfoldActive1[MPRAfoldActive1$Alignment==alig,58]<-MPRAactive1_1[i,33]
		}else{
			MPRAfoldActive1[MPRAfoldActive1$Alignment==alig,57]<-MPRAactive1_1[i,33]
	}}
	MPRAfoldActive1[,57:58]<-replace(MPRAfoldActive1[,57:58],MPRAfoldActive1[,57:58]==0,NA)
	
	### PREPARE DATA FOR DIFFERENTIAL ACTIVITY t TEST ###
	# Rename "_Chimp" fragments to same name as human fragments #
	# NOTE: Species can be identified by "Species" column
	MPRAactivities1<-cbind(MPRAsum1[,1:8],
		Rep12_activity=MPRAsum1$Rep_1_2_cDNA_Tag_Counts-MPRAsum1$Rep_1_2_pDNA_Tag_Counts,
		Rep13_activity=MPRAsum1$Rep_1_3_cDNA_Tag_Counts-MPRAsum1$Rep_1_3_pDNA_Tag_Counts,
		Rep22_activity=MPRAsum1$Rep_2_2_cDNA_Tag_Counts-MPRAsum1$Rep_2_2_pDNA_Tag_Counts,
		Rep23_activity=MPRAsum1$Rep_2_3_cDNA_Tag_Counts-MPRAsum1$Rep_2_3_pDNA_Tag_Counts,
		Lin1_activity=MPRAsum1$Lin1_cDNA_Tag_Counts-MPRAsum1$Lin1_pDNA_Tag_Counts,
		Lin2_activity=MPRAsum1$Lin2_cDNA_Tag_Counts-MPRAsum1$Lin2_pDNA_Tag_Counts)
	MPRAactivities1$Alignment<-sub('_Chimp','',MPRAactivities1$Alignment)
	
	# Create lists by fragment name and species #
	Rep12_nestAct<-MPRAactivities1[MPRAactivities1$Alignment %in% names2List 
		& MPRAactivities1$Alignment %in% MPRAfoldActive1$Alignment,c(2,6,9)] %>% nest(-Alignment, Species)
	Rep13_nestAct<-MPRAactivities1[MPRAactivities1$Alignment %in% names2List 
		& MPRAactivities1$Alignment %in% MPRAfoldActive1$Alignment,c(2,6,10)] %>% nest(-Alignment, Species)
	Rep22_nestAct<-MPRAactivities1[MPRAactivities1$Alignment %in% names2List 
		& MPRAactivities1$Alignment %in% MPRAfoldActive1$Alignment,c(2,6,11)] %>% nest(-Alignment, Species)
	Rep23_nestAct<-MPRAactivities1[MPRAactivities1$Alignment %in% names2List 
		& MPRAactivities1$Alignment %in% MPRAfoldActive1$Alignment,c(2,6,12)] %>% nest(-Alignment, Species)
	Lin1_nestAct<-MPRAactivities1[MPRAactivities1$Alignment %in% names2List 
		& MPRAactivities1$Alignment %in% MPRAfoldActive1$Alignment,c(2,6,13)] %>% nest(-Alignment, Species)
	Lin2_nestAct<-MPRAactivities1[MPRAactivities1$Alignment %in% names2List 
		& MPRAactivities1$Alignment %in% MPRAfoldActive1$Alignment,c(2,6,14)] %>% nest(-Alignment, Species)
	# Assign empty output lists #
	Rep12_tDAActRes<-list(double(length(Rep12_nestAct$Alignment)))
	Rep13_tDAActRes<-list(double(length(Rep12_nestAct$Alignment)))
	Rep22_tDAActRes<-list(double(length(Rep12_nestAct$Alignment)))
	Rep23_tDAActRes<-list(double(length(Rep12_nestAct$Alignment)))
	Lin1_tDAActRes<-list(double(length(Rep12_nestAct$Alignment)))
	Lin2_tDAActRes<-list(double(length(Rep12_nestAct$Alignment)))
	Rep12_tDAActp<-double(length(Rep12_nestAct$Alignment))
	Rep13_tDAActp<-double(length(Rep12_nestAct$Alignment))
	Rep22_tDAActp<-double(length(Rep12_nestAct$Alignment))
	Rep23_tDAActp<-double(length(Rep12_nestAct$Alignment))
	Lin1_tDAActp<-double(length(Rep12_nestAct$Alignment))
	Lin2_tDAActp<-double(length(Rep12_nestAct$Alignment))
	Lin1_tDAActt<-double(length(Rep12_nestAct$Alignment))
	Lin2_tDAActt<-double(length(Rep12_nestAct$Alignment))
	
	# Differential activity t test #
	for(i in 1:length(Rep12_nestAct$Alignment)){
		if(i%%100==0){
			timestamp()
			print(i)
		}
		Rep12_tDAActRes[[i]]<-t.test(Rep12_activity ~ Species,data=Rep12_nestAct$data[[i]])
		Rep13_tDAActRes[[i]]<-t.test(Rep13_activity ~ Species,data=Rep13_nestAct$data[[i]])
		Rep22_tDAActRes[[i]]<-t.test(Rep22_activity ~ Species,data=Rep22_nestAct$data[[i]])
		Rep23_tDAActRes[[i]]<-t.test(Rep23_activity ~ Species,data=Rep23_nestAct$data[[i]])
		Lin1_tDAActRes[[i]]<-t.test(Lin1_activity ~ Species,data=Lin1_nestAct$data[[i]])
		Lin2_tDAActRes[[i]]<-t.test(Lin2_activity ~ Species,data=Lin2_nestAct$data[[i]])
		Rep12_tDAActp[i]<-Rep12_tDAActRes[[i]]$p.value
		Rep13_tDAActp[i]<-Rep13_tDAActRes[[i]]$p.value
		Rep22_tDAActp[i]<-Rep22_tDAActRes[[i]]$p.value
		Rep23_tDAActp[i]<-Rep23_tDAActRes[[i]]$p.value
		Lin1_tDAActp[i]<-Lin1_tDAActRes[[i]]$p.value
		Lin2_tDAActp[i]<-Lin2_tDAActRes[[i]]$p.value
		Lin1_tDAActt[i]<-Lin1_tDAActRes[[i]]$statistic
		Lin2_tDAActt[i]<-Lin2_tDAActRes[[i]]$statistic
	}
	# Fragment names #
	names(Rep12_tDAActp)<-Rep12_nestAct$Alignment
	names(Rep13_tDAActp)<-Rep12_nestAct$Alignment
	names(Rep22_tDAActp)<-Rep12_nestAct$Alignment
	names(Rep23_tDAActp)<-Rep12_nestAct$Alignment
	names(Lin1_tDAActp)<-Rep12_nestAct$Alignment
	names(Lin2_tDAActp)<-Rep12_nestAct$Alignment
	names(Lin1_tDAActt)<-Rep12_nestAct$Alignment
	names(Lin2_tDAActt)<-Rep12_nestAct$Alignment
	
	# Add test results to folded fragment table #
	MPRAfold1<-merge(MPRAfold1,Rep12_tDAActp,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Rep12_ttest_DA_P"
	MPRAfold1<-merge(MPRAfold1,Rep13_tDAActp,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Rep13_ttest_DA_P"
	MPRAfold1<-merge(MPRAfold1,Rep22_tDAActp,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Rep22_ttest_DA_P"
	MPRAfold1<-merge(MPRAfold1,Rep23_tDAActp,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Rep23_ttest_DA_P"
	MPRAfold1<-merge(MPRAfold1,Lin1_tDAActp,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Lin1_ttest_DA_P"
	MPRAfold1<-merge(MPRAfold1,Lin2_tDAActp,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Lin2_ttest_DA_P"
	MPRAfold1<-merge(MPRAfold1,Lin1_tDAActt,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Lin1_ttest_DA_t"
	MPRAfold1<-merge(MPRAfold1,Lin2_tDAActt,by.x="Alignment",by.y=0,all=T)
	colnames(MPRAfold1)[ncol(MPRAfold1)]<-"Lin2_ttest_DA_t"

	# Number of fragments tested in differential activity test #
	length(Rep12_nestAct$Alignment) # 3219
	
	rm(list=c("alig","i","Lin1_nestAct","Lin1_tDAActp","Lin1_tDAActRes","Lin1_tDAActt","Lin2_nestAct","Lin2_tDAActp","Lin2_tDAActRes",
		"Lin2_tDAActt","Rep12_nestAct","Rep12_tDAActp","Rep12_tDAActRes","Rep13_nestAct","Rep13_tDAActp","Rep13_tDAActRes","Rep22_nestAct",
		"Rep22_tDAActp","Rep22_tDAActRes","Rep23_nestAct","Rep23_tDAActp","Rep23_tDAActRes"))
	write.table(MPRAfold1,"data/MPRAfold1.tsv",quote=F,sep="\t",row.names=F)
}

### COLLECT DIFFERENTIALLY ACTIVE FRAGMENT TABLE ###
# As above, permissive screen to retrieve candidates: 
# adjusted P value <0.1 for two or more replicates OR both lineage OR lineage Mann-Whitney U test
if(file.exists("data/diffActive1_1.tsv")){ # If available load; else create
	MPRAda1_1<-read.delim("data/diffActive1_1.tsv",stringsAsFactors=F)
}else{
	bothDA<-MPRAfold1[p.adjust(MPRAfold1$Lin1_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Lin2_ttest_DA_P,method="BH")<.1 
		& !is.na(MPRAfold1$Lin1_ttest_DA_P),]
	bothDA<-cbind(bothDA,Lins=rep("Lin1.2",nrow(bothDA)))
	bothDA$Lins<-as.character(bothDA$Lins)
	allDA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep12_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep13_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")>.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep22_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep23_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")>.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep12.13_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")>.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep12.22_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep12.23_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")>.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep13.22_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")>.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep13.23_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")>.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")>.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	no_Rep22.23_DA<-MPRAfold1[p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.1 & p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.1 
		& p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")>.1 & p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")>.1 & !is.na(MPRAfold1$Rep12_ttest_DA_P),]
	# Mann-Whitney U test results #
	MWU_da<-read.delim("data/MWU_da_fragments.lst",stringsAsFactors=F,header=F)
	MWU_da<-as.vector(t(MWU_da))
	MWU_da<-MPRAfold1[MPRAfold1$Alignment %in% MWU_da,]
	MWU_da<-cbind(MWU_da,MWUda=rep("MWU",nrow(MWU_da)))
	MWU_da$MWUda<-as.character(MWU_da$MWUda)
	# Collect all differentially active fragments in one table #
	MPRAda1_1<-rbind(cbind(allDA,differential=rep("allReps",nrow(allDA))),cbind(no_Rep12_DA,differential=rep("Rep13.22.23",nrow(no_Rep12_DA))),
		cbind(no_Rep13_DA,differential=rep("Rep12.22.23",nrow(no_Rep13_DA))),cbind(no_Rep22_DA,differential=rep("Rep12.13.23",nrow(no_Rep22_DA))),
		cbind(no_Rep23_DA,differential=rep("Rep12.13.22",nrow(no_Rep23_DA))),cbind(no_Rep12.13_DA,differential=rep("Rep22.23",nrow(no_Rep12.13_DA))),
		cbind(no_Rep12.22_DA,differential=rep("Rep13.23",nrow(no_Rep12.22_DA))),cbind(no_Rep12.23_DA,differential=rep("Rep13.22",nrow(no_Rep12.23_DA))),
		cbind(no_Rep13.22_DA,differential=rep("Rep12.23",nrow(no_Rep13.22_DA))),cbind(no_Rep13.23_DA,differential=rep("Rep12.22",nrow(no_Rep13.23_DA))),
		cbind(no_Rep22.23_DA,differential=rep("Rep12.13",nrow(no_Rep22.23_DA))))
	MPRAda1_1$differential<-as.character(MPRAda1_1$differential)
	MPRAda1_1<-merge(x=MPRAda1_1,y=bothDA,all=T)
	MPRAda1_1<-merge(MPRAda1_1,MWU_da,all=T)
	daDummy<-double(nrow(MPRAda1_1))
	for(i in 1:nrow(MPRAda1_1)){
		daDummy[i]<-paste5(MPRAda1_1$differential[i],MPRAda1_1$Lins[i],MPRAda1_1$MWUda[i],sep=".",na.rm=T)
	}
	MPRAda1_1<-cbind(MPRAda1_1[,1:64],differential=daDummy)
	MPRAda1_1$differential<-as.character(MPRAda1_1$differential)
	write.table(MPRAda1_1,"data/diffActive1_1.tsv",sep="\t",quote=F,row.names=F)
	rm(list=c("allDA","bothDA","daDummy","i","MWU_da","no_Rep12_DA","no_Rep13_DA","no_Rep22_DA","no_Rep23_DA","no_Rep12.13_DA","no_Rep12.22_DA","no_Rep12.23_DA","no_Rep13.22_DA",
		"no_Rep13.23_DA","no_Rep22.23_DA"))
}
MPRAda1_1<-cbind(MPRAda1_1,diff=apply(MPRAda1_1[,19:22],1,mean)-apply(MPRAda1_1[,44:47],1,mean))

### MORE CONSERVATIVE FRAGMENT SETS FOR REPORTED STATISTICS ###
# Active fragments have adjusted P value <0.05 and cDNA > pDNA in all replicates
if(file.exists("data/MPRAactive1_2.tsv")){ # If available load; else create
	MPRAactive1_2<-read.delim("data/MPRAactive1_2.tsv",stringsAsFactors=F)
}else{
	a<-nrow(MPRAfrag1)
	MPRAactiveDummy<-data.frame(Rep12=double(a),Rep13=double(a),Rep22=double(a),Rep23=double(a))
	MPRAactiveDummy$Rep12<-replace(MPRAactiveDummy$Rep12,p.adjust(MPRAfrag1$Rep12_ttest_P,method="BH")<.05 & !is.na(MPRAfrag1$Rep12_ttest_P),T)
	MPRAactiveDummy$Rep13<-replace(MPRAactiveDummy$Rep13,p.adjust(MPRAfrag1$Rep13_ttest_P,method="BH")<.05 & !is.na(MPRAfrag1$Rep13_ttest_P),T)
	MPRAactiveDummy$Rep22<-replace(MPRAactiveDummy$Rep22,p.adjust(MPRAfrag1$Rep22_ttest_P,method="BH")<.05 & !is.na(MPRAfrag1$Rep22_ttest_P),T)
	MPRAactiveDummy$Rep23<-replace(MPRAactiveDummy$Rep23,p.adjust(MPRAfrag1$Rep23_ttest_P,method="BH")<.05 & !is.na(MPRAfrag1$Rep23_ttest_P),T)
	MPRAactive1_2<-MPRAfrag1[rowSums(MPRAactiveDummy)>1,]
	lact<-MPRAfrag1[MPRAfrag1$Rep12_cDNAmedian > MPRAfrag1$Rep12_pDNAmedian & MPRAfrag1$Rep13_cDNAmedian > MPRAfrag1$Rep13_pDNAmedian 
		& MPRAfrag1$Rep22_cDNAmedian > MPRAfrag1$Rep22_pDNAmedian & MPRAfrag1$Rep23_cDNAmedian > MPRAfrag1$Rep23_pDNAmedian 
		& !is.na(MPRAfrag1$Rep23_ttest_P) & MPRAfrag1$Alignment %in% MPRAactive1_1$Alignment,]
	MPRAactive1_2<-MPRAactive1_2[MPRAactive1_2$Alignment %in% lact$Alignment,]
	write.table(MPRAactive1_2,"data/MPRAactive1_2.tsv",quote=F,row.names=F,sep="\t")
	MPRAdataDummy<-unique(MPRAdata1[,c(2,4:6)])
	bedAct<-cbind(MPRAactive1_2[,c(1:4,6,7)])
	for(i in 1:nrow(bedAct)){
		if(grepl("Chimp",bedAct$Alignment[i])==T){
			bedAct$Chr[i]<-MPRAdataDummy$Chr[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
			bedAct$Start[i]<-MPRAdataDummy$Start[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
			bedAct$Stop[i]<-MPRAdataDummy$Stop[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
	}}
	bedAct<-unique(bedAct[,-1])
	bedAct$Start<-bedAct$Start-1
	write.table(bedAct,"data/MPRAactive1_2.bed",row.names=F,col.names=F,quote=F,sep="\t")
	rm(list=c("a","bedAct","i","lact","MPRAactiveDummy"))
}

# Differentially active fragments agree in direction and have an average log2 fold change >0.2
if(file.exists("data/diffActive1_2.tsv")){ # If available load; else create
	MPRAda1_2<-read.delim("data/diffActive1_2.tsv",stringsAsFactors=F)
}else{
	a<-nrow(MPRAfold1)
	MPRAdaDummy<-data.frame(Rep12=double(a),Rep13=double(a),Rep22=double(a),Rep23=double(a))
	MPRAdaDummy$Rep12<-replace(MPRAdaDummy$Rep12,p.adjust(MPRAfold1$Rep12_ttest_DA_P,method="BH")<.05 & !is.na(MPRAfold1$Rep12_ttest_DA_P),T)
	MPRAdaDummy$Rep13<-replace(MPRAdaDummy$Rep13,p.adjust(MPRAfold1$Rep13_ttest_DA_P,method="BH")<.05 & !is.na(MPRAfold1$Rep13_ttest_DA_P),T)
	MPRAdaDummy$Rep22<-replace(MPRAdaDummy$Rep22,p.adjust(MPRAfold1$Rep22_ttest_DA_P,method="BH")<.05 & !is.na(MPRAfold1$Rep22_ttest_DA_P),T)
	MPRAdaDummy$Rep23<-replace(MPRAdaDummy$Rep23,p.adjust(MPRAfold1$Rep23_ttest_DA_P,method="BH")<.05 & !is.na(MPRAfold1$Rep23_ttest_DA_P),T)
	MPRAda1_2<-MPRAfold1[rowSums(MPRAdaDummy)>1,]

	MPRAbias1<-cbind(MPRAda1_1[,1:6],
		Rep12_actDiff=MPRAda1_1$Rep12_actMedian_hs-MPRAda1_1$Rep12_actMedian_pt,Rep13_actDiff=MPRAda1_1$Rep13_actMedian_hs-MPRAda1_1$Rep13_actMedian_pt,
		Rep22_actDiff=MPRAda1_1$Rep22_actMedian_hs-MPRAda1_1$Rep22_actMedian_pt,Rep23_actDiff=MPRAda1_1$Rep23_actMedian_hs-MPRAda1_1$Rep23_actMedian_pt,
		Lin1_actDiff=MPRAda1_1$Lin1_actMedian_hs-MPRAda1_1$Lin1_actMedian_pt,Lin2_actDiff=MPRAda1_1$Lin2_actMedian_hs-MPRAda1_1$Lin2_actMedian_pt,
		MPRAda1_1[,57:65],direction=rep(NA,nrow(MPRAda1_1)))
	rownames(MPRAbias1)<-NULL
	MPRAbias1$direction<-replace(MPRAbias1$direction,MPRAbias1$Rep12_actDiff >0 & MPRAbias1$Rep13_actDiff >0 & MPRAbias1$Rep22_actDiff >0 
			& MPRAbias1$Rep23_actDiff >0 & MPRAbias1$Lin1_actDiff >0 & MPRAbias1$Lin2_actDiff >0,'hs')
	MPRAbias1$direction<-replace(MPRAbias1$direction,MPRAbias1$Rep12_actDiff <0 & MPRAbias1$Rep13_actDiff <0 & MPRAbias1$Rep22_actDiff <0 
		& MPRAbias1$Rep23_actDiff <0 & MPRAbias1$Lin1_actDiff <0 & MPRAbias1$Lin2_actDiff <0,'pt')
	lfc.lst<-MPRAbias1$Alignment[MPRAbias1$Rep12_actDiff>0 & MPRAbias1$Rep13_actDiff>0 & MPRAbias1$Rep22_actDiff>0 & MPRAbias1$Rep23_actDiff>0 
			& apply(MPRAbias1[,7:10],1,mean) >.2]
	lfc.lst<-c(lfc.lst,MPRAbias1$Alignment[MPRAbias1$Rep12_actDiff<0 & MPRAbias1$Rep13_actDiff<0 & MPRAbias1$Rep22_actDiff<0 & MPRAbias1$Rep23_actDiff<0 
			& apply(MPRAbias1[,7:10],1,mean) < -.2])
	MPRAda1_2<-MPRAda1_2[MPRAda1_2$Alignment %in% lfc.lst,]
	write.table(MPRAda1_2,"data/diffActive1_2.tsv",quote=F,row.names=F,sep="\t")
	bedtable<-cbind(MPRAda1_2[,c(2:6)],diff=apply(MPRAda1_2[,23:24],1,mean)-apply(MPRAda1_2[,48:49],1,mean))
	bedtable$Start<-bedtable$Start-1
	write.table(bedtable,"data/diffActive1_2.bed",sep="\t",quote=F,row.names=F,col.names=F)
	rm(list=c("a","bedtable","lfc.lst","MPRAbias1","MPRAdaDummy"))
}
MPRAda1_2<-cbind(MPRAda1_2,diff=apply(MPRAda1_2[,19:22],1,mean)-apply(MPRAda1_2[,44:47],1,mean))

# Histograms to visualize log fold change cutoff choice
hist(apply(MPRAfrag1[MPRAfrag1$noTags>11,20:23],1,mean),
	breaks=seq(-.9,5.5,.05),las=1,ylab="# MPRA fragments",xlab=expression(paste(log[2]," activity")),main="")
hist(apply(MPRAfrag1[MPRAfrag1$noTags>11 & MPRAfrag1$Alignment %in% MPRAactive1_1$Alignment,20:23],1,mean),breaks=seq(-.9,5.5,.05),add=T,col=grey(.6))
hist(apply(MPRAfrag1[MPRAfrag1$noTags>11 & MPRAfrag1$Alignment %in% MPRAactive1_1$Alignment & MPRAfrag1$Alignment %out% MPRAactive1_2$Alignment,20:23],1,mean),
	breaks=seq(-.9,5.5,.05),add=T,col=2)

hist(abs(apply(MPRAfold1[!is.na(MPRAfold1$Rep23_ttest_DA_P),19:22],1,mean)-apply(MPRAfold1[!is.na(MPRAfold1$Rep12_ttest_DA_P),44:47],1,mean)),
	breaks=seq(0,3.8,.05),las=1,ylab="# MPRA fragment pairs",xlab="Human - chimpanzee activity",main="")
hist(abs(apply(MPRAfold1[!is.na(MPRAfold1$Rep23_ttest_DA_P) & MPRAfold1$Alignment %in% MPRAda1_1$Alignment,19:22],1,mean)
		-apply(MPRAfold1[!is.na(MPRAfold1$Rep12_ttest_DA_P) & MPRAfold1$Alignment %in% MPRAda1_1$Alignment,44:47],1,mean)),breaks=seq(0,3.8,.05),add=T,col=grey(.6))
hist(abs(apply(MPRAfold1[!is.na(MPRAfold1$Rep23_ttest_DA_P) & MPRAfold1$Alignment %in% MPRAda1_1$Alignment & MPRAfold1$Alignment %out% MPRAda1_2$Alignment,19:22],1,mean)
		-apply(MPRAfold1[!is.na(MPRAfold1$Rep12_ttest_DA_P) & MPRAfold1$Alignment %in% MPRAda1_1$Alignment & MPRAfold1$Alignment %out% MPRAda1_2$Alignment,44:47],1,mean)),
	breaks=seq(0,3.8,.05),add=T,col=2)

### FRAGMENT SUBSET BED FILES ###
MPRAdataDummy<-unique(MPRAdata1[,c(2,4:6)])
bedAct<-cbind(MPRAactive1_1[,c(1:4,6,7)])
chimpBed<-bedAct[grepl("_Chimp",bedAct$Alignment),]
chimpBed$Start<-chimpBed$Start-1
write.table(chimpBed[,-1],"data/MPRAactive1_1_Chimp.bed",sep="\t",row.names=F,col.names=F,quote=F) # Chimp-active fragments
humanBed<-bedAct[grepl("_Chimp",bedAct$Alignment)==F,]
humanBed$Start<-humanBed$Start-1
write.table(humanBed[,-1],"data/MPRAactive1_1_Human.bed",sep="\t",row.names=F,col.names=F,quote=F) # Human-active fragments
chimpSpec<-chimpBed[sub("_Chimp","",chimpBed$Alignment) %out% humanBed$Alignment,]
write.table(chimpSpec[,-1],"data/MPRAactive1_1_ChimpSpec.bed",sep="\t",row.names=F,col.names=F,quote=F) # Fragments with chimp-specific activity in panTro2
humanSpec<-humanBed[humanBed$Alignment %out% sub("_Chimp","",chimpBed$Alignment),]
write.table(humanSpec[,-1],"data/MPRAactive1_1_HumanSpec.bed",sep="\t",row.names=F,col.names=F,quote=F) # Fragments with human-specific activity
for(i in 1:nrow(chimpBed)){
	chimpBed$Chr[i]<-MPRAdataDummy$Chr[MPRAdataDummy$Alignment==sub("_Chimp",'',chimpBed$Alignment[i])]
	chimpBed$Start[i]<-MPRAdataDummy$Start[MPRAdataDummy$Alignment==sub("_Chimp",'',chimpBed$Alignment[i])]
	chimpBed$Stop[i]<-MPRAdataDummy$Stop[MPRAdataDummy$Alignment==sub("_Chimp",'',chimpBed$Alignment[i])]
}
chimpBed$Start<-chimpBed$Start-1
write.table(chimpBed[,-1],"data/MPRAactive1_1_panTro2tohg19.bed",quote=F,row.names=F,col.names=F,sep="\t") # Chimp-active fragments in hg19
for(i in 1:nrow(bedAct)){
	if(grepl("Chimp",bedAct$Alignment[i])==T){
		bedAct$Chr[i]<-MPRAdataDummy$Chr[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
		bedAct$Start[i]<-MPRAdataDummy$Start[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
		bedAct$Stop[i]<-MPRAdataDummy$Stop[MPRAdataDummy$Alignment==sub("_Chimp",'',bedAct$Alignment[i])]
}}
bedBoth<-cbind(bedAct[,2:4],Alignment=bedAct$Alignment,stringsAsFactors=F)
bedBoth$Start<-bedBoth$Start-1
write.table(bedBoth,"data/MPRAactive1_1_both2hg19.bed",quote=F,row.names=F,col.names=F,sep="\t") # Active fragments in hg19 (present 2x if active in both)
bedAct<-unique(bedAct[,-1])
bedAct$Start<-bedAct$Start-1
write.table(bedAct,"data/MPRAactive1_1.bed",sep="\t",quote=F,row.names=F,col.names=F) # Active fragments in hg19 (present 1x if active in both)

# Measured fragment BEDs
humanMeas<-MPRAfrag1[MPRAfrag1$noTags>11 & MPRAfrag1$Species=="hg19",c(2:4)]
humanMeas$Start<-humanMeas$Start-1
write.table(humanMeas,"data/human-measured-frags.bed",quote=F,row.names=F,col.names=F,sep="\t") # Measured human fragments
chimpMeas<-MPRAfrag1[MPRAfrag1$noTags>11 & MPRAfrag1$Species=="PanTro2",c(1:4)]
chimpMeas$Start<-chimpMeas$Start-1
write.table(chimpMeas[,-1],"data/chimp-measured-frags.bed",quote=F,row.names=F,col.names=F,sep="\t") # Measured chimpanzee fragments in panTro2
for(i in 1:nrow(chimpMeas)){
	chimpMeas$Chr[i]<-MPRAdataDummy$Chr[MPRAdataDummy$Alignment==sub("_Chimp",'',chimpMeas$Alignment[i])]
	chimpMeas$Start[i]<-MPRAdataDummy$Start[MPRAdataDummy$Alignment==sub("_Chimp",'',chimpMeas$Alignment[i])]
	chimpMeas$Stop[i]<-MPRAdataDummy$Stop[MPRAdataDummy$Alignment==sub("_Chimp",'',chimpMeas$Alignment[i])]
}
chimpMeas$Start<-chimpMeas$Start-1
write.table(chimpMeas[,-1],"data/panTro2tohg19-measured-frags.bed",quote=F,row.names=F,col.names=F,sep="\t") # Measured chimpanzee fragments in hg19
rm(list=c("bedAct","bedBoth","chimpBed","chimpMeas","chimpSpec","humanBed","humanMeas","humanSpec","i","MPRAdataDummy"))

### OUTPUT FRAGMENT TABLES FOR THE SECOND MPRA ###
# 2 - Differentially active fragments #
DAdata<-MPRAfrag1[MPRAfrag1$Alignment %in% paste(MPRAda1_1$Alignment,"Chimp",sep="_") | MPRAfrag1$Alignment %in% MPRAda1_1$Alignment,1:7]
DAdata<-cbind(DAdata,Ortholog_Seq=double(nrow(DAdata)))
DAdata$Ortholog_Seq<-replace(DAdata$Ortholog_Seq,DAdata$Species=="hg19",paste(DAdata$Alignment[DAdata$Species=="hg19"],"Chimp",sep="_"))
DAdata$Ortholog_Seq<-replace(DAdata$Ortholog_Seq,DAdata$Species=="PanTro2",sub("_Chimp","",DAdata$Alignment[DAdata$Species=="PanTro2"]))
DAdata<-DAdata[,c(1,8,2:7)]
write.table(DAdata,"data/DataForRd2_2_DAfrags.tsv",sep="\t",quote=F,row.names=F)
# 1 - Active fragments #
actDummy<-unique(sub("_Chimp","",MPRAactive1_1$Alignment))
actDummy<-actDummy[actDummy %out% DAdata$Alignment]
actDummy<-c(actDummy,paste(actDummy,"Chimp",sep="_"))
names(actDummy)<-actDummy
actData<-merge(MPRAdata1[,2:9],actDummy,by.x="Alignment",by.y=0)
actData<-unique(actData[,-9])
write.table(actData,"data/DataForRd2_1_ActFrags.tsv",sep="\t",quote=F,row.names=F)
rm(list=c("actData","actDummy","DAdata"))

q(save="no")
