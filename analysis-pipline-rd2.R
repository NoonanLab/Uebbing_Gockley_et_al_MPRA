### FUNCTIONS AND PACKAGES ###
library(tidyr)
library(dplyr)
library(biomaRt)
cohen.d<-function(x,y){
	mean1<-mean(x)
	mean2<-mean(y)
	sd1<-sd(x)
	sd2<-sd(y)
	n1<-length(x)
	n2<-length(y)
	mean_diff<-mean1-mean2
	pooled_sd_n<-((n1-1)*(sd1^2)) + ((n1-1)*(sd2^2))
	pooled_sd_d<-((n1+n2)-2)
	pooled_sd<-sqrt(pooled_sd_n/pooled_sd_d)
	cohens_d<-mean_diff/pooled_sd
	return(cohens_d)
}
mean_diff<-function(x,y){
	mean1<-mean(x)
	mean2<-mean(y)
	mean_diff<-mean1-mean2
	return(mean_diff)
}

### INPUT DATA ###
if(file.exists("data/MPRAnorm2.tsv.gz")){ # If file available, load; else create
	MPRAnorm2<-read.delim(gzfile("data/MPRAnorm2.tsv.gz"),stringsAsFactors=F)
}else{
	table<-read.delim(gzfile("data/MPRAcounts.tsv.gz"),stringsAsFactors=F)
	### NORMALIZATION AND FILTERS ###
	lib.sizes<-colSums(table[,11:14])
	MPRAnorm2<-table
	for(i in 1:4){
		MPRAnorm2[,(10+i)]<-log2((MPRAnorm2[,(10+i)]+.1)*1e6/lib.sizes[i]) # CPM
	}
	MPRAnorm2<-MPRAnorm2[rowSums(MPRAnorm2[,11:12])/2 >= -5.25,] # cutoff for Poisson error (only pDNA!)
	write.table(MPRAnorm2,gzfile("data/MPRAnorm2.tsv.gz"),row.names=F,sep="\t",quote=F)
	rm(list=c("i","lib.sizes","table"))
}

if(file.exists("data/MPRAfrag2.tsv")){ # If file available, load; else create
	MPRAfrag2<-read.delim("data/MPRAfrag2.tsv",stringsAsFactors=F)
}else{
	Rep1_cDNAdata<-MPRAnorm2[,c(2,13)] %>% nest(-Alignment) # group by fragment name
	lens<-double(length(Rep1_cDNAdata$data))
	names(lens)<-Rep1_cDNAdata$Alignment
	for(i in 1:length(lens)){
		lens[i]<-dim(Rep1_cDNAdata$data[[i]])[1]
	}
	Rep1_cDNAsortedVals<-unlist(Rep1_cDNAdata$data)
	names(Rep1_cDNAsortedVals)<-rep(Rep1_cDNAdata$Alignment,times=lens)
	Rep1_cDNAsortedVals12<-Rep1_cDNAsortedVals[names(Rep1_cDNAsortedVals) %in% names(lens[lens>11])]
	Rep1_cDNAmedian12<-tapply(Rep1_cDNAsortedVals12,names(Rep1_cDNAsortedVals12),median)
	
	Rep1_pDNAdata<-MPRAnorm2[,c(2,11)] %>% nest(-Alignment)
	Rep1_pDNAsortedVals<-unlist(Rep1_pDNAdata$data)
	names(Rep1_pDNAsortedVals)<-rep(Rep1_pDNAdata$Alignment,times=lens)
	Rep1_pDNAsortedVals12<-Rep1_pDNAsortedVals[names(Rep1_pDNAsortedVals) %in% names(lens[lens>11])]
	Rep1_pDNAmedian12<-tapply(Rep1_pDNAsortedVals12,names(Rep1_pDNAsortedVals12),median)
	Rep1_actMedian12<-tapply(Rep1_cDNAsortedVals12-Rep1_pDNAsortedVals12,names(Rep1_pDNAsortedVals12),median)
	
	Rep2_cDNAdata<-MPRAnorm2[,c(2,14)] %>% nest(-Alignment)
	Rep2_cDNAsortedVals<-unlist(Rep2_cDNAdata$data)
	names(Rep2_cDNAsortedVals)<-rep(Rep2_cDNAdata$Alignment,times=lens)
	Rep2_cDNAsortedVals12<-Rep2_cDNAsortedVals[names(Rep2_cDNAsortedVals) %in% names(lens[lens>11])]
	Rep2_cDNAmedian12<-tapply(Rep2_cDNAsortedVals12,names(Rep2_cDNAsortedVals12),median)
	
	Rep2_pDNAdata<-MPRAnorm2[,c(2,12)] %>% nest(-Alignment)
	Rep2_pDNAsortedVals<-unlist(Rep2_pDNAdata$data)
	names(Rep2_pDNAsortedVals)<-rep(Rep2_pDNAdata$Alignment,times=lens)
	Rep2_pDNAsortedVals12<-Rep2_pDNAsortedVals[names(Rep2_pDNAsortedVals) %in% names(lens[lens>11])]
	Rep2_pDNAmedian12<-tapply(Rep2_pDNAsortedVals12,names(Rep2_pDNAsortedVals12),median)
	Rep2_actMedian12<-tapply(Rep2_cDNAsortedVals12-Rep2_pDNAsortedVals12,names(Rep2_pDNAsortedVals12),median)
	
	MPRAfrag2<-unique(MPRAnorm2[,2:8])
	MPRAfrag2<-merge(MPRAfrag2,Rep1_cDNAmedian12,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep1_cDNAmedian"
	MPRAfrag2<-merge(MPRAfrag2,Rep2_cDNAmedian12,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep2_cDNAmedian"
	MPRAfrag2<-merge(MPRAfrag2,Rep1_pDNAmedian12,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep1_pDNAmedian"
	MPRAfrag2<-merge(MPRAfrag2,Rep2_pDNAmedian12,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep2_pDNAmedian"
	MPRAfrag2<-merge(MPRAfrag2,Rep1_actMedian12,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep1_actMedian"
	MPRAfrag2<-merge(MPRAfrag2,Rep2_actMedian12,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep2_actMedian"
	
	## Acitivity t tests of all fragments
	Rep1_activitySortedVals12<-Rep1_cDNAsortedVals12-Rep1_pDNAsortedVals12
	sds<-double(length(unique(names(Rep1_activitySortedVals12))))
	names(sds)<-unique(names(Rep1_activitySortedVals12))
	sds<-tapply(Rep1_activitySortedVals12,names(Rep1_activitySortedVals12),sd)
	Rep1_activitySortedVals12sd<-Rep1_activitySortedVals12[names(Rep1_activitySortedVals12) %in% names(sds[sds!=0])]
	
	Rep1_tcvpRes<-double(length(unique(names(Rep1_activitySortedVals12sd))))
	Rep1_tcvpRes<-tapply(Rep1_activitySortedVals12sd-median(Rep1_activitySortedVals12),
		names(Rep1_activitySortedVals12sd),t.test,alternative="greater")
	Rep1_NEGtcvpRes<-tapply(Rep1_activitySortedVals12sd-median(Rep1_activitySortedVals12[grep("NegCNTRL",names(Rep1_activitySortedVals12))]),
		names(Rep1_activitySortedVals12sd),t.test,alternative="greater")
	Rep1_tcvpP<-Rep1_NEGtcvpP<-double(length(Rep1_tcvpRes))
	for(i in 1:length(Rep1_tcvpP)){
		Rep1_tcvpP[i]<-Rep1_tcvpRes[[i]]$p.value
		Rep1_NEGtcvpP[i]<-Rep1_NEGtcvpRes[[i]]$p.value
	}
	names(Rep1_tcvpP)<-names(Rep1_NEGtcvpP)<-names(Rep1_tcvpRes)
	
	Rep2_activitySortedVals12<-Rep2_cDNAsortedVals12-Rep2_pDNAsortedVals12
	sds<-double(length(unique(names(Rep2_activitySortedVals12))))
	names(sds)<-unique(names(Rep2_activitySortedVals12))
	sds<-tapply(Rep2_activitySortedVals12,names(Rep2_activitySortedVals12),sd)
	Rep2_activitySortedVals12sd<-Rep2_activitySortedVals12[names(Rep2_activitySortedVals12) %in% names(sds[sds!=0])]
	
	Rep2_tcvpRes<-double(length(unique(names(Rep2_activitySortedVals12sd))))
	Rep2_tcvpRes<-tapply(Rep2_activitySortedVals12sd-median(Rep2_activitySortedVals12),
		names(Rep2_activitySortedVals12sd),t.test,alternative="greater")
	Rep2_NEGtcvpRes<-tapply(Rep2_activitySortedVals12sd-median(Rep2_activitySortedVals12[grep("NegCNTRL",names(Rep2_activitySortedVals12))]),
		names(Rep2_activitySortedVals12sd),t.test,alternative="greater")
	Rep2_tcvpP<-Rep2_NEGtcvpP<-double(length(Rep2_tcvpRes))
	for(i in 1:length(Rep2_tcvpP)){
		Rep2_tcvpP[i]<-Rep2_tcvpRes[[i]]$p.value
		Rep2_NEGtcvpP[i]<-Rep2_NEGtcvpRes[[i]]$p.value
	}
	names(Rep2_tcvpP)<-names(Rep2_NEGtcvpP)<-names(Rep2_tcvpRes)
	
	MPRAfrag2<-merge(MPRAfrag2,Rep1_tcvpP,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep1_ttest_P"
	MPRAfrag2<-merge(MPRAfrag2,Rep2_tcvpP,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep2_ttest_P"
	MPRAfrag2<-merge(MPRAfrag2,Rep1_NEGtcvpP,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep1_ttestVneg_P"
	MPRAfrag2<-merge(MPRAfrag2,Rep2_NEGtcvpP,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"Rep2_ttestVneg_P"
	
	MPRAfrag2<-merge(MPRAfrag2,lens,by.x="Alignment",by.y=0)
	colnames(MPRAfrag2)[ncol(MPRAfrag2)]<-"n_tags"
	MPRAactive2<-MPRAfrag2[p.adjust(MPRAfrag2$Rep1_ttestVneg_P,method="BH")<.05 & p.adjust(MPRAfrag2$Rep2_ttestVneg_P,method="BH")<.05,]
	write.table(MPRAfrag2,"data/MPRAfrag2.tsv",row.names=F,sep="\t",quote=F)
	rm(list=c("i","lens","Rep1_activitySortedVals12","Rep1_activitySortedVals12sd","Rep1_actMedian12","Rep1_cDNAdata","Rep1_cDNAmedian12",
		"Rep1_cDNAsortedVals","Rep1_cDNAsortedVals12","Rep1_NEGtcvpP","Rep1_NEGtcvpRes","Rep1_pDNAdata","Rep1_pDNAmedian12","Rep1_pDNAsortedVals",
		"Rep1_pDNAsortedVals12","Rep1_tcvpP","Rep1_tcvpRes","Rep2_activitySortedVals12","Rep2_activitySortedVals12sd","Rep2_actMedian12",
		"Rep2_cDNAdata","Rep2_cDNAmedian12","Rep2_cDNAsortedVals","Rep2_cDNAsortedVals12","Rep2_NEGtcvpP","Rep2_NEGtcvpRes","Rep2_pDNAdata",
		"Rep2_pDNAmedian12","Rep2_pDNAsortedVals","Rep2_pDNAsortedVals12","Rep2_tcvpP","Rep2_tcvpRes","sds"))
}

## COMPARISONS BTW MPRA RD 1 & 2
MPRAfrag1<-read.delim("data/MPRAfrag1.tsv",stringsAsFactors=F)
MPRAfrag2copy<-MPRAfrag2
MPRAfrag2copy$Alignment<-sub("_hsSUB_.*_Chimp_C+_B","_hsSUB_Chimp",MPRAfrag2copy$Alignment)
MPRAfrag2copy$Alignment<-sub("_hsSUB_.*_Human_H+_.*","_hsSUB",MPRAfrag2copy$Alignment)
MPRAfrag2copy$Alignment<-sub("_hsSUB_.*_Human_C+_NB","_hsSUB_Chimp",MPRAfrag2copy$Alignment)
MPRAfrag12<-merge(MPRAfrag1[,c(1,20:23)],MPRAfrag2copy[,c(1,12,13)])

# Correlations
cor.test(MPRAfrag12$Rep12_actMedian,MPRAfrag12$Rep13_actMedian) # p-value < 2.2e-16, r = 0.9134224
cor.test(MPRAfrag12$Rep12_actMedian,MPRAfrag12$Rep22_actMedian) # p-value < 2.2e-16, r = 0.9051204
cor.test(MPRAfrag12$Rep12_actMedian,MPRAfrag12$Rep23_actMedian) # p-value < 2.2e-16, r = 0.9066708
cor.test(MPRAfrag12$Rep13_actMedian,MPRAfrag12$Rep22_actMedian) # p-value < 2.2e-16, r = 0.9056954
cor.test(MPRAfrag12$Rep13_actMedian,MPRAfrag12$Rep23_actMedian) # p-value < 2.2e-16, r = 0.9066865
cor.test(MPRAfrag12$Rep22_actMedian,MPRAfrag12$Rep23_actMedian) # p-value < 2.2e-16, r = 0.9117116

cor.test(MPRAfrag12$Rep1_actMedian,MPRAfrag12$Rep2_actMedian) # p-value < 2.2e-16, r = 0.9446098

cor.test(MPRAfrag12$Rep12_actMedian,MPRAfrag12$Rep1_actMedian) # p-value < 2.2e-16, r = 0.9079
cor.test(MPRAfrag12$Rep13_actMedian,MPRAfrag12$Rep1_actMedian) # p-value < 2.2e-16, r = 0.913927
cor.test(MPRAfrag12$Rep22_actMedian,MPRAfrag12$Rep1_actMedian) # p-value < 2.2e-16, r = 0.9189342
cor.test(MPRAfrag12$Rep23_actMedian,MPRAfrag12$Rep1_actMedian) # p-value < 2.2e-16, r = 0.9004216
cor.test(MPRAfrag12$Rep12_actMedian,MPRAfrag12$Rep2_actMedian) # p-value < 2.2e-16, r = 0.8848698
cor.test(MPRAfrag12$Rep13_actMedian,MPRAfrag12$Rep2_actMedian) # p-value < 2.2e-16, r = 0.8924979
cor.test(MPRAfrag12$Rep22_actMedian,MPRAfrag12$Rep2_actMedian) # p-value < 2.2e-16, r = 0.904391
cor.test(MPRAfrag12$Rep23_actMedian,MPRAfrag12$Rep2_actMedian) # p-value < 2.2e-16, r = 0.8941134

# No. of active fragments
MPRAactive1_2<-read.delim("data/MPRAactive1_2.tsv",stringsAsFactors=F)
nrow(MPRAfrag2copy[MPRAfrag2copy$Alignment %in% MPRAactive1_2$Alignment,]) # 3948
nrow(MPRAfrag2copy[(p.adjust(MPRAfrag2copy$Rep1_ttestVneg_P,method="BH")<.05 & p.adjust(MPRAfrag2copy$Rep2_ttestVneg_P,method="BH")<.05) 
		& MPRAfrag2copy$Alignment %in% MPRAactive1_2$Alignment,]) # 2481
rm(MPRAfrag2copy)

## GENERATE OBJECTS CONTAINING DIFFERENT FRAGMENT TYPES
dTable2<-cbind(MPRAnorm2,act1=MPRAnorm2$cDNA1-MPRAnorm2$pDNA1,act2=MPRAnorm2$cDNA2-MPRAnorm2$pDNA2)
Add<-dTable2[dTable2$Subtype=="Add",]
NonAdd<-dTable2[dTable2$Subtype=="NonAdd",]
Add_list<-Add[,c(3,6,7,15,16)] %>% nest(-Pos)

# Generate folded (human vs chimpanzee allele) fragment table
if(file.exists("data/MPRAfold2.tsv")){ # If file available, load; else create
	MPRAfrag2<-read.delim("data/MPRAfold2.tsv",stringsAsFactors=F)
}else{
	MPRAfold2<-MPRAfrag2[grepl("_Human_H+_",MPRAfrag2$Alignment),c(2:4,12,13)]
	colnames(MPRAfold2)[4:5]<-paste0(colnames(MPRAfold2)[4:5],"_hs")
	MPRAfold2<-merge(MPRAfold2,MPRAfrag2[grepl("_Chimp_C+_",MPRAfrag2$Alignment) | (grepl("_Human_C+_",MPRAfrag2$Alignment) & MPRAfrag2$BG=="NB"),
		c(2:4,12,13)])
	colnames(MPRAfold2)[6:7]<-paste0(colnames(MPRAfold2)[6:7],"_pt")
	write.table(MPRAfold2,"data/MPRAfold2.tsv",quote=F,row.names=F,sep="\t")
}

# Test replication of differential activity using folded table
if(file.exists("data/MPRAfoldActive2.tsv")){ # If file available, load; else create
	MPRAfrag2<-read.delim("data/MPRAfoldActive2.tsv",stringsAsFactors=F)
}else{
	MPRAfoldActive2<-MPRAfold2[MPRAfold2$Pos %in% unique(MPRAactive2$Pos),]
	## Validation differential activity test of combinatorial, validation and tile frags
	MPRAactivities2<-cbind(MPRAnorm2[,c(1:7,11:14)],
		Rep1_activity=MPRAnorm2$cDNA1-MPRAnorm2$pDNA1,Rep2_activity=MPRAnorm2$cDNA2-MPRAnorm2$pDNA2)
	MPRAvalActivities<-MPRAactivities2[grep("_Human_H+_|_Chimp_C+_|_Human_C+_NB",MPRAactivities2$Alignment),]
	MPRAvalActivities$Species<-replace(MPRAvalActivities$Species,grep("_Human_C+_NB",MPRAvalActivities$Alignment),"Chimp")
	
	dummy<-unique(MPRAvalActivities[,c(3,6)])
	names2List<-dummy$Pos[duplicated(dummy$Pos)]
	Rep1_nestAct<-MPRAvalActivities[MPRAvalActivities$Pos %in% MPRAactive2$Pos 
		& MPRAvalActivities$Pos %in% names2List,c(3,6,12)] %>% nest(-Pos, Species)
	Rep2_nestAct<-MPRAvalActivities[MPRAvalActivities$Pos %in% MPRAactive2$Pos 
		& MPRAvalActivities$Pos %in% names2List,c(3,6,13)] %>% nest(-Pos, Species)
	Rep1_tDAActRes<-list(double(length(Rep1_nestAct$Pos)))
	Rep2_tDAActRes<-list(double(length(Rep2_nestAct$Pos)))
	Rep1_tDAActp<-Rep1_tDAActt<-rep(NA,length(Rep1_nestAct$Pos))
	Rep2_tDAActp<-Rep2_tDAActt<-rep(NA,length(Rep2_nestAct$Pos))
	for(i in 1:length(Rep1_nestAct$Pos)){
		if(i%%100==0){
			print(i)
			print(Sys.time())
		}
		if(any(table(Rep1_nestAct$data[[i]]$Species)==1)){
			Rep1_tDAActRes[[i]]<-Rep2_tDAActRes[[i]]<-NA
		}else{
		Rep1_tDAActRes[[i]]<-t.test(Rep1_activity ~ Species,data=Rep1_nestAct$data[[i]])
		Rep2_tDAActRes[[i]]<-t.test(Rep2_activity ~ Species,data=Rep2_nestAct$data[[i]])
		Rep1_tDAActp[i]<-Rep1_tDAActRes[[i]]$p.value
		Rep2_tDAActp[i]<-Rep2_tDAActRes[[i]]$p.value
		Rep1_tDAActt[i]<-Rep1_tDAActRes[[i]]$statistic
		Rep2_tDAActt[i]<-Rep2_tDAActRes[[i]]$statistic
	}}
	names(Rep1_tDAActp)<-names(Rep2_tDAActp)<-names(Rep1_tDAActt)<-names(Rep2_tDAActt)<-Rep1_nestAct$Pos
	MPRAfoldActive2<-merge(MPRAfoldActive2,Rep1_tDAActt,by.x="Pos",by.y=0)
	colnames(MPRAfoldActive2)[ncol(MPRAfoldActive2)]<-"Rep1_tDAActt"
	MPRAfoldActive2<-merge(MPRAfoldActive2,Rep1_tDAActp,by.x="Pos",by.y=0)
	colnames(MPRAfoldActive2)[ncol(MPRAfoldActive2)]<-"Rep1_tDAActp"
	MPRAfoldActive2<-merge(MPRAfoldActive2,Rep2_tDAActt,by.x="Pos",by.y=0)
	colnames(MPRAfoldActive2)[ncol(MPRAfoldActive2)]<-"Rep2_tDAActt"
	MPRAfoldActive2<-merge(MPRAfoldActive2,Rep2_tDAActp,by.x="Pos",by.y=0)
	colnames(MPRAfoldActive2)[ncol(MPRAfoldActive2)]<-"Rep2_tDAActp"
	write.table(MPRAfoldActive2,"data/MPRAfoldActive2.tsv",quote=F,sep="\t",row.names=F)
	rm(list=c("dummy","i","names2List","Rep1_nestAct","Rep1_tDAActp","Rep1_tDAActRes","Rep1_tDAActt",
		"Rep2_nestAct","Rep2_tDAActp","Rep2_tDAActRes","Rep2_tDAActt"))
}

# Comparisons of numbers of differentially active fragments btw. MPRA rounds
nrow(MPRAfold2[MPRAfold2$Pos %in% sub("_hsSUB.*","",MPRAactive1_2$Alignment) & MPRAfold2$Subtype=="Add",]) # 840
MPRAfoldRep2<-MPRAfoldActive2[MPRAfoldActive2$Subtype=="Add" & MPRAfoldActive2$Pos %in% sub("_hsSUB.*","",MPRAactive1_2$Alignment),]
nrow(MPRAfoldRep2) # 747
nrow(MPRAfoldRep2[p.adjust(MPRAfoldRep2$Rep1_tDAActp,method="BH")<.05 & p.adjust(MPRAfoldRep2$Rep2_tDAActp,method="BH")<.05,]) # 411
nrow(MPRAfoldRep2[p.adjust(MPRAfoldRep2$Rep1_tDAActp,method="BH")<.05 & p.adjust(MPRAfoldRep2$Rep2_tDAActp,method="BH")<.05 & 
	MPRAfoldRep2$Rep1_actMedian_hs>MPRAfoldRep2$Rep1_actMedian_pt & MPRAfoldRep2$Rep2_actMedian_hs>MPRAfoldRep2$Rep2_actMedian_pt & 
	MPRAfoldRep2$Rep1_actMedian_hs-MPRAfoldRep2$Rep1_actMedian_pt+MPRAfoldRep2$Rep2_actMedian_hs-MPRAfoldRep2$Rep2_actMedian_pt >.1,]) # 215
nrow(MPRAfoldRep2[p.adjust(MPRAfoldRep2$Rep1_tDAActp,method="BH")<.05 & p.adjust(MPRAfoldRep2$Rep2_tDAActp,method="BH")<.05 & 
	MPRAfoldRep2$Rep1_actMedian_hs<MPRAfoldRep2$Rep1_actMedian_pt & MPRAfoldRep2$Rep2_actMedian_hs<MPRAfoldRep2$Rep2_actMedian_pt & 
	MPRAfoldRep2$Rep1_actMedian_hs-MPRAfoldRep2$Rep1_actMedian_pt+MPRAfoldRep2$Rep2_actMedian_hs-MPRAfoldRep2$Rep2_actMedian_pt < -.1,]) # 196

## ANOVA TEST OF hSUB EFFECTS
if(file.exists("data/iEffectsizes.tsv")){ # If effect size tables available, load; else create
	aEffectSizes<-read.delim("data/aEffectsizes.tsv",stringsAsFactors=F)
	iEffectSizes<-read.delim("data/iEffectsizes.tsv",stringsAsFactors=F)
}else{
	test1_list<-test2_list<-n_list<-Add_list
	for(i in 1:length(Add_list$data)){
		if(i%%100==0){
			print(i)
			print(Sys.time())
		}
		skip<-F
		n_list$data[[i]]<-n<-nchar(Add_list$data[[i]]$Allele[1])
		names<-double(n)
		for(m in 1:n){
			names[m]<-paste0("hSub",m)
		}
		dd<-Add_list$data[[i]] %>% separate(Allele,into=names,sep=1:n)
		chk<-dd[,head(colnames(dd),n+1)]
		drop<-NULL
		for(o in names(chk)){
			if(length(unique(chk[[o]]))<2 | nrow(chk)==length(unique(chk[[o]]))){
				drop<-c(drop,o)
		}}
		chk<-within(chk,rm(list=drop))
		n<-ncol(chk)+1
		if(ncol(chk)<1){
			skip<-T
			test1_list$data[[i]]<-NA
			test2_list$data[[i]]<-NA
		}
		if(skip){next
		}else{
			dd1<-cbind(chk,act1=dd$act1)
			dd2<-cbind(chk,act2=dd$act2)
			for(p in (n-1):1){
				if(p>1){
					test1_list$data[[i]]<-anova(lm(paste0("act1 ~ .^",p),data=dd1))
				}else{
					test1_list$data[[i]]<-anova(lm(act1 ~ .,data=dd1))
					break
				}
				if(any(test1_list$data[[i]]$'Pr(>F)'[grepl(paste0("(.+:){",p-1,"}"),rownames(test1_list$data[[i]]))]<.05)){
					break
			}}
			for(p in (n-1):1){
				if(p>1){
					test2_list$data[[i]]<-anova(lm(paste0("act2 ~ .^",p),data=dd2))
				}else{
					test2_list$data[[i]]<-anova(lm(act2 ~ .,data=dd2))
					break
				}
				if(any(test2_list$data[[i]]$'Pr(>F)'[grepl(paste0("(.+:){",p-1,"}"),rownames(test2_list$data[[i]]))]<.05)){
					break
	}}}}
	
	# Parsing test results
	test1sum<-test2sum<-0
	for(i in 1:length(test1_list$data)){
		if(!is.na(test1_list$data[[i]][1,1])){
			test1sum<-test1sum+nrow(test1_list$data[[i]])-1
			test2sum<-test2sum+nrow(test2_list$data[[i]])-1
	}}
	test1results<-data.frame(Pos=double(test1sum),Effect=double(test1sum),Df=double(test1sum),
		'Sum Sq'=double(test1sum),'Mean Sq'=double(test1sum),'F value'=double(test1sum),'Pr(>F)'=double(test1sum))
	test2results<-data.frame(Pos=double(test2sum),Effect=double(test2sum),Df=double(test2sum),
		'Sum Sq'=double(test2sum),'Mean Sq'=double(test2sum),'F value'=double(test2sum),'Pr(>F)'=double(test2sum))
	cnt1<-cnt2<-0
	for(i in 1:length(test1_list$data)){
		if(!is.na(test1_list$data[[i]][1,1])){
			n<-nrow(test1_list$data[[i]])-1
			test1results[(cnt1+1):(cnt1+n),1]<-Add_list$Pos[[i]]
			test1results[(cnt1+1):(cnt1+n),2]<-rownames(test1_list$data[[i]])[1:n]
			test1results[(cnt1+1):(cnt1+n),3:7]<-test1_list$data[[i]][1:n,]
			cnt1<-cnt1+n
			n<-nrow(test2_list$data[[i]])-1
			test2results[(cnt2+1):(cnt2+n),1]<-Add_list$Pos[[i]]
			test2results[(cnt2+1):(cnt2+n),2]<-rownames(test2_list$data[[i]])[1:n]
			test2results[(cnt2+1):(cnt2+n),3:7]<-test2_list$data[[i]][1:n,]
			cnt2<-cnt2+n
	}}
	testResults<-merge(test1results[,c(1,2,6,7)],test2results[,c(1,2,6,7)],by=c("Pos","Effect"),all=T)
	colnames(testResults)[3:6]<-c("Rep1_F","Rep1_P","Rep2_F","Rep2_P")
	sign<-data.frame(Pos=n_list$Pos,n_eVar=unlist(n_list$data))
	sign$Pos<-as.character(sign$Pos)
	signRes<-merge(testResults[p.adjust(testResults$Rep1_P,method="BH")<.05 & !is.na(testResults$Rep1_P),1:4],
		testResults[p.adjust(testResults$Rep2_P,method="BH")<.05 & !is.na(testResults$Rep2_P),c(1,2,5,6)],all=T)
	sign<-merge(sign,signRes,all=T)
	rm(signRes)
	dSign<-data.frame(Pos=n_list$Pos,n_eVar=unlist(n_list$data))
	dSign$Pos<-as.character(dSign$Pos)
	dSign<-merge(dSign,testResults[p.adjust(testResults$Rep1_P,method="BH")<.05 & p.adjust(testResults$Rep2_P,method="BH")<.05 & 
		!is.na(testResults$Rep1_P) & !is.na(testResults$Rep2_P),],all.x=T)
	
	testResultsEffects<-testResultsP<-as.character(rep(NA,length(unique(dSign$Pos))))
	names(testResultsEffects)<-names(testResultsP)<-unique(dSign$Pos)
	for(i in 1:nrow(dSign)){
		if(!is.na(dSign[i,3])){
			testResultsEffects[dSign[i,1]]<-paste(testResultsEffects[dSign[i,1]],dSign[i,3],sep=",")
			testResultsP[dSign[i,1]]<-paste(testResultsP[dSign[i,1]],as.character(dSign[i,5]),sep=",")
			testResultsP[dSign[i,1]]<-paste(testResultsP[dSign[i,1]],as.character(dSign[i,7]),sep=",")
	}}
	testResultsEffects<-sub("NA,","",testResultsEffects)
	testResultsEffects<-testResultsEffects[which(testResultsEffects!="NA")]
	testResultsP<-sub("NA,","",testResultsP)
	testResultsP<-testResultsP[which(testResultsP!="NA")]
	CombiData<-merge(MPRAfrag2[,c(2,5,6,12,13)],testResultsEffects,by.x="Pos",by.y=0)
	colnames(CombiData)[ncol(CombiData)]<-"Significant_effects"
	CombiData<-merge(CombiData,testResultsP,by.x="Pos",by.y=0)
	colnames(CombiData)[ncol(CombiData)]<-"Effects_P"
	rm(list=c("chk","cnt1","cnt2","dd","dd1","dd2","drop","i","m","n","names","o","p","skip"))
	
	# Effect size extraction
	dSignNA<-dSign[!is.na(dSign$Effect),] %>% nest(-Pos)
	iCnt<-aCnt<-0
	for(i in 1:length(dSignNA$Pos)){
		tFrag<-dTable2[dTable2$Pos==dSignNA$Pos[[i]],c(6:7,15:16)]
		if(any(grepl(":",dSignNA$data[[i]]$Effect))){
			dSignNA$data[[i]]<-cbind(dSignNA$data[[i]],meanDiff1=as.double(rep(NA,nrow(dSignNA$data[[i]]))),
				meanDiff2=as.double(rep(NA,nrow(dSignNA$data[[i]]))),effectSize1=as.double(rep(NA,nrow(dSignNA$data[[i]]))),
				effectSize2=as.double(rep(NA,nrow(dSignNA$data[[i]]))),humanEffectSize1=as.double(rep(NA,nrow(dSignNA$data[[i]]))),
				humanEffectSize2=as.double(rep(NA,nrow(dSignNA$data[[i]]))))
			iEffect<-dSignNA$data[[i]]$Effect[grep(":",dSignNA$data[[i]]$Effect)]
			for(j in iEffect){
				if(grepl("^Species:hSub[1-7]$",j)){
					iCnt<-iCnt+1
					tp<-substr(j,13,13)
					msp<-c(mean_diff(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Chimp",3:4])),
						mean_diff(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Chimp",3:4])))
					mhS<-c(mean_diff(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Human",3:4])),
						mean_diff(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Chimp",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Chimp",3:4])))
					dSignNA$data[[i]]$meanDiff1[dSignNA$data[[i]]$Effect==j]<-msp[which.max(abs(msp))]
					dSignNA$data[[i]]$meanDiff2[dSignNA$data[[i]]$Effect==j]<-mhS[which.max(abs(mhS))]
					sp<-c(cohen.d(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Chimp",3:4])),
						cohen.d(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Chimp",3:4])))
					hS<-c(cohen.d(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Human",3:4])),
						cohen.d(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Chimp",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Chimp",3:4])))
					dSignNA$data[[i]]$effectSize1[dSignNA$data[[i]]$Effect==j]<-sp[which.max(abs(sp))]
					dSignNA$data[[i]]$effectSize2[dSignNA$data[[i]]$Effect==j]<-hS[which.max(abs(hS))]
					sp<-cohen.d(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Chimp",3:4]))
					hS<-cohen.d(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H" & tFrag$Species=="Human",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C" & tFrag$Species=="Human",3:4]))
					dSignNA$data[[i]]$humanEffectSize1[dSignNA$data[[i]]$Effect==j]<-sp
					dSignNA$data[[i]]$humanEffectSize2[dSignNA$data[[i]]$Effect==j]<-hS
				}else if(grepl("^hSub[1-6]:hSub[2-7]$",j)){
					iCnt<-iCnt+1
					tp1<-substr(j,5,5)
					tp2<-substr(j,11,11)
					mh1<-c(mean_diff(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="H",3:4])),
						mean_diff(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="C",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="C",3:4])))
					mh2<-c(mean_diff(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="C",3:4])),
						mean_diff(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="C",3:4])))
					dSignNA$data[[i]]$meanDiff1[dSignNA$data[[i]]$Effect==j]<-mh1[which.max(abs(mh1))]
					dSignNA$data[[i]]$meanDiff2[dSignNA$data[[i]]$Effect==j]<-mh2[which.max(abs(mh2))]
					h1<-c(cohen.d(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="H",3:4])),
						cohen.d(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="C",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="C",3:4])))
					h2<-c(cohen.d(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="C",3:4])),
						cohen.d(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="C",3:4])))
					dSignNA$data[[i]]$effectSize1[dSignNA$data[[i]]$Effect==j]<-h1[which.max(abs(h1))]
					dSignNA$data[[i]]$effectSize2[dSignNA$data[[i]]$Effect==j]<-h2[which.max(abs(h2))]
					h1<-cohen.d(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="C" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]))
					h2<-cohen.d(unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="H",3:4]),
							unlist(tFrag[substr(tFrag$Allele,tp1,tp1)=="H" & substr(tFrag$Allele,tp2,tp2)=="C",3:4]))
					dSignNA$data[[i]]$humanEffectSize1[dSignNA$data[[i]]$Effect==j]<-h1
					dSignNA$data[[i]]$humanEffectSize2[dSignNA$data[[i]]$Effect==j]<-h2
		}}}else{
			dSignNA$data[[i]]<-cbind(dSignNA$data[[i]],meanDiff=as.double(rep(NA,nrow(dSignNA$data[[i]]))),
				effectSize=as.double(rep(NA,nrow(dSignNA$data[[i]]))))
			for(j in dSignNA$data[[i]]$Effect){
				aCnt<-aCnt+1
				if(j=="Species"){
					dSignNA$data[[i]]$meanDiff[dSignNA$data[[i]]$Effect==j]<-mean_diff(unlist(tFrag[tFrag$Species=="Human",3:4]),
						unlist(tFrag[tFrag$Species=="Chimp",3:4]))
					dSignNA$data[[i]]$effectSize[dSignNA$data[[i]]$Effect==j]<-cohen.d(unlist(tFrag[tFrag$Species=="Human",3:4]),
						unlist(tFrag[tFrag$Species=="Chimp",3:4]))
				}else{
					tp<-substr(j,5,5)
					dSignNA$data[[i]]$meanDiff[dSignNA$data[[i]]$Effect==j]<-mean_diff(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H",3:4]),
						unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C",3:4]))
					dSignNA$data[[i]]$effectSize[dSignNA$data[[i]]$Effect==j]<-cohen.d(unlist(tFrag[substr(tFrag$Allele,tp,tp)=="H",3:4]),
						unlist(tFrag[substr(tFrag$Allele,tp,tp)=="C",3:4]))
	}}}}
	iEffectSizes<-data.frame(Pos=double(iCnt),Effect=double(iCnt),meanDiff1=double(iCnt),meanDiff2=double(iCnt),effectSize1=double(iCnt),
		effectSize2=double(iCnt),humanEffectSize1=double(iCnt),humanEffectSize2=double(iCnt))
	aEffectSizes<-data.frame(Pos=double(aCnt),Effect=double(aCnt),meanDiff=double(aCnt),effectSize=double(aCnt))
	iCnt<-aCnt<-n<-m<-0
	for(i in 1:length(dSignNA$Pos)){
		if(any(grepl(":",dSignNA$data[[i]]$Effect))){
			if(any(grepl("^Species:hSub[1-7]$|^hSub[1-6]:hSub[2-7]",dSignNA$data[[i]]$Effect))){
				tmp<-dSignNA$data[[i]][grepl("^Species:hSub[1-7]$|^hSub[1-6]:hSub[2-7]$",dSignNA$data[[i]]$Effect),]
				m<-nrow(tmp)
				iEffectSizes[(iCnt+1):(iCnt+m),1]<-dSignNA$Pos[[i]]
				iEffectSizes[(iCnt+1):(iCnt+m),2:8]<-tmp[,c(2,7:12)]
				iCnt<-iCnt+m
		}}else{
			n<-nrow(dSignNA$data[[i]])
			aEffectSizes[(aCnt+1):(aCnt+n),1]<-dSignNA$Pos[[i]]
			aEffectSizes[(aCnt+1):(aCnt+n),2:4]<-dSignNA$data[[i]][,c(2,7,8)]
			aCnt<-aCnt+n
	}}
	write.table(aEffectSizes,"data/aEffectsizes.tsv",quote=F,row.names=F,sep="\t")
	write.table(iEffectSizes,"data/iEffectsizes.tsv",quote=F,row.names=F,sep="\t")
	rm(list=c("aCnt","h1","h2","hS","i","iCnt","iEffect","j","m","mh1","mh2","mhS","msp","n","n_list","sp","tFrag","tmp","tp","tp1","tp2"))
}

pdf("plots/effectSizes.pdf")
par(mfcol=c(2,2),las=1)
hist(abs(aEffectSizes$effectSize[aEffectSizes$Effect=="Species"]),xlim=c(0,5.5),breaks=seq(0,5.5,.1),main='Additive background effects')
hist(abs(iEffectSizes$effectSize1[grep("Species:",iEffectSizes$Effect)]),xlim=c(0,5.5),breaks=seq(0,5.5,.1),main='Interactive background effects')
hist(abs(aEffectSizes$effectSize[aEffectSizes$Effect!="Species"]),xlim=c(0,5.5),breaks=seq(0,5.5,.1),main='Additive eVar effects')
hist(abs(c(iEffectSizes$effectSize1[grep("^eVar[1-6]:eVar[2-7]$",iEffectSizes$Effect)],iEffectSizes$effectSize2)),xlim=c(0,5.5),breaks=seq(0,5.5,.1),
	main='Interactive eVar effects')
dev.off()

if(file.exists("data/hSub-summary.tsv")){ # If effect size tables available, load; else create
	sumTable<-read.delim("data/hSub-summary.tsv",stringsAsFactors=F)
}else{
	MPRAfold1<-read.delim("data/MPRAfold1.tsv",stringsAsFactors=F)
	colnames(MPRAfold1)[5:6]<-c("hSub_BPChange","hSub_Pos")
	hSs<-MPRAfold1[,1:6] %>% nest(-c(Alignment,hSub_Pos,hSub_BPChange))
	hSs$hSub_Pos<-strsplit(hSs$hSub_Pos,";")
	hSs$hSub_BPChange<-strsplit(hSs$hSub_BPChange,";")
	n<-unlist(length(unlist(hSs$hSub_Pos)))
	dummyBED<-data.frame(Pos=character(n),Chr=character(n),Start=integer(n),Stop=integer(n),Effect=character(n),hSub_Pos=integer(n),
		hSub_BPChange=integer(n),stringsAsFactors=F)
	cnt<-0
	for(i in 1:length(hSs$data)){
		n<-length(hSs$hSub_Pos[[i]])
		dummy<-data.frame(Pos=as.character(rep(sub("_hsSUB","",hSs$Alignment[[i]]),n)),Chr=as.character(rep(hSs$data[[i]]$Chr,n)),
			Start=rep(hSs$data[[i]]$Start,n),Stop=rep(hSs$data[[i]]$Stop,n),Effect=as.character(paste0("hSub",1:length(hSs$hSub_Pos[[i]]))),
			hSub_Pos=as.integer(hSs$hSub_Pos[[i]])+hSs$data[[i]]$Start,hSub_BPChange=hSs$hSub_BPChange[[i]],stringsAsFactors=F)
		dummyBED[(cnt+1):(cnt+n),]<-dummy
		cnt<-cnt+n
	}
	# hSub effects summary table
	sumH2<-aEffectSizes[aEffectSizes$Effect!="Species",c(1:4,2)]
	colnames(sumH2)[5]<-"Effects"
	sumH3<-iEffectSizes[grepl("Species:",iEffectSizes$Effect),c(1,2,4,6,2)]
	sumH3$Effect<-sub("Species:","",sumH3$Effect)
	names(sumH3)[3:5]<-c("meanDiff","effectSize","Effects")
	sumH2<-rbind(sumH2,sumH3)
	sumH3<-iEffectSizes[grepl("Species:",iEffectSizes$Effect)==F,c(1:3,5,2)]
	sumH3$Effect<-sub(":hSub[1-7]$","",sumH3$Effect)
	names(sumH3)[3:5]<-c("meanDiff","effectSize","Effects")
	sumH2<-rbind(sumH2,sumH3)
	sumH3<-iEffectSizes[grepl("Species:",iEffectSizes$Effect)==F,c(1,2,4,6,2)]
	sumH3$Effect<-sub("^hSub[1-7]:hSub","hSub",sumH3$Effect)
	names(sumH3)[3:5]<-c("meanDiff","effectSize","Effects")
	sumH2<-rbind(sumH2,sumH3)
	sumTAB<-merge(dummyBED,sumH2)
	hSnFrag<-read.delim("data/hSubs-vs-combFrags.BED",stringsAsFactors=F,header=F)[,c(1,3,4)] # "190124_eVars-vs-combFrags.BED"
	colnames(hSnFrag)<-c("Chr","hSub_Pos","nFrag")
	sumTAB<-merge(sumTAB,hSnFrag)
	sumLst<-sumTAB[,c(1,2,7:11)] %>% nest(-c(Chr,hSub_Pos,hSub_BPChange,nFrag))
	n<-nrow(sumLst)
	sumTable<-data.frame(Chr=sumLst$Chr,Pos=sumLst$hSub_Pos,baseChange=sumLst$hSub_BPChange,nFrag=sumLst$nFrag,max_add_mean_diff=rep(NA,n),
		average_add_mean_diff=rep(NA,n),max_add_effect_size=rep(NA,n),average_add_effect_size=rep(NA,n),interactions=NA,
		max_intera_mean_diff=rep(NA,n),average_intera_mean_diff=rep(NA,n),max_intera_effect_size=rep(NA,n),average_intera_effect_size=rep(NA,n),
		stringsAsFactors=F)
	for(i in 1:n){
		addTmp<-sumLst$data[[i]][grepl(":",sumLst$data[[i]]$Effects)==F,]
		intTmp<-sumLst$data[[i]][grepl(":",sumLst$data[[i]]$Effects),]
		if(nrow(addTmp)>0){
			sumTable$max_add_mean_diff[i]<-addTmp$meanDiff[which.max(abs(addTmp$meanDiff))]
			sumTable$average_add_mean_diff[i]<-mean(addTmp$meanDiff)
			sumTable$max_add_effect_size[i]<-addTmp$effectSize[which.max(abs(addTmp$effectSize))]
			sumTable$average_add_effect_size[i]<-mean(addTmp$effectSize)
		}
		if(nrow(intTmp)>0){
			sumTable$max_intera_mean_diff[i]<-intTmp$meanDiff[which.max(abs(intTmp$meanDiff))]
			sumTable$average_intera_mean_diff[i]<-mean(intTmp$meanDiff)
			sumTable$max_intera_effect_size[i]<-intTmp$effectSize[which.max(abs(intTmp$effectSize))]
			sumTable$average_intera_effect_size[i]<-mean(intTmp$effectSize)
			nr<-nrow(intTmp[grepl("Species",intTmp$Effects),])
			if(nrow(intTmp)==nr){
				sumTable$interactions[i]<-"Background"
			}else if(nr==0){
				sumTable$interactions[i]<-"hSub"
			}else{
				sumTable$interactions[i]<-"Background,hSub"
	}}}
	hSubFeat<-read.delim("data/hSubs-vs-GainEnh-HAR-HACNS.tsv",stringsAsFactors=F)
	sumTable<-merge(sumTable,hSubFeat)
	TS3_Reilly<-read.delim("data/TableS3_geneID.tsv",stringsAsFactors=F)
	TS3_Reilly$external_gene_name<-toupper(TS3_Reilly$external_gene_name)
	tfbs<-unique(read.delim("data/hSubs-TFBS_intersect_symbol_expressed.BED",stringsAsFactors=F,header=F)[,c(1,4:6,8)])
	colnames(tfbs)<-c("Chr","TFBS","hg_score","pt_score","Pos")
	tfbs$TFBS<-toupper(tfbs$TFBS)
	tfbs<-merge(tfbs,TS3_Reilly[,c(2,6)],by.x="TFBS",by.y="external_gene_name",all.x=T)
	tfbsLst<-tfbs %>% nest(-c(Chr,Pos))
	tfbsTable<-data.frame(Chr=tfbsLst$Chr,Pos=tfbsLst$Pos,TFBS_hg_uniq=NA,hg_uniq_WGCNA=NA,TFBS_hg_stronger=NA,hg_stronger_WGCNA=NA,TFBS_pt_stronger=NA,
		pt_stronger_WGCNA=NA,TFBS_pt_uniq=NA,pt_uniq_WGCNA=NA,stringsAsFactors=F)
	for(i in 1:nrow(tfbsLst)){
		tmp<-tfbsLst$data[[i]]
		tfbsTable$TFBS_hg_uniq[i]<-paste(tmp$TFBS[!is.na(tmp$hg_score) & is.na(tmp$pt_score) & tmp$hg_score>43],collapse=",")
		tfbsTable$hg_uniq_WGCNA[i]<-paste(tmp$ModuleAssignment[!is.na(tmp$hg_score) & is.na(tmp$pt_score) & tmp$hg_score>43],collapse=",")
		tfbsTable$TFBS_hg_stronger[i]<-paste(tmp$TFBS[(!is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$hg_score>tmp$pt_score) |
			(!is.na(tmp$hg_score) & is.na(tmp$pt_score) & tmp$hg_score<43)],collapse=",")
		tfbsTable$hg_stronger_WGCNA[i]<-paste(tmp$ModuleAssignment[(!is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$hg_score>tmp$pt_score) |
			(!is.na(tmp$hg_score) & is.na(tmp$pt_score) & tmp$hg_score<43)],collapse=",")
		tfbsTable$TFBS_pt_stronger[i]<-paste(tmp$TFBS[(!is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$hg_score<tmp$pt_score) |
			(is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$pt_score<43)],collapse=",")
		tfbsTable$pt_stronger_WGCNA[i]<-paste(tmp$ModuleAssignment[(!is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$hg_score<tmp$pt_score) |
			(is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$pt_score<43)],collapse=",")
		tfbsTable$TFBS_pt_uniq[i]<-paste(tmp$TFBS[is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$pt_score>43],collapse=",")
		tfbsTable$pt_uniq_WGCNA[i]<-paste(tmp$ModuleAssignment[is.na(tmp$hg_score) & !is.na(tmp$pt_score) & tmp$pt_score>43],collapse=",")
	}
	sumTable<-merge(sumTable,tfbsTable,all.x=T)
	sumTable[,15:22]<-replace(sumTable[,15:22],is.na(sumTable[,15:22]),"")
	hSgeneTAD<-read.delim("data/hSubs-genes-in-TAD.tsv",stringsAsFactors=F,header=F)[,c(1,2,4)]
	colnames(hSgeneTAD)<-c("Chr","Pos","gene")
	hSgeneTAD<-hSgeneTAD[hSgeneTAD$gene!=".",]
	hSgeneTAD$gene<-toupper(hSgeneTAD$gene)
	hSgeneTAD<-merge(hSgeneTAD,TS3_Reilly[,c(2,6)],by.x="gene",by.y="external_gene_name",all.x=T)
	hSgeneLst<-hSgeneTAD %>% nest(-c(Chr,Pos))
	hSgeneTable<-data.frame(Chr=hSgeneLst$Chr,Pos=hSgeneLst$Pos,expr_genes_in_TAD=NA,expr_genes_WGCNA=NA)
	for(i in 1:nrow(hSgeneLst)){
		hSgeneTable$expr_genes_in_TAD[i]<-paste(hSgeneLst$data[[i]]$gene,collapse=",")
		hSgeneTable$expr_genes_WGCNA[i]<-paste(hSgeneLst$data[[i]]$ModuleAssignment,collapse=",")
	}
	sumTable<-merge(sumTable,hSgeneTable,all.x=T)
	phyP<-read.delim("data/hSubs-hg19phyloP.tsv",stringsAsFactors=F,header=F)
	colnames(phyP)<-c("Chr","Pos","phyloP")
	sumTable<-merge(sumTable,phyP,all.x=T)
	HARhic<-read.delim("data/Won_et_al_HAR_targets_named.tsv",stringsAsFactors=F)
	HIC<-read.delim("data/GZ_GE-geneID.tsv",stringsAsFactors=F)
	colnames(HIC)<-c("Enhancer","hic.GZ.promoter")
	HIC<-rbind(HIC[HIC$Enhancer!="HACNS_711",],HARhic[HARhic$Enhancer!="HGE_1901",-(2:4)])
	sumTable<-merge(sumTable,HIC,all.x=T)
	rm(list=c("addTmp","cnt","dummy","dummyBED","HARhic","HIC","hSgeneLst","hSgeneTable","hSgeneTAD","hSs","hSubFeat","i","intTmp","n","nr","sumH2",
		"sumH3","sumLst","sumTAB","tfbs","tfbsLst","tfbsTable"))
	
	# Some BED file output #
	ahSub<-sumTable[!is.na(sumTable$average_add_effect_size),]
	ihSub<-sumTable[is.na(sumTable$average_add_effect_size),]
	hsSub<-ahSub[ahSub$average_add_effect_size>0,2:3]
	hsSub<-rbind(hsSub,ihSub[ihSub$average_intera_effect_size>0,2:3])
	hsSub$End<-hsSub$Pos
	hsSub$Pos<-hsSub$Pos-1
	write.table(hsSub,"data/hsBiased-hSubs.bed",quote=F,row.names=F,col.names=F,sep="\t") # Human-biased hSubs
	ptSub<-ahSub[ahSub$average_add_effect_size<0,2:3]
	ptSub<-rbind(ptSub,ihSub[ihSub$average_intera_effect_size<0,2:3])
	ptSub$End<-ptSub$Pos
	ptSub$Pos<-ptSub$Pos-1
	write.table(ptSub,"data/ptBiased-hSubs.bed",quote=F,row.names=F,col.names=F,sep="\t") # Chimpanzee-biased hSubs
	rm(list=c("ahSub","hsSub","ihSub","ptSub"))
	
	ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
	transTab<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=ensembl)
	
	slation<-double(nrow(sumTable))
	modAssign<-double(nrow(sumTable))
	for(i in 1:nrow(sumTable)){
		if(sumTable$hic.GZ.promoter[i]!="" & !is.na(sumTable$hic.GZ.promoter[i])){
			tmp<-unlist(strsplit(sumTable$hic.GZ.promoter[i],","))
			tmp2<-double(length(tmp))
			amp2<-double(length(tmp))
			for(j in 1:length(tmp)){
				tmp3<-transTab$external_gene_name[transTab$ensembl_gene_id==tmp[j]]
				amp3<-TS3_Reilly$ModuleAssignment[TS3_Reilly$Gene==tmp[j]]
				if(length(tmp3)==1){
					tmp2[j]<-tmp3
				}else{
					tmp2[j]<-NA
				}
				if(length(amp3)==1){
					amp2[j]<-amp3
				}
				else{
					amp2[j]<-NA
			}}
			slation[i]<-paste(tmp2,collapse=",")
			modAssign[i]<-paste(amp2,collapse=",")
	}}
	slation<-replace(slation,slation==0,NA)
	modAssign<-replace(modAssign,modAssign==0,NA)
	sumTable$hic.GZ.promoter<-slation
	sumTable<-cbind(sumTable,hic.GZ.promoter_WGCNA=modAssign,stringsAsFactors=F)
	PsychHIC<-unique(read.delim("data/hiccups_loopedPromoters_nonsubsampled.txt",stringsAsFactors=F,header=F)[,26:27])
	colnames(PsychHIC)<-c("Enhancer","hic.NPC")
	PsychHIC<-data.frame(PsychHIC %>% group_by(Enhancer) %>% summarize(hic.NPC=paste(hic.NPC,collapse=",")))
	sumTable<-merge(sumTable,PsychHIC,all.x=T)
	write.table(sumTable,"data/hSub-summary.tsv",quote=F,row.names=F,sep="\t")
	rm(list=c("amp2","amp3","ensembl","i","j","modAssign","phyP","PsychHIC","slation","tmp","tmp2","tmp3","transTab","TS3_Reilly"))
}

q(save="no")
