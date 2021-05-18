library(tidyr)

## Load and prepare our data for comparisons
# If any files are missing, create by running 'analysis-pipeline-rd1.R'
MPRAdata1<-read.delim("data/COUNTS_ANNOTATED.txt.gz",stringsAsFactors=F) # raw data
MPRAdata1<-MPRAdata1[!is.na(MPRAdata1$Tag) & rowSums(MPRAdata1[,16:19])>0,]

MPRAfrag1<-read.delim("data/MPRAfrag1.tsv",stringsAsFactors=F) # Summarized by fragment
MPRAmeas1<-MPRAfrag1[MPRAfrag1$noTags>11,] # measured fragments (>11 barcodes)

MPRAfold1<-read.delim("data/MPRAfold1.tsv",stringsAsFactors=F) # "Folded" data summarized by fragment
MPRAactive1_1<-read.delim("data/MPRAactive1_1.tsv",stringsAsFactors=F)
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
MPRAsum1<-read.delim("data/MPRAsum1.tsv",stringsAsFactors=F)
MPRAactivities1<-cbind(MPRAsum1[,1:8],
	Rep12_activity=MPRAsum1$Rep_1_2_cDNA_Tag_Counts-MPRAsum1$Rep_1_2_pDNA_Tag_Counts,
	Rep13_activity=MPRAsum1$Rep_1_3_cDNA_Tag_Counts-MPRAsum1$Rep_1_3_pDNA_Tag_Counts,
	Rep22_activity=MPRAsum1$Rep_2_2_cDNA_Tag_Counts-MPRAsum1$Rep_2_2_pDNA_Tag_Counts,
	Rep23_activity=MPRAsum1$Rep_2_3_cDNA_Tag_Counts-MPRAsum1$Rep_2_3_pDNA_Tag_Counts,
	Lin1_activity=MPRAsum1$Lin1_cDNA_Tag_Counts-MPRAsum1$Lin1_pDNA_Tag_Counts,
	Lin2_activity=MPRAsum1$Lin2_cDNA_Tag_Counts-MPRAsum1$Lin2_pDNA_Tag_Counts)
MPRAactivities1$Alignment<-sub('_Chimp','',MPRAactivities1$Alignment)
load("data/fragnamelist.RData")
Rep12_nestAct<-MPRAactivities1[MPRAactivities1$Alignment %in% names2List 
	& MPRAactivities1$Alignment %in% MPRAfoldActive1$Alignment,c(2,6,9)] %>% nest(-Alignment, Species)

MPRAfold1$Rep12_actDiff<-MPRAfold1$Rep12_actMedian_hs-MPRAfold1$Rep12_actMedian_pt
MPRAfold1$Rep13_actDiff<-MPRAfold1$Rep13_actMedian_hs-MPRAfold1$Rep13_actMedian_pt
MPRAfold1$Rep22_actDiff<-MPRAfold1$Rep22_actMedian_hs-MPRAfold1$Rep22_actMedian_pt
MPRAfold1$Rep23_actDiff<-MPRAfold1$Rep23_actMedian_hs-MPRAfold1$Rep23_actMedian_pt
MPRAfold1$direction<-NULL
MPRAfold1$direction<-replace(MPRAfold1$direction,MPRAfold1$Rep12_actDiff >0 
	& MPRAfold1$Rep13_actDiff >0 & MPRAfold1$Rep22_actDiff >0 & MPRAfold1$Rep23_actDiff >0,'hs')
MPRAfold1$direction<-replace(MPRAfold1$direction,MPRAfold1$Rep12_actDiff <0 
	& MPRAfold1$Rep13_actDiff <0 & MPRAfold1$Rep22_actDiff <0 & MPRAfold1$Rep23_actDiff <0,'pt')

MPRAactive1_2<-read.delim("data/MPRAactive1_2.tsv",stringsAsFactors=F)
nrow(MPRAactive1_2)
# Active fragment pairs in our original results (human, chimp or both): 3202
nrow(MPRAactive1_2[grepl("_Chimp",MPRAactive1_2$Alignment)==F,])
# Human active fragment pairs in our original results: 1550
nrow(MPRAactive1_2[grepl("_Chimp",MPRAactive1_2$Alignment),])
# Chimp active fragment pairs in our original results: 1652

MPRAda1_2<-read.delim("data/diffActive1_2.tsv",stringsAsFactors=F)
nrow(MPRAda1_2)
# DA fragment pairs in our original results: 673

## MPRAnalyze ##
# The method was run in 62 batches, the results of which are loaded and compared below.

 # Activity test #
human<-read.delim("data/MPRAna_human.tsv",stringsAsFactors=F)
chimp<-read.delim("data/MPRAna_chimp.tsv",stringsAsFactors=F)

human_sig<-human[p.adjust(human$pval.mad,method="BH")<.05,]
# Active human fragments: 3891
chimp_sig<-chimp[p.adjust(chimp$pval.mad,method="BH")<.05,]
# Active chimp fragments: 1557
# Total active fragments: 5448
nrow(human_sig[rownames(human_sig) %in% MPRAactive1_2$Alignment,])
# Active human fragments that were also active in our results: 1033
nrow(chimp_sig[paste(rownames(chimp_sig),"Chimp",sep="_") %in% MPRAactive1_2$Alignment,])
# Active chimp fragments that were also active in our results: 25
MPRAnalyze_act_vector<-c(rownames(human_sig),paste(rownames(chimp_sig),"Chimp",sep="_"))
# Active fragments: 5448
length(MPRAnalyze_act_vector[MPRAnalyze_act_vector %in% MPRAactive1_2$Alignment])
# Active fragments that were also active in our results: 1058

active<-merge(human,chimp,by=0)
colnames(active)<-sub(".x","_hs",colnames(active))
colnames(active)<-sub(".y","_pt",colnames(active))
active_sig<-active[p.adjust(active$pval.mad_hs,method="BH")<.05 | p.adjust(active$pval.mad_pt,method="BH")<.05,]
# Active in either allele (not counted twice if active in both): 5103
nrow(active[p.adjust(active$pval.mad_hs,method="BH")<.05 & p.adjust(active$pval.mad_pt,method="BH")<.05,])
# Active in both alleles: 345

 # Differential activity test #
da<-read.delim("data/MPRAna_da.tsv",stringsAsFactors=F)

active_da<-da[rownames(da) %in% Rep12_nestAct$Alignment,]
MPRAnalyze_da_vector<-rownames(active_da)[p.adjust(active_da$pval,method="BH")<.05]
# DA fragment pairs that we also tested for differential activity: 2381

da_sig<-da[p.adjust(da$pval,method="BH")<.05,]
# DA fragments: 6949
nrow(da_sig[abs(da_sig$logFC)>.2,])
# DA fragment pairs with a logFC > 0.2: 2963
nrow(da_sig[rownames(da_sig) %in% active_sig$Row.names,])
# DA among active (according to MPRAnalyze) fragments: 2275

nrow(da_sig[rownames(da_sig) %in% MPRAda1_2$Alignment,])
# DA fragment pairs that were also DA in our results: 451

nrow(da_sig)-sum(table(MPRAfold1$direction[MPRAfold1$Alignment %in% rownames(da_sig)]))
# DA fragment pairs that disagree in the direction of bias: 4478
nrow(da_sig[rownames(da_sig) %in% active_sig$Row.names,])-sum(table(MPRAfold1$direction[MPRAfold1$Alignment %in% rownames(da_sig[rownames(da_sig) %in% active_sig$Row.names,])]))
# DA fragment pairs (among active pairs) that disagree in the direction of bias: 1064
MPRAfold1_dummy<-MPRAfold1[MPRAfold1$Alignment %in% rownames(da_sig[abs(da_sig$logFC > .2) & rownames(da_sig) %in% Rep12_nestAct$Alignment,]),]
nrow(MPRAfold1_dummy[!is.na(MPRAfold1_dummy$direction),])
# DA fragment pairs that agree in direction and show logFC > 0.2: 845

## mpralm ##
library(mpra)

 # Differential activity test #
MPRAsummed<-aggregate(MPRAdata1[,12:19],by=list(Category=MPRAdata1$Alignment),FUN=sum)
colnames(MPRAsummed)[1]<-"Alignment"
MPRAsummed$Species<-"hg19"
MPRAsummed$Species<-replace(MPRAsummed$Species,grepl("_Chimp",MPRAsummed$Alignment),"PanTro2")
MPRAsummedFold<-MPRAsummed
MPRAsummedFold$Alignment<-sub("_Chimp","",MPRAsummedFold$Alignment)
MPRAsummedFold<-merge(MPRAsummedFold[MPRAsummedFold$Species=='hg19',-10],
	MPRAsummedFold[MPRAsummedFold$Species=='PanTro2',1:9],by='Alignment')
colnames(MPRAsummedFold)<-sub('.x','_hs',colnames(MPRAsummedFold),fixed=T)
colnames(MPRAsummedFold)<-sub('.y','_pt',colnames(MPRAsummedFold),fixed=T)

rna<-as.matrix(cbind(MPRAsummedFold[,c(2:5,10:13)]))
dna<-as.matrix(MPRAsummedFold[,c(6:9,14:17)])
MPRAset<-MPRASet(DNA=dna,RNA=rna,barcode=NULL,eid=MPRAsummedFold$Alignment,eseq=NULL)
rownames(MPRAset)<-MPRAsummedFold$Alignment

design<-data.frame(intcpt=1,Chimp=grepl("_pt",colnames(MPRAset)))
block_vector<-rep(1:4,2)
mpralm_allele_fit<-mpralm(object=MPRAset,design=design,aggregate="none",normalize=T,
	block=block_vector,model_type="corr_groups",plot=T)
toptab_allele<-topTable(mpralm_allele_fit,coef=2,number=Inf)
nrow(toptab_allele[p.adjust(toptab_allele$P.Value,method="BH")<.05,])
# DA fragment pairs from whole dataset: 5268

active_toptab<-toptab_allele[rownames(toptab_allele) %in% Rep12_nestAct$Alignment,]
mpralm_da_vector<-rownames(active_toptab)[p.adjust(active_toptab$P.Value,method="BH")<.05]
# DA fragment pairs that we also tested for differential activity: 1764
nrow(MPRAda1_2[MPRAda1_2$Alignment %in% rownames(active_toptab[p.adjust(active_toptab$P.Value,method="BH")<.05,]),])
# DA fragment pairs from our results that are also DA in mpralm: 672

MPRAfold1_dummy<-MPRAfold1[MPRAfold1$Alignment %in% rownames(active_toptab[p.adjust(active_toptab$P.Value,method="BH")<.05 & abs(active_toptab$logFC > .2),]),]
nrow(MPRAfold1_dummy[!is.na(MPRAfold1_dummy$direction),])
# DA fragment pairs that agree in direction and show logFC > 0.2: 1135

nrow(toptab_allele[p.adjust(toptab_allele$P.Value,method="BH")<.05,])-
	sum(table(MPRAfold1$direction[MPRAfold1$Alignment %in% rownames(toptab_allele[p.adjust(toptab_allele$P.Value,method="BH")<.05,])]))
# DA fragment pairs that disagree in the direction of bias: 1457
nrow(active_toptab[p.adjust(active_toptab$P.Value,method="BH")<.05,])-
	sum(table(MPRAfold1$direction[MPRAfold1$Alignment %in% rownames(active_toptab[p.adjust(active_toptab$P.Value,method="BH")<.05,])]))
# DA fragment pairs that disagree in the direction of bias: 267

## Ryu et al. ##

 # Activity test #
Ryu_act_vector<-MPRAmeas1$Alignment[MPRAmeas1$Rep12_actMedian > quantile(MPRAmeas1$Rep12_actMedian,.75) & 
	MPRAmeas1$Rep13_actMedian > quantile(MPRAmeas1$Rep13_actMedian,.75) & 
	MPRAmeas1$Rep22_actMedian > quantile(MPRAmeas1$Rep22_actMedian,.75) & 
	MPRAmeas1$Rep23_actMedian > quantile(MPRAmeas1$Rep23_actMedian,.75)]
length(Ryu_act_vector)
# Active fragments: 7684
length(Ryu_act_vector[grepl("Chimp",Ryu_act_vector)==F])
# Active human fragments: 3697
length(Ryu_act_vector[grepl("Chimp",Ryu_act_vector)])
# Active chimp fragments: 3987
length(unique(sub("_Chimp","",Ryu_act_vector)))
# Active in either allele (not counted twice if active in both): 5717
nrow(MPRAactive1_2[MPRAactive1_2$Alignment %in% Ryu_act_vector[grepl("Chimp",Ryu_act_vector)==F],])
# Active human fragments that were also active in our results: 1529
nrow(MPRAactive1_2[MPRAactive1_2$Alignment %in% Ryu_act_vector[grepl("Chimp",Ryu_act_vector)],])
# Active chimp fragments that were also active in our results: 1626

ryuList<-read.delim("~/project/MPRA/190424_Ryu_etal/TableS3.tsv",stringsAsFactors=F)
ryu2<-ryuList[grepl("2x",ryuList$HAR),]
actHAR<-read.delim("~/project/MPRA//features/2xHAR_active.tsv",stringsAsFactors=F)
ryuAct<-merge(ryu2,actHAR)
map<-read.delim("~/project/MPRA/features/MPRA_HAR_intersect.bed",stringsAsFactors=F,header=F)

length(unique(map$V8[map$V8 %in% ryuAct$HAR]))
# Active HARs accoridng to us and Ryu: 51
length(unique(map$V8[ map$V4 %in% sub("_Chimp","",Ryu_act_vector) & map$V8 %in% ryuAct$HAR]))
# Active HARs according to us using Ryu's approach and Ryu: 51

## Venn overlaps ##

 # Activity test results #
length(MPRAnalyze_act_vector[MPRAnalyze_act_vector %in% Ryu_act_vector]) # 1865
nrow(MPRAactive1_2[MPRAactive1_2$Alignment %in% MPRAnalyze_act_vector,]) # 1058
nrow(MPRAactive1_2[MPRAactive1_2$Alignment %in% Ryu_act_vector,]) # 3155
nrow(MPRAactive1_2[MPRAactive1_2$Alignment %in% Ryu_act_vector & 
	MPRAactive1_2$Alignment %in% MPRAnalyze_act_vector,]) # 1043

 # Differential activity test results #
length(mpralm_da_vector[mpralm_da_vector %in% MPRAnalyze_da_vector]) # 1239
nrow(MPRAda1_2[MPRAda1_2$Alignment %in% MPRAnalyze_da_vector,]) # 492
nrow(MPRAda1_2[MPRAda1_2$Alignment %in% mpralm_da_vector,]) # 672
nrow(MPRAda1_2[MPRAda1_2$Alignment %in% MPRAnalyze_da_vector & 
	MPRAda1_2$Alignment %in% mpralm_da_vector,]) # 491

q(save="no")
