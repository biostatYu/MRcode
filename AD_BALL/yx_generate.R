library(data.table)
exp = c("AD_Asthma_adultonset_2019AJHG_maf_GRCh37",
	"AD_Asthma_childonset_2019AJHG_maf_GRCh37",
	"AD_Crohn_disease_2017LangeNG",
	"AD_Graves_disease_2021Sakaue",
	"AD_Hashimoto_Disease_2021Sakaue",
	"AD_Hypothyroidism_2021Sakaue",
	"AD_Inflammatory_bowel_disease_2017LangeNG",
	"AD_Multiple_Sclerosis_2018_qc",
	"AD_PBC_2021Cordell",
	"AD_Primary_sclerosing_cholangitis2016",
	"AD_Rheumatoid_Arthritis_European_2022_NG",
	"AD_Rheumatoid_Arthritis_pos_European_2022_NG",
	"AD_Systemic_Lupus_Erythematosus_2015_NG",
	"AD_systemic_sclerosis",
	"AD_T1D_2020_Diabetes",
	"AD_Ulcerative_Colitis_2017LangeNG",
	"AD_vitiligo_NG_2016")
RES = NULL
for(i in 1:length(exp)){
	setwd("E:\\data\\AD")
	dd0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))[,1:10]
	print(head(dd0))
	print(dim(dd0))
	colnames(dd0) = c("SNP","CHR","POS","INC_ALLELE","DEC_ALLELE","BETA","SE","Z","P","N")
	setwd("E:\\data\\AD\\index_snp")
	SNP = as.character(fread(paste(exp[i],"_5E-8_0.001_10000_clump.txt",sep=""),header=TRUE)$SNP)
	dd = dd0[which(dd0$SNP%in%SNP),]
	RES = rbind(RES,data.frame(Exp = exp[i],dd))
	print(i)
}

setwd("E:\\data\\AD\\index_snp")
write.table(RES,"index_SNP_AD.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)



RES = NULL
exp = c("Blood_B_ALL_2019NC")
for(i in 1:length(exp)){
	setwd("E:\\data\\blood")
	dd0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))
	setwd("E:\\data\\blood\\index_snp")
	SNP = as.character(fread(paste(exp[i],"_5E-8_0.001_10000_clump.txt",sep=""),header=TRUE)$SNP)
	dd = dd0[which(dd0$SNP%in%SNP),]
	RES = rbind(RES,data.frame(Exp = exp[i],dd))
	print(i)
}

setwd("E:\\data\\blood\\index_snp")
write.table(RES,"index_SNP_B_ALL.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)





setwd("E:\\data\\AD\\index_snp")
index_snp = data.frame(read.table("index_SNP_AD.txt",head=T))
setwd("E:\\data\\blood")
mx0 = data.frame(fread(paste("Blood_B_ALL_2019NC.txt.gz",sep = ""),sep='\t',head=T))
index = which(duplicated(mx0$SNP))
if(length(index)>0){mx0 = mx0[-index,]}
yx0 = merge(index_snp,mx0,by="SNP",all.x=TRUE)
setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Results\\MR")
write.table(yx0,paste("yx_AD_B_ALL.txt",sep=""),quote=FALSE,sep="\t",col.names=T,row.names=FALSE)



setwd("E:\\data\\blood\\index_snp")
index_snp = data.frame(read.table("index_SNP_B_ALL.txt",head=T))
library(data.table)
exp = c("AD_Asthma_adultonset_2019AJHG_maf_GRCh37",
	"AD_Rheumatoid_Arthritis_European_2022_NG")
RES = NULL
for(i in 1:length(exp)){
	setwd("E:\\data\\AD")
	dd0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))[,1:10]
	colnames(dd0) = c("SNP","CHR","POS","INC_ALLELE","DEC_ALLELE","BETA","SE","Z","P","N")
	yx0 = merge(index_snp, dd0,by="SNP",all.x=TRUE)
	RES = rbind(RES,data.frame(Out = exp[i],yx0))
	print(i)
}
setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Results\\MR")
write.table(RES,paste("yx_B_ALL_AD.txt",sep=""),quote=FALSE,sep="\t",col.names=T,row.names=FALSE)

