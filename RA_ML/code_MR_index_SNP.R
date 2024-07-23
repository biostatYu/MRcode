library(data.table)
exp = c("AD_Rheumatoid_Arthritis_European_2022_NG",
	"AD_Rheumatoid_Arthritis_pos_European_2022_NG")
RES = NULL
for(i in 1:length(exp)){
	try({
	setwd("/home/data/summarydata")
	dd0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))[,1:10]
	setwd("/home/data/summarydata/index_snp")
	SNP = as.character(fread(paste(exp[i],"_5E-8_0.001_10000_clump.txt",sep=""),header=TRUE)$SNP)
	dd = dd0[which(dd0$SNP%in%SNP),]
	RES = rbind(RES,data.frame(Exp = exp[i],dd))
	print(i)
	})
}

setwd("/home/yxh/data/Summarydata/A_result_MR/")
write.table(RES,"index_snp_RA.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)




source("/home/yxh/software/MRfunction_noraps.R")
setwd("/home/yxh/data/Summarydata/A_result_MR/")
index_snp = data.frame(read.table("index_snp_RA.txt",head=T))
library(data.table)
exp = c("Blood_Pheweb_EUR_Malignant_lymphoma",
	"Blood_finnUKBB_C3_CLL_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_C3_DLBCL_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_FOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_HODGKIN_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_TNK_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_NONFOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_NONHODGKIN_NAS_EXALLC_GRCh37_qc")
RES = NULL
for(i in 1:length(exp)){
	setwd("/home/data/summarydata")
	dd0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))[,1:10]
	dd1 = dd0[which(dd0$SNP%in%index_snp$SNP),]
	yx0 = MR_process_impute(index_snp, dd1, as.character(index_snp$SNP))$yx
	RES = rbind(RES,data.frame(Out = exp[i],yx0))
	print(i)
}

setwd("/home/yxh/data/Summarydata/A_result_MR/RA_ML")
write.table(RES,paste("yx_RA_ML.txt",sep=""),quote=FALSE,sep="\t",col.names=T,row.names=FALSE)







library(data.table)
exp = c("Blood_Pheweb_EUR_Malignant_lymphoma",
	"Blood_finnUKBB_C3_CLL_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_C3_DLBCL_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_FOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_HODGKIN_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_TNK_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_NONFOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_NONHODGKIN_NAS_EXALLC_GRCh37_qc")
RES = NULL
for(i in 1:length(exp)){
	try({
	setwd("/home/data/summarydata")
	dd0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))[,1:10]
	setwd("/home/data/summarydata/index_snp")
	SNP = as.character(fread(paste(exp[i],"_5E-8_0.001_10000_clump.txt",sep=""),header=TRUE)$SNP)
	dd = dd0[which(dd0$SNP%in%SNP),]
	RES = rbind(RES,data.frame(Exp = exp[i],dd))
	print(i)
	})
}

setwd("/home/yxh/data/Summarydata/A_result_MR/")
write.table(RES,"index_snp_ML.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)




source("/home/yxh/software/MRfunction_noraps.R")
setwd("/home/yxh/data/Summarydata/A_result_MR/")
index_snp = data.frame(read.table("index_snp_ML.txt",head=T))
library(data.table)
exp = c("AD_Rheumatoid_Arthritis_European_2022_NG",
	"AD_Rheumatoid_Arthritis_pos_European_2022_NG")
RES = NULL
for(i in 1:length(exp)){
	setwd("/home/data/summarydata")
	dd0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))[,1:10]
	dd1 = dd0[which(dd0$SNP%in%index_snp$SNP),]
	yx0 = MR_process_impute(index_snp, dd1, as.character(index_snp$SNP))$yx
	RES = rbind(RES,data.frame(Out = exp[i],yx0))
	print(i)
}

setwd("/home/yxh/data/Summarydata/A_result_MR/RA_ML")
write.table(RES,paste("yx_ML_RA.txt",sep=""),quote=FALSE,sep="\t",col.names=T,row.names=FALSE)

