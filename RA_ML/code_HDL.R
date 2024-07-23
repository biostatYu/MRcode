args=as.numeric(commandArgs(TRUE))
REP=args[1]
library(HDL)
library(data.table)
LD.path <- "/home/data/ref/HDL"
ax = c(	"AD_Rheumatoid_Arthritis_pos_European_2022_NG")
DD = c(	"Blood_Pheweb_EUR_Malignant_lymphoma","Blood_finnUKBB_C3_CLL_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_C3_DLBCL_EXALLC_GRCh37_qc","Blood_finnUKBB_CD2_FOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_HODGKIN_LYMPHOMA_EXALLC_GRCh37_qc","Blood_finnUKBB_CD2_NONFOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
	"Blood_finnUKBB_CD2_NONHODGKIN_NAS_EXALLC_GRCh37_qc","Blood_finnUKBB_CD2_TNK_LYMPHOMA_EXALLC_GRCh37_qc")
tmp = data.frame(x1 = rep(1:length(ax), each=length(DD)),x2 = rep(1:length(DD), times=length(ax)))
for(ii in REP){
	i = tmp[ii,1];
	j = tmp[ii,2];
	print("GO")
	data.1 = data.frame(fread(paste("/home/data/summarydata/",ax[i],".txt.gz",sep=""), head=T))
	data.2 = data.frame(fread(paste("/home/data/summarydata/",DD[j],".txt.gz",sep=""),header =TRUE))

	index = which(data.1$CHR == 6 & data.1$POS >= 28477797 & data.1$POS <= 33448354)
	if (length(index)>0) {data.1.1=data.1[-index,]}					
	index = which(data.2$CHR == 6 & data.2$POS >= 28477797 & data.2$POS <= 33448354)
	if (length(index)>0) {data.2.1=data.2[-index,]}
	data.1.1 = data.1.1[,c("SNP","INC_ALLELE","DEC_ALLELE","N","BETA","SE","Z")]
	data.2.1 = data.2.1[,c("SNP","INC_ALLELE","DEC_ALLELE","N","BETA","SE","Z")]
	colnames(data.1.1) = c("SNP","A1","A2","N","b","se","Z")
	colnames(data.2.1) = c("SNP","A1","A2","N","b","se","Z")
	print("Read done")
	print("Going")
	res.HDL.1 <- HDL.rg(data.1.1, data.2.1, LD.path)
	print(res.HDL.1)
	RES.1 = c(res.HDL.1$rg, res.HDL.1$rg.se, res.HDL.1$P)
	setwd("/home/yxh/data/Summarydata/A_result_LDSC/HDL/RA_ML")
	fwrite(data.frame(RES.1),paste("RA_ML_",ax[i],"_",DD[j],".txt",sep=""),row.names= F,col.names= T,sep="\t",quote = F)
	print("noMHC done")
}
