args=as.numeric(commandArgs(TRUE))
REP=args[1]
library(HDL)
library(data.table)
LD.path <- "/data/ref/HDL"
ax = c( "AD_Asthma_2022CellGenom_maf_GRCh37_qc",
	"AD_Asthma_adultonset_2019AJHG_maf_GRCh37",
	"AD_Asthma_childonset_2019AJHG_maf_GRCh37",
	"AD_Crohn_disease_2017LangeNG",
	"AD_Graves_disease_2021Sakaue",
	"AD_Hashimoto_Disease_2021Sakaue",
	"AD_Hypothyroidism_2021Sakaue",
	"AD_Inflammatory_bowel_disease_2017LangeNG",
	"AD_Multiple_Sclerosis_2018_qc",
	"AD_PBC_2021Cordell",
	"AD_Primary_sclerosing_cholangitis2016",
	"AD_Rheumatoid_Arthritis_European_2014_nature",
	"AD_Rheumatoid_Arthritis_European_2022_NG",
	"AD_Rheumatoid_Arthritis_pos_European_2022_NG",
	"AD_Systemic_Lupus_Erythematosus_2015_NG",
	"AD_systemic_sclerosis",
	"AD_T1D_2020_Diabetes",
	"AD_Ulcerative_Colitis_2017LangeNG",
	"AD_vitiligo_NG_2016")
DD = c("Blood_B_ALL_2019NC")
tmp = data.frame(x1 = rep(1:length(ax), each=length(DD)),x2 = rep(1:length(DD), times=length(ax)))
for(ii in REP){
	i = tmp[ii,1];
	j = tmp[ii,2];
	print("GO")
	data.1 = data.frame(fread(paste("/data/Summarydata/autoimmune/",ax[i],".txt.gz",sep=""), head=T))
	data.2 = data.frame(fread(paste("/data/Summarydata/blood/data_qc/",DD[j],".txt.gz",sep=""),header =TRUE))

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
	setwd("/data/Summarydata/ldsc/HDL_results/AD_BALL")
	fwrite(data.frame(RES.1),paste("CZ_",ax[i],"_",DD[j],".txt",sep=""),row.names= F,col.names= T,sep="\t",quote = F)


}
