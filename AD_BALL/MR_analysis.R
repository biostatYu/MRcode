
x.1 = c("AD_Asthma_adultonset_2019AJHG_maf_GRCh37",
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

library(data.table)
library(MendelianRandomization)
library(progress)
setwd("E:\\Cloud_disk\\program")
source("MRfunction.R")
pb <- progress_bar$new(total = 100)
RES = NULL
setwd("G:\\data\\summarydata")
DD0 = data.frame(fread(paste("Blood_B_ALL_2019NC.txt.gz",sep=""),header=TRUE))[,1:10]

m=length(x.1)
pb <- progress_bar$new(total = m)
rex = matrix(NA, m, 54)
for (j in 1:m){
try({
	setwd("G:\\data\\summarydata")
	ax0 = data.frame(fread(paste(x.1[j],".txt.gz",sep=""),header=TRUE))[,1:10]
	setwd("G:\\data\\summarydata\\C_data_autoimmune\\index_snp")
	snp = as.character(fread(paste(x.1[j],"_5E-8_0.001_10000_clump.txt",sep=""),header=TRUE)$SNP)
	yx0 = MR_process_impute(ax0,DD0,snp)
	yx = yx0$yx

	index4 = which(yx$CHR.x == 6 & yx$POS.x >= 28477797 & yx$POS.x <= 33448354)
	if (length(index4)>0) {yx=yx[-index4,]}

	yx=na.omit(data.frame(yx$BETA.x,yx$SE.x,yx$BETA.y,yx$SE.y,yx$P.y))
	colnames(yx)=c("BETA.x","SE.x","BETA.y","SE.y","P.y")
	yx=as.data.frame(yx)
	
	fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
	cooksd <- cooks.distance(fit)
	sample_size <- length(yx$BETA.y)
	index = which(cooksd > (4/sample_size))
	yx$SNP[index]
	if(length(index)>0){yx = yx[-index,]}

	index2=which(yx$P.y<=(0.05/length(yx$P.y)))
	if (length(index2)>0) {yx=yx[-index2,]}
	index3 = which(abs(yx$BETA.y)>100)
	if (length(index3)>0) {yx=yx[-index3,]}
	index4 = which(yx$BETA.y == 0)
	if (length(index4)>0) {yx=yx[-index4,]}
	n = dim(yx)[1]
	N = mean(ax0$N)
	R2 = yx$BETA.x*yx$BETA.x/(yx$BETA.x*yx$BETA.x+yx$SE.x*yx$SE.x*N)
	F = R2*(N-2)/(1-R2)
	Res = MR_sum(yx,result="OR")$result
	rex[j,] = c(j,sum(R2),min(F),n,Res[1,1:6],Res[2,1:6],Res[5,1:6],Res[6,1:6],Res[3,1:6],Res[4,1:6],Res[7,1:6],Res[8,1:6],Res[1,7:8])
	pb$tick()
	Sys.sleep(1/m)
	})
}
rexx = data.frame(x.1, rex)
colnames(rexx) = c("Exp","ID","PVE","F_min","n_snp",
		"IVW_fix","SE","OR","Down","Up","Pval",
		"IVW_random","SE","OR","Down","Up","Pval",
		"weighted_Mode","SE","OR","Down","Up","Pval",
		"weighted_Median","SE","OR","Down","Up","Pval",
		"Egger_slope","SE","OR","Down","Up","Pval",
		"Egger_intercept","SE","OR","Down","Up","Pval",
		"DIVW","SE","OR","Down","Up","Pval",
		"MR-RAPS","SE","OR","Down","Up","Pval","Heter.est","Heter.Pval")
setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Results\\MR")
write.table(rexx,paste("Res_sensitive_AD_B_ALL.txt",sep=""),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)





x.1 = c("AD_Asthma_adultonset_2019AJHG_maf_GRCh37", "AD_Rheumatoid_Arthritis_European_2022_NG")

library(data.table)
library(MendelianRandomization)
library(progress)
setwd("E:\\Cloud_disk\\program")
source("MRfunction.R")
pb <- progress_bar$new(total = 100)
RES = NULL
setwd("G:\\data\\summarydata")
ax0 = data.frame(fread(paste("Blood_B_ALL_2019NC.txt.gz",sep=""),header=TRUE))[,1:10]
setwd("G:\\data\\summarydata\\C_data_blood\\index_snp")
snp = as.character(fread(paste("Blood_B_ALL_2019NC_5E-8_0.001_10000_clump.txt",sep=""),header=TRUE)$SNP)

m=length(x.1)
pb <- progress_bar$new(total = m)
rex = matrix(NA, m, 54)
for (j in 1:m){
try({
	setwd("G:\\data\\summarydata")
	DD0 = data.frame(fread(paste(x.1[j],".txt.gz",sep=""),header=TRUE))[,1:10]
	yx0 = MR_process_impute(ax0,DD0,snp)
	yx = yx0$yx

	index4 = which(yx$CHR.x == 6 & yx$POS.x >= 28477797 & yx$POS.x <= 33448354)
	if (length(index4)>0) {yx=yx[-index4,]}

	yx=na.omit(data.frame(yx$BETA.x,yx$SE.x,yx$BETA.y,yx$SE.y,yx$P.y))
	colnames(yx)=c("BETA.x","SE.x","BETA.y","SE.y","P.y")
	yx=as.data.frame(yx)
	
	fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
	cooksd <- cooks.distance(fit)
	sample_size <- length(yx$BETA.y)
	index = which(cooksd > (4/sample_size))
	yx$SNP[index]
	if(length(index)>0){yx = yx[-index,]}

	index2=which(yx$P.y<=(0.05/length(yx$P.y)))
	if (length(index2)>0) {yx=yx[-index2,]}
	index3 = which(abs(yx$BETA.y)>100)
	if (length(index3)>0) {yx=yx[-index3,]}
	index4 = which(yx$BETA.y == 0)
	if (length(index4)>0) {yx=yx[-index4,]}
	n = dim(yx)[1]
	N = mean(ax0$N)
	R2 = yx$BETA.x*yx$BETA.x/(yx$BETA.x*yx$BETA.x+yx$SE.x*yx$SE.x*N)
	F = R2*(N-2)/(1-R2)
	Res = MR_sum(yx,result="OR")$result
	rex[j,] = c(j,sum(R2),min(F),n,Res[1,1:6],Res[2,1:6],Res[5,1:6],Res[6,1:6],Res[3,1:6],Res[4,1:6],Res[7,1:6],Res[8,1:6],Res[1,7:8])
	pb$tick()
	Sys.sleep(1/m)
	})
}
rexx = data.frame(x.1, rex)
colnames(rexx) = c("Exp","ID","PVE","F_min","n_snp",
		"IVW_fix","SE","OR","Down","Up","Pval",
		"IVW_random","SE","OR","Down","Up","Pval",
		"weighted_Mode","SE","OR","Down","Up","Pval",
		"weighted_Median","SE","OR","Down","Up","Pval",
		"Egger_slope","SE","OR","Down","Up","Pval",
		"Egger_intercept","SE","OR","Down","Up","Pval",
		"DIVW","SE","OR","Down","Up","Pval",
		"MR-RAPS","SE","OR","Down","Up","Pval","Heter.est","Heter.Pval")
setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Results\\MR")
write.table(rexx,paste("Res_sensitive_B_ALL_AD.txt",sep=""),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)


