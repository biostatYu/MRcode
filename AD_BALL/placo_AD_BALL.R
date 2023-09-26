args=as.numeric(commandArgs(TRUE))
REP=args[1]
library(data.table)
library(pbapply)
source("/software/PLACO.R")
ax = c("AD_Asthma_adultonset_2019AJHG_maf_GRCh37",
	"AD_Crohn_disease_2017LangeNG",
	"AD_Inflammatory_bowel_disease_2017LangeNG",
	"AD_PBC_2021Cordell",
	"AD_Multiple_Sclerosis_2018_qc",
	"AD_Rheumatoid_Arthritis_European_2022_NG",
	"AD_Hypothyroidism_2021Sakaue")
DD = c("Blood_B_ALL_2019NC")
tmp = data.frame(x1 = rep(1:length(ax), each=length(DD)),x2 = rep(1:length(DD), times=length(ax)))
for(ii in REP){
	i = tmp[ii,1];
	j = tmp[ii,2];
	print("GO")
	covid.data = data.frame(fread(paste("/data/Summarydata/autoimmune/",ax[i],".txt.gz",sep=""), head=T))
	m.data = data.frame(fread(paste("/data/Summarydata/blood/data_qc/",DD[j],".txt.gz",sep=""),header =TRUE))
	print("Read done")
	yx = merge(covid.data,m.data,by = "SNP")
	print("Merge done")
	index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
	if (length(index1)<length(yx$INC_ALLELE.x))
	{yx$Z.y[-index1]=yx$Z.y[-index1]*(-1)}
	print("Step 1")
	yx$Z2.x = yx$Z.x^2
	yx$Z2.y = yx$Z.y^2
	ind.x = which(yx$Z2.x>80)
	if(length(ind.x)>0){yx = yx[-ind.x,]}
	ind.y = which(yx$Z2.y>80)	
	if(length(ind.y)>0){yx = yx[-ind.y,]}
	Z.matrix = as.matrix(yx[,c("Z.x","Z.y")])
	P.matrix = as.matrix(yx[,c("P.x","P.y")])
	colnames(Z.matrix) <- paste("Z",1:2,sep="")
	colnames(P.matrix) <- paste("P",1:2,sep="")
	print("Step 2")
	VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)
	p = dim(yx)[1]
	print("Going")
	out <- pbsapply(1:p, function(i)placo(Z=Z.matrix[i,], VarZ=VarZ))
	PLACO = t(out)
	PLACO.1 = data.frame(yx,PLACO)
	PLACO.1$N_com = (PLACO.1$N.x+PLACO.1$N.y)/2
	PLACO.2 = PLACO.1[,c("SNP","CHR.x","POS.x","INC_ALLELE.x","DEC_ALLELE.x","P.x","P.y","T.placo","p.placo","N_com")]
	colnames(PLACO.2) = c("SNP","CHR","POS","INC_ALLELE","DEC_ALLELE","P.x","P.y","T.placo","P.placo","N")
	setwd("/data/Summarydata/placo/results")
	fwrite(PLACO.2,paste("AD_BALL_",ax[i],"_",DD[j],".txt.gz",sep=""),row.names= F,col.names= T,sep="\t",quote = F)
	print("Done")
}
