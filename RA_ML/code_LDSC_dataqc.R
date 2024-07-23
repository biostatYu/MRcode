args=as.numeric(commandArgs(TRUE))
REP=args[1]
library(data.table)
setwd("/home/data/ref/EUR")
bim = data.frame(fread((paste("EUR_1000Genome_phase3_all.bim",sep="")),head=F,sep="\t"))
dx = c("RA_ML_AD_Rheumatoid_Arthritis_pos_European_2022_NG_Blood_C3_DLBCL_EXALLC_meta_out_finnUKBB_GRCh37_qc",
	"RA_ML_AD_Rheumatoid_Arthritis_pos_European_2022_NG_Blood_CD2_HODGKIN_LYMPHOMA_EXALLC_meta_out_finnUKBB_GRCh37_qc",
	"RA_ML_AD_Rheumatoid_Arthritis_pos_European_2022_NG_Blood_CD2_NONFOLLICULAR_LYMPHOMA_EXALLC_meta_out_finnUKBB_GRCh37_qc",
	"RA_ML_AD_Rheumatoid_Arthritis_pos_European_2022_NG_Blood_CD2_NONHODGKIN_NAS_EXALLC_meta_out_finnUKBB_GRCh37_qc",		
	"RA_ML_AD_Rheumatoid_Arthritis_pos_European_2022_NG_Blood_Pheweb_EUR_Malignant_lymphoma")
m=length(dx)
for (j in REP){
	setwd("/home/yxh/data/Summarydata/A_result_PLACO/results")
	ax=data.frame(fread(paste(dx[j],".txt.gz",sep=""), head=T,sep = "\t", fill=TRUE))
	attach(ax)
	zx = data.frame(SNP,INC_ALLELE,DEC_ALLELE,T.placo,P.placo,N,CHR,POS)
	colnames(zx)=c("SNP","INC_ALLELE","DEC_ALLELE","Z","P","N","CHR","POS")
	ax = zx
	index4 = which(ax$CHR == 6 & ax$POS >= 28477797 & ax$POS <= 33448354)
	if (length(index4)>0) {ax=ax[-index4,]}
	a1=paste(ax$INC_ALLELE, ax$DEC_ALLELE, sep="")
	index1=which((a1=="AT") | (a1=="TA") | (a1=="CG") | (a1=="GC") | (nchar(a1)!=2))
	if (length(index1)>0) {ax=ax[-index1,]}
	ax=ax[which(substr(ax$SNP, 1, 2)=="rs"),]
	maxx=merge(ax, bim, by.x="SNP", by.y="V2")
	#rm(ax)
	a1=paste(maxx$INC_ALLELE, maxx$DEC_ALLELE, sep="")
	a2=paste(maxx$DEC_ALLELE, maxx$INC_ALLELE, sep="")
	a3=paste(maxx$V5, maxx$V6, sep="")
	index2=which((a1==a3) | (a2==a3))
	maxx=maxx[index2, ]
	a1=paste(maxx$INC_ALLELE, maxx$DEC_ALLELE, sep="")
	a2=paste(maxx$DEC_ALLELE, maxx$INC_ALLELE, sep="")
	a3=paste(maxx$V5, maxx$V6, sep="")
	index3=which(a2==a3)
	#maxx$BETA[index3]=-maxx$BETA[index3]
	maxx$Z[index3] =-maxx$Z[index3]

	maxx=as.data.frame(cbind(as.character(maxx$SNP), maxx$V5, maxx$V6, maxx$Z, maxx$N))
	colnames(maxx)=c("SNP","A1","A2","Z","N")
	maxx=na.omit(maxx)
	setwd("/home/yxh/data/Summarydata/A_result_LDSC/PLACO")
	fwrite(data.frame(maxx),paste("ldsc_",dx[j],".txt.gz",sep=""),quote=F,sep='\t',col.names=T,row.names=F)
	print(j)
}
