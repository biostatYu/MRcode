args=as.numeric(commandArgs(TRUE))
REP=args[1]
library(data.table)
library(hyprcoloc)
out = paste("GCST",c(90001391:90002121),sep="")
loci.all = data.frame(fread("/home/yxh/data/Summarydata/A_result_HyprColoc/data/HyPrColoc_GenomicRiskLoci_RA_ML.csv"))
RES = NULL
for(i in REP){
	immunecell.tt = data.frame(fread(paste("/home/data/summarydata/C_data_Immune_cell_v1/data/",out[i],"_buildGRCh37_qc2.txt.gz",sep="")))[,c("SNP","A1","A2","BETA","SE")]
	print("step 1")
	tt = merge(loci.all,immunecell.tt,by="SNP",all.x=TRUE)
	loci.name = names(table(tt$loci.tt.uniqID.z.))
	print("step 2")
	for(j in 1:length(loci.name)){
		try({
		tt1 = na.omit(tt[which(tt$loci.tt.uniqID.z.%in%loci.name[j]),])
		ind = which(duplicated(tt1$SNP))
		if(length(ind)>0){tt1 = tt1[-ind,]}
		bete.HyPrColoc = as.matrix(tt1[,c("BETA.x","BETA.y","BETA")])
		se.HyPrColoc = as.matrix(tt1[,c("SE.x","SE.y","SE")])
		traits = c(tt1$Trait.1[1],tt1$Trait.2[1],out[i])
		colnames(bete.HyPrColoc) = traits; colnames(se.HyPrColoc) = traits
		rsid = tt1$SNP
		rownames(bete.HyPrColoc) = rsid; rownames(se.HyPrColoc) = rsid
		if(dim(bete.HyPrColoc)[1]>3){res = hyprcoloc(bete.HyPrColoc, se.HyPrColoc, trait.names = traits, snp.id = rsid)
		setwd("/home/yxh/data/Summarydata/A_result_HyprColoc/results/RA_ML")
		write.table(data.frame(res$results),paste(out[i],"_",loci.name[j],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
		}
		print(j)
		})
	}
	print(i)
	print("Done")
}