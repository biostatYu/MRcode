
#clumped
	bash
	for risk in PrimaryBiliaryCirrhosis2015 Primary_sclerosing_cholangitis2016 Inflammatory_bowel_disease_2017LangeNG Crohn_disease_2017LangeNG Ulcerative_Colitis_2017LangeNG Systemic_Lupus_Erythematosus_2015_NG Multiple_Sclerosis_2018_qc Rheumatoid_Arthritis_European_2014_nature T1D_2020_Diabetes; do
	plink=/home/yxh/software/plink/plink
	ref=/data/yxh/data/ref/EUR
	out=/data/yxh/data/summarydata/autoimmune/index_snp
	cd /data/yxh/data/summarydata/autoimmune
	${plink} --bfile ${ref}/EUR_1000Genome_phase3_all --clump-snp-field SNP --clump-field P --clump ${risk}.txt --clump-p1 5E-8 --clump-r2 0.01 --clump-kb 500
	mv plink.clumped ${out}/plink_clump_${risk}_5E8.txt
	rm plink.log
	rm plink.nosex
	done



#------------------------------------------------------------DATA--------------------------------------------------------#
	library(data.table)
	library(MendelianRandomization)
	library(progress)

	out = c("Crohn_disease",
		"Inflammatory_bowel_disease",
		"Multiple_Sclerosis",
		"Rheumatoid_Arthritis",
		"Systemic_Lupus_Erythematosus",
		"T1D",
		"Ulcerative_Colitis",
		"PrimaryBiliaryCirrhosis",
		"Primary_sclerosing_cholangitis2016")
	for(i in 1:length(out)){
		setwd("/data/yxh/data/summarydata/autoimmune")
		mx0 = data.frame(fread(paste(out[i],".txt",sep = ""),sep='\t',head=T))
		setwd("/data/yxh/data/summarydata/metabolite/indexsnp")
		index_snp = data.frame(read.table("metabolites_index_SNP_impute.txt",head=T))
		x.1 = names(table(index_snp$Metabolites))
		m=length(x.1)
		rex = matrix(NA, m, 6)
		yx0 = merge(index_snp,mx0,by="SNP",all.x=TRUE)
		setwd("/data/yxh/data/summarydata/metabolite/indexsnp/autoimmune")
		write.table(yx0,paste("yx_",out[i],"_metabolites.txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
		print(i)
	}




#------------------------------------------------------------IVW-----------------------------------------------------------#

	library(data.table)
	library(MendelianRandomization)
	library(progress)

	out = c("Crohn_disease",
		"Inflammatory_bowel_disease",
		"Multiple_Sclerosis",
		"Rheumatoid_Arthritis",
		"Systemic_Lupus_Erythematosus",
		"T1D",
		"Ulcerative_Colitis",
		"PrimaryBiliaryCirrhosis",
		"Primary_sclerosing_cholangitis2016")
		for(i in 1:length(out)){
			setwd("/data/yxh/data/summarydata/metabolite/indexsnp/autoimmune")
			yx0 = read.table(paste("yx_",out[i],"_metabolites.txt",sep=""),header=TRUE)
			x.1 = names(table(yx0$Metabolites))
			m=length(x.1)
			pb <- progress_bar$new(total = m)
			rex = matrix(NA, m, 6)
			for (j in 1:m){
				try({
					yx = yx0[which(yx0$Metabolites==x.1[j]),]
					index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
					if (length(index1)<length(yx$INC_ALLELE.x))
					{yx$BETA.y[-index1]=yx$BETA.y[-index1]*(-1)}
					index4 = which(yx$CHR.x == 6 & yx$POS.x >= 28477797 & yx$POS.x <= 33448354)
					if (length(index4)>0) {yx=yx[-index4,]}
					yx=na.omit(data.frame(yx$BETA.x,yx$SE.x,yx$BETA.y,yx$SE.y,yx$P.y))
					colnames(yx)=c("BETA.x","SE.x","BETA.y","SE.y","P.y")
					yx=as.data.frame(yx)
					index2=which(yx$P.y<=(0.05/length(yx$P.y)))
					if (length(index2)>0) {yx=yx[-index2,]}
					index3 = which(abs(yx$BETA.y)>100)
					if (length(index3)>0) {yx=yx[-index3,]}
					index4 = which(yx$BETA.y == 0)
					if (length(index4)>0) {yx=yx[-index4,]}
					fit1 = mr_ivw(mr_input(bx = yx$BETA.x, bxse = yx$SE.x, by = yx$BETA.y, byse = yx$SE.y),model="fixed")
					fit2 = mr_ivw(mr_input(bx = yx$BETA.x, bxse = yx$SE.x, by = yx$BETA.y, byse = yx$SE.y),model="random")
					if(fit1$Heter.Stat[2]<0.05){fit = fit2}else{fit = fit1}
					res = c(dim(yx)[1],fit$Estimate,fit$StdError,fit$Pvalue,fit$Heter.Stat)
					rex[j,] = res
					#setwd("/data/yxh/data/summarydata/metabolite/result/autoimmune/single")
					#write.table(res,paste(out[i],"_Metabolites_",x.1[j],".txt",sep = ""),quote=F,sep="\t",col.names=F,row.names=F)
					rm(fit1);rm(fit2);rm(fit);rm(res)
					pb$tick()
					Sys.sleep(1/m)
				})
			}
					rexx = cbind(as.character(x.1),rex)
					colnames(rexx) = c("Metabolites","N","beta","se","pvalue","Heter.Stat","Heter.p")
					setwd("/data/yxh/data/summarydata/metabolite/result/autoimmune")
					write.table(rexx,paste("Metabolites_",out[i],"_IVW.txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
					print(i)
		}



#------------------------------------------------------------sensitive analysis-----------------------------------------------------------#

	library(data.table)
	setwd("/home/yxh/program/MR")
	source("MRfunction.R")
	library(progress)
	library(MRPRESSO)

	out = c("Crohn_disease",
		"Inflammatory_bowel_disease",
		"Multiple_Sclerosis",
		"Rheumatoid_Arthritis",
		"Systemic_Lupus_Erythematosus",
		"T1D",
		"Ulcerative_Colitis",
		"PrimaryBiliaryCirrhosis",
		"Primary_sclerosing_cholangitis2016")
		for(i in 1:length(out)){
			setwd("/data/yxh/data/summarydata/metabolite/indexsnp/autoimmune")
			yx0 = read.table(paste("yx_",out[i],"_metabolites.txt",sep=""),header=TRUE)
			x.1 = names(table(yx0$Metabolites))
			m=length(x.1)
			pb <- progress_bar$new(total = m)
			rex = matrix(NA, m, 36)
			for (j in 1:m){
				try({
					yx = yx0[which(yx0$Metabolites==x.1[j]),]
					index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
					if (length(index1)<length(yx$INC_ALLELE.x))
					{yx$BETA.y[-index1]=yx$BETA.y[-index1]*(-1)}
					index4 = which(yx$CHR.x == 6 & yx$POS.x >= 28477797 & yx$POS.x <= 33448354)
					if (length(index4)>0) {yx=yx[-index4,]}
					yx=na.omit(data.frame(yx$BETA.x,yx$SE.x,yx$BETA.y,yx$SE.y,yx$P.y))
					colnames(yx)=c("BETA.x","SE.x","BETA.y","SE.y","P.y")
					yx=as.data.frame(yx)
					index2=which(yx$P.y<=(0.05/length(yx$P.y)))
					if (length(index2)>0) {yx=yx[-index2,]}
					index3 = which(abs(yx$BETA.y)>100)
					if (length(index3)>0) {yx=yx[-index3,]}
					index4 = which(yx$BETA.y == 0)
					if (length(index4)>0) {yx=yx[-index4,]}
					n=dim(yx)[1]
					R2 = yx$BETA.x*yx$BETA.x/(yx$BETA.x*yx$BETA.x+yx$SE.x*yx$SE.x*7824)
					F = R2*(7824-2)/(1-R2)
					Res = MR_sum(yx,result="OR")$result
					rex[j,] = c(j,sum(R2),min(F),n,Res[1,1:6],Res[2,1:6],Res[5,1:6],Res[3,1:6],Res[4,1:6],Res[1,7:8])
					rm(Res)
					pb$tick()
					Sys.sleep(1/m)
				})
			}
					rexx = cbind(as.character(x.1),rex)
					colnames(rexx) = c("Metabolites","ID","PVE","F_min","n_snp",
							"IVW_fix","SE","OR","Down","Up","Pval",
							"IVW_random","SE","OR","Down","Up","Pval",
							"weighted_Median","SE","OR","Down","Up","Pval",					
							"Egger_slope","SE","OR","Down","Up","Pval",
							"Egger_intercept","SE","OR","Down","Up","Pval","Heter.est","Heter.Pval")
					rm(rex)
					setwd("/data/yxh/data/summarydata/metabolite/result/autoimmune")
					write.table(rexx,paste("Metabolites_",out[i],"_sensitive.txt",sep=""),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
					print(i)
		}

#------------------------------------------------------------MR-PRESSO-----------------------------------------------------------#

	library(data.table)
	setwd("/home/yxh/program/MR")
	source("MRfunction.R")
	library(progress)
	library(MRPRESSO)

	out = c("Crohn_disease",
		"Inflammatory_bowel_disease",
		"Multiple_Sclerosis",
		"Rheumatoid_Arthritis",
		"Systemic_Lupus_Erythematosus",
		"T1D",
		"Ulcerative_Colitis",
		"PrimaryBiliaryCirrhosis",
		"Primary_sclerosing_cholangitis2016")

		for(i in 4:length(out)){
			setwd("/data/yxh/data/summarydata/metabolite/indexsnp/autoimmune")
			yx0 = read.table(paste("yx_",out[i],"_metabolites.txt",sep=""),header=TRUE)
			x.1 = names(table(yx0$Metabolites))
			m=length(x.1)
			pb <- progress_bar$new(total = m)
			rex = matrix(NA, m, 42)
			for (j in 1:m){
				try({
					yx = yx0[which(yx0$Metabolites==x.1[j]),]
					index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
					if (length(index1)<length(yx$INC_ALLELE.x))
					{yx$BETA.y[-index1]=yx$BETA.y[-index1]*(-1)}
					index4 = which(yx$CHR.x == 6 & yx$POS.x >= 28477797 & yx$POS.x <= 33448354)
					if (length(index4)>0) {yx=yx[-index4,]}
					yx=na.omit(data.frame(yx$BETA.x,yx$SE.x,yx$BETA.y,yx$SE.y,yx$P.y))
					colnames(yx)=c("BETA.x","SE.x","BETA.y","SE.y","P.y")
					yx=as.data.frame(yx)
					index2=which(yx$P.y<=(0.05/length(yx$P.y)))
					if (length(index2)>0) {yx=yx[-index2,]}
					index3 = which(abs(yx$BETA.y)>100)
					if (length(index3)>0) {yx=yx[-index3,]}
					index4 = which(yx$BETA.y == 0)
					if (length(index4)>0) {yx=yx[-index4,]}
					
					fit=mr_presso(BetaOutcome="BETA.y",BetaExposure="BETA.x",SdOutcome="SE.y",SdExposure="SE.x",
							OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=yx,NbDistribution=1000,SignifThreshold=0.05)
					res.mrpresso = data.frame(fit$`Main MR results`)
					rs1 = fit$`MR-PRESSO results`$`Global Test`$RSSobs
					rs2 = fit$`MR-PRESSO results`$`Global Test`$Pvalue
					global = c(rs1,rs2)
					ouliers = fit$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
					res.mrpresso$Low = res.mrpresso$Causal.Estimate-1.96*res.mrpresso$Sd;
					res.mrpresso$Up = res.mrpresso$Causal.Estimate+1.96*res.mrpresso$Sd;

					raw = c(res.mrpresso$Causal.Estimate[1],res.mrpresso$Sd[1],res.mrpresso$Low[1],res.mrpresso$Up[1],res.mrpresso$P.value[1])
					adj = c(res.mrpresso$Causal.Estimate[2],res.mrpresso$Sd[2],res.mrpresso$Low[2],res.mrpresso$Up[2],res.mrpresso$P.value[2])
					if(length(ouliers)==0){
						n = dim(yx); res = c(j,n,global,raw,adj, rep(NA,27)); rex[j,] = res	
					}else if(ouliers=="No significant outliers"){
						n = dim(yx); res = c(j,n,global,raw,adj, rep(NA,27)); rex[j,] = res	
					}else{
						yx = yx[-ouliers,]
						n = dim(yx)
						Res = MR_sum(yx,result="BETA")$result
						res = c(j,n,global,raw,adj,Res[1,1:5],Res[2,1:5],Res[6,1:5],Res[3,1:5],Res[4,1:5],Res[1,6:7])
						rex[j,] = res
						}
					setwd("/data/yxh/data/summarydata/metabolite/result/autoimmune/single")
					write.table(res,paste(out[i],"_gut_mrpresso_",x.1[j],".txt",sep = ""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
					pb$tick()
					Sys.sleep(1/m)
					#print(j)
					})
				}
					rexx = cbind(as.character(x.1),rex)
					colnames(rexx) = c("Gut","ID","N","global.test","global.pval",
							"mrpresso","SE","Down","Up","Pval",
							"mrpresso_out_adjust","SE","Down","Up","Pval",
							"IVW_fix_out_adjust","SE","Down","Up","Pval",
							"IVW_random_out_adjust","SE","Down","Up","Pval",
							"weighted_Median_out_adjust","SE","Down","Up","Pval",					
							"Egger_slope_out_adjust","SE","Down","Up","Pval",
							"Egger_intercept_out_adjust","SE","Down","Up","Pval","Heter.est","Heter.Pval")
					rm(rex)
					setwd("/data/yxh/data/summarydata/metabolite/result/autoimmune")
					write.table(rexx,paste("Metabolites_",out[i],"_mrpresso.txt",sep=""),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
					print(i)
			}



#--------------------------OP to RA------------------------------#	
	meta = c("90001480","90001481","90001482","90001483","90001485",
		"90001487","90001493","90001498","90001499","90001500",
		"90001508","90001510","90001511","90001518","90001535",
		"90001538","90001542","90001544","90001546","90001563",
		"90001658","90001659","90001670","90001671","90001672",
		"90001677","90001841","90001845","90001937","90001980",
		"90001981","90001985","90001989","90002044","90002073","90002103")
	library(data.table)
	setwd("D:\\work\\Program\\MR")
	source("MRfunction.R")
	setwd("H:\\data\\summarydata")
	DD0=data.frame(fread(paste("Bone_mineral_density_all_GEFOS_2018AJHG.txt",sep=""),head=T,sep="\t"))
	setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\index_snp")
	snp = as.character(read.table(paste("BMD_TotalBody_all.txt",sep=""),head=T)$SNP)
	RES=NULL
	for(i in 1:36){
		setwd("H:\\data\\summarydata\\immue\\data_qc")
		ax0 = data.frame(fread(paste("GCST",meta[i],"_buildGRCh37_qc.txt.gz",sep = ""),head=T,sep='\t',fill= T))
		ax = ax0[which(ax0$SNP%in%snp),] 
		zx = with(ax,data.frame(SNP,CHR,POS,toupper(A1),toupper(A2),BETA,SE,BETA/SE,P))
		colnames(zx)=c("SNP","CHR","POS","INC_ALLELE","DEC_ALLELE","BETA","SE","Z","P")
		
		data = MR_process(zx,DD0,snp)	
		yx = data$yx;
		m = data$m		
		n = dim(yx)[1]

		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		#plot(cooksd, pch=20, cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
	#	"rs2476601" "rs9268839" 
		if(length(index)>0){yx = yx[-index,]}
		Res = MR_sum(yx,result = "BETA")
		res = Res$result
		RES = rbind(RES,cbind(meta[i],res))
		print(i)
	}
		setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\reverseMR")
		write.csv(RES,paste("result_reverse_TB-BMD.csv",sep=""))


	meta = c("90002073","90001472","90001464","90001459","90001473")
	library(data.table)
	setwd("D:\\work\\Program\\MR")
	source("MRfunction.R")
	setwd("H:\\data\\summarydata")
	DD0=data.frame(fread(paste("BMD_FemoralNeck.txt",sep=""),head=T,sep="\t"))
	setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\index_snp")
	snp = as.character(read.table(paste("BMD_FemoralNeck.txt",sep=""),head=T)$SNP)
	RES=NULL
	for(i in 1:5){
		setwd("H:\\data\\summarydata\\immue\\data_qc")
		ax0 = data.frame(fread(paste("GCST",meta[i],"_buildGRCh37_qc.txt.gz",sep = ""),head=T,sep='\t',fill= T))
		ax = ax0[which(ax0$SNP%in%snp),] 
		zx = with(ax,data.frame(SNP,CHR,POS,toupper(A1),toupper(A2),BETA,SE,BETA/SE,P))
		colnames(zx)=c("SNP","CHR","POS","INC_ALLELE","DEC_ALLELE","BETA","SE","Z","P")
		
		data = MR_process(zx,DD0,snp)	
		yx = data$yx;
		m = data$m		
		n = dim(yx)[1]

		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		#plot(cooksd, pch=20, cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
	#	"rs2476601" "rs9268839" 
		if(length(index)>0){yx = yx[-index,]}
		Res = MR_sum(yx,result = "BETA")
		res = Res$result
		RES = rbind(RES,cbind(meta[i],res))
		print(i)
	}
		setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\reverseMR")
		write.csv(RES,paste("result_reverse_FN-BMD.csv",sep=""))



	meta = c("90001672","90001489","90001423","90001495","90002074","90002071")
	library(data.table)
	setwd("D:\\work\\Program\\MR")
	source("MRfunction.R")
	setwd("H:\\data\\summarydata")
	DD0=data.frame(fread(paste("BMD_LumbarSpine.txt",sep=""),head=T,sep="\t"))
	setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\index_snp")
	snp = as.character(read.table(paste("BMD_LumbarSpine.txt",sep=""),head=T)$SNP)
	RES=NULL
	for(i in 1:6){
		setwd("H:\\data\\summarydata\\immue\\data_qc")
		ax0 = data.frame(fread(paste("GCST",meta[i],"_buildGRCh37_qc.txt.gz",sep = ""),head=T,sep='\t',fill= T))
		ax = ax0[which(ax0$SNP%in%snp),] 
		zx = with(ax,data.frame(SNP,CHR,POS,toupper(A1),toupper(A2),BETA,SE,BETA/SE,P))
		colnames(zx)=c("SNP","CHR","POS","INC_ALLELE","DEC_ALLELE","BETA","SE","Z","P")
		
		data = MR_process(zx,DD0,snp)	
		yx = data$yx;
		m = data$m		
		n = dim(yx)[1]

		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		#plot(cooksd, pch=20, cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
	#	"rs2476601" "rs9268839" 
		if(length(index)>0){yx = yx[-index,]}
		Res = MR_sum(yx,result = "BETA")
		res = Res$result
		RES = rbind(RES,cbind(meta[i],res))
		print(i)
	}
		setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\reverseMR")
		write.csv(RES,paste("result_reverse_LS-BMD.csv",sep=""))



	meta = c("90001980","90001585","90001462","90001458","90001771")
	library(data.table)
	setwd("D:\\work\\Program\\MR")
	source("MRfunction.R")
	setwd("H:\\data\\summarydata")
	DD0=data.frame(fread(paste("BMD_Forearm.txt",sep=""),head=T,sep="\t"))
	setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\index_snp")
	snp = as.character(read.table(paste("BMD_Forearm.txt",sep=""),head=T)$SNP)
	RES=NULL
	for(i in 1:5){
		setwd("H:\\data\\summarydata\\immue\\data_qc")
		ax0 = data.frame(fread(paste("GCST",meta[i],"_buildGRCh37_qc.txt.gz",sep = ""),head=T,sep='\t',fill= T))
		ax = ax0[which(ax0$SNP%in%snp),] 
		zx = with(ax,data.frame(SNP,CHR,POS,toupper(A1),toupper(A2),BETA,SE,BETA/SE,P))
		colnames(zx)=c("SNP","CHR","POS","INC_ALLELE","DEC_ALLELE","BETA","SE","Z","P")
		
		data = MR_process(zx,DD0,snp)	
		yx = data$yx;
		m = data$m		
		n = dim(yx)[1]

		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		#plot(cooksd, pch=20, cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
	#	"rs2476601" "rs9268839" 
		if(length(index)>0){yx = yx[-index,]}
		Res = MR_sum(yx,result = "BETA")
		res = Res$result
		RES = rbind(RES,cbind(meta[i],res))
		print(i)
	}
		setwd("E:\\Cloud_disk\\program\\MR\\MR_immu\\reverseMR")
		write.csv(RES,paste("result_reverse_FA-BMD.csv",sep=""))



	library(data.table)
	setwd("E:\\Cloud_disk\\program")
	source("MRfunction.R")
	library(progress)
	library(MRPRESSO)
	out = c("Crohn_disease",
		"Inflammatory_bowel_disease",
		"Multiple_Sclerosis",
		"Rheumatoid_Arthritis",
		"Systemic_Lupus_Erythematosus",
		"T1D",
		"Ulcerative_Colitis",
		"PrimaryBiliaryCirrhosis",
		"Primary_sclerosing_cholangitis2016")
		for(i in 1:9){
			setwd("E:\\Cloud_disk\\program\\MR\\MR_meta_immu\\index_snp")
			yx0 = read.table(paste("yx_",out[i],"_metabolites.txt",sep=""),header = TRUE)						
			x.1 = read.table(paste("sensitive_immune_id_",out[i],".txt",sep=""),header =FALSE)[,1]
			m=length(x.1)
			pb <- progress_bar$new(total = m)
			rex = matrix(NA, m, 41)
				for (j in 1:m){
					try({
						yx = yx0[which(yx0$Metabolites==x.1[j]),]
						index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
						if (length(index1)<length(yx$INC_ALLELE.x))
						{yx$BETA.y[-index1]=yx$BETA.y[-index1]*(-1)}
						index4 = which(yx$CHR.x == 6 & yx$POS.x >= 28477797 & yx$POS.x <= 33448354)
						if (length(index4)>0) {yx=yx[-index4,]}
						yx=na.omit(data.frame(yx$BETA.x,yx$SE.x,yx$BETA.y,yx$SE.y,yx$P.y))
						colnames(yx)=c("BETA.x","SE.x","BETA.y","SE.y","P.y")
						yx=as.data.frame(yx)
						index2=which(yx$P.y<=(0.05/length(yx$P.y)))
						#if (length(index2)>0) {yx=yx[-index2,]}
						index3 = which(abs(yx$BETA.y)>100)
						if (length(index3)>0) {yx=yx[-index3,]}
						index4 = which(yx$BETA.y == 0)
						if (length(index4)>0) {yx=yx[-index4,]}						
						fit=mr_presso(BetaOutcome="BETA.y",BetaExposure="BETA.x",SdOutcome="SE.y",SdExposure="SE.x",
								OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=yx,NbDistribution=1000,SignifThreshold=0.05)
						res.mrpresso = data.frame(fit$`Main MR results`)		
						ouliers = fit$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
						res.mrpresso$Low = res.mrpresso$Causal.Estimate-1.96*res.mrpresso$Sd;
						res.mrpresso$Up = res.mrpresso$Causal.Estimate+1.96*res.mrpresso$Sd;
						global = c(fit$`MR-PRESSO results`$`Global Test`$RSSobs,fit$`MR-PRESSO results`$`Global Test`$Pvalue)
						raw = c(res.mrpresso$Causal.Estimate[1],res.mrpresso$Sd[1],res.mrpresso$Low[1],res.mrpresso$Up[1],res.mrpresso$P.value[1])
						adj = c(res.mrpresso$Causal.Estimate[2],res.mrpresso$Sd[2],res.mrpresso$Low[2],res.mrpresso$Up[2],res.mrpresso$P.value[2])
						if(length(ouliers)==0){
							n = dim(yx)[1]; res = c(j,n,global,raw,adj, rep(NA,27)); rex[j,] = res	
						}else if(ouliers=="No significant outliers"){
							n = dim(yx)[1]; res = c(j,n,global,raw,adj, rep(NA,27)); rex[j,] = res	
						}else{
							yx = yx[-ouliers,]
							n = dim(yx)[1]
							Res = MR_sum(yx,result="BETA")$result
							res = c(j,n,global,raw,adj,Res[1,1:5],Res[2,1:5],Res[6,1:5],Res[3,1:5],Res[4,1:5],Res[1,6:7])
							rex[j,] = res
							}
						rm(fit)
						setwd("E:\\Cloud_disk\\program\\MR\\MR_meta_immu\\result\\single")
						write.table(res,paste(out[i],"_mrpresso_",x.1[j],".txt",sep = ""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
						pb$tick()
						Sys.sleep(1/m)
						})
					}
						rexx = cbind(as.character(x.1),rex)
						colnames(rexx) = c("Immune","ID","N","global.test","global.test.p",
								"mrpresso","SE","Down","Up","Pval",
								"mrpresso_out_adjust","SE","Down","Up","Pval",
								"IVW_fix_out_adjust","SE","Down","Up","Pval",
								"IVW_random_out_adjust","SE","Down","Up","Pval",
								"weighted_Median_out_adjust","SE","Down","Up","Pval",					
								"Egger_slope_out_adjust","SE","Down","Up","Pval",
								"Egger_intercept_out_adjust","SE","Down","Up","Pval","Heter.est","Heter.Pval")
						rm(rex)
						setwd("E:\\Cloud_disk\\program\\MR\\MR_meta_immu\\result")
						write.table(rexx,paste("mrpresso_",out[i],".txt",sep=""),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
						print(i)
				}


