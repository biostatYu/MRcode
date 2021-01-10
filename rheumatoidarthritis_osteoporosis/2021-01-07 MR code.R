

#-------------------------------------------------------MR Analysis--------------------------------------------------------------#

	#-------------------------BMD-----------------------#
		source("MRfunction.R")
		library(ggplot2)

		yx = read.table("index_snp_BMD_Rheumatoid_Arthritis.txt",header = T)
		dim(yx)
		index.2 = which(as.numeric(yx$P.y)<(0.05/dim(yx)[1]))
		yx$SNP[index.2]
		yx = yx[-index.2,]
	
		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		plot(cooksd, pch=20, cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
		yx = yx[-index,]

		#BW related snp removed
		bw_snp = read.table("RA_related.txt", header = T, stringsAsFactors = F)[,1]
		yx = yx[-which(yx$SNP%in%bw_snp),]
		n = dim(yx)[1]
		Res = MR_sum(yx,result = "BETA")
		res = Res$result

		library(MRPRESSO)
		fit=mr_presso(BetaOutcome="BETA.y",BetaExposure="BETA.x",SdOutcome="SE.y",SdExposure="SE.x",
			OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=yx,NbDistribution=10000,SignifThreshold=0.05)

	#-------------------------Fracture-----------------------#
		yx = read.table("index_snp_Fractures_Rheumatoid_Arthritis.txt",header = T)
		dim(yx)
		index.2 = which(as.numeric(yx$P.y)<(0.05/dim(yx)[1]))
		yx$SNP[index.2]
		
		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		plot(cooksd, pch=20, cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
		yx = yx[-index,]

		#BW related snp removed
		bw_snp = read.table("RA_related.txt", header = T, stringsAsFactors = F)[,1]
		yx = yx[-which(yx$SNP%in%bw_snp),]

		n = dim(yx)[1]
		Res = MR_sum(yx,result = "OR")
		res = Res$result

		library(MRPRESSO)
		fit=mr_presso(BetaOutcome="BETA.y",BetaExposure="BETA.x",SdOutcome="SE.y",SdExposure="SE.x",
			OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=yx,NbDistribution=10000,SignifThreshold=0.05)


	#-------------------------osteoporosis Asian-----------------------#
		yx = read.table("index_snp_Rheumatoid_Arthritis_osteoporosis_Asian.txt",header = T)
		dim(yx)
		index.2 = which(as.numeric(yx$P.y)<(0.05/dim(yx)[1]))
		yx$SNP[index.2]
	
		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		plot(cooksd, pch=20, cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
		yx = yx[-index,]

		#BW related snp removed
		bw_snp = read.table("RA_related.txt", header = T, stringsAsFactors = F)[,1]
		yx = yx[-which(yx$SNP%in%bw_snp),]
		n = dim(yx)[1]
		Res = MR_sum(yx,result = "OR")
		res = Res$result

		library(MRPRESSO)
		fit=mr_presso(BetaOutcome="BETA.y",BetaExposure="BETA.x",SdOutcome="SE.y",SdExposure="SE.x",
			OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=yx,NbDistribution=10000,SignifThreshold=0.05)

	#-------------------------Multivariable MR-----------------------#
		library(data.table)
		library(MendelianRandomization)
		m=c("Asthma","Crohn_disease","Inflammatory_bowel_disease","MultipleSclerosi","systemic_lupus_erythematosus","T1D","Ulcerative_colitis")
		library(data.table)
		mm = length(m)
		rex = matrix(NA, mm, 5)
		rex.2 = matrix(NA, mm, 5)
		yx0 = data.frame(fread(paste("index_snp_BMD_Rheumatoid_Arthritis.txt",sep = ""),sep='\t',head=T))
		snp=yx0$SNP
		for (i in 1:mm){
			x0 = data.frame(fread(paste(m[i],".txt",sep=""),sep='\t',head=T))
			ax1=x0[which(x0$SNP%in%snp),]
			yx1=merge(yx0,ax1,by="SNP",all.x=T)
			index1=which(toupper(yx1$INC_ALLELE.x)==toupper(yx1$INC_ALLELE))
			if (length(index1)<length(yx1$INC_ALLELE.x)){yx1$BETA[-index1]=yx1$BETA[-index1]*(-1)}
			yx=na.omit(yx1)
			fit = mr_mvivw(mr_mvinput(bx = cbind(yx$BETA.x, yx$BETA), bxse = cbind(yx$SE.x, yx$SE),
			   by = yx$BETA.y, byse = yx$SE.y))
			res = cbind(fit$Estimate,fit$StdError,fit$CILower,fit$CIUpper,fit$Pvalue)
			rex[i,] = res[1,]
			rex.2[i,] = res[2,]
			print(i)
		}
		colnames(rex) = c("BETA","SE","Low","Up","P")
		colnames(rex.2) = c("BETA","SE","Low","Up","P")
		write.table(cbind(m,rex),paste("multivariable_RA_main.txt",sep = ""),quote=F,sep="\t",col.names=F,row.names=F)
		write.table(cbind(m,rex.2),paste("multivariable_RA_covarites.txt",sep = ""),quote=F,sep="\t",col.names=F,row.names=F)


#--------------------------reverse------------------------------#	
		library(data.table)
		DD0=data.frame(fread(paste("BMD.txt",sep=""),head=T,sep="\t"))
		ax0 = data.frame(fread(paste("Rheumatoid_Arthritis.txt",sep = ""),head=T,sep='\t'))

		snp = as.character(read.table(paste("snp_BMD_Biobank2.txt",sep=""),head=T)[,1])
		data = MR_process(ax0,DD0,snp)	
		yx = data$yx;
		m = data$m		
		n = dim(yx)[1]

		fit=lm(yx$BETA.y ~ yx$BETA.x -1, weights = 1/yx$SE.y**2)
		cooksd <- cooks.distance(fit)
		sample_size <- length(yx$BETA.y)
		index = which(cooksd > (4/sample_size))
		yx$SNP[index]
		yx = yx[-index,]
		Res = MR_sum(yx)
		res = Res$result


#-------------------------------------UKB-------------------------------------#
	#	phenotype
		library(data.table)
		race = fam1[,which(colnames(fam1)=="X21000.0.0")] # 1001;1002;1003;1
		age = fam1[,which(colnames(fam1)=="X21003.0.0")]
		sex = fam1[,which(colnames(fam1)=="X31.0.0")] # gender; male = 1; female = 0

		height  = fam1[,which(colnames(fam1)=="X50.0.0")]  # standing height
		weight  = fam1[,which(colnames(fam1)=="X21002.0.0")]  # weight
		eid = fam1[,which(colnames(fam1)=="eid")] 
		BMD.T1  = fam1[,which(colnames(fam1)=="X77.0.0")] # Heel bone ultrasound T-score, manual entry
		BMD.T2  = fam1[,which(colnames(fam1)=="X78.0.0")] # Heel bone mineral density (BMD) T-score, automated
		BMD.T3  = fam1[,which(colnames(fam1)=="X4106.0.0")] # 	Heel bone mineral density (BMD) T-score, automated (left)
		BMD1  = fam1[,which(colnames(fam1)=="X3084.0.0")] # Heel bone mineral density (BMD), manual entry
		BMD2  = fam1[,which(colnames(fam1)=="X3148.0.0")] # Heel bone mineral density (BMD)
		BMD3  = fam1[,which(colnames(fam1)=="X4105.0.0")] # Heel bone mineral density (BMD) (left)
		BMD = cbind(BMD1,BMD2,BMD3)
		BMD.T = cbind(BMD.T1,BMD.T2,BMD.T3)

		Summary = function(x){
			if(is.na(x[3])){
				if(is.na(x[2])){
					if(is.na(x[1])){tmp = NA;
					}else(tmp = x[1])
				}else(tmp = x[2])
			}else{tmp = x[3]}
		}
		bmd = apply(BMD,1,Summary)
		bmd.T = apply(BMD.T,1,Summary)

		Pheno = cbind(eid,bmd.T,bmd,bw,age,sex,race,height,weight)

	#genotype
		gx = fread(paste("1_snp.raw",sep = ""), header = T)[,-(2:6)]
		for(i in 2:22){
		try({
		Gx = fread(paste(i,"_snp.raw",sep = ""), header = T)[,-(1:6)]
		gx = cbind(gx, Gx)
		print(i)
		})
		}
		fwrite(data.frame(gx),"Genotype_ra.txt", row.names = F, col.names = T,quote = F, sep = "\t")
		
		# merge 33IDÓë48ID
		gx = data.frame(fread("Genotype_ra.txt", header = T))
		temp = data.frame(fread("chr_1.fam", header = F))
		data = merge(temp, gx, by.x="V1", by.y="FID", all.x = T)
		#data = merge(temp, gx, by.x="X1", by.y="FID")
		GE = data[,-(1:6)]
		p=dim(GE)[2]
		infor=matrix(NA,p,2)
		for (ii in 1:p){
			GE[,ii]=as.numeric(paste(GE[,ii]))
			index=which(is.na(GE[,ii]))
			infor[ii,1] = length(index)
			if (length(index)>0){GE[index,ii] = mean(GE[-index,ii])}
			infor[ii,2] = sum(GE[,ii]==0)
			print(ii)
			}
		index=which(infor[,1]>(dim(GE)[1]*0.1))
		if(length(index)>0){GE=GE[,-index]}

	#---------------------------------------analysis---------------------------------------#
		library(data.table)
		cx = data.frame(fread("Phenotype.txt", header = T))
		gx = data.frame(fread("Genotype_qc.txt", header = T))
		temp = colnames(gx)	
		temp1 = strsplit(temp,"_")
		rsid = c(); A1 = c()
		for(i in 1:length(temp1)){rsid = c(rsid,temp1[[i]][1]); A1 = c(A1,temp1[[i]][2])}
		A1 = substr(A1,1,1)

		index.snp0 = read.table("yx_RA.txt",header = T)
		index.snp0 = na.omit(index.snp0)
		index = which(as.character(index.snp0$SNP)%in%rsid)
		index_snp1 = index.snp0[index,]

		GRS_ind = NULL
		for(j in 1:dim(index_snp1)[1]){
			RS = index_snp1$SNP[j]
			index_RS = which(rsid == RS)
			geno_temp = gx[,index_RS]
			if(index_snp1$INC_ALLELE[j] == A1[index_RS]){
			grs = geno_temp*index_snp1$BETA[j]
			}else{
			grs = (2-geno_temp)*index_snp1$BETA[j]
			}
			GRS_ind = cbind(GRS_ind,grs)
			print(j)
		}
		
		GRS = apply(GRS_ind,1,sum)
		cx_GRS = cbind(cx,GRS)

#-------------------------------------------------------------------------GRS---------------------------------------------------------------------------------------#
	
			#---------------------------------------Male + Female---------------------------------------#
				library(data.table)
				library(nnet)
				cx_GRS = data.frame(fread("GRS_RA.txt", header = T))
				index = which(cx_GRS$race%in%c("1001","1002","1003","1"))
				cx_GRS = cx_GRS[index,]
				index = union(which(cx_GRS$bmd>mean(na.omit(cx_GRS$bmd))+4.5*sd(na.omit(cx_GRS$bmd))),which(cx_GRS$bmd<mean(na.omit(cx_GRS$bmd))-4.5*sd(na.omit(cx_GRS$bmd))))
				cx_GRS = cx_GRS[-index,]

				vari = c("Fracture","bmd","GRS","age","sex","smoke","drink","height","weight")
				cx1 = cx_GRS[,vari]
				cx1[,-c(1,5:7)] = data.frame(apply(cx1[,-c(1,5:7)],2,scale))
				fit.1 = lm(bmd~GRS+age+height+weight+sex, data = cx1)
				res.1.1 = summary(fit.1)$coefficients

				fit.1 = glm(Fracture~GRS+age+height+weight+sex, data = cx1, family=binomial(link='logit'))
				res.1.2 = summary(fit.1)$coefficients
				#---------------------------------------Male---------------------------------------#
				index = which(cx1$sex == 0)
				cx2 = cx1[index,]
				fit.2 = lm(bmd~GRS+age+height+weight,  data = cx2)
				res.2.1 = summary(fit.2)$coefficients
	
				fit.1 = glm(Fracture~GRS+age+height+weight, data = cx2, family=binomial(link='logit'))
				res.2.2 = summary(fit.1)$coefficients
			
			#---------------------------------------Female---------------------------------------#
				
				index = which(cx1$sex == 1)
				cx3 = cx1[index,]	 
				fit.2 = lm(bmd~GRS+age+height+weight, data = cx3)
				res.3.1 = summary(fit.2)$coefficients

				fit.1 = glm(Fracture~GRS+age+height+weight, data = cx3, family=binomial(link='logit'))
				res.3.2 = summary(fit.1)$coefficients

				res = rbind(res.1.1, res.1.2, res.2.1, res.2.2, res.3.1, res.3.2)
				colnames(res) = c("BETA","SE","t","P")



#-------------------------------------------------------------------------mediation---------------------------------------------------------------------------------------#
			
			library(mediation)
			library(data.table)
			cx_GRS = data.frame(fread("GRS_RA.txt", header = T))
			index = which(cx_GRS$race%in%c("1001","1002","1003","1"))
			cx_GRS = cx_GRS[index,]
			index = union(which(cx_GRS$bmd>mean(na.omit(cx_GRS$bmd))+4.5*sd(na.omit(cx_GRS$bmd))),which(cx_GRS$bmd<mean(na.omit(cx_GRS$bmd))-4.5*sd(na.omit(cx_GRS$bmd))))
			cx_GRS = cx_GRS[-index,]

			vari = c("Fracture","bmd","GRS","age","sex","smoke","drink","height","weight")
			cx1 = cx_GRS[,vari]
			cx1[,-c(1,5:7)] = data.frame(apply(cx1[,-c(1,5:7)],2,scale))
			cx1 = na.omit(cx1)
			#---------------------------------------bmd---------------------------------------#

			fit.a = lm(bmd ~ GRS + age + height + weight, data = cx1)
			fit.c = glm(Fracture ~ GRS + bmd + age + height + weight, data = cx1, family=binomial)
			contcont <- mediate(fit.a, fit.c, treat = "GRS", mediator = "bmd",sims = 1000)
			tmp = summary(contcont)
			tmp
			res = rbind(c(tmp$d0,tmp$d0.ci,tmp$d0.p),
					c(tmp$z0,tmp$z0.ci,tmp$z0.p),
					c(tmp$n0,tmp$n0.ci,tmp$n0.p),
					c(tmp$tau.coef,tmp$tau.ci,tmp$tau.p))

			#--------male--------#
			cx.GRS.male = cx2.2[which(cx2.2$sex == 1),]
			fit.a = lm(bmd ~ GRS + age + smoke + drink + height + weight, data = cx.GRS.male)
			fit.c = glm(Fracture ~ GRS + bmd + age + smoke + drink + height + weight, data = cx.GRS.male, family=binomial)
			contcont <- mediate(fit.a, fit.c, treat = "GRS", mediator = "bmd",sims = 1000)
			tmp = summary(contcont)
			res = rbind(c(tmp$d0,tmp$d0.ci,tmp$d0.p),
					c(tmp$z0,tmp$z0.ci,tmp$z0.p),
					c(tmp$n0,tmp$n0.ci,tmp$n0.p),
					c(tmp$tau.coef,tmp$tau.ci,tmp$tau.p))

			#--------female--------#
			cx.GRS.female = cx2.2[which(cx2.2$sex == 0),]
			fit.a = lm(bmd ~ GRS + age + smoke + drink + height + weight, data = cx.GRS.female)
			fit.c = glm(Fracture ~ GRS + bmd + age + smoke + drink + height + weight, data = cx.GRS.female, family=binomial)
			contcont <- mediate(fit.a, fit.c, treat = "GRS", mediator = "bmd")
			tmp = summary(contcont)
			res = rbind(c(tmp$d0,tmp$d0.ci,tmp$d0.p),
					c(tmp$z0,tmp$z0.ci,tmp$z0.p),
					c(tmp$n0,tmp$n0.ci,tmp$n0.p),
					c(tmp$tau.coef,tmp$tau.ci,tmp$tau.p))



#-------------------------------------data process-------------------------------------#
	library(MendelianRandomization)
	#DD0 ½á¾Ö
	#ax0 ±©Â¶
	#snp index_exposure(vector)

	MR_process = function(DD0,ax0,snp){
		if(length(snp)>1){
			ax = ax0[which(ax0$SNP%in%snp),]
			DD = DD0[which(DD0$SNP%in%snp),]
			}else if(length(snp)==1){
				snp = ax0$SNP
				ax = ax0
				DD = DD0[which(DD0$SNP%in%snp),]
			}
		yx=merge(ax,DD,by="SNP")
		index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
		if (length(index1)<length(yx$INC_ALLELE.x)){
			yx$BETA.y[-index1]=yx$BETA.y[-index1]*(-1)
			}
		#index = which(yx$P.y<(0.05/length(yx$P.y)))
		m = dim(yx)[1]
		#if (length(index)>0){yx=yx[-index,]}
		n = dim(yx)[1]
		return(list(yx = yx, m = m, n=n))
		}


	MR_sum = function(yx,result="OR"){
		MRinput = mr_input(bx = yx$BETA.x ,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y)
		IVWObject1 = mr_ivw(MRinput,model="fixed")
		IVWObject2 = mr_ivw(MRinput,model="random")
		EggerObject <- mr_egger(MRinput,robust = FALSE,penalized = FALSE,correl = FALSE,distribution = "t-dist",alpha = 0.05)
		MaxlikObject <- mr_maxlik(MRinput,correl = FALSE,distribution = "normal",alpha = 0.05)
		weighted_MedianObject <- mr_median(MRinput, weighting = "weighted",distribution = "normal", alpha = 0.05, iterations = 10000,seed = 314159265)
		simple_MedianObject <- mr_median(MRinput, weighting = "simple",distribution = "normal", alpha = 0.05, iterations = 10000,seed = 314159265)
		#if(IVWObject1$Heter.Stat[2]<0.05){IVWObject = IVWObject2}else{IVWObject = IVWObject1}
		a = c(IVWObject1$Estimate,IVWObject1$StdError,IVWObject1$CILower,IVWObject1$CIUpper,IVWObject1$Pvalue,IVWObject1$Heter.Stat)
		f = c(IVWObject2$Estimate,IVWObject2$StdError,IVWObject2$CILower,IVWObject2$CIUpper,IVWObject2$Pvalue,IVWObject2$Heter.Stat)
		#a = c(IVWObject$Estimate,IVWObject$CILower,IVWObject$CIUpper,IVWObject$Pvalue,IVWObject$Heter.Stat)
		b = c(EggerObject$Estimate,EggerObject$StdError.Est,EggerObject$CILower.Est,EggerObject$CIUpper.Est,EggerObject$Pvalue.Est,EggerObject$Heter.Stat)
		c = c(EggerObject$Intercept,EggerObject$StdError.Int,EggerObject$CILower.Int,EggerObject$CIUpper.Int,EggerObject$Pvalue.Int)
		d = c(MaxlikObject$Estimate,MaxlikObject$StdError,MaxlikObject$CILower,MaxlikObject$CIUpper,MaxlikObject$Pvalue,MaxlikObject$Heter.Stat)
		e = c(weighted_MedianObject$Estimate,weighted_MedianObject$StdError,weighted_MedianObject$CILower,weighted_MedianObject$CIUpper,weighted_MedianObject$Pvalue)
		g = c(simple_MedianObject$Estimate,simple_MedianObject$StdError,simple_MedianObject$CILower,simple_MedianObject$CIUpper,simple_MedianObject$Pvalue)
		if(result == "OR"){
			temp_a = c(a[1:2],exp(a[c(1,3:4)]),a[5:7]);
			temp_b = c(b[1:2],exp(b[c(1,3:4)]),b[5:7]);
			temp_c = c(c[1:2],c[c(1,3:4)],c[5],NA,NA) ;
			temp_d = c(d[1:2],exp(d[c(1,3:4)]),d[5:7]);
			temp_e = c(e[1:2],exp(e[c(1,3:4)]),e[5],NA,NA);
			temp_f = c(f[1:2],exp(f[c(1,3:4)]),f[5:7]);
			temp_g = c(g[1:2],exp(g[c(1,3:4)]),g[5],NA,NA)
			res = rbind(temp_a,temp_f,temp_b,temp_c,temp_d,temp_e,temp_g)
			rownames(res) = c("IVW_fix","IVW_random","Egger_slope","Egger_intercept","Maxlik","weighted_Median","simple_Median")		
			colnames(res) = c("BETA","SE","OR","Down","Up","Pval","Heter.est","Heter.Pval")
		}else if(result == "BETA"){
			temp_a = c(a[1:4],a[5:7]);
			temp_f = c(f[1:4],f[5:7]);
			temp_b = c(b[1:4],b[5:7]);
			temp_c = c(c[1:4],c[5],NA,NA);
			temp_d = c(d[1:4],d[5:7]);
			temp_e = c(e[1:4],e[5],NA,NA);
			temp_g = c(g[1:4],g[5],NA,NA)
			res = rbind(temp_a,temp_f,temp_b,temp_c,temp_d,temp_e,temp_g)
			rownames(res) = c("IVW_fix","IVW_random","Egger_slope","Egger_intercept","Maxlik","weighted_Median","simple_Median")
			colnames(res) = c("BETA","SE","Down","Up","Pval","Heter.est","Heter.Pval")
			}
		return(list(IVW_fix = temp_a, IVW_random = temp_f, Egger_est = temp_b, Egger_int =temp_c, Maxlik = temp_d, weighted_Median = temp_e, simple_Median = temp_g, result = res))
	}

