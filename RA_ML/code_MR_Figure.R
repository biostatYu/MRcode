#---------------
# RA_ML
#---------------

	library(data.table)
	library(MendelianRandomization)
	library(progress)
	setwd("E:\\Cloud_disk\\program")
	source("MRfunction.R")
	pb <- progress_bar$new(total = 100)
	RES = NULL
	setwd("F:\\summarydata")
	DD0 = data.frame(fread(paste("Blood_Pheweb_EUR_Malignant_lymphoma.txt.gz",sep=""),header=TRUE))[,1:10]

	x.1 = c("AD_Rheumatoid_Arthritis_pos_European_2022_NG")
	m=length(x.1)
	pb <- progress_bar$new(total = m)
	rex = matrix(NA, m, 54)
	j=1
	setwd("F:\\summarydata")
	ax0 = data.frame(fread(paste(x.1[j],".txt.gz",sep=""),header=TRUE))[,1:10]
	setwd("F:\\summarydata\\index_snp")
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

	#			      BETA          SE        OR        Down         Up        Pval Heter.est  Heter.Pval
	#	IVW_fix         0.07348641 0.026998092 1.0762539  1.02078417 1.13473790 0.006490557  99.55418 0.007577511
	#	IVW_random      0.07348641 0.032666934 1.0762539  1.00950529 1.14741596 0.024476533  99.55418 0.007577511
	#	Egger_slope     0.03752836 0.076419669 1.0382414  0.89136085 1.20932537 0.624972898  99.15230 0.006519697
	#	Egger_intercept 0.00511850 0.009822193 0.0051185 -0.01448668 0.02472368 0.604004328        NA          NA
	#	weighted_mode   0.06230490 0.054289072 1.0642868  0.95685853 1.18377624 0.251112688        NA          NA
	#	weighted_Median 0.08172315 0.049069621 1.0851553  0.98565239 1.19470327 0.095822382        NA          NA
	#	DIVW            0.07486847 0.032795230 1.0777424  1.01064728 1.14929182 0.022435827        NA          NA
	#	MR-RAPS         0.07556980 0.027416593 1.0784985  1.02207337 1.13803868 0.005844992        NA          NA

	library(MRPRESSO)
	fit=mr_presso(BetaOutcome="BETA.y",BetaExposure="BETA.x",SdOutcome="SE.y",SdExposure="SE.x",
	OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=yx,NbDistribution=1000,SignifThreshold=0.05) 
	#	$`Main MR results`
	#	  Exposure       MR Analysis Causal Estimate         Sd   T-stat    P-value
	#	1   BETA.x               Raw      0.07995005 0.03247974 2.461536 0.01694035
	#	2   BETA.x Outlier-corrected              NA         NA       NA         NA
	#
	#	$`MR-PRESSO results`
	#	$`MR-PRESSO results`$`Global Test`
	#	$`MR-PRESSO results`$`Global Test`$RSSobs
	#	[1] 62.45819
	#
	#	$`MR-PRESSO results`$`Global Test`$Pvalue
	#	[1] 0.317


	MRinput = mr_input(bx = yx$BETA.x ,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y)
	IVWObject = mr_ivw(MRinput,model="fixed")
	EggerObject <- mr_egger(MRinput,robust = FALSE,penalized = FALSE,correl = FALSE,distribution = "t-dist",alpha = 0.05)
	weighted_MedianObject <- mr_median(MRinput, weighting = "weighted",distribution = "normal", alpha = 0.05, iterations = 10000,seed = 314159265)

##### a
	a = c(IVWObject$Estimate,0)
	b = c(EggerObject$Estimate,EggerObject$Intercept)
	#d = c(weighted_MedianObject$Estimate,0)
	m = c("IVW","MR-Egger")
	c = data.frame(m,rbind(a,b));colnames(c) = c("Method","Beta","Int")

	pa1 = ggplot(data = yx, aes(x = BETA.x, y = BETA.y)) +
	    geom_point(colour = "black", alpha = 0.5, size = 2) +
		theme_bw() + #去掉背景色
		theme(panel.background = element_rect(fill = "white"),
		panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
		panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
		axis.text = element_text(size = 15),
		axis.title = element_text(size = 15),
		axis.text.x  = element_text(margin = margin(b = 7),size = 15),
		axis.text.y  = element_text(margin = margin(l = 3),size = 15),
		legend.title=element_text(size=12),
		 legend.text=element_text(size=12)) +
	    xlab(paste("Rheumatoid arthritis")) +
	    ylab(paste("Malignant lymphoma"))  + 
	    geom_hline(yintercept=0, linetype = "dashed",size = 0.3,colour = "gray50") + 
	    geom_vline(xintercept=0, linetype = "dashed",size = 0.3,colour = "gray50") +
#	    scale_x_continuous(limits = c(-0.35,0.3))+
#	    scale_y_continuous(limits = c(-0.2,0.43))+
	    geom_errorbar(aes(ymin = BETA.y - qnorm(0.975)*SE.y, ymax = BETA.y + qnorm(0.975)*SE.y), alpha = 0.3, colour = "black") +
	    geom_errorbarh(aes(xmin = BETA.x - qnorm(0.975)*SE.x, xmax = BETA.x + qnorm(0.975)*SE.x), alpha = 0.3, colour = "black") +
	    geom_abline(data = c, 
		      aes(intercept = Int, slope = Beta, linetype = Method, color = Method), 
		      show.legend = T, size = 1)+
	    labs(tag = '(A)')
##### b
	
	res2 = matrix(NA,dim(yx)[1],2)
	for( z in 1:dim(yx)[1]){
	  temp = yx[z,]
	  MRinput2 = mr_input(bx = temp$BETA.x ,bxse = temp$SE.x,by = temp$BETA.y,byse = temp$SE.y)
	  IVWObject2 = mr_ivw(MRinput2,model="fixed")
	  res2[z,] = c(IVWObject2$Estimate,IVWObject2$StdError)
	}
	res2 = data.frame(res2)
	colnames(res2) = c("beta","se")
	pb1 = ggplot(data = res2, aes(x = beta, y = 1/se)) +
		geom_point(colour = "black", alpha = 0.4,size = 2) +
		theme_bw() + #去掉背景色
		theme(panel.background = element_rect(fill = "white"),
		panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
		panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
		axis.text = element_text(size = 15),
		axis.title = element_text(size = 15),
		axis.text.x  = element_text(margin = margin(b = 7),size = 15),
		axis.text.y  = element_text(margin = margin(l = 3),size = 15),
		legend.title=element_text(size=12),
		 legend.text=element_text(size=12)) +
	  #scale_x_continuous(limits = c(-0.005,0.005))+
	  #scale_y_continuous(limits = c(0,15))+
	  xlab(paste("Individual causal effect")) +
	  ylab(paste("1/(SE of casual effect)"))  + 
	  geom_errorbarh(aes(xmin = beta - qnorm(0.975)*se, xmax = beta + qnorm(0.975)*se), alpha = 0.5) +
	  geom_vline(xintercept=IVWObject$Estimate,size = 1, linetype = "dashed",colour = "red",alpha = 0.8) +
	    labs(tag = '(B)')

	library(patchwork)
	pp1 = pa1+pb1+plot_layout(ncol=2)
	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\RA_ML\\mr")
	ggsave("Figure scatter and funnel plot.pdf",pp1,width = 12,height = 5,dpi=300)


#---------------
# RA_ML
#---------------

	library(data.table)
	library(MendelianRandomization)
	library(progress)
	setwd("E:\\Cloud_disk\\program")
	source("MRfunction.R")
	exp = c("AD_Rheumatoid_Arthritis_pos_European_2022_NG")
	out = c("Blood_Pheweb_EUR_Malignant_lymphoma",
		"Blood_finnUKBB_C3_CLL_EXALLC_GRCh37_qc",
		"Blood_finnUKBB_C3_DLBCL_EXALLC_GRCh37_qc",
		"Blood_finnUKBB_CD2_FOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
		"Blood_finnUKBB_CD2_HODGKIN_LYMPHOMA_EXALLC_GRCh37_qc",
		"Blood_finnUKBB_CD2_TNK_LYMPHOMA_EXALLC_GRCh37_qc",
		"Blood_finnUKBB_CD2_NONFOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc",
		"Blood_finnUKBB_CD2_NONHODGKIN_NAS_EXALLC_GRCh37_qc")
	exp.name = c("RA")
	out.name = c("ML","CLL","DLBCL","FL","HL","TNKL","NFL","NHL")
	RES = NULL
	pica = list()
	picb = list()
	for(i in 1:length(exp)){
	for(j in 1:length(out)){

		setwd("F:\\summarydata")
		DD0 = data.frame(fread(paste(out[j],".txt.gz",sep=""),header=TRUE))[,1:10]
		m=length(exp)
		pb <- progress_bar$new(total = m)
		rex = matrix(NA, m, 54)
		setwd("F:\\summarydata")
		ax0 = data.frame(fread(paste(exp[i],".txt.gz",sep=""),header=TRUE))[,1:10]
		setwd("F:\\summarydata\\index_snp")
		snp = as.character(fread(paste(exp[i],"_5E-8_0.001_10000_clump.txt",sep=""),header=TRUE)$SNP)
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
		RES = rbind(RES,data.frame(Exp=exp[i],Out=out[j],Res))

		library(MRPRESSO)
		fit=mr_presso(BetaOutcome="BETA.y",BetaExposure="BETA.x",SdOutcome="SE.y",SdExposure="SE.x",
			OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=yx,NbDistribution=5000,SignifThreshold=0.05) 
		print(fit)


		MRinput = mr_input(bx = yx$BETA.x ,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y)
		IVWObject = mr_ivw(MRinput,model="fixed")
		EggerObject <- mr_egger(MRinput,robust = FALSE,penalized = FALSE,correl = FALSE,distribution = "t-dist",alpha = 0.05)
		weighted_MedianObject <- mr_median(MRinput, weighting = "weighted",distribution = "normal", alpha = 0.05, iterations = 10000,seed = 314159265)
		indd = which(yx$SE.y>0.5); if(length(indd)>0){yx = yx[-indd,]}
	##### a
		a = c(IVWObject$Estimate,0)
		b = c(EggerObject$Estimate,EggerObject$Intercept)
		#d = c(weighted_MedianObject$Estimate,0)
		m = c("IVW","MR-Egger")
		c = data.frame(m,rbind(a,b));colnames(c) = c("Method","Beta","Int")

		pica[[paste(exp[i],out[j],sep="_")]] = ggplot(data = yx, aes(x = BETA.x, y = BETA.y)) +
		    geom_point(colour = "black", alpha = 0.5, size = 2) +
			theme_bw() + #去掉背景色
			theme(panel.background = element_rect(fill = "white"),
			panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
			panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
			axis.text = element_text(size = 15),
			axis.title = element_text(size = 15),
			axis.text.x  = element_text(margin = margin(b = 7),size = 15),
			axis.text.y  = element_text(margin = margin(l = 3),size = 15),
			legend.title=element_text(size=12),
			 legend.text=element_text(size=12)) +
		    xlab(paste(exp.name[i])) +
		    ylab(paste(out.name[j]))  + 
		    geom_hline(yintercept=0, linetype = "dashed",size = 0.3,colour = "gray50") + 
		    geom_vline(xintercept=0, linetype = "dashed",size = 0.3,colour = "gray50") +
	#	    scale_x_continuous(limits = c(-0.35,0.3))+
	#	    scale_y_continuous(limits = c(-0.2,0.43))+
		    geom_errorbar(aes(ymin = BETA.y - qnorm(0.975)*SE.y, ymax = BETA.y + qnorm(0.975)*SE.y), alpha = 0.3, colour = "black") +
		    geom_errorbarh(aes(xmin = BETA.x - qnorm(0.975)*SE.x, xmax = BETA.x + qnorm(0.975)*SE.x), alpha = 0.3, colour = "black") +
		    geom_abline(data = c, 
			      aes(intercept = Int, slope = Beta, linetype = Method, color = Method), 
			      show.legend = T, size = 1)+
		    labs(tag = '(A)')
	##### b
		
		res2 = matrix(NA,dim(yx)[1],2)
		for( z in 1:dim(yx)[1]){
		  temp = yx[z,]
		  MRinput2 = mr_input(bx = temp$BETA.x ,bxse = temp$SE.x,by = temp$BETA.y,byse = temp$SE.y)
		  IVWObject2 = mr_ivw(MRinput2,model="fixed")
		  res2[z,] = c(IVWObject2$Estimate,IVWObject2$StdError)
		}
		res2 = data.frame(res2)
		colnames(res2) = c("beta","se")
		picb[[paste(exp[i],out[j],sep="_")]] = ggplot(data = res2, aes(x = beta, y = 1/se)) +
			geom_point(colour = "black", alpha = 0.4,size = 2) +
			theme_bw() + #去掉背景色
			theme(panel.background = element_rect(fill = "white"),
			panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
			panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
			axis.text = element_text(size = 15),
			axis.title = element_text(size = 15),
			axis.text.x  = element_text(margin = margin(b = 7),size = 15),
			axis.text.y  = element_text(margin = margin(l = 3),size = 15),
			legend.title=element_text(size=12),
			 legend.text=element_text(size=12)) +
		  #scale_x_continuous(limits = c(-0.005,0.005))+
		  #scale_y_continuous(limits = c(0,15))+
		  xlab(paste("Individual causal effect")) +
		  ylab(paste("1/(SE of casual effect)"))  + 
		  geom_errorbarh(aes(xmin = beta - qnorm(0.975)*se, xmax = beta + qnorm(0.975)*se), alpha = 0.5) +
		  geom_vline(xintercept=IVWObject$Estimate,size = 1, linetype = "dashed",colour = "red",alpha = 0.8) +
		    labs(tag = '(B)')
	}
	}
	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\RA_ML\\mr")
	write.csv(RES,"result_MR.csv")

	library(patchwork)
	pp1 = pica[[paste(exp[1],out[1],sep="_")]]+picb[[paste(exp[1],out[1],sep="_")]]+plot_layout(ncol=2)
	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\RA_ML\\mr")
	ggsave("Figure scatter and funnel plot 1.pdf",pp1,width = 12,height = 5,dpi=300)

	pp2 =	pica[[paste(exp[1],out[1],sep="_")]]	+	
		pica[[paste(exp[1],out[2],sep="_")]]	+	pica[[paste(exp[1],out[3],sep="_")]]+
		pica[[paste(exp[1],out[4],sep="_")]]	+	pica[[paste(exp[1],out[5],sep="_")]]+
		pica[[paste(exp[1],out[6],sep="_")]]	+	pica[[paste(exp[1],out[7],sep="_")]]+
		pica[[paste(exp[1],out[8],sep="_")]]	+	plot_layout(ncol=4)
	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\RA_ML\\mr")
	ggsave("Figure scatter and funnel plot 2.pdf",pp2,width = 24,height = 8,dpi=300)

	pp3 =	picb[[paste(exp[1],out[1],sep="_")]]	+	
		picb[[paste(exp[1],out[2],sep="_")]]	+	picb[[paste(exp[1],out[3],sep="_")]]+
		picb[[paste(exp[1],out[4],sep="_")]]	+	picb[[paste(exp[1],out[5],sep="_")]]+
		picb[[paste(exp[1],out[6],sep="_")]]	+	picb[[paste(exp[1],out[7],sep="_")]]+
		picb[[paste(exp[1],out[8],sep="_")]]	+	plot_layout(ncol=4)
	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\RA_ML\\mr")
	ggsave("Figure scatter and funnel plot 3.pdf",pp3,width = 18,height = 8,dpi=300)
