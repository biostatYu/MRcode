library(ggplot2)
library(reshape)
library(grDevices)
library(cowplot)
library(fdrtool)
library(RColorBrewer)
library(dplyr)
setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Results\\")
cc = read.table("result_all_gene_combine.txt",header = T,sep="\t")
dd = cc[which(cc$P<0.05),c(1,33:50)]
paste(paste(colnames(dd)[3:19],"=sum(",colnames(dd)[3:19],")",sep=""),collapse = ",")
dd <- dd %>%
	group_by(SYMBOL) %>%
	summarize(sig.placo=sum(sig.placo),
		positionmapping=sum(positionmapping),
		eQTL_B_cell_naive=sum(eQTL_B_cell_naive),
		eQTL_CEDAR_B.cell_CD19=sum(eQTL_CEDAR_B.cell_CD19),
		eQTL_Cells_EBV.transformed_lymphocytes=sum(eQTL_Cells_EBV.transformed_lymphocytes),
		eQTL_eQTLGen_cis_eQTLs=sum(eQTL_eQTLGen_cis_eQTLs),
		eQTL_eQTLGen_trans_eQTLs=sum(eQTL_eQTLGen_trans_eQTLs),
		eQTL_Spleen=sum(eQTL_Spleen),eQTL_Whole_Blood=sum(eQTL_Whole_Blood),
		SMR_AD_Cells_EBV.transformed_lymphocytes=sum(SMR_AD_Cells_EBV.transformed_lymphocytes),
		SMR_AD_eQTLgen=sum(SMR_AD_eQTLgen),
		SMR_AD_Spleen=sum(SMR_AD_Spleen),
		SMR_AD_Whole_Blood=sum(SMR_AD_Whole_Blood),
		SMR_BALL_Cells_EBV.transformed_lymphocytes=sum(SMR_BALL_Cells_EBV.transformed_lymphocytes),
		SMR_BALL_eQTLgen=sum(SMR_BALL_eQTLgen),
		SMR_BALL_Spleen=sum(SMR_BALL_Spleen),
		SMR_BALL_Whole_Blood=sum(SMR_BALL_Whole_Blood))
colnames(dd)[2] = c("Loci_signal")
#dd$Trait = factor(dd$Trait,levels = c("AD_Asthma_adultonset_2019AJHG_maf_GRCh37",
#					"AD_Hypothyroidism_2021Sakaue",
#					"AD_PBC_2021Cordell",
#					"AD_Inflammatory_bowel_disease_2017LangeNG",
#					"AD_Crohn_disease_2017LangeNG",					
#					"AD_Rheumatoid_Arthritis_European_2022_NG",
#					"AD_Multiple_Sclerosis_2018_qc"), labels = c("B-ALL&AOA","B-ALL&HT","B-ALL&PBC","B-ALL&IBD","B-ALL&CD","B-ALL&RA","B-ALL&MS"))

#dd$Loci_signal = factor(dd$Loci_signal)

ee = data.frame(dd)
c0 = dd$Loci_signal; c0[which(c0>1)] = 1
c1 = apply(ee[,c(4:10)],1,sum); c1[which(c1>1)] = 1
c2 = apply(ee[,c(11:18)],1,sum); c2[which(c2>1)] = 1
c3 = c0 + c1 + c2
ind = which(c3 == 3)
ee = ee[ind,-c(3,8,12,16)] 
#ee[,-1][which(ee[,-1]>1)] = 1
ff = melt(ee, id.var = c("SYMBOL"))
ff$variable = factor(ff$variable,levels = colnames(dd)[-c(1)])
ee$SYMBOL = factor(ee$SYMBOL,levels = ee$SYMBOL)

label_text = ee$SYMBOL
theme_change <- theme(
	  plot.background = element_blank(),
	  panel.grid.minor = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.background = element_blank(),
	  panel.border = element_blank(),
	  axis.line = element_blank(),
	  axis.ticks = element_blank(),
	  axis.text.x  = element_text(margin = margin(l = 3),size = 13,colour="black",angle = 45, hjust = 1),
	  axis.text.y  = element_text(margin = margin(l = 3),size = 13,colour="black"),
	  axis.title.x = element_blank(),
	  axis.title.y = element_blank(),
	  legend.title = element_text(size = 13, face = "bold"),
	  legend.text = element_text(size = 13),
	   legend.key.size=unit(0.7,'cm')
)

## output the graphics
F1 = ggplot(ff, aes(y = SYMBOL, x = variable, fill = factor(value))) +
	geom_tile(color = "grey", size = 0.5) +
	scale_fill_brewer(labels = c(0,1,2,3,4,6), palette="PuBu") +
	theme_change +
	scale_y_discrete(limits = rev(unique(label_text)), position = "left") +
	guides(fill = guide_legend(title="")) #+theme(legend.position="top")



setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
ggsave("Figure gene heatmap.pdf",F1,width = 9,height = 13,dpi=300)






