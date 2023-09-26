	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
	RES = read.table(paste("MAGMAKEGG.txt",sep=""),header = T,sep="\t")[,c("group","gene.set","beta","se","p")]
	colnames(RES) = c("exp","Pathway","BETA","SE","P")

	library(ggplot2)
	library(RColorBrewer)
	dd = RES
	dd$exp=factor(dd$exp,levels = c("AOA","HT","PBC","IBD","CD","RA","MS"), 
			labels = c("B-ALL&AOA","B-ALL&HT","B-ALL&PBC","B-ALL&IBD","B-ALL&CD","B-ALL&RA","B-ALL&MS"))
	dd$Pathway = factor(dd$Pathway)
	dd$logP=(-1)*log10(dd$P)
	p = ggplot(data=dd,
		   aes(x = Pathway,y = logP))+
	  geom_bar(stat="identity",color="grey",fill="#277DA1")+
	  xlab('Pathways')+ ylab("-log10(P value)")+
	  facet_grid(exp~.,scales = "free", space = "free") +
	  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "blue") +
	  geom_hline(yintercept = 5.490689, linetype = "dashed", color = "red") +
	  #scale_y_continuous(limits = c(0,1.5))+
	  theme_bw() + 
	  theme(legend.position = "None") +
	  theme(panel.background = element_rect(fill = "white"),
		panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
		panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
		axis.text = element_text(size = 15),
		axis.title = element_text(size = 15),
		strip.text = element_text(size=15),
		strip.text.x = element_text(margin = margin(b = 7),size = 15),#x label
		axis.text.x = element_text(margin = margin(l = 3),size = 15),
		#axis.text.y = element_blank(),
		#axis.ticks.y = element_blank(),
		legend.text = element_text(size = 13),
		axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +  
		coord_flip()+ guides(shape = guide_legend(order = 0))
	p

	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
	ggsave("Figure Pathways.pdf",p,width = 15,height = 15,dpi=300)





	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
	RES = read.table(paste("MAGMATissue.txt",sep=""),header = T,sep="\t")[,c("Group","Tissues","BETA","SE","P")]
	colnames(RES) = c("exp","Tissues","BETA","SE","P")

	library(ggplot2)
	library(RColorBrewer)
	dd = RES
	dd$exp=factor(dd$exp,levels = c("B-ALL&AOA","B-ALL&HT","B-ALL&PBC","B-ALL&IBD","B-ALL&CD","B-ALL&RA","B-ALL&MS"), 
			labels = c("B-ALL&AOA","B-ALL&HT","B-ALL&PBC","B-ALL&IBD","B-ALL&CD","B-ALL&RA","B-ALL&MS"))
	dd$Tissues = factor(dd$Tissues)
	dd$logP=(-1)*log10(dd$P)
	p = ggplot(data=dd,
		   aes(x = Tissues,y = logP,fill=logP))+
	  geom_bar(stat="identity",color="grey",fill="#277DA1")+
	  xlab('Tissues')+ ylab("-log10(P value)")+
	  facet_grid(exp~.,scales = "free", space = "free") +
	  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "blue") +
	  geom_hline(yintercept = 3.033424, linetype = "dashed", color = "red") +
	  #scale_fill_gradient(low = "lightblue", high = "darkblue") + 
	  #scale_y_continuous(limits = c(0,1.5))+
	  theme_bw() + 
	  theme(legend.position = "None") +
	  theme(panel.background = element_rect(fill = "white"),
		panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
		panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
		axis.text = element_text(size = 15),
		axis.title = element_text(size = 15),
		strip.text = element_text(size=15),
		strip.text.x = element_text(margin = margin(b = 7),size = 15),#x label
		axis.text.x = element_text(margin = margin(l = 3),size = 15),
		#axis.text.y = element_blank(),
		#axis.ticks.y = element_blank(),
		legend.text = element_text(size = 13),
		axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +  
		coord_flip()+ guides(shape = guide_legend(order = 0))
	p

	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
	ggsave("Figure Tissues.pdf",p,width = 12,height = 13,dpi=300)






	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
	RES = read.table(paste("PPI.txt",sep=""),header = T,sep="\t")[,c("Group","Description","Log10.P.")]
	colnames(RES) = c("exp","Pathways","P")

	library(ggplot2)
	library(RColorBrewer)
	dd = RES
	dd$exp=factor(dd$exp)
	dd$Pathways = factor(dd$Pathways,levels = rev(dd$Pathways))
	dd$logP=(-1)*dd$P
	p = ggplot(data=dd,
		   aes(x = Pathways,y = logP,fill=exp))+
	  geom_bar(stat="identity",color="grey")+
	  xlab('Tissues')+ ylab("-log10(P value)")+
	  #facet_grid(exp~.,scales = "free", space = "free") +
	  #geom_hline(yintercept = 1.30103, linetype = "dashed", color = "blue") +
	  #geom_hline(yintercept = 3.033424, linetype = "dashed", color = "red") +
	  #scale_fill_gradient(low = "lightblue", high = "darkblue") + 
	  #scale_y_continuous(limits = c(0,1.5))+
	  theme_bw() + 
	  #theme(legend.position = "None") +
	  theme(panel.background = element_rect(fill = "white"),
		panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
		panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
		axis.text = element_text(size = 15),
		axis.title = element_text(size = 15),
		strip.text = element_text(size=15),
		strip.text.x = element_text(margin = margin(b = 7),size = 15),#x label
		axis.text.x = element_text(margin = margin(l = 3),size = 15),
		#axis.text.y = element_blank(),
		#axis.ticks.y = element_blank(),
		legend.text = element_text(size = 13),
		axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +  
		coord_flip()+ guides(shape = guide_legend(order = 0))
	p

	setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
	ggsave("Figure PPI 2.pdf",p,width = 12,height = 5,dpi=300)


