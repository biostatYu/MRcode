library(MendelianRandomization)
library(ggplot2)
library(patchwork)
setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
dd=data.frame(fread("Figure MR.txt",head=T,sep='\t'))
dd$MR=factor(dd$MR,levels = c(1),labels=c("B-ALL as outcome"))
dd$Sig=factor(dd$Sig,levels = c(1,0),labels=c("Sig","Nosig"))
dd$Traits = factor(dd$Traits, levels=rev(c("Adult-onset asthma (AOA)",
	"Childhood-onset asthma (COA)",
	"Grave's disease (GD)",
	"Hashimoto's disease (HD)",
	"Hypothyroidism (HT)",
	"Primary biliary cholangitis (PBC)",
	"Primary sclerosing cholangitis (PSC)",
	"Inflammatory bowel disease (IBD)",
	"Crohn's disease (CD)",
	"Ulcerative Colitis (UC)",
	"Rheumatoid arthritis (RA)",
	"Multiple sclerosis (MS)",
	"Systemic sclerosis (SS)",
	"Systemic lupus erythematosus (SLE)",
	"Type 1 diabetes (T1D)",
	"Vitiligo")))
limits <- aes(ymax = dd$Up, ymin=dd$Down)

p = ggplot(data=dd,
           aes(x = Traits,y = OR, ymin = Down, ymax = Up, color = Sig))+
  geom_pointrange(aes(color=Sig), linewidth = 0.8, alpha = 1)+
  geom_hline(yintercept = 1, linetype=2)+
  scale_colour_manual(values = c("Sig" = "#CC0000", "Nosig" = "#0099FF"))+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=Down, ymax=Up),width=0.3,cex=1)+ 
  facet_grid(.~MR, scales="free") +
  #scale_y_continuous(limits = c(0,1.5))+
  theme_bw() + 
  theme(legend.position = "none",
	panel.background = element_rect(fill = "white"),
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
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  labs(x = "", y = "95% confidence intervals")  +  coord_flip()+ guides(shape = guide_legend(order = 0))
p

setwd("E:\\Cloud_disk\\program\\Hematology\\Pleiotropy\\AD_BALL\\Figure")
ggsave("Figure MR.pdf",p,width = 7,height = 6,dpi=300)


