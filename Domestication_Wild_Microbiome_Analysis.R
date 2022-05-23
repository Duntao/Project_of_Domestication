#Core code for Project of domestication.
#Pls don not hestitate to emial me if your want to conduct deeply communication with us
library(dplyr)
library(ggplot2)
library(tidyverse)
library(igraph)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(psych)
library(corrplot)
library(pheatmap)
library(tidyjson)
library(jsonlite)
library(tidyr)
library(clusterProfiler)
library(readxl)
library(data.table)
library(rio)
library(edgeR)
library(aplot)

#cap
p1 <- plot_ordination(physeq_bac,cap_1,color="type",shape="type")+
  geom_point(size=3)+
  scale_shape_manual(values=c(15,16))+
  scale_colour_brewer(palette="Set2")+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle(bquote(paste(.(round(var_tab_soil16s["constrained","proportion"],2)*100),
                       "% of total variation (CI=",
                       .(round(cca_ci(cap_1)[1],2)*100),"%, ",
                       .(round(cca_ci(cap_1)[2],2)*100),"%); ", "p=0.001" )))+
  theme_bw()+theme(plot.title = element_text(size = 10, face = "bold"),
      legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())

p1

p2 <- plot_ordination(physeq_bac,cap_2,color="compartment",shape="compartment")+
  geom_point(size=3)+
  scale_shape_manual(values=c(15,16))+
  scale_colour_brewer(palette="Set2")+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle(bquote(paste(.(round(var_tab_soil16s["constrained","proportion"],2)*100),"% of total variation (CI=",
                       .(round(cca_ci(cap_2)[1],2)*100),"%, ",
                       .(round(cca_ci(cap_2)[2],2)*100),"%); ", "p=0.001" )))+
  theme_bw()+theme(plot.title = element_text(size = 10, face = "bold"),
                   legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())


p2

p3 <- plot_ordination(physeq_bac,cap_3,color="geno_type",shape="geno_type")+
  geom_point(size=3)+
  scale_shape_manual(values=c(16,16,16,16,16,16))+
  scale_colour_npg()+theme_bw()+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle(bquote(paste(.(round(var_tab_soil16s["constrained","proportion"],2)*100),"% of total variation (CI=",
                       .(round(cca_ci(cap_3)[1],2)*100),"%, ",
                       .(round(cca_ci(cap_3)[2],2)*100),"%); ", "p=0.001" )))+
  theme(plot.title = element_text(size = 10, face = "bold"))
p3

#NMDS
class(forPlot3)
forPlot3$type<-factor(forPlot3$type,levels=c("DT1","DT2","DT3","WT1","WT2","WT3"))
pnm_dtwt_16s<-ggplot(data = forPlot3)+
  labs(x="MDS1",y="MDS2")+
  theme(legend.title=element_blank(),legend.text=element_text(size=10),legend.key.height=unit(0.6,"cm"))+
  geom_point(aes(x=MDS1,y=MDS2,colour=type),size=3)+scale_colour_brewer(palette="Set3")+
  annotate("text",x=-Inf,y=Inf,hjust=-0.2,vjust=2,label="Abundant (Stress: 0.153)")
pnm_dtwt_16s
ggsave("pnm_dtwt_16s.PDF",pnm_dtwt_16s,width = 5,height = 3 )

#the other part of code upon your reasonable request.
