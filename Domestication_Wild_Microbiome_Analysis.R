#Core code for Project of domestication
#Pls don not hestitate to email me if your want to conduct deeply communication with us
#The loading R packages upon yourself system status and dadaset
library(circlize);library(cowplot)
library(grid);library(ggplotify)
library(tidyverse);library(dplyr)
library(vegan);library(ggplot2)
library(multcomp);library(agricolae)
library(devtools);library(ellipse)
library(TSA);library(reshape2)
library(ggsci);library(gridExtra)
library(ggpubr);library(tidyr)
library(plyr);library(dplyr)
library(grid);library(Cairo)
library(ellipse);library(gplots)
library(RColorBrewer);library(png)
library(betapart);library(picante)
library(ggpubr);library(EcolUtils)
library(spaa);library(EcoSimR)
library(VennDiagram);library(parallel)
library(snowfall);library(tidyverse)
library(nycflights13);library(Hmisc)
library(picante);library(phyloseq)
library(ape);library(caper);library(metagenomeSeq)
library(data.table)
library(qvalue);library(igraph)
library(VennDiagram)
#devtools::install_github("Hy4m/linkET")
library(linkET);library(dplyr)
library(patchwork);library(lme4)
library(lmerTest);library(sampling)
library(caret);library(robustbase)
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

#adonis
adonis_result_all <- adonis(formula=(t(df_all)~type*geno_type), df_group_all, distance = 'bray', permutations = 999)
adonis_result_all

##glm
asv_mod_16s<-manyglm(asv_abund_16s~type*compartment*geno_type,data=df_group2,family="negative.binomial")
asv_sum_16s<-summary.manyglm(asv_mod_16s,test="wald",nBoot = 100, p.uni = 'adjusted')
asv_sum_16s
anova.manyglm(asv_mod_16s,nBoot=100)

asv_mod_its<-manyglm(asv_abund_its~type*compartment,data=df_group2,family="negative.binomial")
asv_sum_its<-summary.manyglm(asv_mod_its,test="wald",nBoot = 100, p.uni = 'adjusted')
asv_sum_its
anova.manyglm(asv_mod_its,nBoot=100)

#glm
m1.lme4 = lmer(shannon ~ type+compartment+geno_type+(1|block),data = alp_16s)
summary(m1.lme4)
anova(m1.lme4)

m2.lme4 = lmer(eveness2 ~ type+compartment+(1|block),data = df_alp_group_its)
summary(m2.lme4)
anova(m2.lme4)

#edgeR
##ggplot
rich_all<-rbind(Dome_rich,wild_rich,ns_rich)
View(rich_all)
rich_all$type<-factor(rich_all$type,levels=c("Domestication","Wild","ns"))
p_wd_16s<-ggplot(rich_all,aes(x=logCPM,y=logFC,colour=type))+geom_point(size=0.8)+
  scale_colour_manual(values=c("#2c7fb8","#fec44f","#bdbdbd"))+
  labs(y="Log2 (Fold change)",x="Log2 (Count per million)")+
  ggtitle("Bacteria")+
  theme(legend.title=element_blank(),legend.position="none")
p_wd_16s

rich_all<-rbind(Bulk_rich,Rhizosphere_rich,ns_rich)
View(rich_all)
rich_all$type<-factor(rich_all$type,levels=c("Bulk","Rhizosphere","ns"))
pw_br_16s<-ggplot(rich_all,aes(x=logCPM,y=logFC,colour=type))+geom_point(size=0.8)+
  scale_colour_manual(values=c("#2ca25f","#d95f0e","#bdbdbd"))+
  labs(y="Log2 (Fold change)",x="Log2 (Count per million)")+
  ggtitle("Bacteria")+
  theme(legend.title=element_blank(),legend.position="none")
pw_br_16s

p_final<-ggarrange(pb_wd_16s,pd_br_16s,pt_wd_16s,pw_br_16s,ncol=2,nrow=2,labels=c("a","b","c","d"))
p_final
p_final<-ggarrange(p_wd_16s,p_br_16s,ncol=1,nrow=2,labels=c("a","b"))
ggsave("p_final_all.PDF",p_final,width = 6.5,height = 6 )

#joyplot
p1<-ggplot(gather_rf3, aes(x = id, y = genus_name,height=relative,fill=type,color=relative)) + 
  ggridges::geom_density_ridges2(stat="identity",color="#FFFFFF",alpha=3/4,scale=2.5,size=.1)+
  scale_fill_manual(values=c("#2c7fb8","#fec44f"))+
  scale_x_continuous(breaks=seq(0,60,4))+
  theme(legend.position="top", legend.spacing.x = unit(0.2, 'cm'),
        legend.text=element_text(size=12),legend.title=element_blank(),
        axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )
p1
ggsave("randomforest_bacteria_rhizosphere.pdf",p1,width=4,height=3.5)


#network analysis
#The code for network construct can be found here:https://github.com/scwatts/FastSpar

#node feature
###each node toplogical features
node_top<-function(g){
  betweeness.centrality<-betweenness(g,v=V(g),directed = FALSE, weights = NA,nobigint = TRUE, normalized = FALSE)
  closeness.centrality<-closeness(g,vids=V(g),weights = NA, normalized = FALSE)
  node.transitivity<-transitivity(g,type = c("local"), vids = NULL,weights = NULL)
  ec<-eigen_centrality(g)
  node.degree<-degree(g,v=V(g),mode="all")
  node.topology<-data.frame(betweeness.centrality,closeness.centrality,node.transitivity,ec,node.degree)
  return(node.topology)
}

compare_means(Betweeness_centrality ~ Type, node_ddr_wir, method = "kruskal.test")
compare_means(Closeness_centrality ~ Type, node_ddr_wir, method = "kruskal.test")
compare_means(Node_transitivity ~ Type, node_ddr_wir, method = "kruskal.test")
compare_means(Node_degree ~ Type, node_ddr_wir, method = "kruskal.test")
node_ddr_wir$type<-factor(node_ddr_wir$Type,levels=c("Bacteria","Fungi"))
p1<-ggplot(node_ddr_wir,aes(x=Type,y=Betweeness_centrality,fill=Type))+geom_boxplot()+
  ylab("Betweeness centrality")+geom_jitter(width=0.3,height=0.2,size=0.4)+
  scale_fill_manual(values=c("#3182BD","#99dbc9"))+
  ggtitle("Kruskal.test: P=0.088")+
  theme(
    legend.title=element_blank(),legend.position="none",
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=10,face="bold",colour="black"),
    axis.text.x=element_text(size=8,face="bold",colour="black"),
    axis.text.y=element_text(size=8))
ggsave("p1.pdf",p1,width=2.5,height=2.5)

p2<-ggplot(node_ddr_wir,aes(x=Type,y=Closeness_centrality,fill=Type))+geom_boxplot()+
  ylab("Closeness centrality")+geom_jitter(width=0.3,height=0.2,size=0.4)+
  ggtitle("Kruskal.test: P<0.001")+
  scale_fill_manual(values=c("#3182BD","#99dbc9"))+
  theme(
    legend.title=element_blank(),legend.position="none",
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=10,face="bold",colour="black"),
    axis.text.x=element_text(size=8,face="bold",colour="black"),
    axis.text.y=element_text(size=8))
ggsave("p2.pdf",p2,width=2.5,height=2.5)

p3<-ggplot(node_ddr_wir,aes(x=Type,y=Node_transitivity,fill=Type))+geom_boxplot()+
  ylab("Eigenvector centrality")+geom_jitter(width=0.3,height=0.2,size=0.4)+
  ggtitle("Kruskal.test: P=0.56")+
  scale_fill_manual(values=c("#3182BD","#99dbc9"))+
  theme(
    legend.title=element_blank(),legend.position="none",
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=10,face="bold",colour="black"),
    axis.text.x=element_text(size=8,face="bold",colour="black"),
    axis.text.y=element_text(size=8))
ggsave("p3.pdf",p3,width=2.5,height=2.5)

p4<-ggplot(node_ddr_wir,aes(x=Type,y=Node_degree,fill=Type))+geom_boxplot()+
  ylab("Node degree")+geom_jitter(width=0.3,height=0.2,size=0.4)+
  ggtitle("Kruskal.test: P<0.001")+
  scale_fill_manual(values=c("#3182BD","#99dbc9"))+
  theme(
    legend.title=element_blank(),legend.position="none",
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=10,face="bold",colour="black"),
    axis.text.x=element_text(size=8,face="bold",colour="black"),
    axis.text.y=element_text(size=8))
p4
ggsave("p4.pdf",p4,width=2.5,height=2.5)


p5<-grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
p5
ggsave("p5_wt.pdf",p5,width=4.5,height=4)
dev.off()

##metabolism analysis#######################################
##aplot for ggtree_heatmap
tree3<-aa
str(tree3)
View(tree3$tip.label)
tree_tl<-as.data.frame(tree3$tip.label)
colnames(tree_tl)<-"id"
View(tree_tl)
View(otu_dt4_16s_re)
class(otu_dt4_16s_re)
otu_dt4_16s_re$id<-rownames(otu_dt4_16s_re)
tree_tl2<-dplyr::left_join(tree_tl,otu_dt4_16s_re,by="id")

View(tree_tl2)
tree3$tip.label<-as.vector(tree_tl2$Genus)
str(tree3)
##ggtree
ggtree_plot<-ggtree::ggtree(tree3,aes(color=group),branch.length='none')%<+%dd2+geom_tippoint(aes(shape=Phylum),size=3)+
  geom_tiplab(aes(label=Genus),offset=0.5,size=2,align=TRUE)+
  scale_shape_manual(values=c(15:19))
ggtree_plot

ggtree_plot2<-ggtree::ggtree(tree3,aes(color=group),branch.length='none')
ggtree_plot2

##ggplot-heatmap
View(otu_dt4_16s_re)
df_otu<-corr_dt_16s_mp
View(df_otu)
df_corr<-t(df_otu)
rownames(df_corr)<-rownames(otu_dt4_16s_re)
View(df_corr);max(df_corr);min(df_corr)
class(df_corr)
df_corr<-as.data.frame(df_corr)
df_corr$ID<-rownames(df_corr)
View(df_corr)
df_corr_gat<-gather(df_corr,metabolite,value,-"ID")
View(df_corr_gat)


df_corr2<-t(df_otu)
View(df_corr2)
df_corr2<-as.data.frame(df_corr2)
df_corr2$ID<-rownames(df_corr2)
View(df_corr2)
df_corr_gat2<-gather(df_corr2,metabolite,value,-"ID")
View(df_corr_gat2)
df_corr_gat$metabolite<-factor(df_corr_gat$metabolite,levels=colnames(df_corr))

g_map<-ggplot(df_corr_gat2,aes(metabolite,ID))+
  geom_tile(aes(fill=value),colour="white")+
  scale_fill_gradient2(name="value",mid="white",
                       low="#009e73",high="#56b4e9")+
  scale_y_discrete(position = "right")+
  theme_bw()+
  theme(axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(hjust=0,angle=-90)
        )
g_map
g_map %>% insert_left(ggtree_plot2,width=0.28)

View(df_otu)
df_corr3<-df_otu
df_corr3[df_corr3<0]<-0
df_corr3[df_corr3>0]<-1
class(df_corr3)

#positive
df_corr3<-as.data.frame(df_corr3)
cor_p<-as.data.frame(apply(df_corr3,1,sum))
colnames(cor_p)<-"Positive"
View(cor_p)

#negative
df_corr4<-df_otu
df_corr4[df_corr4>0]<-0
df_corr4[df_corr4<0]<-1
class(df_corr4)
df_corr4<-as.data.frame(df_corr4)
cor_n<-as.data.frame(apply(df_corr4,1,sum))
colnames(cor_n)<-"Negative"
View(cor_n)
cor_pn<-cbind(cor_p,cor_n)
View(cor_pn)
cor_pn$profile<-rownames(cor_pn)
cor_pn_gg<-gather(cor_pn,pn,Count,-profile)
View(cor_pn_gg)
p_pn<-ggplot(cor_pn_gg,aes(x=profile,y=Count,fill=pn))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#FA9FB5","#67a9cf"))+
  scale_y_continuous(position="right")+
  theme_bw()+
  theme(legend.title=element_blank(),
    axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
    axis.title.y=element_text())
p_pn

p_tree_final<-g_map %>% insert_left(ggtree_plot2,width=0.28)%>%
  insert_top(p_pn,height=0.18)
p_tree_final
ggsave("p_tree_dt_16s_final.pdf",p_tree_final,height=6,width=11)

#kegg
p_c<-pheatmap(df_c5,cluster_row=T,cluster_col=F,gaps_row=c(5),colorRampPalette(c("black","white","firebrick3"))(80),
              gaps_col=c(9),legend_breaks = -3:3,cutree_rows = 3,annotation_row=df_c4[3],annotation_col=annotation_col[2])
p_c
ggsave("p_c.pdf",p_c,width=9,height=5)

p_n<-pheatmap(df_n5,cluster_row=T,cluster_col=F, gaps_col=c(9),cutree_rows = 3,annotation_row=df_n4[3],
              colorRampPalette(c("black","white","firebrick3"))(56),annotation_col=annotation_col[2])
p_n
ggsave("p_n.pdf",p_n,width=9,height=7)

p_p<-pheatmap(df_p5,cluster_row=T,cluster_col=F,colorRampPalette(c("black","white","firebrick3"))(56),
              gaps_col=c(9),cutree_rows = 3,annotation_row=df_p4[3],
              annotation_col=annotation_col[2])
p_p
ggsave("p_p.pdf",p_p,width=9,height=7)

#associations between asv and plant atrributes
corr_meta_core<-function(df1,df2){
  corr_1<-rcorr(as.matrix(df1),as.matrix(df2),type = "spearman")
  #str(corr_1)
  corr_1r<-corr_1$r[1:9,-c(1:9)]
  corr_1p<-corr_1$P[1:9,-c(1:9)]
  corr_1r[abs(corr_1r)<0.4]<-0
  #View(corr_1r)
  corr_1p[corr_1p<0.05]<-1
  corr_1p[corr_1p<1]<-0
  #View(corr_1p)
  corr_1_rp<-corr_1r*corr_1p
  #View(corr_1_rp)
  return(corr_1_rp)
}
##dt_16s
corr_dt_16s_mp<-corr_meta_core(pcr[,1:9],pcr[,10:19])
write.csv(corr_dt_16s_mp,"corr_dt_16s_mp.csv")
View(corr_dt_16s_mp)
p1<-pheatmap(t(corr_dt_16s_mp),cluster_rows=F,
             gaps_col=c(6),cluster_cols=F)
p1
ggsave("p_corr_dt_16s_mp.pdf",p1,width=8.5,height=5)


##ggplot
p2<-ggplot(pcr,aes(x=BOTU4467,y=RL))+geom_point()+
    geom_smooth(method="lm",level=0.9,fill="green")+
    geom_smooth(method="lm",level=0.7,fill="blue")
p2

p3<-ggplot(pcr,aes(x=BOTU7109,y=RL))+geom_point()+
  geom_smooth(method="lm")
p3

p4<-ggplot(pcr,aes(x=BOTU2326,y=RL))+geom_point()+
  geom_smooth(method="lm")
p4

p5<-ggplot(pcr,aes(x=FOTU2059,y=RL))+geom_point()+
  geom_smooth(method="lm")
p5

p6<-ggplot(pcr,aes(x=BOTU7968,y=RL))+geom_point()+
  geom_smooth(method="lm")
p6

p7<-ggarrange(p2,p3,p4,p5,p6,ncol=3,nrow=2,
              labels=c("a","b","c","d","e"))
p7
ggsave("p7.pdf",p7,width=9.5,height=7.7)

##Other codes for this study can be found in the corresponding R packages##
