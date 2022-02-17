## Title: Data Analysis code for eDNA Wetmarket
## Author: Shelby E. McIlroy
## Date: 17 Feb 2022

#Load Libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(BiodiversityR)

##Load data
OTUmat<-read.csv(file.choose(),header=T) #use file "data/OTUmatrix.csv"
METAmat<-read.csv(file.choose(),header=T) #use file "data/METAmatrix.csv"
TAXAmat<-read.csv(file.choose(),header=T) #use filt "data/TAXAmatrix.csv"

##Format data for phyloseq
OTU_rows<-OTUmat$Sample
OTUmat1<-OTUmat[,-1]
rownames(OTUmat1)<-OTU_rows
OTUmat2<-as.matrix(OTUmat1)
OTU<-otu_table(OTUmat2,taxa_are_rows=T)

TAXA_rows<-TAXAmat$X
TAXAmat1<-TAXAmat[,-1]
rownames(TAXAmat1)<-TAXA_rows
TAX<-tax_table(as.matrix(TAXAmat1))

META_rows<-METAmat$X
METAmat1<-METAmat[,-1]
rownames(METAmat1)<-META_rows
METAmat1$Day<-as.character(METAmat1$Day)
META<-sample_data(METAmat1)

##Make a phyloseq object
ps<-phyloseq(OTU,META,TAX)
ps<-prune_taxa(taxa_sums(ps)>0,ps)

##Species Accumulation Curves
ps_spec<-t(as(otu_table(ps),"matrix"))
ps_spec<-as.data.frame(ps_spec)
ps_samp<-as(sample_data(ps),"data.frame")

Accum.1 <- accumcomp(ps_spec, y=ps_samp, factor='Method', 
                     method='exact', conditioned=FALSE, plotit=FALSE)
Accum.1

accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=5)
head(accum.long1)

##Plot Species Accumulation Curve
ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  geom_line(aes(colour=Grouping), size=1.25) +
  geom_point(data=subset(accum.long1, labelit==TRUE), 
             aes(colour=Grouping), size=0.5) +
  geom_ribbon(aes(fill=Grouping), alpha=0.2, show.legend=FALSE) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.155,0.84),legend.background = element_rect(fill='transparent'))+
  scale_colour_manual(name = element_blank(),
                      labels = c("eDNA Combined","Filtration", "Precipitation", "Visual"),
                      values = c("darkorchid4","#FF0000", "#00A08A", "#F2AD00"))  +
  scale_fill_manual(values = c("darkorchid4","#FF0000", "#00A08A", "#F2AD00"))+
  labs(x = "Samples", y = "Taxa")


##Principle Component Analysis

ps_single<-subset_samples(ps,Method!="eD") ##remove "combined eDNA" data
ps_f<-tax_glom(ps_single,"Family")
ps_f<-prune_taxa(taxa_sums(ps_f)>0,ps_f)

ps_f_num<-ps_f
x<-1:62
ps_f_num@tax_table@.Data[,5]<-x
Familynames<-ps_f@tax_table@.Data[,5]

j_dist_f<-phyloseq::distance(ps_f_num,method = "jaccard")
j_ord_f<-ordinate(ps_f_num,method="PCoA",distance=j_dist_f)
pf<-plot_ordination(ps_f_num,j_ord_f,type="biplot",color="Method",shape="Location")

pf + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_colour_manual(name = "Extraction Method",
                      labels = c("Taxa", "Filtration", "Precipitation", "Visual"),
                      values = c("white", "#FF0000", "#00A08A", "#F2AD00"))   +
  scale_shape_manual(name = "Location",
                     labels = c("Taxa", "Kennedy Town", "Shek Tong Tsui", "Sai Ying Pun"),
                     values = c(20, 15, 16, 17))+
  scale_size_manual(name = "Type",
                    labels = c("Data","Taxa"),
                    values = c(3,0)) +
  geom_text(aes(label = Family), size = 2.5,check_overlap = TRUE,color="black",fontface=2)  




