#          ######### Script Statistics ########


###        I Packages and data importation ####


library(ggplot2)
theme_set(theme_minimal()) 
library(data.table)
library(vegan)
library(devtools)
library(mixOmics)
library(FactoMineR)
library(pastecs)
source("http://peterhaschke.com/Code/multiplot.R")


# Euka :
euka=read.table("/home/faure/Documents/Stage/OK_2_finaux/alleuka.OK_2_calcizoo", sep="\t", header = FALSE)
names(euka)=c("PFT", "ID_barcode", "Taxogroup", "Station", "Abundance")
summary(euka)

# Proka :
proka=read.table("/home/faure/Documents/Stage/OK_2_finaux/allproca.OK_2", sep="\t", header = FALSE)
names(proka)=c("PFT", "ID", "Taxogroup", "Station", "Abundance")
summary(proka)

# PFT
# READ ME : The first two parts of this script lead to the creation of the tabpft table, if you
# possess only the files euka and proka, then run everything from part I, if you already possess
# the file tabpft_91clean.csv, then just launch the next line and go to part III.
tabpft=read.table("/home/faure/Documents/Stage/OK_2_finaux/tabpft_91clean.csv", header = TRUE)

# Environmental variables :
envi=read.table("/home/faure/Documents/Stage/OK_2_finaux/tabenvi.csv", header=TRUE)

# Geographical correspondance :
corres = read.table("/home/faure/Documents/Stage/OK_2_finaux/corres_geo.csv", header=TRUE)

###        II Preliminary work on data  #####


# We want to create contingency tables :

# Euka :

tabeuka <- xtabs(euka$Abundance ~ euka$PFT + euka$Station)
attr(tabeuka, "class") <- NULL
attr(tabeuka, "call") <- NULL
tabeuka=as.data.frame(tabeuka)
# We know some stations have more than four filters due to bugs, we get rid of bugged samples :
tabeuka=tabeuka[,-which(names(tabeuka)=="TARA_124_SRF_0.22-3")]
tabeuka=tabeuka[,-which(names(tabeuka)=="TARA_124_SRF_3-20")]
tabeuka=tabeuka[,-which(names(tabeuka)=="TARA_100_DCM_0.8-3")]
tabeuka=tabeuka[,-which(names(tabeuka)=="TARA_047_SRF_0.8-20")]

# We now want to keep only samples for which 4 filters have been used
for (i in c(1:length(names(tabeuka)))) {
  names(tabeuka)[i]=substr(names(tabeuka)[i],1,12)
}
uniqc=rle(names(tabeuka))
station_bug=uniqc$values[which(uniqc$lengths!=4)]
# We obtain the list of stations where less than 4 filters have been used :
station_bug
# We sum abundances on each samples :
tabeuka=sapply(unique(names(tabeuka)), 
               function(x) rowSums( tabeuka[ , grep(x, names(tabeuka)), drop=FALSE]) )
tabeuka=as.data.frame(tabeuka)
# We get rid of bugged stations :
tabeuka=tabeuka[,-which(names(tabeuka) %in% station_bug)]


# Proka : 

tabproka <- xtabs(proka$Abundance ~ proka$PFT + proka$Station)
attr(tabproka, "class") <- NULL
attr(tabproka, "call") <- NULL
tabproka=as.data.frame(tabproka)


# Unification of the two :

tabeuka.glob=tabeuka[,names(tabeuka) %in% names(tabproka)]
tabproka.glob=tabproka[,names(tabproka) %in% names(tabeuka)]
tabpft=rbind(tabeuka.glob,tabproka.glob)
rownames(tabpft)
# We rename Calcifiers to phyto calcifiers for more precision, and we pool all pico nano autotrophs
rownames(tabpft)[1]="Phyto_Calcifiers"
tabpft[22,]=tabpft[12,]+tabpft[21,]
tabpft=tabpft[-c(12,21),]
rownames(tabpft)[20]="Pico_nano_phyto"
head(tabpft)
# We want now to have stations in lines
tabpft=as.data.frame(t(tabpft))
# And we want to add two columns, one for the depth and one for the station
tabpft[,21]=substr(rownames(tabpft),10,12)
tabpft[,22]=substr(rownames(tabpft),6,8)
names(tabpft)[21:22]=c("Depth","Station")
tabpft[,21]=as.factor(tabpft[,21])
tabpft[,22]=as.factor(tabpft[,22])
summary(tabpft)

# We can now save the file created :
#write.table(tabpft, file="/home/faure/Documents/Stage/OK_2_finaux/tabpft_91clean.csv", sep="\t", row.names=TRUE, col.names=TRUE)



###        III Exploration graphs  #####

#       Eukaryotes histogram of barcode quantity by PFT

eukahist <- within(euka, PFT <- factor(PFT,levels=c("Meso_macro_zoo","Proto_zoo","Pico_nano_zoo","Meso_macro_mixo",
                                                    "Proto_mixo","Pico_nano_mixo","Phyto_silicifier", "Zoo_silicifier",
                                                    "Para_gen", "Para_meta", "Para_algal", "Pico_nano_phyto", "Mixed_phyto", 
                                                    "DMS_prod", "Calcifiers","Zoo_Calcifiers",  "others_noPFT")))
eukahist[,6]=rep(0,nrow(eukahist))
colnames(eukahist)[6]="Groupe"
index=levels(eukahist$PFT)
values=c("Zoo", "Zoo", "Zoo", "Mixo", "Mixo","Mixo","Silicificateurs", "Silicificateurs","Para", "Para", "Para", "Phyto",  "Phyto","DMS", "Calcificateurs","Calcificateurs", "Others")
eukahist$Groupe = as.factor(values[match(eukahist$PFT,index)])

eukahist <- within(eukahist, Groupe <- factor(Groupe, levels=names(sort(table(Groupe),  decreasing=TRUE))))
g4 = ggplot(data = eukahist, aes(x=Groupe, fill=PFT)) + geom_bar(stat="count") +
  labs(title="Nombre de barcodes par PFT",subtitle="Chez les eucaryotes") +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5)) +
  scale_fill_manual(values=c("steelblue2", "steelblue", "steelblue4", "yellow2", "gold", "yellow3", "green4", 
                             "darkorange3", "violetred", "violetred1","violetred4","springgreen",
                             "springgreen3", "yellowgreen", "firebrick1", "firebrick3",  "grey74"))
g4

#       Ekaryotes barplot of abundance

eukahist = xtabs(Abundance~PFT+Groupe, data=eukahist)
attr(eukahist, "class") <- NULL
attr(eukahist, "call") <- NULL
eukahist=as.data.frame(eukahist)
head(eukahist)
eukahist[,9]=c("Zoo", "Zoo", "Zoo", "Mixo", "Mixo","Mixo","Silicificateurs", "Silicificateurs","Para", "Para", "Para", "Phyto",  "Phyto","DMS", "Calcificateurs", "Calcificateurs", "Others")
for (i in 1:17) {
  eukahist[i,10]=sum(eukahist[i,1:8])
}
eukahist=eukahist[,-c(1:8)]
eukahist[,3]=as.factor(rownames(eukahist))
names(eukahist)=c("Groupe","Abundance", "PFT")
eukahist = within(eukahist, Groupe <- factor(Groupe,levels=c("Zoo", "Mixo", "Silicificateurs", 
                                                             "Para", "Phyto",  "DMS", "Calcificateurs",
                                                             "Others")))
eukahist <- within(eukahist, PFT <- factor(PFT,levels=c("Meso_macro_zoo","Proto_zoo","Pico_nano_zoo","Meso_macro_mixo",
                                                        "Proto_mixo","Pico_nano_mixo","Phyto_silicifier", "Zoo_silicifier",
                                                        "Para_gen", "Para_meta", "Para_algal", "Pico_nano_phyto", "Mixed_phyto", 
                                                        "DMS_prod", "Calcifiers","Zoo_Calcifiers","others_noPFT")))
g5 = ggplot(data = eukahist, aes(x=Groupe, y=Abundance, fill=PFT)) + geom_bar(stat="identity") +
  labs(title="Abondance des barcodes par PFT",subtitle="Chez les eucaryotes") +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5)) +
  scale_fill_manual(values=c("steelblue2", "steelblue", "steelblue4", "yellow2", "gold", "yellow3", "green4", 
                             "darkorange3", "violetred", "violetred1","violetred4","springgreen",
                             "springgreen3", "yellowgreen", "firebrick1", "firebrick3", "grey74"))
g5


#      Prokaryotes histogram of barcode Abundance by PFT

g4 = ggplot(proka, aes(x=PFT, y=log(Abundance), fill=PFT)) + geom_bar(stat="identity") +
  scale_fill_manual(values=c("beige", "bisque2", "Royalblue4", "springgreen")) +
  labs(title="Abondance par PFT",subtitle="Chez les procaryotes") +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5)) 
g4

##              Plots on the 91 Samples

#We need to create a new table with summed abundance for each PFT, regardless of stations
summary(tabpft)
tabpfthist=rep(0,20)
for (i in 1:20) {
  tabpfthist[i]=sum(tabpft[,i])
  names(tabpfthist)[i]=names(tabpft)[i]
}
tabpfthist=as.data.frame(tabpfthist)
names(tabpfthist)="Abundance"
# Now we want to add PFTs and Groupes
tabpfthist[,2]=rownames(tabpfthist)
tabpfthist[,3]=c("Calcificateurs", "DMS","Mixo","Zoo","Phyto","Others","Para","Para","Para",
                 "Silicificateurs", "Mixo", "Zoo", "Mixo","Zoo","Calcificateurs", "Silicificateurs",
                 "Proca","Proca","Proca","Phyto")
names(tabpfthist)[2:3]=c("PFT","Groupe")

tabpfthist = within(tabpfthist, Groupe <- factor(Groupe,levels=c("Zoo", "Mixo", "Silicificateurs", 
                                                             "Para", "Phyto",  "DMS", "Calcificateurs","Proca",
                                                             "Others")))
tabpfthist <- within(tabpfthist, PFT <- factor(PFT,levels=c("Meso_macro_zoo","Proto_zoo","Pico_nano_zoo","Meso_macro_mixo",
                                                        "Proto_mixo","Pico_nano_mixo","Phyto_silicifier", "Zoo_silicifier",
                                                        "Para_gen", "Para_meta", "Para_algal", "Pico_nano_phyto", "Mixed_phyto", 
                                                        "DMS_prod", "Phyto_Calcifiers","Zoo_Calcifiers","N2_auto","N2_hetero",
                                                        "Pico_nano_hetero","others_noPFT")))

g4 = ggplot(tabpfthist, aes(x=Groupe, y=Abundance, fill=PFT)) + geom_bar(stat="identity") +
  labs(title="Abondance des barcodes par PFT",subtitle="Sur nos 91 couples station/profondeur finaux") +
  scale_fill_manual(values=c("steelblue2", "steelblue", "steelblue4", "yellow2", "gold", "yellow3", "green4", 
                             "darkorange3", "violetred", "violetred1","violetred4","springgreen",
                             "springgreen3", "yellowgreen", "firebrick1", "firebrick3", "beige", "bisque2", 
                             "Royalblue4","grey74")) +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5)) 
g4

# For MAPS see corresponding script.

###        IV PCA on PFTs  #####

# First we add geographical info to use it in our analyses :
corres2=corres[corres$StationID %in% rownames(tabpft),]
tabpft=merge(x=tabpft, y=corres2, by.x="row.names", by.y="StationID")
head(tabpft)
rownames(tabpft)=tabpft[,1]
tabpft=tabpft[,-1]
# Number of samples in each Ocean region :
ggplot(data=tabpft, aes(x=oceanSeaRegion)) + geom_bar(stat="count")

# We have three stations where depth is indicate as MIX instead of DCM
tabpft[which(tabpft$Depth=="MIX"),21]=c("DCM","DCM","DCM")

# PCA :
# We use the hellinger transformation :
tabpft.hel=decostand(tabpft[,1:20], method="hellinger")
acp=rda(tabpft.hel)
summary(acp) # 52% on the first two axis 

ev=as.data.frame(acp$CA$eig[1:9])
ggplot(data=ev, aes(x=factor(rownames(ev)), y=ev[,1])) + geom_bar(stat="identity") + 
  geom_hline(yintercept=mean(ev[,1]), linetype = 'dashed', color="red") +
  labs(y="Eigen Value", x="Axes de l'analyse de correspondance")
# The four first axis are significant

# We need individuals coordinates
coord.acp.sta = as.data.frame(acp$CA$u[,1:4])
# and variable coordinates
coord.acp.pft = as.data.frame(acp$CA$v[,1:4])
# to build the biplot :
# first, we reorganize factors in order to have a more interesting legend :
tabpft <- within(tabpft,oceanSeaRegion <- factor(oceanSeaRegion,levels=c("NAO","SAO","MS","RS","IO","NPO","SPO","SO")))
# then we build the graph :
g1 = ggplot() + geom_point(data=coord.acp.sta, aes(x=PC1, y=PC2, col=tabpft[,23], shape=tabpft$Depth)) + 
  geom_vline(xintercept = 0, linetype='dotted') +
  geom_hline(yintercept = 0, linetype='dotted') +
  labs(x="PC1 (30.68%)", y="PC2 (22.13%)", title="PCA individuals plot", subtitle = "Post correction filtrages") +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5))
g2 = ggplot() + geom_segment(data=coord.acp.pft, aes(xend=PC1, yend=PC2),x=0,y=0, color="red",  arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(data=coord.acp.pft, aes(label=rownames(coord.acp.pft), x=PC1, y=PC2), size=3) +
  labs(x="PC1 (30.68%)", y="PC2 (22.13%)", title="PCA species plot", subtitle = "Post correction filtrages") +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5)) +
  geom_vline(xintercept = 0, linetype='dotted') +
  geom_hline(yintercept = 0, linetype='dotted') +
  xlim(-0.9,0.9) +
  ylim(-0.9,0.9)
g3 = ggplot() + geom_point(data=coord.acp.sta, aes(x=PC3, y=PC4, col=tabpft[,23], shape=tabpft$Depth)) + 
  geom_vline(xintercept = 0, linetype='dotted') +
  geom_hline(yintercept = 0, linetype='dotted') +
  labs(x="PC3 (14.25%)", y="PC4 (10.69%)", title="PCA individuals plot", subtitle = "Post correction filtrages") +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5))
g4 = ggplot() + geom_segment(data=coord.acp.pft, aes(xend=PC3, yend=PC4),x=0,y=0, color="red",  arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(data=coord.acp.pft, aes(label=rownames(coord.acp.pft), x=PC3, y=PC4), size=3) +
  labs(x="PC3 (14.25%)", y="PC4 (10.69%)", title="PCA species plot", subtitle = "Post correction filtrages") +
  theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5)) +
  geom_vline(xintercept = 0, linetype='dotted') +
  geom_hline(yintercept = 0, linetype='dotted') +
  xlim(-0.9,0.9) +
  ylim(-0.9,0.9)

multiplot(g1,g3,g2,g4, cols=2)
# Some work is then achieved on inkscape in order to obtain the final figure.

###        IV PCA on Environmental variables  #####

# Once again we want to add geographical locations inorder to use them in our analysis
envitot=merge(envi,corres, by.x="row.names", by.y=2)
rownames(envitot)=envitot[,1]
envitot=envitot[,-1]
summary(envitot)
# We get rid of Arctique Ocean samples which aren't in our samples
envi=envitot[envitot$oceanSeaRegion!="AO",]

# We want to transform the dataset so that it can be used for a PCA.
# We begin by getting rid of constant and only NA variables :
envi=envi[,-c(20,23:24)]
# Then we replace NA values by the mean of the corresponding column :
for (i in 1:ncol(envi[,-51])) {
  envi[is.na(envi[,i]),i]=mean(envi[,i], na.rm=TRUE)
}
rm(i)
summary(envi) # No more NAs

# The data set is huge, we want to simplify it by making a selection of variables, 
# we use Escoufier vectors:
envi.sel=escouf(envi[,-51])
plot(envi.sel)
# We choose a 0.9 treshold
envi=extract(envi.sel,level=0.9)
# We add lat, lon, and ocean Sea Region
envi[,length(names(envi))+1]=envitot[envitot$oceanSeaRegion!="AO",]$latitude
names(envi)[length(names(envi))]="latitude"
envi[,length(names(envi))+1]=envitot[envitot$oceanSeaRegion!="AO",]$longitude
names(envi)[length(names(envi))]="longitude"
envi[,length(names(envi))+1]=envitot[envitot$oceanSeaRegion!="AO",]$oceanSeaRegion
names(envi)[length(names(envi))]="oceanSeaRegion"
summary(envi)

# Now we reorganize factors like for PFTs
envi <- within(envi,oceanSeaRegion <- factor(oceanSeaRegion,levels=c("NAO","SAO","MS","RS","IO","NPO","SPO","SO")))

acp.sel=PCA(envi, quali.sup=length(names(envi)), quanti.sup=c(length(names(envi))-1,length(names(envi))-2))
summary(acp.sel)
# Not too bad, 50% on the two first axis, >70% on the first four
g1 = ggplot() + geom_segment(data=as.data.frame(acp.sel$var$coord[,1:4]), aes(x=0,y=0,xend=Dim.1,yend=Dim.2), color="red") +
  geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x="Dim 1 (28.010%)", y="Dim 2 (22.758%)") +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_text(data=as.data.frame(acp.sel$var$coord[,1:4]), aes(label=rownames(as.data.frame(acp.sel$var$coord[,1:4])),x=Dim.1,y=Dim.2)) +
  theme_minimal()
g2 = ggplot() + geom_point(data=as.data.frame(acp.sel$ind$coord[,1:4]), aes(x=Dim.1,y=Dim.2, col=envi$oceanSeaRegion)) +
  scale_color_discrete(name="Zone océanique") +
  labs(x="Dim 1 (28.010%)", y="Dim 2 (22.758%)") +
  geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  geom_text(data=as.data.frame(acp.sel$quali.sup$coord), aes(x=Dim.1,y=Dim.2, label=rownames(as.data.frame(acp.sel$quali.sup$coord))), size=2.1) +
  theme_minimal()
g3 = ggplot() + geom_segment(data=as.data.frame(acp.sel$var$coord[,1:4]), aes(x=0,y=0,xend=Dim.3,yend=Dim.4), color="red") +
  geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x="Dim 3 (12.015%)", y="Dim 4 (8.000%)") +
  xlim(-1,1) +
  ylim(-1,1) +
  geom_text(data=as.data.frame(acp.sel$var$coord[,1:4]), aes(label=rownames(as.data.frame(acp.sel$var$coord[,1:4])),x=Dim.3,y=Dim.4)) +
  theme_minimal()
g4 = ggplot() + geom_point(data=as.data.frame(acp.sel$ind$coord[,1:4]), aes(x=Dim.3,y=Dim.4, col=envi$oceanSeaRegion)) +
  scale_color_discrete(name="Zone océanique") +
  labs(x="Dim 3 (12.015%)", y="Dim 4 (8.000%)") +
  geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  geom_text(data=as.data.frame(acp.sel$quali.sup$coord), aes(x=Dim.3,y=Dim.4, label=rownames(as.data.frame(acp.sel$quali.sup$coord))), size=2.1) +
  theme_minimal()

multiplot(g1,g3,g2,g4, cols=2)
# Some work is then achieved on inkscape in order to obtain the final figure.

###        V RDA  #####

# We want to conduct our RDA on the complete environmental dataset
envi=read.table("/home/faure/Documents/Stage/OK_2_finaux/tabenvi.csv", header=TRUE)
envi=envi[rownames(envi) %in% rownames(tabpft),]
envi=merge(corres,envi,by.x=2,by.y="row.names")
rownames(envi)=envi[,1]
envi=envi[,-1]
summary(envi) 
# PAR in situ only has NAs so we get rid of it
envi=envi[,-which(names(envi)=="PAR.in.situ")]
# In order to be able to do the analysis, we have to replace NAs with mean values of  the variable. 
for (i in 1:ncol(envi)) {
  envi[is.na(envi[,i]),i]=mean(envi[,i], na.rm=TRUE)
}
rm(i)
summary(envi) 
# We have constant values for Sea Ice variables, because no more Arctic stations are in the dataset
# So we get rid of them.
envi=envi[,-c(23,24)]

rda.pft=rda(tabpft.hel~.,envi, na.action = na.exclude)
summary(rda.pft)
RsquareAdj(rda.pft)
# Adj Rsq = 36.8%, not too bad
anova(rda.pft, step=1000, by="axis") # 10 axes are significant
# We use the ordistep function to select variables through permutation tests :
set.seed(1500)
rda.step.forward = ordistep(rda(tabpft.hel~1,data=envi, na.action=na.exclude), scope=formula(rda.pft), direction="both", pstep=10000)
formula(rda.step.forward)

rda.pft=rda(formula=formula(rda.step.forward), data=envi, na.action = na.exclude)
summary(rda.pft)
RsquareAdj(rda.pft)
#Adj Rsq is now down to 30.1%
anova(rda.pft, step=1000, by="axis") # 5 axes are significant

# Triplot :

# Sites with site names
station <- scores(rda.pft, scaling=2)$sites
g1 <- ggplot() + geom_point(aes(station[,1], station[,2], shape=tabpft[,21]))+
  xlab("RDA1 (16.7%)") +  ylab("RDA2 (14.9%)")

# Adding red segments for the species scores and add their names
species <- scores(rda.pft, scaling=2)$species
g2 <- g1  + geom_segment(aes(xend = species[,1], yend = species[,2]),x=0,y=0,size = 0.5, color = 'red') +
  geom_text(aes(species[,1]*1.1, species[,2]*1.1, label=rownames(species)), color="red")

# Adding centroids for the Oceanic zone factor
zone <- as.data.frame(rda.pft$CCA$centroids[,1:2])
rownames(zone)=substr(rownames(zone),15,17)
zone[,3]=rownames(zone)
zone = within(zone,V3 <- factor(V3,levels=c("NAO","SAO","MS","RS","IO","NPO","SPO","SO")))
g3 <- g2 + geom_point(data=as.data.frame(zone), aes(x=RDA1,y=RDA2, col=zone$V3), shape=19, size=3) +
  scale_color_discrete(name="Zone Océanique")

g3 <- g3 +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(title="Redundancy Analysis") +
  theme(plot.title=element_text(hjust=0.5))
g3

# Adding envi variables 
envi.scores <- scores(rda.pft, choices = 1:2, display="bp", scaling=2)[-c(3:9),] # On retire les zones océnaiques et on ne conserve que les deux premiers axes
g4 <- g3 + geom_segment(data=as.data.frame(envi.scores),aes(xend = RDA1*1.5, yend = RDA2*1.5),x=0,y=0,size = 0.5, linetype="F1",color = 'blue',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(envi.scores[,1]*1.7, envi.scores[,2]*1.7, label=rownames(envi.scores)), color="blue")
g4










