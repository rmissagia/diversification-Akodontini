#Article: Missagia et al. Decoupled patterns of diversity and disparity characterize an ecologically specialized lineage of neotropical cricetids 
#Authors: Rafaela V. Missagia1, Daniel M. Casali, Bruce D. Patterson, Fernando A. Perini
#Description: R code to reproduce analyses performed in the study.

#read packages
library(geomorph)
library(phytools)
library(geiger)
library(evobiR)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(scales)
library(PhylogeneticEM)
library(RevGadgets)
library(openxlsx)
library(tidyr)
library(reshape)
library(BAMMtools)
library(coda)

#read_tree - Chronogram from Maestri et al. 2017, with 226 tips dropped. Ultrametric rounding correction performed with phytool::force.ultrametric
tr <- read.tree("AKO_TREE.nwk")

###################################################################################################################################################
########################################################GEOMETRIC_MORPHOMETRICS####################################################################
###################################################################################################################################################

#loading TPS files, defining speciesID
ventral_prov <- readland.tps(file = "ventral.tps", specID = "imageID", negNA = TRUE)
ventral <- estimate.missing(ventral_prov, method = "TPS")
lateral_prov <- readland.tps(file = "lateral.tps", specID = "imageID", negNA = TRUE)
lateral <- estimate.missing(lateral_prov, method = "TPS")
mandible_prov <- readland.tps(file = "mandible.tps", specID = "imageID", negNA = TRUE)
mandible <- estimate.missing(mandible_prov, method = "TPS")
landmarks<-list(ventral,lateral,mandible)

#lists defining groups (spp)
ventral_list <- read.csv("ventral_list.csv", header = T, sep = ";")
lateral_list <- read.csv("lateral_list.csv", header = T, sep = ";")
mandible_list <- read.csv("mandible_list.csv", header = T, sep = ";")
groups<-list(ventral_list,lateral_list,mandible_list)

#GPA
GPA<-vector(mode="list",length=length(landmarks))
lands2d<-vector(mode="list",length=length(landmarks))
GPA2<-vector(mode="list",length=length(landmarks))
for (i in 1:length(groups))
{
	GPA[[i]] <- gpagen(landmarks[[i]])
	lands2d[[i]] <- two.d.array(GPA[[i]]$coords)
	GPA2[[i]] <- arrayspecs(lands2d[[i]], length(landmarks[[i]][,,1][,1]), 2)
}

#symmetric component - ventral view only
ind<-as.factor(groups[[1]][,2])#list of individuals for symmetric component extraction
pairs<-as.matrix(read.csv("ventral_lmpairs.csv", header=T, sep = ","))#pairs of landmarks on each side of the skull
gdf <- geomorph.data.frame(shape = GPA2[[1]], ind = ind)
sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, land.pairs = pairs, data = gdf, RRPP = TRUE, iter = 149)
sym2 <- two.d.array(sym$symm.shape)
sym2

#logsize
logCS_species_avg<-vector(mode="list",length=length(groups))
logCS_species_sd<-vector(mode="list",length=length(groups))
for (i in 1:length(groups))
{
	spp<-as.factor(groups[[i]][,1])#list of specific groupings
	Csize2d <- log(GPA[[i]]$Csize)
	logCS_species_avg[[i]] <- aggregate(Csize2d, by = list(spp), FUN = mean)
	colnames(logCS_species_avg[[i]])<-c("Taxon","Avg.")
	logCS_species_sd[[i]] <- aggregate(Csize2d, by = list(spp), FUN = sd)
	colnames(logCS_species_sd[[i]])<-c("Taxon","SD")
}
logCS_species<-list(logCS_species_avg,logCS_species_sd)
names(logCS_species)<-c("logCS_avg","logCS_sd")
saveRDS(logCS_species,"logCS_species.rds")

#taking the average value of landmarks by species
colnames(sym2)<-colnames(lands2d[[1]])
landmarks2<-list(sym2,lands2d[[2]],lands2d[[3]])
sppmean<-vector(mode="list",length=length(landmarks2))
sppmean3d<-vector(mode="list",length=length(landmarks2))
for (i in 1:length(landmarks2))
{	
	spp<-as.factor(groups[[i]][,1])#list of specific groupings
	sppmean[[i]] <-aggregate(landmarks2[[i]], by = list(spp), FUN=mean)[,2:(length(landmarks2[[i]][1,])+1)]
	rownames(sppmean[[i]])<-levels(spp)
	sppmean3d[[i]] <-arrayspecs(sppmean[[i]], length(landmarks2[[i]][1,])/2, 2, sep=".")
}

#PCA (species)
pca<-vector(mode="list",length=length(sppmean3d))
pcs<-vector(mode="list",length=length(sppmean3d))
for (i in 1:length(sppmean3d))
{
	pca[[i]] <- gm.prcomp(sppmean3d[[i]],tol=0)
	pcs[[i]] <- pca[[i]]$x
}
pcs
saveRDS(pcs,"PCs_species.rds")

#loading outlines for plotting shape changes on each component
lands_ventral<-readland.tps(file = "ventral_out.TPS")
lands_lateral<-readland.tps(file = "lateral_out.TPS")
lands_mandible<-readland.tps(file = "mandible_out.TPS")	
lands<-list(lands_ventral,lands_lateral,lands_mandible)
out_ventral<-warpRefOutline(file = "ventral_out.csv", lands[[1]][,,1], mshape(sppmean3d[[1]]))
out_lateral<-warpRefOutline(file = "lateral_out.csv", lands[[2]][,,1], mshape(sppmean3d[[2]]))
out_mandible<-warpRefOutline(file = "mandible_out.csv", lands[[3]][,,1], mshape(sppmean3d[[3]]))
out<-list(out_ventral,out_lateral,out_mandible)

##plotting shape changes by outlines
titles<-c("Ventral","Lateral","Mandible")
pdf("Ako_outlines.pdf")
for (i in 1:length(sppmean3d))
{
	ref <- mshape(sppmean3d[[i]])	
	
	plotRefToTarget(ref,pca[[i]]$shapes$shapes.comp1$min,outline=out[[i]]$outline)
	title(titles[[i]],"pc1 min")
	plotRefToTarget(ref,pca[[i]]$shapes$shapes.comp1$max,outline=out[[i]]$outline)
	title(titles[[i]],"pc1 max")	
	plotRefToTarget(ref,pca[[i]]$shapes$shapes.comp2$min,outline=out[[i]]$outline)
	title(titles[[i]],"pc2 min")
	plotRefToTarget(ref,pca[[i]]$shapes$shapes.comp2$max,outline=out[[i]]$outline)
	title(titles[[i]],"pc2 max")
}
dev.off()

###################################################################################################################################################
###########################################################PHYLO SIGNAL############################################################################
###################################################################################################################################################

##Size correlations and getting the average size
logCS_species<-readRDS("logCS_species.rds")
logCS_species_avg<-logCS_species$logCS_avg
VL<-cor.test(logCS_species_avg[[1]]$Avg.,logCS_species_avg[[2]]$Avg.)
VM<-cor.test(logCS_species_avg[[1]]$Avg.,logCS_species_avg[[3]]$Avg.)
LM<-cor.test(logCS_species_avg[[2]]$Avg.,logCS_species_avg[[3]]$Avg.)
size_cor<-list(VL,VM,LM)
cor_table<-as.data.frame(matrix(ncol=4,nrow=length(size_cor)))
rownames(cor_table)<-c("Skull - ventral x Skull - lateral","Skull - ventral x Mandible","Skull - lateral x Mandible")
colnames(cor_table)<-c("test statistic","df","p-value","correlation")
for (i in 1:length(size_cor))
{
	cor_table[i,1]<-size_cor[[i]]$statistic
	cor_table[i,2]<-size_cor[[i]]$parameter
	cor_table[i,3]<-size_cor[[i]]$p.value
	cor_table[i,4]<-size_cor[[i]]$estimate
}
cor_table
write.csv(cor_table,"size_correlations_table.csv")

size_avg<-as.data.frame(matrix(ncol=2,nrow=length(tr$tip.label)))
colnames(size_avg)<-c("Taxon","Size")
for (i in 1:length(tr$tip.label))
{
	size_avg[i,1]<-levels(logCS_species_avg[[1]]$Taxon)[i]
	size_avg[i,2]<-mean(logCS_species_avg[[1]]$Avg.[i],logCS_species_avg[[2]]$Avg.[i],logCS_species_avg[[3]]$Avg.[i])
}
size_avg
size_matrix<-as.matrix(size_avg$Size)
rownames(size_matrix)<-size_avg$Taxon
colnames(size_matrix)<-"logCS"
size_matrix
write.csv(size_matrix,"avg_size_table.csv")

to_sig<-list(sppmean3d[[1]],sppmean3d[[2]],sppmean3d[[3]],size_matrix)
names(to_sig)<-c("Ventral","Lateral","Mandible","Size")
kmult_results<-vector(mode="list",length=length(to_sig))
kmult_table<-as.data.frame(matrix(ncol=3,nrow=length(to_sig)))
rownames(kmult_table)<-c("Skull - ventral","Skull - lateral","Mandible","Size")
colnames(kmult_table)<-c("Phylogenetic signal","p-value","Effect size")
for (i in 1:length(to_sig))
{
	kmult_results[[i]]<-physignal(to_sig[[i]], tr, iter = 999, seed = NULL, print.progress = T)
	kmult_table[i,1]<-kmult_results[[i]]$phy.signal
	kmult_table[i,2]<-kmult_results[[i]]$pvalue
	kmult_table[i,3]<-kmult_results[[i]]$Z
}
kmult_table
write.csv(kmult_table,"kmult_table.csv")

###################################################################################################################################################
##############################################################DISPARITY############################################################################
###################################################################################################################################################

##Allometry

allometry<-vector(mode="list",length=length(sppmean3d))
for (i in 1:length(sppmean3d))
{	
	pgls_data <-geomorph.data.frame(SHAPE = sppmean3d[[i]], SIZE = size_matrix[,1], phy = tr)
	allometry[[i]] <- procD.pgls(SHAPE~SIZE, data = pgls_data, phy, iter = 999, print.progress = TRUE)
}
summary(allometry[[1]])
summary(allometry[[2]])
summary(allometry[[3]])

wb<-createWorkbook()
addWorksheet(wb, "Ventral")
addWorksheet(wb, "Lateral")
addWorksheet(wb, "Mandible")
writeData(wb, 1, allometry[[1]]$aov.table)
writeData(wb, 2, allometry[[2]]$aov.table)
writeData(wb, 3, allometry[[3]]$aov.table)
saveWorkbook(wb, "allometry_aov_tables.xlsx", overwrite=TRUE)

##Clade disparities

cladeNames<-data.frame(spp=rownames(sppmean[[1]]),clade=c(rep("Clade B",30),rep("Clade A",4),rep("Clade B",3),rep("Clade A",2),rep("Clade B",15),rep("Clade A",2),rep("Clade B",3)))
cladeDisp<-vector(mode="list",length=length(allometry))
for (i in 1:length(allometry))
{		
	cladeDisp[[i]] <- morphol.disparity(f=allometry[[i]], data=allometry[[i]], transform=F, groups=cladeNames$clade)
}
summary(cladeDisp[[1]])
summary(cladeDisp[[2]])
summary(cladeDisp[[3]])

##Size difference
size_diff<-morphol.disparity(size_matrix[,1]~1,groups=cladeNames$clade)
summary(size_diff)

clade_comp<-list(cladeDisp[[1]],cladeDisp[[2]],cladeDisp[[3]],size_diff)
clade_disp_table<-as.data.frame(matrix(ncol=4,nrow=length(clade_comp)))
rownames(clade_disp_table)<-c("Skull - ventral","Skull - lateral","Mandible","Size")
colnames(clade_disp_table)<-c("Procrustes variance Clade A","Procrustes variance Clade B","Pairwise absolute distance","p-value")
for (i in 1:length(clade_comp))
{
	clade_disp_table[i,1]<-clade_comp[[i]]$Procrustes.var[[1]]
	clade_disp_table[i,2]<-clade_comp[[i]]$Procrustes.var[[2]]
	clade_disp_table[i,3]<-clade_comp[[i]]$PV.dist[[2]]
	clade_disp_table[i,4]<-clade_comp[[i]]$PV.dist.Pval[[2]]
}
clade_disp_table	
write.csv(clade_disp_table,"clade_disp_table.csv")

##Partial disparities

PD<-vector(mode="list",length=length(landmarks2))
for (i in 1:length(landmarks2))
{
	spp <- as.factor(groups[[i]][,1])#list of specific groupings
	gdf <- geomorph.data.frame(landmarks2[[i]], species = spp)
	partialDisp <- morphol.disparity(landmarks2[[i]]~1, group = spp, data=gdf, partial = T)
	PD[[i]] <- partialDisp$Procrustes.var
}
names(PD)<-titles
PD
write.csv(PD,"part_disp_table.csv")

##Node disparities
#get the first 8 PCs (90% of the variance) and size data

disp_node_data<-list(pca[[1]]$x[,1:8],pca[[2]]$x[,1:8],pca[[3]]$x[,1:8],size_matrix)
disp_node_results<-vector(mode="list",length=length(disp_node_data))
for (i in 1:length(disp_node_data))
{
	disp_node_results[[i]]<-disparity(tr, disp_node_data[[i]])
}
disp_node_results

#table
init<-length(tr$tip.label)+1
end<-length(tr$tip.label)+tr$Nnode
disp_node_table<-data.frame(disp_node_results[[1]],disp_node_results[[2]],disp_node_results[[3]],disp_node_results[[4]])
colnames(disp_node_table)<-c("Ventral","Lateral","Mandible","Size")
rownames(disp_node_table)<-c(init:end)
disp_node_table
write.csv(disp_node_table,"node_disp_table.csv")

#plot
tr <- read.tree("AKO_TREE.nwk")#read tree
par(mfrow=c(1,4))
labels<-as.data.frame(tr$tip.label)
colnames(labels)<-"taxa"
sep<-separate(labels,taxa,c("G","S"),sep="_")
unt<-unite(sep,taxa,c("G","S"),sep=" ")
tr$tip.label<-unt$taxa
color<-c(rep("green",51),rep("#E31C79",8),"black",rep("green",50),rep("#E31C79",7))
d<- data.frame(node=c(1:end),color=color)
title_list<-c("A) Skull (ventral view)","B) Skull (lateral view)","C) Mandible","D) Size")

data1<-data.frame(c(rep(0,length(tr$tip.label)),disp_node_table[,1]))
colnames(data1)<-c("disp")
data2<-data.frame(c(rep(0,length(tr$tip.label)),disp_node_table[,2]))
colnames(data2)<-c("disp")
data3<-data.frame(c(rep(0,length(tr$tip.label)),disp_node_table[,3]))
colnames(data3)<-c("disp")
data4<-data.frame(c(rep(0,length(tr$tip.label)),disp_node_table[,4]))
colnames(data4)<-c("disp")

fig1<-revts(ggtree(tr, size=0.8, layout = "rectangular", ladderize=TRUE, right=TRUE)%<+% d + aes(color=I(color))) + geom_nodepoint(aes(size=data1$disp)) + geom_tiplab(fontface="italic", size=2, color="black")+
xlim_tree(3) + scale_x_continuous(breaks=seq(-8,0,1), labels=abs(seq(-8,0,1))) + theme_tree2()+ labs (size="Disparity", title=title_list[1])+
theme(legend.title=element_text(size=10), legend.position="bottom", legend.text=element_text(size=8),axis.text.x=element_text(size=8))
fig2<-revts(ggtree(tr, size=0.8, layout = "rectangular", ladderize=TRUE, right=TRUE)%<+% d + aes(color=I(color))) + geom_nodepoint(aes(size=data2$disp)) + geom_tiplab(fontface="italic", size=2, color="black")+
xlim_tree(3) + scale_x_continuous(breaks=seq(-8,0,1), labels=abs(seq(-8,0,1))) + theme_tree2()+ labs (size="", title=title_list[2])+
theme(legend.title=element_text(size=10), legend.position="bottom", legend.text=element_text(size=8),axis.text.x=element_text(size=8))
fig3<-revts(ggtree(tr, size=0.8, layout = "rectangular", ladderize=TRUE, right=TRUE)%<+% d + aes(color=I(color))) + geom_nodepoint(aes(size=data3$disp)) + geom_tiplab(fontface="italic", size=2, color="black")+
xlim_tree(3) + scale_x_continuous(breaks=seq(-8,0,1), labels=abs(seq(-8,0,1))) + theme_tree2()+ labs (size="", title=title_list[3])+
theme(legend.title=element_text(size=10), legend.position="bottom", legend.text=element_text(size=8),axis.text.x=element_text(size=8))
fig4<-revts(ggtree(tr, size=0.8, layout = "rectangular", ladderize=TRUE, right=TRUE)%<+% d + aes(color=I(color))) + geom_nodepoint(aes(size=data4$disp)) + geom_tiplab(fontface="italic", size=2, color="black")+
xlim_tree(3) + scale_x_continuous(breaks=seq(-8,0,1), labels=abs(seq(-8,0,1))) + theme_tree2()+ labs (size="", title=title_list[4])+
theme(legend.title=element_text(size=10), legend.position="bottom", legend.text=element_text(size=8),axis.text.x=element_text(size=8))

fig<-ggarrange(fig1,fig2,fig3,fig4, ncol = 4, nrow = 1, common.legend = FALSE, legend="bottom")
fig
ggsave("Node_disparities.svg",fig,width=14,height=6)

####################################################################################################################################################
##############################################################MORPHO_RATES##########################################################################
####################################################################################################################################################

###Phylogenetic EM - obtaining reasonable values for the number of rate shifts to be added as a pior in RevBayes script
pcs<-readRDS("PCs_species.rds")
tr <- read.tree("AKO_TREE.nwk")#read tree
Skull_ventral<-t(ReorderData(tr, pcs[[1]][,1:8]))
Skull_lateral<-t(ReorderData(tr, pcs[[2]][,1:8]))
Mandible<-t(ReorderData(tr, pcs[[3]][,1:8]))
size<-read.csv("avg_size_table.csv",row.names=1)
Size<-t(ReorderData(tr, size))

skull_vent<-PhyloEM(tr, Skull_ventral, process = "scOU", Ncores = 6, methods.segmentation="best_single_move")
saveRDS(skull_vent,"skull_vent.rds")
skull_lat<-PhyloEM(tr, Skull_lateral, process = "scOU", Ncores = 6, methods.segmentation="best_single_move")
saveRDS(skull_lat,"skull_lat.rds")
mandible<-PhyloEM(tr, Mandible, process = "scOU", Ncores = 6, methods.segmentation="best_single_move")
saveRDS(mandible,"mandible.rds")
size<-PhyloEM(tr, Size, process = "OU", Ncores = 6, methods.segmentation="best_single_move")
saveRDS(size,"size.rds")

pdf("PhylogeneticEM_plots.pdf",width=15, height=5)
plot(skull_lat)
plot(skull_vent)
plot(mandible)
plot(size)
dev.off()

### processing RevBayes outputs
#read RevBayes trees
skull_ventral<-readTrees("skull_ventral_MAP.tre")
skull_lateral<-readTrees("skull_lateral_MAP.tre")
mandible<-readTrees("jaw_map.tre")
size<-readTrees("size_MAP.tre")
all_trees<-list(skull_ventral,skull_lateral,mandible,size)
title_list<-c("Ventral","Lateral","Mandible","Size")

#make rate plots
plots<-vector(mode="list",length=length(title_list))
for (i in 1:length(title_list))
{
	plots[[i]]<-plotTree(all_trees[[i]], ladderize=F, color_branch_by="branch_rates") + ggtitle(title_list[[i]])
}
pdf("RevBayes_rates.pdf")
plots
dev.off()

##get species rates from trees
labels<-data.frame(c(1:length(skull_lateral[[1]][[1]]@phylo$tip.label)),skull_lateral[[1]][[1]]@phylo$tip.label)
colnames(labels)<-c("Number","Taxon")
all_rates<-vector(mode="list",length=length(all_trees))
for (i in 1:length(all_trees))
{
	rates_data<-as.data.frame(all_trees[[i]][[1]][[1]]@data)
	rates_data$node<-as.numeric(rates_data$node)
	rates_data$branch_rates<-as.numeric(rates_data$branch_rates)
	x1<-rates_data[order(rates_data$node),]
	x2<-x1[is.na(x1$posterior),][c(2,6)]
	colnames(x2)<-c(title_list[i],"Number")
	rates_data<-merge(labels,x2,by="Number")
	all_rates[[i]]<-rates_data
}
names(all_rates)<-title_list

###join rates
skull_rates<-merge(all_rates[[1]][,2:3],all_rates[[2]][,2:3],by="Taxon")
other_rates<-merge(all_rates[[3]][,2:3],all_rates[[4]][,2:3],by="Taxon")
ako_rates<-merge(skull_rates,other_rates,by="Taxon")
ako_rates
write.csv(rates_table,"rates_table.csv")

cladeNames<-data.frame(Taxon=ako_rates$Taxon,clade=c(rep("Clade B",30),rep("Clade A",4),rep("Clade B",3),rep("Clade A",2),rep("Clade B",15),rep("Clade A",2),rep("Clade B",3)))
ako_rates2<-merge(ako_rates,cladeNames,by="Taxon")
ako_rates3<-ako_rates2[,-1]
rownames(ako_rates3)<-ako_rates2$Taxon

###test rate difference between clades
tree<-all_trees[[1]][[1]][[1]]@phylo
vrates<-procD.pgls(Ventral~Size+clade, data = ako_rates3, phy=tree, iter = 999, print.progress = TRUE)
lrates<-procD.pgls(Lateral~Size+clade, data = ako_rates3, phy=tree, iter = 999, print.progress = TRUE)
mrates<-procD.pgls(Mandible~Size+clade, data = ako_rates3, phy=tree, iter = 999, print.progress = TRUE)
srates<-procD.pgls(Size~clade, data = ako_rates3, phy=tree, iter = 999, print.progress = TRUE)
summary(vrates)
summary(lrates)
summary(mrates)
summary(srates)

wb<-createWorkbook()
addWorksheet(wb, "Ventral")
addWorksheet(wb, "Lateral")
addWorksheet(wb, "Mandible")
addWorksheet(wb, "Size")
writeData(wb, 1, vrates$aov.table, rowNames = TRUE)
writeData(wb, 2, lrates$aov.table, rowNames = TRUE)
writeData(wb, 3, mrates$aov.table, rowNames = TRUE)
writeData(wb, 4, srates$aov.table, rowNames = TRUE)
saveWorkbook(wb, "clade_rates_aov_tables.xlsx", overwrite = TRUE)

###############################################################################################################################################
###############################################################PHYLOMORPHOSPACE################################################################
###############################################################################################################################################

#get rates
all_rates2<-vector(mode="list",length=length(all_trees))
for (i in 1:length(all_trees))
{
	rates_data2<-as.data.frame(all_trees[[i]][[1]][[1]]@data)
	rates_data2$node<-as.numeric(rates_data2$node)
	rates_data2$branch_rates<-as.numeric(rates_data2$branch_rates)
	x<-rates_data2[order(rates_data2$node),]
	rates_data2<-x[c(2,6)]
	all_rates2[[i]]<-rates_data2
}
names(all_rates2)<-title_list
all_rates2

#create continuous color list
pcs<-readRDS("PCs_species.rds")
control_list<-vector(mode="list",length=length(pcs))
for (i in 1:length(pcs))
{
	rates<-all_rates2[[i]]
	ord<-rates[order(rates$branch_rates),]
	ord$label<-NA
	ord$label[1]<-1
	for (j in 2:(length(ord$branch_rates)-1))
	{
		if (ord$branch_rates[j]==ord$branch_rates[j-1]){ord$label[j]<-ord$label[j-1]} 
		else {ord$label[j]<-ord$label[j-1]+1}
	}
	palette<- colorRampPalette(c("blue","red"))(length(unique(ord$label))-1)
	cols<-data.frame(c(1:(length(unique(ord$label))-1)),palette)
	colnames(cols)<-c("label","hex")
	cols_joint<-merge(ord,cols,by="label")
	cols_joint[117,]<-c(31,0,60,"black")
	cols_joint$node<-as.numeric(cols_joint$node)
	cols_final<-cols_joint[order(cols_joint$node),]
	cols_final
	col<-cols_final$hex
	names(col)<-cols_final$node
	control_list[[i]]<-list(col.edge=col,col.node=col)
}
names(control_list)<-title_list[1:3]

#phylomorphospace

tree<-all_trees[[1]][[1]][[1]]@phylo

svg("skull_ventral_pms.svg")
colors<-control_list[[1]]$col.node
plot(NA,xlim=range(pcs[[1]][,1]),ylim=range(pcs[[1]][,2]),asp=1,axes=FALSE,xlab="PC1", ylab="PC2")
axis(1,at=seq(-0.06,0.06,by=0.02),cex.axis=1)
axis(2,at=seq(-0.07,0.07,by=0.02),cex.axis=1)
phylomorphospace(tree, pcs[[1]][,1:2], label=c("off"), node.size=c(0,1.2), control=control_list[[1]],add=T)
	add.color.bar(0.045,sort(colors[colors!="black"]),title="Rates",lims=range(all_rates2[[1]]$branch_rates[-60]),digits=10,prompt=FALSE,
	x=0.01*par()$usr[1],y=0.8*par()$usr[3])
dev.off()

svg("skull_lateral_pms.svg")
colors<-control_list[[2]]$col.node
plot(NA,xlim=range(pcs[[2]][,1]),ylim=range(pcs[[2]][,2]),asp=1,axes=FALSE,xlab="PC1", ylab="PC2")
axis(1,at=seq(-0.08,0.08,by=0.02),cex.axis=1)
axis(2,at=seq(-0.06,0.06,by=0.02),cex.axis=1)
phylomorphospace(tree, pcs[[2]][,1:2], label=c("off"), node.size=c(0,1.2), control=control_list[[2]],add=T)
	add.color.bar(0.045,sort(colors[colors!="black"]),title="Rates",lims=range(all_rates2[[2]]$branch_rates[-60]),digits=10,prompt=FALSE,
	x=0.01*par()$usr[1],y=0.8*par()$usr[3])
dev.off()

svg("mandible_pms.svg")
colors<-control_list[[3]]$col.node
plot(NA,xlim=range(pcs[[3]][,1]),ylim=range(pcs[[3]][,2]),asp=1,axes=FALSE,xlab="PC1", ylab="PC2")
axis(1,at=seq(-0.09,0.09,by=0.02),cex.axis=1)
axis(2,at=seq(-0.07,0.07,by=0.02),cex.axis=1)
phylomorphospace(tree, pcs[[3]][,1:2], label=c("off"), node.size=c(0,1.2), control=control_list[[3]],add=T)
	add.color.bar(0.045,sort(colors[colors!="black"]),title="Rates",lims=range(all_rates2[[3]]$branch_rates[-60]),digits=10,prompt=FALSE,
	x=0.01*par()$usr[1],y=0.8*par()$usr[3])
dev.off()

############################################################################################################################################
#######################################################LINEAGE DIVERSIFICATION##############################################################
############################################################################################################################################

tr <- read.tree("AKO_TREE.nwk")#read tree

##Read bammdata - BAMM outputs from two parallel runs
edata1 <- getEventData(tr, eventdata = "event_data_1.txt", burnin=0.25)
edata2 <- getEventData(tr, eventdata = "event_data_2.txt", burnin=0.25)
mcmcout1 <- read.csv("mcmc_out_1.txt", header=T)
mcmcout2 <- read.csv("mcmc_out_2.txt", header=T)

##Assessing MCMC convergence
#plot the log-likelihood trace of your MCMC output file
par(mfrow=c(2,1))
plot(mcmcout1$logLik ~ mcmcout1$generation)
plot(mcmcout2$logLik ~ mcmcout2$generation)

#burn-in
burnstart <- floor(0.25 * nrow(mcmcout1))
postburn1 <- mcmcout1[burnstart:nrow(mcmcout1),]
postburn2 <- mcmcout2[burnstart:nrow(mcmcout2),]

#assess ESS
effectiveSize(postburn1$N_shifts)
effectiveSize(postburn1$logLik)
effectiveSize(postburn2$N_shifts)
effectiveSize(postburn2$logLik)

#How many rate shifts?
shift_probs_1 <- summary(edata1)
shift_probs_2 <- summary(edata2)

#Single best shift configuration
best1 <- getBestShiftConfiguration(edata1, expectedNumberOfShifts=1)
best2 <- getBestShiftConfiguration(edata2, expectedNumberOfShifts=1)
par(mfrow=c(1,1))
plot(best1, lwd = 2)
addBAMMshifts(best1, cex=2.5)
plot(best2, lwd = 2)
addBAMMshifts(best2, cex=2.5)

#Bayesian credible sets of shift configurations
pdf("BAMM_credible_shift.pdf")
css1 <- credibleShiftSet(edata1, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css1$number.distinct
summary(css1)
plot(css1)
title("run1",line=1)
css2 <- credibleShiftSet(edata2, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css2$number.distinct
summary(css2)
plot(css2)
title("run2",line=1)
dev.off()

#Summary_table
wb<-createWorkbook()
addWorksheet(wb, "shifts_run1")
addWorksheet(wb, "shifts_run2")
addWorksheet(wb, "credible_shift_set_run1")
addWorksheet(wb, "credible_shift_set_run2")
writeData(wb, 1, shift_probs_1, rowNames = TRUE)
writeData(wb, 2, shift_probs_2, rowNames = TRUE)
writeData(wb, 3, summary(css1), rowNames = TRUE)
writeData(wb, 4, summary(css2), rowNames = TRUE)
saveWorkbook(wb, "BAMM_summary_tables.xlsx", overwrite = TRUE)

#Obtaining tip lambda, mu and diversification rates
lambda1<-edata1$meanTipLambda
lambda2<-edata2$meanTipLambda
lambda_final<-data.frame(Taxon=tr$tip.label,lambda=rowMeans(cbind(lambda1,lambda2)))
mu1<-edata1$meanTipMu
mu2<-edata2$meanTipMu
mu_final<-data.frame(Taxon=tr$tip.label,mu=rowMeans(cbind(mu1,mu2)))
BAMM_DR<-merge(lambda_final, mu_final, by="Taxon")
BAMM_DR
write.csv(BAMM_DR,"BAMM_tip_table.csv")

#Test diversification rate differences between clades
ako_rates<-read.csv("rates_table.csv")
cladeNames<-data.frame(Taxon=ako_rates$Taxon,clade=c(rep("Clade B",30),rep("Clade A",4),rep("Clade B",3),rep("Clade A",2),rep("Clade B",15),rep("Clade A",2),rep("Clade B",3)))
AkodonNames<-data.frame(Taxon=ako_rates$Taxon,ako=c(rep("Akodon",30),rep("others",29)))
allNames<-merge(cladeNames,AkodonNames,by="Taxon")
BAMM_full<-merge(BAMM_DR,allNames,by="Taxon")
BAMM_full2<-BAMM_full[,-1]
rownames(BAMM_full2)<-BAMM_full$Taxon
BAMM_full2

speciation_clades<-procD.pgls(lambda~clade, data = BAMM_full2, phy=tr, iter = 999, print.progress = TRUE)
extinction_clades<-procD.pgls(mu~clade, data = BAMM_full2, phy=tr, iter = 999, print.progress = TRUE)
speciation_ako<-procD.pgls(lambda~ako, data = BAMM_full2, phy=tr, iter = 999, print.progress = TRUE)
extinction_ako<-procD.pgls(mu~ako, data = BAMM_full2, phy=tr, iter = 999, print.progress = TRUE)
summary(speciation_clades)
summary(extinction_clades)
summary(speciation_ako)
summary(extinction_ako)

wb<-createWorkbook()
addWorksheet(wb, "speciation_clades")
addWorksheet(wb, "extinction_clades")
addWorksheet(wb, "speciation_ako")
addWorksheet(wb, "extinction_ako")
writeData(wb, 1, speciation_clades$aov.table, rowNames = TRUE)
writeData(wb, 2, extinction_clades$aov.table, rowNames = TRUE)
writeData(wb, 3, speciation_ako$aov.table, rowNames = TRUE)
writeData(wb, 4, extinction_ako$aov.table, rowNames = TRUE)
saveWorkbook(wb, "BAMM_clade_aov_tables.xlsx", overwrite = TRUE)

############################################################################################################################################
##############################################################SUMMARY VARIABLES#############################################################
############################################################################################################################################

##RAW_DATA
BAMM_DR<-read.csv("BAMM_tip_table.csv")
PD<-read.csv("part_disp_table.csv")
ako_rates<-read.csv("rates_table.csv")
DV<-as.data.frame(PD$Ventral)
DL<-as.data.frame(PD$Lateral)
DM<-as.data.frame(PD$Mandible)
SIZE<-read.csv("avg_size_table.csv",row.names=1)
data<-data.frame(ako_rates,DV,DL,DM,SIZE,BAMM_DR$lambda,BAMM_DR$mu)
rownames(data)<-NULL
colnames(data)<-c("N","Taxon","RV","RL","RM","RS","DV","DL","DM","SZ","SP","EX")
data<-data[,-1]
data
write.csv(data,"ako_data_RAW.csv")

#Z-SCORES
data_new<-data[,-1]
dd<-dim(data_new)
data2<-as.data.frame(matrix(nrow=dd[1],ncol=dd[2]))
rownames(data2)<-data[,1]
colnames(data2)<-c("RV","RL","RM","RS","DV","DL","DM","SZ","SP","EX")
for (i in 1:dd[2])
{
	data2[,i]<-(data_new[,i] - mean(data_new[,i]))/(sd(data_new[,i]))
}
data2
write.csv(data2,"ako_data_ZSCORE.csv")

#HEATMAP
tr<-read.tree("AKO_TREE.nwk")
labels<-as.data.frame(tr$tip.label)
colnames(labels)<-"taxa"
sep<-separate(labels,taxa,c("G","S"),sep="_")
unt<-unite(sep,taxa,c("G","S"),sep=" ")
tr$tip.label<-unt$taxa
d<- data.frame(node=c(1:117),color=c(rep("green",51),rep("#E31C79",8),"black",rep("green",50),rep("#E31C79",7)))
tree<-revts(ggtree(tr, layout = "rectangular", ladderize=TRUE, size=1, right=TRUE)%<+% d + aes(color=I(color),legend=T)) +
geom_tiplab(size=3.5, color="black", fontface = 3) + xlim_tree(12) +
scale_x_continuous(breaks=seq(-8,0,1), labels=abs(seq(-8,0,1))) + theme_tree2()
tree

rownames(data2)<-sort(unt$taxa)
map<-gheatmap(tree, data2, offset=4, width=1, colnames=FALSE, legend_title="SD units",low = "white", high = "black",
colnames_position = "bottom", colnames_angle = 0, hjust = 0.5, colnames_offset_y=0) +
scale_x_ggtree(breaks=seq(-8,0,1), labels=abs(seq(-8,0,1))) + theme(legend.position="right", axis.text = element_text(size = 12),
legend.title = element_text(size=12),legend.text = element_text(size=13))
map

cladeNames<-data.frame(Taxon=data$Taxon,Clade=c(rep("Clade B",30),rep("Clade A",4),rep("Clade B",3),rep("Clade A",2),rep("Clade B",15),rep("Clade A",2),rep("Clade B",3)))
f_data<-merge(data,cladeNames,by="Taxon")
g_data<-gather(f_data,key="Metric",value = "Value",c("RV","RL","RM","RS","DV","DL","DM","SZ","SP","EX"))
	
cols<-c("#E31C79","green")
g_data$Clade<-factor(g_data$Clade,levels=c("Clade A","Clade B"))
g_data$Metric<-factor(g_data$Metric,levels=c("RV","RL","RM","RS","DV","DL","DM","SZ","SP","EX"))
clade_plot<-ggplot(g_data, aes(x=Clade, y=Value, colour=Clade, fill=Clade),alpha=0.9) +
geom_boxplot(alpha=0.5, show.legend=T) + labs(x="", y="") +
scale_fill_manual(values=cols) +
scale_colour_manual(values=cols) +
#stat_summary(fun=ci, geom="line", linewidth=0.6) +
#stat_summary(fun=median, geom="point", size=2.5) +
theme_bw() +
theme(legend.title=element_blank(), legend.position="right", legend.text=element_text(size=12), plot.title=element_text(size=14,face= "bold",hjust=0.5),
axis.text.x=element_blank(), axis.text.y=element_text(size=10), axis.ticks.x=element_blank(), axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_text(size = 12)) +
scale_y_continuous(labels=label_number(accuracy=0.00001 ,decimal.mark="."))
full_clade_plot<-clade_plot+facet_wrap(vars(Metric), ncol = 10, nrow = 1, scales="free")
full_clade_plot

comp_plot<-ggarrange(map,full_clade_plot, ncol = 1, nrow = 2, common.legend = F,
font.label=list(size = 16, face="bold"), align="hv", labels=c("A","B"), heights = c(6,1))
comp_plot
ggsave("Heatmap.svg",comp_plot,width=14,height=10)

############################################################################################################################################
##################################################PHENOTYPE x DIVERSIFICATION###############################################################
############################################################################################################################################

###BAMMtools
edata1 <- getEventData(tr, eventdata = "event_data_1.txt", burnin=0.25)
edata2 <- getEventData(tr, eventdata = "event_data_2.txt", burnin=0.25)
edata<-list(edata1,edata2)
data<-read.csv("ako_data_ZSCORE.csv")
data_BAMM_test<-data[,2:9]
taxa_names<-data[,1]
prov<-as.data.frame(matrix(nrow=length(edata),ncol=4))
BAMM_list<-list(prov,prov)
BAMM_final<-vector(mode="list",length=length(BAMM_list))
for (n in 1:length(BAMM_list))
{
	BAMM_summary<-BAMM_list[[n]]
	
	#Speciation
	BAMM_lambda_results<-vector(mode="list",length=length(data_BAMM_test))
	for (i in 1:length(data_BAMM_test))
	{
		mod_data<-data_BAMM_test[,i]
		names(mod_data)<-taxa_names
		BAMM_lambda_results[[i]]<-traitDependentBAMM(edata[[n]], mod_data, reps = 1000, rate="speciation")
	}
	names(BAMM_lambda_results)<-names(data_BAMM_test)
	
	#Extinction
	BAMM_mu_results<-vector(mode="list",length=length(data_BAMM_test))
	for (i in 1:length(data_BAMM_test))
	{
		mod_data<-data_BAMM_test[,i]
		names(mod_data)<-taxa_names
		BAMM_mu_results[[i]]<-traitDependentBAMM(edata[[n]], mod_data, reps = 1000, rate="extinction")
	}
	names(BAMM_mu_results)<-names(data_BAMM_test)
	
	#Summarizing results
	BAMM_all_results<-c(BAMM_lambda_results,BAMM_mu_results)
	
	labs<-c("Skull ventral - rates","Skull lateral - rates","Mandible - rates","Size - rates",
	"Skull ventral - disparity","Skull lateral - disparity","Mandible - disparity","Size - logCS")
	labels<-c(labs,labs)
	colnames(BAMM_summary)<-c("Phenotype","Rate","Average correlation","p-value")
	for (i in 1:length(labels))
	{
		BAMM_summary[i,1]<-labels[i]
		BAMM_summary[i,2]<-BAMM_all_results[[i]]$rate
		BAMM_summary[i,3]<-BAMM_all_results[[i]]$estimate
		BAMM_summary[i,4]<-BAMM_all_results[[i]]$p.value
	}
	BAMM_final[[n]]<-BAMM_summary
}
BAMM_to_export<-cbind(BAMM_final[[1]],BAMM_final[[2]][,3:4])
colnames(BAMM_to_export)<-c("Phenotype","Rate","Avg. cor. 1","p-value 1","Avg. cor. 2","p-value 2")
BAMM_to_export$`p-adj. 1`<-p.adjust(BAMM_to_export$`p-value 1`)
BAMM_to_export$`p-adj. 2`<-p.adjust(BAMM_to_export$`p-value 2`)
BAMM_to_export<-BAMM_to_export[,c(1:4,7,5:6,8)]
BAMM_to_export
write.csv(BAMM_to_export,"BAMM_correlations_table.csv")
saveRDS(BAMM_final,"BAMM_correlation_results.rds")

###PGLS
tr <- read.tree("AKO_TREE.nwk")#read tree
data<-read.csv("ako_data_ZSCORE.csv")#Z-score data
data2<-data[,2:11]
rownames(data2)<-data[,1]
data2

#Speciation
RV_SP<-procD.pgls(RV~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RL_SP<-procD.pgls(RL~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RM_SP<-procD.pgls(RM~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RS_SP<-procD.pgls(RS~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DV_SP<-procD.pgls(DV~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DL_SP<-procD.pgls(DL~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DM_SP<-procD.pgls(DM~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)
SZ_SP<-procD.pgls(SZ~SP, data = data2, phy=tr, iter = 999, print.progress = TRUE)

#Extinction
RV_EX<-procD.pgls(RV~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RL_EX<-procD.pgls(RL~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RM_EX<-procD.pgls(RM~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RS_EX<-procD.pgls(RS~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DV_EX<-procD.pgls(DV~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DL_EX<-procD.pgls(DL~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DM_EX<-procD.pgls(DM~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)
SZ_EX<-procD.pgls(SZ~EX, data = data2, phy=tr, iter = 999, print.progress = TRUE)

#Phenotype
RV_DV<-procD.pgls(RV~DV, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RL_DL<-procD.pgls(RL~DL, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RM_DM<-procD.pgls(RM~DM, data = data2, phy=tr, iter = 999, print.progress = TRUE)
RS_SZ<-procD.pgls(RS~SZ, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DV_RV<-procD.pgls(DV~RV, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DL_RL<-procD.pgls(DL~RL, data = data2, phy=tr, iter = 999, print.progress = TRUE)
DM_RM<-procD.pgls(DM~RM, data = data2, phy=tr, iter = 999, print.progress = TRUE)
SZ_RS<-procD.pgls(SZ~RS, data = data2, phy=tr, iter = 999, print.progress = TRUE)

PGLS_results<-list(RV_SP, RL_SP, RM_SP, RS_SP, DV_SP, DL_SP, DM_SP, SZ_SP, RV_EX, RL_EX, RM_EX, RS_EX, DV_EX, DL_EX, DM_EX, SZ_EX, RV_DV, RL_DL, RM_DM, RS_SZ, DV_RV, DL_RL, DM_RM, SZ_RS)
names(PGLS_results)<-c("RV_SP","RL_SP","RM_SP","RS_SP","DV_SP","DL_SP","DM_SP","SZ_SP","RV_EX","RL_EX","RM_EX","RS_EX","DV_EX","DL_EX","DM_EX","SZ_EX","RV_DV","RL_DL","RM_DM","RS_SZ","DV_RV","DL_RL","DM_RM","SZ_RS")
saveRDS(PGLS_results,"PGLS_results.rds")

wb<-createWorkbook()
addWorksheet(wb, "RV_SP")
addWorksheet(wb, "RL_SP")
addWorksheet(wb, "RM_SP")
addWorksheet(wb, "RS_SP")
addWorksheet(wb, "DV_SP")
addWorksheet(wb, "DL_SP")
addWorksheet(wb, "DM_SP")
addWorksheet(wb, "SZ_SP")
addWorksheet(wb, "RV_EX")
addWorksheet(wb, "RL_EX")
addWorksheet(wb, "RM_EX")
addWorksheet(wb, "RS_EX")
addWorksheet(wb, "DV_EX")
addWorksheet(wb, "DL_EX")
addWorksheet(wb, "DM_EX")
addWorksheet(wb, "SZ_EX")
addWorksheet(wb, "RV_DV")
addWorksheet(wb, "RL_DL")
addWorksheet(wb, "RM_DM")
addWorksheet(wb, "RS_SZ")
addWorksheet(wb, "DV_RV")
addWorksheet(wb, "DL_RL")
addWorksheet(wb, "DM_RM")
addWorksheet(wb, "SZ_RS")
writeData(wb, 1,  RV_SP$aov.table, rowNames = TRUE)
writeData(wb, 2,  RL_SP$aov.table, rowNames = TRUE)
writeData(wb, 3,  RM_SP$aov.table, rowNames = TRUE)
writeData(wb, 4,  RS_SP$aov.table, rowNames = TRUE)
writeData(wb, 5,  DV_SP$aov.table, rowNames = TRUE)
writeData(wb, 6,  DL_SP$aov.table, rowNames = TRUE)
writeData(wb, 7,  DM_SP$aov.table, rowNames = TRUE)
writeData(wb, 8,  SZ_SP$aov.table, rowNames = TRUE)
writeData(wb, 9,  RV_EX$aov.table, rowNames = TRUE)
writeData(wb, 10, RL_EX$aov.table, rowNames = TRUE)
writeData(wb, 11, RM_EX$aov.table, rowNames = TRUE)
writeData(wb, 12, RS_EX$aov.table, rowNames = TRUE)
writeData(wb, 13, DV_EX$aov.table, rowNames = TRUE)
writeData(wb, 14, DL_EX$aov.table, rowNames = TRUE)
writeData(wb, 15, DM_EX$aov.table, rowNames = TRUE)
writeData(wb, 16, SZ_EX$aov.table, rowNames = TRUE)
writeData(wb, 17, RV_DV$aov.table, rowNames = TRUE)
writeData(wb, 18, RL_DL$aov.table, rowNames = TRUE)
writeData(wb, 19, RM_DM$aov.table, rowNames = TRUE)
writeData(wb, 20, RS_SZ$aov.table, rowNames = TRUE)
writeData(wb, 21, DV_RV$aov.table, rowNames = TRUE)
writeData(wb, 22, DL_RL$aov.table, rowNames = TRUE)
writeData(wb, 23, DM_RM$aov.table, rowNames = TRUE)
writeData(wb, 24, SZ_RS$aov.table, rowNames = TRUE)
saveWorkbook(wb, "PGLS_aov_tables.xlsx", overwrite = TRUE)

PGLS_coeff<-as.data.frame(matrix(nrow=length(PGLS_results),ncol=5))
rownames(PGLS_coeff)<-names(PGLS_results)
colnames(PGLS_coeff)<-c("Intercept","Slope","F-statistic","Effect size","p-value")
for (i in 1:length(PGLS_results))
{
	PGLS_coeff[i,1]<-round(PGLS_results[[i]]$pgls.coefficients[,1][[1]],3)
	PGLS_coeff[i,2]<-round(PGLS_results[[i]]$pgls.coefficients[,1][[2]],3)
	PGLS_coeff[i,3]<-round(PGLS_results[[i]]$aov.table$F[[1]],3)
	PGLS_coeff[i,4]<-round(PGLS_results[[i]]$aov.table$Z[[1]],3)
	PGLS_coeff[i,5]<-round(PGLS_results[[i]]$aov.table$`Pr(>F)`[[1]],3)
}
PGLS_coeff
PGLS_coeff$`p-adj.`<-p.adjust(PGLS_coeff$`p-value`)
rownames(PGLS_coeff)<-c("Skull ventral rates x Speciation rates","Skull lateral rates x Speciation rates",
"Mandible rates x Speciation rates","Size rates x Speciation rates","Skull ventral disparity x Speciation rates",
"Skull lateral disparity x Speciation rates","Mandible disparity x Speciation rates","Size x Speciation rates",
"Skull ventral rates x Extinction rates","Skull lateral rates x Extinction rates","Mandible rates x Extinction rates",
"Size rates x Extinction rates","Skull ventral disparity x Extinction rates","Skull lateral disparity x Extinction rates",
"Mandible disparity x Extinction rates","Size x Extinction rates","Skull ventral rates x Skull ventral disparity",
"Skull lateral rates x Skull lateral disparity","Mandible rates x Mandible disparity","Size rates x Size",
"Skull ventral disparity x Skull ventral rates","Skull lateral disparity x Skull lateral rates","Mandible disparity x Mandible rates","Size x Size rates")
PGLS_coeff
write.csv(PGLS_coeff,"PGLS_summary_table.csv")

#Plot PGLS regressions
models<-data.frame(model=names(PGLS_results))
models2<-separate(models,col=model,into=c("y","x"),sep="_")
lab1<-unique(models2$y)
lab2<-unique(models2$x)
index<-data.frame(model=unique(c(lab1,lab2)),label=c("Skull ventral - rates","Skull lateral - rates","Mandible - rates","Size - rates",
"Skull ventral - disparity","Skull lateral - disparity","Mandible - disparity","Size - logCS","Speciation rates","Extinction rates"))
models3<-models2
models3$D<-index$label[match(models3$y,index$model,)]
models3$I<-index$label[match(models3$x,index$model,)]

plot_list<-vector(mode="list",length=length(PGLS_results))
for (i in 1:length(PGLS_results))
{
	to_plot<-data.frame(x=PGLS_results[[i]]$X[,2],y=PGLS_results[[i]]$Y[,1])
	plot_list[[i]]<-ggplot(data = to_plot, aes(x, y)) + geom_point() + labs(x=models3$I[i],y=models3$D[i]) +
	geom_abline(intercept=PGLS_coeff$Intercept[i], slope=PGLS_coeff$Slope[i],col="red") +
	theme_bw() + theme(legend.title=element_blank(), legend.position="right", legend.text=element_text(size=12), 
	axis.text=element_text(size=10), axis.ticks.x=element_blank(), axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),
	panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}

PGLS_plot<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],
plot_list[[18]],plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]]
,ncol = 8, nrow = 4, common.legend = TRUE, legend="bottom")
 ggsave("PGLS_plot.svg",PGLS_plot,width=16,height=8)