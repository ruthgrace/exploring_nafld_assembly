#!/usr/bin/env Rscript
library(zCompositions)
library(compositions)
library(ALDEx2)
library(stringr)

d <- read.table("all_counts_with_seq_length_zeros_removed.txt",sep="\t",quote="",comment.char="",header=TRUE)

colnames(d)[c(3:22)] <- str_extract(colnames(d)[c(3:22)],"^[A-Z]+_[0-9]+")

carbs <- d[which(d$subsys1 == "Carbohydrates"),]
lipids <- d[which(d$subsys1 == "Fatty Acids, Lipids, and Isoprenoids"),]
amino <- d[which(d$subsys1 == "Amino Acids and Derivatives"),]

carbs.samples <- carbs[,c(3:(ncol(carbs)-4))]
lipids.samples <- lipids[,c(3:(ncol(lipids)-4))]
amino.samples <- amino[,c(3:(ncol(amino)-4))]
d.samples <- d[,c(3:(ncol(d)-4))]

# aggregate by subsys4
grouping <- carbs$subsys4
carbs.pca.data <- aggregate(carbs.samples,list(grouping),sum)
rownames(carbs.pca.data) <- carbs.pca.data$Group.1
carbs.pca.data <- carbs.pca.data[,c(2:ncol(carbs.pca.data))]

grouping <- lipids$subsys4
lipids.pca.data <- aggregate(lipids.samples,list(grouping),sum)
rownames(lipids.pca.data) <- lipids.pca.data$Group.1
lipids.pca.data <- lipids.pca.data[,c(2:ncol(lipids.pca.data))]

grouping <- amino$subsys4
amino.pca.data <- aggregate(amino.samples,list(grouping),sum)
rownames(amino.pca.data) <- amino.pca.data$Group.1
amino.pca.data <- amino.pca.data[,c(2:ncol(amino.pca.data))]

grouping <- d$subsys4
d.pca.data <- aggregate(d.samples,list(grouping),sum)
rownames(d.pca.data) <- d.pca.data$Group.1
d.pca.data <- d.pca.data[,c(2:ncol(d.pca.data))]

sparseness <- apply(d.pca.data,1,function(x) { return(length(which(x > 0))); } )
d.filter.pca.data  <- d.pca.data[which(sparseness > 2),]

groups <- colnames(carbs.pca.data)
groups[grepl("^H",groups)] <- "nash"
groups[grepl("^C",groups)] <- "healthy"

# adjust zeros
carbs.pca.adj.zero <- t(cmultRepl(t(carbs.pca.data),method="CZM"))
lipids.pca.adj.zero <- t(cmultRepl(t(lipids.pca.data),method="CZM"))
amino.pca.adj.zero <- t(cmultRepl(t(amino.pca.data),method="CZM"))
d.pca.adj.zero <- t(cmultRepl(t(d.pca.data),method="CZM"))
d.filter.pca.adj.zero <- t(cmultRepl(t(d.filter.pca.data),method="CZM"))

carbs.pca.adj.zero <- carbs.pca.adj.zero[order(apply(carbs.pca.adj.zero,1,sum),decreasing=TRUE),]
lipids.pca.adj.zero <- lipids.pca.adj.zero[order(apply(lipids.pca.adj.zero,1,sum),decreasing=TRUE),]
amino.pca.adj.zero <- amino.pca.adj.zero[order(apply(amino.pca.adj.zero,1,sum),decreasing=TRUE),]
d.pca.adj.zero <- d.pca.adj.zero[order(apply(d.pca.adj.zero,1,sum),decreasing=TRUE),]
d.filter.pca.adj.zero <- d.filter.pca.adj.zero[order(apply(d.filter.pca.adj.zero,1,sum),decreasing=TRUE),]

carbs.names <- rownames(carbs.pca.adj.zero)
lipids.names <- rownames(lipids.pca.adj.zero)
amino.names <- rownames(amino.pca.adj.zero)
d.names <- rownames(d.pca.adj.zero)
d.filter.names <- rownames(d.filter.pca.adj.zero)

carbs.prop <- apply(carbs.pca.adj.zero,2,function(x){x/sum(x)})
lipids.prop <- apply(lipids.pca.adj.zero,2,function(x){x/sum(x)})
amino.prop <- apply(amino.pca.adj.zero,2,function(x){x/sum(x)})
d.prop <- apply(d.pca.adj.zero,2,function(x){x/sum(x)})
d.filter.prop <- apply(d.filter.pca.adj.zero,2,function(x){x/sum(x)})

carbs.clr <- t(apply(carbs.prop,2,function(x){log(x) - mean(log(x))}))
lipids.clr <- t(apply(lipids.prop,2,function(x){log(x) - mean(log(x))}))
amino.clr <- t(apply(amino.prop,2,function(x){log(x) - mean(log(x))}))
d.clr <- t(apply(d.prop,2,function(x){log(x) - mean(log(x))}))
d.filter.clr <- t(apply(d.filter.prop,2,function(x){log(x) - mean(log(x))}))

carbs.pcx <- prcomp(carbs.clr)
lipids.pcx <- prcomp(lipids.clr)
amino.pcx <- prcomp(amino.clr)
d.pcx <- prcomp(d.clr)
d.filter.pcx <- prcomp(d.filter.clr)

conds <- data.frame(groups)
colnames(conds) <- "cond"

palette=palette(c(rgb(1,0,0,0.6), rgb(0,0,1,0.3), rgb(0,1,1,0.6)))
dev.off()

pdf("biplots.pdf")

layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
par(mgp=c(2,0.5,0))
# make a covariance biplot of the data with compositions function
coloredBiplot(carbs.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(carbs.pcx$sdev[1]^2)/mvar(carbs.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(carbs.pcx$sdev[2]^2)/mvar(carbs.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Carbohydrate functions biplot")
barplot(carbs.pcx$sdev^2/mvar(carbs.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

myrgb <- col2rgb("black")
mycolor <- rgb(myrgb[1], myrgb[2], myrgb[3], 0.3)

#pca only, no biplot
mylabels <- str_extract(rownames(carbs.pcx$x),"^[A-Z]+_[0-9]+")
layout(matrix(c(1),1,1,byrow=T),widths=c(7),heights=c(7))
plot(carbs.pcx$x[,1],carbs.pcx$x[,2],col="white",xlim=c(min(carbs.pcx$x[,1])-10,max(carbs.pcx$x[,1])+10),ylim=c(min(carbs.pcx$x[,2]) - 10, max(carbs.pcx$x[,2]) + 10),
xlab=paste("PC1 ", round (sum(carbs.pcx$sdev[1]^2)/mvar(carbs.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(carbs.pcx$sdev[2]^2)/mvar(carbs.clr),3), sep=""),
,main="Principal Components Analysis\nCarbohydrate subset")
text(carbs.pcx$x[,1],carbs.pcx$x[,2],labels = mylabels,col=c(rep("black",10),rep("red",10)))

points <- c(rep("o", length(dimnames(carbs.pcx$rotation)[[1]])))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(carbs.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(carbs.pcx$sdev[1]^2)/mvar(carbs.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(carbs.pcx$sdev[2]^2)/mvar(carbs.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Carbohydrate functions biplot")
barplot(carbs.pcx$sdev^2/mvar(carbs.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot


layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(lipids.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(lipids.pcx$sdev[1]^2)/mvar(lipids.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(lipids.pcx$sdev[2]^2)/mvar(lipids.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Lipid functions biplot")
barplot(lipids.pcx$sdev^2/mvar(lipids.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

mylabels <- str_extract(rownames(lipids.pcx$x),"^[A-Z]+_[0-9]+")
layout(matrix(c(1),1,1,byrow=T),widths=c(7),heights=c(7))
plot(lipids.pcx$x[,1],lipids.pcx$x[,2],col="white",xlim=c(min(lipids.pcx$x[,1])-10,max(lipids.pcx$x[,1])+10),ylim=c(min(lipids.pcx$x[,2]) - 10, max(lipids.pcx$x[,2]) + 10),
xlab=paste("PC1 ", round (sum(lipids.pcx$sdev[1]^2)/mvar(lipids.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(lipids.pcx$sdev[2]^2)/mvar(lipids.clr),3), sep=""),
,main="Principal Components Analysis\nLipid subset")
text(lipids.pcx$x[,1],lipids.pcx$x[,2],labels = mylabels,col=c(rep("black",10),rep("red",10)))

points <- c(rep("o", length(dimnames(lipids.pcx$rotation)[[1]])))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(lipids.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(lipids.pcx$sdev[1]^2)/mvar(lipids.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(lipids.pcx$sdev[2]^2)/mvar(lipids.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Lipid functions biplot")
barplot(lipids.pcx$sdev^2/mvar(lipids.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot


layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(amino.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(amino.pcx$sdev[1]^2)/mvar(amino.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(amino.pcx$sdev[2]^2)/mvar(amino.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Amino acid functions biplot")
barplot(amino.pcx$sdev^2/mvar(amino.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

mylabels <- str_extract(rownames(amino.pcx$x),"^[A-Z]+_[0-9]+")
layout(matrix(c(1),1,1,byrow=T),widths=c(7),heights=c(7))
plot(amino.pcx$x[,1],amino.pcx$x[,2],col="white",xlim=c(min(amino.pcx$x[,1])-10,max(amino.pcx$x[,1])+10),ylim=c(min(amino.pcx$x[,2]) - 10, max(amino.pcx$x[,2]) + 10),
xlab=paste("PC1 ", round (sum(amino.pcx$sdev[1]^2)/mvar(amino.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(amino.pcx$sdev[2]^2)/mvar(amino.clr),3), sep=""),
,main="Principal Components Analysis\nAmino acid subset")
text(amino.pcx$x[,1],amino.pcx$x[,2],labels = mylabels,col=c(rep("black",10),rep("red",10)))

points <- c(rep("o", length(dimnames(amino.pcx$rotation)[[1]])))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(amino.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(amino.pcx$sdev[1]^2)/mvar(amino.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(amino.pcx$sdev[2]^2)/mvar(amino.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Amino acid functions biplot")
barplot(amino.pcx$sdev^2/mvar(amino.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot


layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(d.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.pcx$sdev[1]^2)/mvar(d.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.pcx$sdev[2]^2)/mvar(d.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Principal Components Analysis")
barplot(d.pcx$sdev^2/mvar(d.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

mylabels <- str_extract(rownames(d.pcx$x),"^[A-Z]+_[0-9]+")
layout(matrix(c(1),1,1,byrow=T),widths=c(7),heights=c(7))
plot(d.pcx$x[,1],d.pcx$x[,2],col="white",xlim=c(min(d.pcx$x[,1])-10,max(d.pcx$x[,1])+10),ylim=c(min(d.pcx$x[,2]) - 10, max(d.pcx$x[,2]) + 10),
xlab=paste("PC1 ", round (sum(d.pcx$sdev[1]^2)/mvar(d.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.pcx$sdev[2]^2)/mvar(d.clr),3), sep=""),
,main="Principal Components Analysis")
text(d.pcx$x[,1],d.pcx$x[,2],labels = mylabels,col=c(rep("black",10),rep("red",10)))

points <- c(rep("o", length(dimnames(d.pcx$rotation)[[1]])))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(d.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.pcx$sdev[1]^2)/mvar(d.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.pcx$sdev[2]^2)/mvar(d.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Principal Components Analysis")
barplot(d.pcx$sdev^2/mvar(d.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot


layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(d.filter.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.filter.pcx$sdev[1]^2)/mvar(d.filter.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.filter.pcx$sdev[2]^2)/mvar(d.filter.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Principal Components Analysis\nwith sparsity filter")
barplot(d.filter.pcx$sdev^2/mvar(d.filter.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

mylabels <- str_extract(rownames(d.filter.pcx$x),"^[A-Z]+_[0-9]+")
layout(matrix(c(1),1,1,byrow=T),widths=c(7),heights=c(7))
plot(d.filter.pcx$x[,1],d.filter.pcx$x[,2],col="white",xlim=c(min(d.filter.pcx$x[,1])-10,max(d.filter.pcx$x[,1])+10),ylim=c(min(d.filter.pcx$x[,2]) - 10, max(d.filter.pcx$x[,2]) + 10),
xlab=paste("PC1 ", round (sum(d.filter.pcx$sdev[1]^2)/mvar(d.filter.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.filter.pcx$sdev[2]^2)/mvar(d.filter.clr),3), sep=""),
,main="Principal Components Analysis\nwith sparsity filter")
text(d.filter.pcx$x[,1],d.filter.pcx$x[,2],labels = mylabels,col=c(rep("black",10),rep("red",10)))

points <- c(rep("o", length(dimnames(d.filter.pcx$rotation)[[1]])))
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(d.filter.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.filter.pcx$sdev[1]^2)/mvar(d.filter.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.filter.pcx$sdev[2]^2)/mvar(d.filter.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Principal Components Analysis\nwith sparsity filter")
barplot(d.filter.pcx$sdev^2/mvar(d.filter.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

library(ALDEx2)
# 
# x <- aldex(d.pca.data, groups)
# write.table(x,file="aldex_output_pca_plot_script.txt",sep="\t",quote=FALSE)

x <- read.table(file="aldex_output_pca_plot_script.txt",sep="\t",quote="",comment.char="",header=TRUE,row.names=1)

print(paste("rownames are in the same order in aldex output: ",all.equal(rownames(d.pca.data),rownames(x)),sep=" "))

high.effect <- which(abs(x$effect) > 1)

points <- c(rep("", length(dimnames(d.pcx$rotation)[[1]])))
points[high.effect] <- dimnames(d.pcx$rotation)[[1]][high.effect]
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(d.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.pcx$sdev[1]^2)/mvar(d.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.pcx$sdev[2]^2)/mvar(d.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Principal Components Analysis")
barplot(d.pcx$sdev^2/mvar(d.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot


points <- c(rep("", length(dimnames(carbs.pcx$rotation)[[1]])))
points[which(dimnames(carbs.pcx$rotation)[[1]] %in% rownames(x)[high.effect])] <- dimnames(carbs.pcx$rotation)[[1]][which(dimnames(carbs.pcx$rotation)[[1]] %in% rownames(x)[high.effect])]
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(carbs.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(carbs.pcx$sdev[1]^2)/mvar(carbs.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(carbs.pcx$sdev[2]^2)/mvar(carbs.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Carbohydrate functions biplot")
barplot(carbs.pcx$sdev^2/mvar(carbs.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot


points <- c(rep("", length(dimnames(lipids.pcx$rotation)[[1]])))
lipids.high.effect <- which(x$effect[match(dimnames(lipids.pcx$rotation)[[1]], rownames(x))] > 0.5)
points[lipids.high.effect] <- dimnames(lipids.pcx$rotation)[[1]][lipids.high.effect]
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(lipids.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(lipids.pcx$sdev[1]^2)/mvar(lipids.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(lipids.pcx$sdev[2]^2)/mvar(lipids.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Lipid functions biplot")
barplot(lipids.pcx$sdev^2/mvar(lipids.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

points <- c(rep("", length(dimnames(amino.pcx$rotation)[[1]])))
amino.high.effect <- which(x$effect[match(dimnames(amino.pcx$rotation)[[1]], rownames(x))] > 0.5)
points[amino.high.effect] <- dimnames(amino.pcx$rotation)[[1]][amino.high.effect]
layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
coloredBiplot(amino.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(amino.pcx$sdev[1]^2)/mvar(amino.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(amino.pcx$sdev[2]^2)/mvar(amino.clr),3), sep=""),
xlabs.col=c(rep("black",10),rep("red",10)),
ylabs=points,
expand=0.8,var.axes=FALSE, scale=1, main="Amino acids functions biplot")
barplot(amino.pcx$sdev^2/mvar(amino.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

dev.off()

