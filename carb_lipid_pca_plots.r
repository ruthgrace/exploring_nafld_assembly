#!/usr/bin/env Rscript

d <- read.table("all_counts_with_seq_length_zeros_removed.txt",sep="\t",quote="",comment.char="",header=TRUE)

carbs <- d[which(d$subsys1 == "Carbohydrates"),]
lipids <- d[which(d$subsys1 == "Fatty Acids, Lipids, and Isoprenoids"),]

carbs.samples <- carbs[,c(3:(ncol(carbs)-4))]
lipids.samples <- lipids[,c(3:(ncol(lipids)-4))]

# aggregate by subsys4
grouping <- carbs$subsys4
carbs.pca.data <- aggregate(carbs.samples,list(grouping),sum)
rownames(carbs.pca.data) <- carbs.pca.data$Group.1
carbs.pca.data <- carbs.pca.data[,c(2:ncol(carbs.pca.data))]

grouping <- lipids$subsys4
lipids.pca.data <- aggregate(lipids.samples,list(grouping),sum)
rownames(lipids.pca.data) <- lipids.pca.data$Group.1
lipids.pca.data <- lipids.pca.data[,c(2:ncol(lipids.pca.data))]













library(zCompositions)
library(compositions)
library(ALDEx2)
library(stringr)

groups <- colnames(carbs.pca.data)
groups[grepl("^H",groups)] <- "nash"
groups[grepl("^C",groups)] <- "healthy"

# adjust zeros
carbs.pca.adj.zero <- t(cmultRepl(t(carbs.pca.data),method="CZM"))
lipids.pca.adj.zero <- t(cmultRepl(t(lipids.pca.data),method="CZM"))

carbs.pca.adj.zero <- carbs.pca.adj.zero[order(apply(carbs.pca.adj.zero,1,sum),decreasing=TRUE),]
lipids.pca.adj.zero <- lipids.pca.adj.zero[order(apply(lipids.pca.adj.zero,1,sum),decreasing=TRUE),]

carbs.names <- rownames(carbs.pca.adj.zero)
lipids.names <- rownames(lipids.pca.adj.zero)

carbs.prop <- apply(carbs.pca.adj.zero,2,function(x){x/sum(x)})
lipids.prop <- apply(lipids.pca.adj.zero,2,function(x){x/sum(x)})

carbs.clr <- t(apply(carbs.prop,2,function(x){log(x) - mean(log(x))}))
lipids.clr <- t(apply(lipids.prop,2,function(x){log(x) - mean(log(x))}))

carbs.pcx <- prcomp(carbs.clr)
lipids.pcx <- prcomp(lipids.clr)

conds <- data.frame(groups)
colnames(conds) <- "cond"

palette=palette(c(rgb(1,0,0,0.6), rgb(0,0,1,0.6), rgb(0,1,1,0.6)))

pdf("biplots.pdf")

layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
par(mgp=c(2,0.5,0))
# make a covariance biplot of the data with compositions function
coloredBiplot(carbs.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(carbs.pcx$sdev[1]^2)/mvar(carbs.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(carbs.pcx$sdev[2]^2)/mvar(carbs.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Carbohydrate functions biplot")
barplot(carbs.pcx$sdev^2/mvar(carbs.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(lipids.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(lipids.pcx$sdev[1]^2)/mvar(lipids.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(lipids.pcx$sdev[2]^2)/mvar(lipids.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=FALSE, scale=1, main="Lipid functions biplot")
barplot(lipids.pcx$sdev^2/mvar(lipids.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

dev.off()



