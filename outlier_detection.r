library(zCompositions)
library(compositions)
library(ALDEx2)
library(edgeR)
library(gplots)

d <- read.table("all_counts_with_seq_length_zeros_removed.txt",sep="\t",quote="",comment.char="",header=TRUE)

d.agg <- d
d.agg <- d.agg[,c(3:22)]
d.agg <- aggregate(d.agg,by=list(d$subsys4),sum)

rownames(d.agg) <- d.agg$Group.1
d.agg <- d.agg[,c(2:ncol(d.agg))]

d.agg <- t(d.agg)

# remove rows with 0 counts
d.agg.gt0 <- t(d.agg[,apply(d.agg, 2, sum) > 0])

# estimate 0 values (zCompositions)
d.agg.n0 <- cmultRepl(t(d.agg.gt0), method="CZM", label=0)

# clr transform
d.agg.n0.clr <- t(apply(d.agg.n0, 1, function(x){log2(x) - mean(log2(x))}))

mvar.agg.clr <- mvar(d.agg.n0.clr)
pcx.agg  <- prcomp(d.agg.n0.clr)
par(mfrow=c(1,2))

pdf("outlier_exploration.pdf")

plot(pcx.agg$x[,1], pcx.agg$x[,2], pch=NA,
    xlim=c(min(pcx.agg$x[,1] )-5, max(pcx.agg$x[,1] + 5)),
    xlab=paste("PC1: ", round(sum(pcx.agg$sdev[1]^2)/mvar.agg.clr, 3)),
    ylab=paste("PC2: ", round(sum(pcx.agg$sdev[2]^2)/mvar.agg.clr, 3)),
    main="PCA score plot")
text(pcx.agg$x[,1], pcx.agg$x[,2], labels=rownames(pcx.agg$x), cex=0.7)
#hist(apply(pcx.agg$x,1,function(x){sum(x^2/mvar.agg.clr)}), breaks=100, xlab=NULL)
#which(apply(pcx.agg$x,1,function(x){sum(x^2/mvar.agg.clr)}) > 1)


# density of scores along PC2 for SNF
plot(density(pcx.agg$x[,2][pcx.agg$x[,1] <0]), lty=1, lwd=3,
    xlim=c(-30,70), main="Density plot", xlab="PC2 Value", ylab="Density")
# density for wild type
points(density(pcx.agg$x[,2][pcx.agg$x[,1] >0]), lty=2, lwd=3, type="l")
# supports the idea that +/- 20 is a reasonable cutoff for sample exclusion
# of the first two components.

# check each independently
mvar.s <- mvar(d.agg.n0.clr[1:10,])
pcx.s <- prcomp(d.agg.n0.clr[1:10,])
# biplot(pcx.s, cex=c(0.5,0.01), var.axes=FALSE, scale=0)

mvar.w <- mvar(d.agg.n0.clr[11:20,])
pcx.w <- prcomp(d.agg.n0.clr[11:20,])
# biplot(pcx.w, cex=c(0.5,0.01), var.axes=FALSE, scale=0)

# histogram of the variance of the sample as a function of the total variance
# identifies a similar set of outliers, but some clearly have significant variance
# on other component axes
plot(density(apply(pcx.s$x,1,function(x){sum(x^2/mvar.s)})), main="SNF2",
    xlab="Var Fraction", ylab="Density" )
cut.s <- median(apply(pcx.s$x,1,function(x){sum(x^2/mvar.s)})) +
    2 * IQR(apply(pcx.s$x,1,function(x){sum(x^2/mvar.s)}))
abline(v=cut.s, lty=2)

plot(density(apply(pcx.w$x,1,function(x){sum(x^2/mvar.w)})), main="WT",
    xlab="Var Fraction", ylab="Density" )
cut.w <- median(apply(pcx.w$x,1,function(x){sum(x^2/mvar.w)})) +
    2 * IQR(apply(pcx.w$x,1,function(x){sum(x^2/mvar.w)}))
abline(v=cut.w, lty=2)

# list the outliers
bad.s <- names(which(apply(pcx.s$x,1,function(x){sum(x^2/mvar.s)})>cut.s))
bad.w <- names(which(apply(pcx.w$x,1,function(x){sum(x^2/mvar.w)})>cut.w))



