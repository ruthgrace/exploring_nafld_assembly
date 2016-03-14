#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

print(args)
if (length(args) != 2) {
  print("This script takes two arguments: the count file path, and the file path of the counts with all-zero features removed. The first two columns of the count table should be refseqid and length, and the last four should be subsys4, subsys1, subsys2, and subsys3. In the middle should be all the counts per sample.")
} else {
  d <- read.table(args[1],sep="\t",header=TRUE,quote="",comment.char="",stringsAsFactors=FALSE)
  rowsums <- apply(d,1,function(x) { sum(as.numeric(x[c(3:(length(x)-4))])) } )
  d <- d[which(rowsums != 0),]
  print(paste("removed",length(which(rowsums==0)),"features with zero counts for all samples from",length(rowsums),"total features, leaving",length(which(rowsums!=0)),"features left"))
  write.table(d,file=args[2],sep="\t",quote=FALSE,row.names=FALSE)
}