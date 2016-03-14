#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

print(args)
if (length(args) != 3) {
  print("Please run with 3 arguments, the file path for the two tab-separated count files you want to amalgamate, and the file path for the amalgamated output. The first two columns of your count tables should be refseqid and length, and the last four columns should be subsys4, subsys1, subsys2, and subsys3")
} else {
  counts.a <- read.table(args[1],sep="\t",header=TRUE,quote="",comment.char="",stringsAsFactors=FALSE)
  counts.b <- read.table(args[2],sep="\t",header=TRUE,quote="",comment.char="",stringsAsFactors=FALSE)
  
  samples.a <- colnames(counts.a)[c(3:(ncol(counts.a)-4))]
  samples.b <- colnames(counts.b)[c(3:(ncol(counts.b)-4))]
  a.b <- samples.a[which(!(samples.a %in% samples.b))]
  b.a <- samples.b[which(!(samples.b %in% samples.a))]

  if (length(a.b) > 0) {
    print("Aborted. Samples in first count table but not second:")
    print(a.b)
  }
  if (length(b.a) > 0) {
    print("Aborted. Samples in second count table but not first:")
    print(b.a)
  }

  if (length(a.b) == 0 && length(b.a) == 0) {
    amalgamated <- data.frame(matrix(nrow=(nrow(counts.a)+nrow(counts.b)),ncol=ncol(counts.a)))
    colnames(amalgamated) <- colnames(counts.a)
    amalgamated[c(1:nrow(counts.a)),] <- counts.a
    counts.b[,c(3:(ncol(counts.b)-4))] <- counts.b[,2+(match( colnames(counts.a)[c(3:(ncol(counts.a)-4))], colnames(counts.b)[c(3:(ncol(counts.b)-4))] ))]
    amalgamated[c((nrow(counts.a)+1):nrow(amalgamated)),] <- counts.b
    write.table(amalgamated,file=args[3],sep="\t",quote=FALSE,row.names=FALSE)
  }
}