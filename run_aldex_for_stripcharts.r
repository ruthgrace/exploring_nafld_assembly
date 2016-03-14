# runs ALDEx for stripcharts

options(error=recover)

library(ALDEx2)
source("AitchisonTransform.r")

sampleindexes <- c(3:22)

lengthindex <- 2

# outfolder <- "subsys4_counts"
# outfolder <- "carbohydrates"
 outfolder <- "lipids"

nash <- c("CL_166_BL_R2_R1", "CL_169_BL_R1", "CL_139_BL_2_R1", "CL_173_2_R1", "CL_144_2_R1", "CL_177_R1", "CL_160_R1", "CL_165_R1", "CL_119_R1", "CL_141_BL_R2_R1")

healthy <- c("HLD_100_R1", "HLD_102_R1", "HLD_111_2_R1", "HLD_80_R1", "HLD_85_R1", "HLD_28_R1", "HLD_47_R1", "HLD_72_2_R1", "HLD_112_R1", "HLD_23_R1")

# read in fully annotated counts with hierarchies
# d <- read.table("annotated_counts_with_refseq_length_non_zero_features.txt",sep="\t",header=TRUE,row.names=NULL,quote="",comment.char="")
# d <- read.table("annotated_carbohydrate_counts_with_refseq_length.txt",sep="\t",header=TRUE,row.names=NULL,quote="",comment.char="")
 d <- read.table("annotated_lipid_counts_with_refseq_length.txt",sep="\t",header=TRUE,row.names=NULL,quote="",comment.char="")

originaldata <- d

d <- d[which(apply(d[,sampleindexes],1,sum)!=0),]

#put subsystems in hierarchies with "|||" as separator
subsys <- paste(d$subsys4, d$subsys1, d$subsys2, d$subsys3, sep="|||")

subsys.unique <- unique(subsys)

subsys.nums <- c(1:length(subsys.unique))

d.samples.mat <- as.matrix(d[,sampleindexes])

aggregate.subsys.counts <- function(x,d.samples.mat,subsys,subsys.unique) {
  indices <- which(subsys==subsys.unique[x])
  if (length(indices) == 1) {
    return(d.samples.mat[which(subsys==subsys.unique[x]),])
  }
  else {
    return(apply(d.samples.mat[which(subsys==subsys.unique[x]),],2,sum))
  }
}

d.aggregate <- sapply(subsys.nums,function(x) { return(aggregate.subsys.counts(x,d.samples.mat,subsys,subsys.unique)) } )

d.aggregate <- t(d.aggregate)

rownames(d.aggregate) <- subsys.unique

#get mean lengths per unique subsys
aggregate.subsys.lengths <- function(x,d.lengths,subsys,subsys.unique) {
  indices <- which(subsys==subsys.unique[x])
  if (length(indices) == 1) {
    return(d.lengths[which(subsys==subsys.unique[x])])
  }
  else {
    return(mean(d.lengths[which(subsys==subsys.unique[x])]))
  }
}

d.aggregate.lengths <- sapply(subsys.nums,function(x) { return(aggregate.subsys.lengths(x,d[,lengthindex],subsys,subsys.unique)) } )

# columns are refseq, length, counts for all samples, group, but here we'll just use subsys for both group and refseq
d.transform.in <- data.frame(matrix(NA,nrow=length(subsys.unique),ncol=(length(sampleindexes)+3)))
d.transform.in[,1] <- subsys.unique
d.transform.in[,ncol(d.transform.in)] <- subsys.unique
d.transform.in[,2] <- d.aggregate.lengths
d.transform.in[,3:(3+length(sampleindexes)-1)] <- d.aggregate

colnames(d.transform.in) <- c("subsys","length",colnames(d[,sampleindexes]),"subsys")

write.table(d.transform.in,file=paste(outfolder,"AitchisonTransform_input_for_stripcharts_merged_subsys.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)

d.transformed <- aitchison.transform.reads(filename=paste(outfolder,"AitchisonTransform_input_for_stripcharts_merged_subsys.txt",sep="/"),rounded=TRUE, subjects = 20, firstsubjectindex = 3, lastsubjectindex = 22, groupindex = 23,lengthindex=2,outputfolder="subsys4_counts")

### for some reason the aitchison transform method sometimes refuses to run unless you manually run it line by line. here are the relevant input params from the previous line
# filename = paste(outfolder,"AitchisonTransform_input_for_stripcharts_merged_subsys.txt",sep="/")
# rounded = TRUE
# subjects = 20
# firstsubjectindex = 3
# lastsubjectindex = 22
# groupindex = 23
# lengthindex=2
# outputfolder=outfolder

### the aldex input is already output by aitchison.transform function
# write.table(d.aggregate,file=paste(outfolder, "ALDEx_input_for_stripcharts_merged_subsys.txt",sep="/"),,sep="\t",quote=FALSE)

aldex.data <- d.transformed

aldex.data <- data.frame(aldex.data,check.names=FALSE)
rownames(aldex.data) <- aldex.data$subsys
dropColumns <- c("indices","subsys")
aldex.data <- aldex.data[,which(!(colnames(aldex.data) %in% dropColumns))]

conditions <- colnames(aldex.data)
conditions[which(conditions %in% nash)] <- "nash"
conditions[which(conditions %in% healthy)] <- "healthy"

x <- aldex(aldex.data, conditions, mc.samples=128)

write.table(x,file=paste(outfolder, "ALDEx_output_for_stripcharts_merged_subsys.txt",sep="/"),,sep="\t",quote=FALSE)

x.separate.subsys <- data.frame(matrix(NA,nrow=nrow(x),ncol=(ncol(x)+4)))

x.separate.subsys[,c(1:4)] <- t(data.frame(strsplit(rownames(x),split="|||",fixed=TRUE)))

x.separate.subsys[,c(5:ncol(x.separate.subsys))] <- as.matrix(x)

colnames(x.separate.subsys) <- c("subsys4","subsys1","subsys2","subsys3",colnames(x))

write.table(x.separate.subsys,file=paste(outfolder, "ALDEx_output_for_stripcharts.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)

x.ordered <- x.separate.subsys[order(abs(x.separate.subsys$effect),decreasing=TRUE),]

write.table(x.ordered,file=paste(outfolder, "ALDEx_output_for_stripcharts_ordered_by_effect.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)

pdf(paste(outfolder,"ALDEx_all_hierarchies_output.pdf",sep="/"))

aldex.plot(x)

dev.off()
