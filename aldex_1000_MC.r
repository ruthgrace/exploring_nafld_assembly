#!/usr/bin/env Rscript

nash <- c("CL_166_BL_R2_R1", "CL_169_BL_R1", "CL_139_BL_2_R1", "CL_173_2_R1", "CL_144_2_R1", "CL_177_R1", "CL_160_R1", "CL_165_R1", "CL_119_R1", "CL_141_BL_R2_R1")
healthy <- c("HLD_100_R1", "HLD_102_R1", "HLD_111_2_R1", "HLD_80_R1", "HLD_85_R1", "HLD_28_R1", "HLD_47_R1", "HLD_72_2_R1", "HLD_112_R1", "HLD_23_R1")

aldex.data <- read.table("/Volumes/data/ruth/nafld_assembly/assembly_test_blast/subsys4_counts/AitchisonTransform_input_for_stripcharts_merged_subsys_AitchisonTransformedDataForALDExInput.txt",sep="\t",row.names=1,quote="",comment.char="",header=TRUE)

conditions <- colnames(aldex.data)
conditions[which(conditions %in% nash)] <- "nash"
conditions[which(conditions %in% healthy)] <- "healthy"

x.1000 <- aldex(aldex.data, conditions, mc.samples=1000)

write.table(x.1000,file=paste(outfolder, "ALDEx_output_for_stripcharts_merged_subsys_1000_samples.txt",sep="/"),,sep="\t",quote=FALSE)

x.1000.separate.subsys <- data.frame(matrix(NA,nrow=nrow(x.1000),ncol=(ncol(x.1000)+4)))

x.1000.separate.subsys[,c(1:4)] <- t(data.frame(strsplit(rownames(x.1000),split="|||",fixed=TRUE)))

x.1000.separate.subsys[,c(5:ncol(x.1000.separate.subsys))] <- as.matrix(x.1000)

colnames(x.1000.separate.subsys) <- c("subsys4","subsys1","subsys2","subsys3",colnames(x.1000))

write.table(x.1000.separate.subsys,file=paste(outfolder, "ALDEx_output_for_stripcharts_1000_samples.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)

x.1000.ordered <- x.1000.separate.subsys[order(abs(x.1000.separate.subsys$effect),decreasing=TRUE),]

write.table(x.1000.ordered,file=paste(outfolder, "ALDEx_output_for_stripcharts_ordered_by_effect_1000_samples.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)

pdf(paste(outfolder,"ALDEx_all_hierarchies_output_1000_samples.pdf",sep="/"))

aldex.plot(x.1000)

dev.off()
