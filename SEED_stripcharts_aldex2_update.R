#Modifed Jul 25, 2012 - JM
# CAN source() this

#Update Feb 26 for fixed subsys4 at e-6
#Update Mar 19, 2013 for B. cereus data
#Update 2-Jul-2014: For the new ALDEx2 format

#------------------------------------------------------------------------------------
# New data format - A line for each UNIQUE SEED hierarchy
#------------------------------------------------------------------------------------

make_stripchart <- function(inputfile,outputfolder=".") {

 ### headers for this file are
 # [1] "subsys4"         "subsys1"         "subsys2"         "subsys3"        
 #  [5] "rab.all"         "rab.win.healthy" "rab.win.nash"    "diff.btw"       
 #  [9] "diff.win"        "effect"          "overlap"         "we.ep"          
 # [13] "we.eBH"          "wi.ep"           "wi.eBH"        
 	d <- read.table(file=inputfile, header=TRUE, sep="\t", quote="",comment.char="")	

	#color for non-DE subsys4
	base_col="#00000050"

	#Define your columns of interest
	col<-"we.eBH"
	#diff<-"diff.btw"
	diff<-"effect"
	#d[,col]

	#make the x-axis symmetrical based on the highest value
	xlim <- c(- max(abs(d[,diff])) -0.5, max(abs(d[,diff])) + 0.5)

	#----------------

	#put subsys1, subsys2, subsys3, subsys4 in a list
	levels <- paste( "subsys", 1:4, sep="" )

	#assuming only subsys1 to subsys3
	for (i in 1:3){
		
	#e.g. levels[1] sets g to "subsys1"
	#g<-"subsys1"
		g<-levels[i]
		
	#get all the unique functions for the current level. [[]] removes quotes around "subsys1"
		groups<-unique(d[[g]])

		ylim<-c(length(groups) - (length(groups)+0.5), length(groups) + 0.5)


	#For ALDEx2.0.2.7
		cutoffb <- ((d[,diff] < 0) & (d[,col] <= 0.05))
		cutoffn <- ((d[,diff] > 0) & (d[,col] <= 0.05))
		nocut <- (d[,col] > 0.05)

		bv_sig <- unique( data.frame( group=d[[g]], absolute=d[,diff])[ cutoffb , ] )
		n_sig <- unique( data.frame( group=d[[g]], absolute=d[,diff])[ cutoffn , ] )
		no_sig <- unique( data.frame( group=d[[g]], absolute=d[,diff])[nocut,])
		
	pdf(file=paste(outputfolder,"/","subsys",i,".pdf",sep=""),width=7,height=max(7,length(groups)/4))
	#	png(file=paste("subsys",i,".png",sep=""), width=9, height=(length(groups) / 5), units="in", res=300)
	#height=(length(groups) / 5),
		
		par(mar=c(5,20,0.5,0.5), cex.axis=0.8, cex.lab=0.6, las=1)
	#c(bottom, left, top, right)
		stripchart(absolute ~ group, data=no_sig, method="jitter", jitter=0.25, pch=20, col=base_col, xlim=xlim, cex=0.8, cex.axis=0.7)
		
	#Use add=TRUE to overplot
		stripchart(absolute ~ group, data=n_sig, method="jitter", jitter=0.25, pch=20, col="blue", xlim=xlim, cex=0.8, add=TRUE, cex.axis=0.7)
		stripchart(absolute ~ group, data=bv_sig, method="jitter", jitter=0.25, pch=20, col="red", xlim=xlim, cex=0.8, add=TRUE, cex.axis=0.7)

		abline(v=0, col="black", lty=2)
	#	title(xlab="Median Absolute Difference" ~~Log[2]~~, line=2, cex.lab=0.8)
		title(xlab="Median Effect Size" ~~(Log[2]), line=2, cex.lab=0.8)
		
		#draw horizonal lines
		for (j in 0.5:(length(groups)+0.5)){
			abline(h=j, lty=3, col="grey70")
		}
		
		dev.off()
		
	}
}
	#-----------
	# END new data
	#-----------

make_stripchart("subsys4_counts/ALDEx_output_for_stripcharts.txt","subsys4_counts")
# make_stripchart("carbohydrates/ALDEx_output_for_stripcharts.txt","carbohydrates")
# make_stripchart("lipids/ALDEx_output_for_stripcharts.txt","lipids")
# make_stripchart("subsys4_counts/ALDEx_output_for_stripcharts_ordered_by_effect_carbohydrates_only.txt", "subsys4_counts/carbohydrates")
# make_stripchart("subsys4_counts/ALDEx_output_for_stripcharts_ordered_by_effect_lipids_only.txt", "subsys4_counts/lipids")