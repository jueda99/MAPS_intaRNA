#
#


args = commandArgs(trailingOnly=TRUE)   # get values for legend

# args[1] = padj threshold
# args[2] = l2fc threshold
# args[3] = working directory
# args[4] = ncRNA ID
# args[5] = classification code





############### import libraries ###############

library(circlize)

################################################





############### define working directory ###############


setwd(args[3])
#print(getwd())

########################################################





############### import data ###############

colorListFile = paste0("linkFiles/",args[5], "_",args[4],"_classification.txt")
colorList = read.csv(colorListFile)

############################################




############# circos function ##############

circos_function <- function(labelDf, classDf, rcdDf){
  # initialize circos
  circos.par(gap.degree = 0,
    cell.padding = c(0, 0, 0, 0),
	  start.degree = 90,
	  "canvas.xlim" = c(-1.38, 1.38), "canvas.ylim" = c(-1.38, 1.38) # modisy size of ideogram
    )
  circos.initialize("a", xlim = c(0, 4290252))

  ### create label layer
  circos.labels(labelDf$name,
	  labelDf$position,
	  labelDf$label,
	  side = "outside",
	  cex = 6.5,
	  line_lwd = 4,
	  connection_height = mm_h(40),
	  col = labelDf$color,
	  line_col = labelDf$color
    )

  circos.track(ylim = c(0, 1), bg.border = NA)

  circos.axis(labels.cex = 4,
    major.tick.length = mm_y(12),
	  lwd = 5,
	  direction = "inside",
	  col = "gray28",
	  labels.col = "gray67"
    )
    
    
  circos.genomicLink(classDf,
		rcdDf,
		rou = 0.905,
		border = NA,
		col = classDf$color
    )

}

############################################







############## circos with classification ################

classInputColor = paste0("linkFiles/",args[5],"_",args[4],"_links_classDeseq.txt")
rcdInputColor = paste0("linkFiles/",args[5], "_",args[4],"_links_rcd_classDeseq.txt")
labelInputColor = paste0("linkFiles/",args[5], "_",args[4],"_label_classDeseq.txt")

dataClassColor = read.csv(classInputColor ,sep="\t")
dataRcdColor = read.csv(rcdInputColor ,sep="\t")
dataLabelColor = read.csv(labelInputColor ,sep="\t")

rilName = paste0("plots/",args[5],"_",args[4],"_classCircosR.png")
png(rilName, width = 9100, height = 6300)#11000, height = 9000)

circos_function(dataLabelColor,dataClassColor,dataRcdColor)

### add legend

lgd_points = legend("bottomright", inset=.02, title="classes",
		    colorList$className, fill=colorList$color, cex=8)

### add title
titleTxt = paste0(args[4], " DESeq2 circos (", args[5],")")
indivTitle = title(main = titleTxt, cex.main=20, line = -15)

parameterTxt = paste0("padj < ", args[1],"\nlog2fc >", args[2], "\ntotal: ",length(dataClassColor$chr)," sequences")

mtext(parameterTxt,
      side=1,
      cex = 8,       # text size
      line = -9,     # move up
      at = -1.5)        # move left

circos.clear()

dev.off()


############################################################







############## deseq circos ################

classInput = paste0("linkFiles/",args[5], "_",args[4],"_links_deseq.txt")
rcdInput = paste0("linkFiles/",args[5], "_",args[4],"_links_rcd_deseq.txt")
labelInput = paste0("linkFiles/",args[5], "_",args[4],"_label_deseq.txt")

dataClass = read.csv(classInput ,sep="\t")
dataRcd = read.csv(rcdInput ,sep="\t")
dataLabel = read.csv(labelInput ,sep="\t")

deseqName = paste0("plots/",args[4],"_circosR_deseq.png")
png(deseqName, width = 9100, height = 6300)#11000, height = 9000)


circos_function(dataLabel,dataClass,dataRcd)

### add legend
if (length(unique(dataClass$color))==1){
	lgd_points = legend("bottomright", inset=.02, title="log2fc", c("positive"), fill=c("black"), cex=10)
}else if (length(unique(dataClass$color))==2){
	lgd_points = legend("bottomright", inset=.02, title="log2fc", c("positive","negative"), fill=c("black","blue"), cex=10)
}

### add title
titleTxt = paste0(args[4], " DESeq2 circos")
indivTitle = title(main = titleTxt, cex.main=20, line = -15)

parameterTxt = paste0("padj < ", args[1],"\nlog2fc >", args[2], "\ntotal: ",length(dataClass$chr)," sequences")

mtext(parameterTxt,
      side=1,
      cex = 8,       # text size
      line = -9,     # move up
      at = -1.5)        # move left

circos.clear()

dev.off()





############################################