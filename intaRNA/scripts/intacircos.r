#
#


args = commandArgs(trailingOnly=TRUE)   # get values for legend

# args[1] = interaction energy threshold
# args[2] = working directory
# args[3] = ncRNA ID






############### import libraries ###############

library(circlize)

################################################





############### define working directory ###############


setwd(args[2])
#print(getwd())

########################################################





############### import data ###############




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





############## deseq circos ################

classInput = paste0(args[3],"/linkFiles/",args[3],"_links_inta.txt")
rcdInput = paste0(args[3],"/linkFiles/",args[3],"_links_rcd_inta.txt")
labelInput = paste0(args[3],"/linkFiles/",args[3],"_label_inta.txt")

dataClass = read.csv(classInput ,sep="\t")
dataRcd = read.csv(rcdInput ,sep="\t")
dataLabel = read.csv(labelInput ,sep="\t")

deseqName = paste0("plots/",args[3],"_circosR_inta.png")
png(deseqName, width = 9100, height = 6300)#11000, height = 9000)


circos_function(dataLabel,dataClass,dataRcd)

### add legend
lgd_points = legend("bottomright", inset=.02, title="targets", c("sequence id"), fill=c("black"), cex=10)


### add title
titleTxt = paste0(args[3], " intaRNA circos")
indivTitle = title(main = titleTxt, cex.main=20, line = -15)

parameterTxt = paste0("interaction energy < ", args[1],"kcal/mol\ntotal: ",length(dataClass$chr)," sequences")

mtext(parameterTxt,
      side=1,
      cex = 8,       # text size
      line = -9,     # move up
      at = -1.5)        # move left

circos.clear()

dev.off()





############################################