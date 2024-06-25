

args = commandArgs(trailingOnly=TRUE)   # get values for legend

# 1 = nom des fichiers link 
# 2 = working directory



##################### get arguments ####################

inputName = args[1]
wdPath = args[2]
padj = args[3]
l2fc = args[4]
intaEng = args[5]


########################################################



############### import libraries ###############


library(circlize)



################################################





############### define working directory ###############


setwd(wdPath)
#print(getwd())

########################################################






############### import data ###############


legendList = c('MAPS','IntaRNA','MAPS & IntaRNA')
colorList = c('red','blue','purple')



############################################


############# circos function ##############

circos_function <- function(labelDf, classDf, rcdDf){
  # initialize circos
  circos.par(gap.degree = 0,
    cell.padding = c(0, 0, 0, 0),
	  start.degree = 90,
	  "canvas.xlim" = c(-1.25, 1.25), "canvas.ylim" = c(-1.25, 1.25) # modisy size of ideogram
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





############## rilSeq circos with classification ################

classInput = paste0("linkFiles/",inputName,"_links.txt")
rcdInput = paste0("linkFiles/",inputName,"_links_rcd.txt")
labelInput = paste0("linkFiles/",inputName,"_label.txt")

dataClass = read.csv(classInput ,sep="\t")
dataRcd = read.csv(rcdInput ,sep="\t")
dataLabel = read.csv(labelInput ,sep="\t")

circosName = paste0("plots/",inputName,"_circosR.png")
png(circosName, width = 9100, height = 7000)#11000, height = 9000)


circos_function(dataLabel,dataClass,dataRcd)


### add legend
lgd_points = legend("bottomright", inset=.02, title="classes",
		    legendList, fill=colorList, cex=8)


### add title
titileTxt = paste0(inputName, " - circos plot (MAPS, IntaRNA)")
indivTitle = title(main = titileTxt, cex.main=20, line = -15)

nbMaps = sum(dataLabel$color == "red") + sum(dataLabel$color == "purple")
nbInta = sum(dataLabel$color == "blue") + sum(dataLabel$color == "purple")

parameterTxt =  paste0("MAPS (",nbMaps,") : log2FC >= ", l2fc, ", padj <= ", padj, "\nIntaRNA (",nbInta,") : interaction energy <= ", intaEng, " kcal/mol")

mtext(parameterTxt, side=1, cex = 8, line = -9, at = -1.5, adj = 0)

circos.clear()

dev.off()