"""
usage:
    linkfiles_dri.py maps_deseqFile padj log2fc rilOutput_csv nbInteractionFilter intaRNAoutputFile energyFilter gffInput

    python3 linkfiles_dri.py ../Bureau/cDiff/maps_deseq2/results/rcd21vsctrl.complete.txt 0.05 6 ../Bureau/cDiff/rilSeq/ril_results/RCd21/RCd21vsCD630.csv 100 ../Bureau/cDiff/intaRNA/results/tssOr100/passFilter_list_rcd21.tsv -15 ../Bureau/cDiff/data/CD630_CDS-ncRNA.gff > ../printed_aaaa.txt

todo:
rewrite output path and file name for pipeline
"""
import sys 
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

mapsOutFile = sys.argv[1]
padjFilter = float(sys.argv[2])
l2fcFilter = float(sys.argv[3])
intaOutFile = sys.argv[4]
enrgFilter= float(sys.argv[5])
gffFileInput = sys.argv[6]
ncRNAid = sys.argv[7]

## import info about ncRNA

rcdLabelDict = {}

with open('labelDictFile.txt') as f:
    lines = [line.replace("\n","") for line in f]
for i in lines:
    rcdLabelDict[i.split(":")[0]]=i.replace("[","").replace("]","").split(":")[1].split(sep=",")

with open('cdiff-portal-genes-download-2024-05-16T05-06-25.csv') as f:
    lines = [line.replace("\n","") for line in f]
dictNames = {}
for i in lines:
    if i.startswith("\"Title")==False:
        split=i.split("\"")
        dictNames[split[1]]=split[-4]



## import gff file

gffFile = open(gffFileInput, 'r')
gffDict={}
for line in gffFile:
    if line[0]=="#":
        continue
    else: #if line.split("\t")[2]=="gene":
        lineList = line.split("\t")
        geneId = lineList[-1].split("=")[-1].replace("\"","").replace("\n","")
        gffDict[geneId] = [lineList[3],lineList[4]]


## import DESeq2 file
deseqData = pd.read_csv(mapsOutFile, sep="\t", usecols=['Id','log2FoldChange','pvalue','padj'])
deseqData = deseqData.drop(deseqData[deseqData.log2FoldChange < l2fcFilter].index)
deseqData = deseqData.drop(deseqData[deseqData.padj > padjFilter].index)



## import intaRNA file
intaData = pd.read_csv(intaOutFile, sep="\t")
intaData = intaData.drop(intaData[intaData["Interaction energy (kcal/mol)"] > enrgFilter].index)


# merge inputs, two columns: seqID and a color indicating source of seq
#       blue = ril , red = MAPS, yellow3 = inta
# mix colors if belong in multiple groups


allDF = pd.concat([deseqData["Id"], intaData["Gene name"]], axis=0, ignore_index=True, names="ID") # concatanate seqIDs from three inputs

unqDF = allDF.drop_duplicates(ignore_index=True) # remove duplicate seqID
unqDF = unqDF.to_frame() # convert series to pandas dataframe

unqDF=unqDF.rename(columns={0: "ID"}) # rename ID column

unqDF.insert(1,"color"," ") # add an empty column for color attribution
unqDF.insert(2,"value"," ") # add an empty column for value attribution



# color attribution
for j in range(len(unqDF)):
    for k in deseqData["Id"]:
        if unqDF["ID"][j]==k:
            unqDF.loc[j,"color"]="red"
            
            unqDF.loc[j,"value"]=1
            
    for k in intaData["Gene name"]:
        if unqDF["ID"][j]==k:
            if unqDF["color"][j]==" ":
                unqDF.loc[j,"color"]="blue"
                unqDF.loc[j,"value"]=1
            elif unqDF["color"][j]=="red":
                unqDF.loc[j,"color"]="purple"
                unqDF.loc[j,"value"]=2


# output writing
outdir = "linkFiles"

outputPosFile = open(f"{outdir}/{ncRNAid}_links.txt", 'w')  
outputRcdFile = open(f"{outdir}/{ncRNAid}_links_rcd.txt", 'w')
outputLabelFile = open(f"{outdir}/{ncRNAid}_label.txt", 'w')

#write header / first line
outputPosFile.writelines(f"chr\tstart\tend\tvalue\tcolor\n")
outputRcdFile.writelines(f"chr\tstart\tend\tcolor\n")
outputLabelFile.writelines(f"name\tposition\tlabel\tcolor\n")

# add label for ncRNA
outputLabelFile.writelines(f"{rcdLabelDict[ncRNAid.lower()][0]}\t{rcdLabelDict[ncRNAid.lower()][1]}\t{rcdLabelDict[ncRNAid.lower()][2]}\t{rcdLabelDict[ncRNAid.lower()][3]}\n")


for ind in unqDF.index:
    posMean = round((int(gffDict[unqDF['ID'][ind]][0])+int(gffDict[unqDF['ID'][ind]][1]))/2)
    outputPosFile.writelines(f"a\t{gffDict[unqDF['ID'][ind]][0]}\t{gffDict[unqDF['ID'][ind]][1]}\t{unqDF['value'][ind]}\t{unqDF['color'][ind]}\n")
    outputRcdFile.writelines(f"{rcdLabelDict[ncRNAid.lower()][4]}\t{rcdLabelDict[ncRNAid.lower()][5]}\t{rcdLabelDict[ncRNAid.lower()][6]}\t{rcdLabelDict[ncRNAid.lower()][7]}\n")
    
    # add gene name if exists in database
    if unqDF['ID'][ind] not in dictNames or dictNames[unqDF['ID'][ind]]=="":
        outputLabelFile.writelines(f"a\t{posMean}\t{unqDF['ID'][ind]}\t{unqDF['color'][ind]}\n")
    else:
        outputLabelFile.writelines(f"a\t{posMean}\t{unqDF['ID'][ind]} ({dictNames[unqDF['ID'][ind]]})\t{unqDF['color'][ind]}\n")

