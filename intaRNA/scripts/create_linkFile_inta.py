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


intaOutFile = sys.argv[1]
gffFileInput = sys.argv[2]
ncRNAid = sys.argv[3]
labelFile = sys.argv[4]
symbFile = sys.argv[5]
## import info about ncRNA

rcdLabelDict = {}

with open(labelFile) as f: # '../../cdiffrib/circlize/scripts/labelDictFile.txt'
    lines = [line.replace("\n","") for line in f]
for i in lines:
    rcdLabelDict[i.split(":")[0]]=i.replace("[","").replace("]","").split(":")[1].split(sep=",")

with open(symbFile) as f: #'../../cdiffrib/cdiff-portal-genes-download-2024-05-16T05-06-25.csv'
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
        geneId = lineList[-1].split("=")[-1].replace("\"","").replace("\n","").strip()
        gffDict[geneId] = [lineList[3],lineList[4]]


## import intaRNA file
intaData = pd.read_csv(intaOutFile, sep="\t")



# output writing

outputPosFile = open(f"{ncRNAid}/linkFiles/{ncRNAid}_links_inta.txt", 'w')  
outputRcdFile = open(f"{ncRNAid}/linkFiles/{ncRNAid}_links_rcd_inta.txt", 'w')
outputLabelFile = open(f"{ncRNAid}/linkFiles/{ncRNAid}_label_inta.txt", 'w')

#write header / first line
outputPosFile.writelines(f"chr\tstart\tend\tvalue\tcolor\n")
outputRcdFile.writelines(f"chr\tstart\tend\tcolor\n")
outputLabelFile.writelines(f"name\tposition\tlabel\tcolor\n")

ncrnaIdx = "cd630_"+ncRNAid.lower()

# add label for ncRNA
outputLabelFile.writelines(f"{rcdLabelDict[ncrnaIdx][0]}\t{rcdLabelDict[ncrnaIdx][1]}\t{rcdLabelDict[ncrnaIdx][2]}\tred\n")


for index, column in intaData.iterrows():
    posMean = round((int(gffDict[column['Gene name']][0])+int(gffDict[column['Gene name']][1]))/2)
    outputPosFile.writelines(f"a\t{gffDict[column['Gene name']][0]}\t{gffDict[column['Gene name']][1]}\t{column['Interaction energy (kcal/mol)']}\tblack\n")
    outputRcdFile.writelines(f"{rcdLabelDict[ncrnaIdx][4]}\t{rcdLabelDict[ncrnaIdx][5]}\t{rcdLabelDict[ncrnaIdx][6]}\t{rcdLabelDict[ncrnaIdx][7]}\n")
    
    # add gene name if exists in database
    if column['Gene name'] not in dictNames or dictNames[column['Gene name']]=="":
        outputLabelFile.writelines(f"a\t{posMean}\t{column['Gene name']}\tblack\n")
    else:
        outputLabelFile.writelines(f"a\t{posMean}\t{column['Gene name']} ({dictNames[column['Gene name']]})\tblack\n")

