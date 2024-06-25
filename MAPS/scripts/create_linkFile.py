'''
python3 ~/cdiffrib/circlize/scripts/create_linkFile_Deseq.py generated_files 0.05 6 data/RCd1vsempty.complete_2_replicates.txt data/CD630_MaGe_CTN.gff Ma2 RCd1

'''
import pandas as pd
import sys
import numpy as np
import random

outdir = sys.argv[1]
padjTh = float(sys.argv[2])
l2fcTh = sys.argv[3]
if sys.argv[3] == "none":
    l2fcTh = False
else:
    l2fcTh = float(sys.argv[3])
deseqFile = sys.argv[4] 
gffInput = sys.argv[5]
ncRnaName = sys.argv[6]
fileName = sys.argv[7]



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

'''
rcdLabelDict = {
    rcdX = [labelFile, rcdLink]
}
'''


### data import

deseqData = pd.read_csv(deseqFile, sep="\t", usecols=['Id','log2FoldChange','pvalue','padj']) # padj peut etre vide des fois donc prendre pvalue pour filtre ?
#print(deseqData)

if type(l2fcTh) == float:
    deseqData = deseqData[(deseqData.padj < padjTh) & (deseqData.log2FoldChange > l2fcTh)]
else:
    deseqData = deseqData[(deseqData.padj < padjTh)]

deseqDict = deseqData.set_index('Id').T.to_dict('list')


# color list

colorList = ["#00FF00", "#0000FF", "#FF00FF", "#800000", "#008000", "#000080", 
             "#808000", "#800080", "#008080", "#808080", "#C00000", "#00C000", "#0000C0", "#C0C000", "#C000C0", "#00C0C0",
             "#C0C0C0", "#400000", "#004000", "#000040", "#404000", "#400040", "#004040", "#404040", "#200000", "#002000",
             "#000020", "#202000", "#200020", "#002020", "#202020", "#600000", "#006000", "#000060", "#606000", "#600060", 
             "#006060", "#606060", "#A00000", "#00A000", "#0000A0", "#A0A000", "#A000A0", "#00A0A0", "#A0A0A0", "#E00000", 
             "#00E000", "#0000E0", "#E0E000", "#E000E0", "#00E0E0", "#E0E0E0"]

# dict for class:color 

colorDict = {}
count = 0

# import tss classification

with open('classificationTSSCOPY.tsv') as f:
    lines = [line.replace("\n","") for line in f]
for m in lines:
    lineSplit=m.split(";")
    if lineSplit[0] in deseqDict:
        if "/" in lineSplit[1]:
            splitTssList=lineSplit[1].split("/")
            for n in splitTssList:
                if "?" in n or "*" in n or "Sig" not in n:
                    pass
                else:
                    deseqDict[lineSplit[0]].append(n)
                    if n not in colorDict:
                        colorDict[n] = colorList[count]
                        count+=1

        else:
            if "?" in lineSplit[1] or "*" in lineSplit[1] or lineSplit[1]=="" or "Sig" not in lineSplit[1]:
                    pass
            else:
                deseqDict[lineSplit[0]].append(lineSplit[1])
                if lineSplit[1] not in colorDict:
                        colorDict[lineSplit[1]] = colorList[count]
                        count+=1

## import gff file

gffFile = open(gffInput, 'r')
gffDict={}
for line in gffFile:
    if line[0]=="#":
        continue
    else:#if line.split("\t")[2]=="gene":
        seqName = line.split("\t")[-1].split("\"")[1]
        gffDict[seqName] = [line.split("\t")[3],line.split("\t")[4]]





########################### writing linkfile without classification ###########################


outputPosFile = open(f"{outdir}/{fileName}_{ncRnaName}_links_deseq.txt", 'w')  # fichier pour tracer liens des gènes avec log2fc positif
outputRcdFile = open(f"{outdir}/{fileName}_{ncRnaName}_links_rcd_deseq.txt", 'w')
outputLabelFile = open(f"{outdir}/{fileName}_{ncRnaName}_label_deseq.txt", 'w')

#write header / first line
outputPosFile.writelines(f"chr\tsrart\tend\tvalue\tcolor\n")
outputRcdFile.writelines(f"chr\tsrart\tend\tcolor\n")

outputLabelFile.writelines(f"name\tposition\tlabel\tcolor\n")

tmpList = rcdLabelDict[ncRnaName.lower()]
outputLabelFile.writelines(f"{tmpList[0]}\t{tmpList[1]}\t{tmpList[2]}\tred\n")






if l2fcTh == False:
    for locus in deseqDict:
        posLocus = round((int(gffDict[locus][0])+int(gffDict[locus][1]))/2)
        if deseqDict[locus][0]>0:
            outputPosFile.writelines(f"a\t{gffDict[locus][0]}\t{gffDict[locus][1]}\t{deseqDict[locus][0]}\tblack\n")
            outputRcdFile.writelines(f"{tmpList[4]}\t{tmpList[5]}\t{tmpList[6]}\tblack\n")
            if locus not in dictNames or dictNames[locus]=="":
                outputLabelFile.writelines(f"a\t{posLocus}\t{locus}\tblack\n")
            else:
                outputLabelFile.writelines(f"a\t{posLocus}\t{locus} ({dictNames[locus]})\tblack\n")
        if deseqDict[locus][0]<0:
            outputPosFile.writelines(f"a\t{gffDict[locus][0]}\t{gffDict[locus][1]}\t{deseqDict[locus][0]}\tblue\n")
            outputRcdFile.writelines(f"{tmpList[4]}\t{tmpList[5]}\t{tmpList[6]}\tblue\n")
            if locus not in dictNames or dictNames[locus]=="":
                outputLabelFile.writelines(f"a\t{posLocus}\t{locus}\tblue\n")
            else:
                outputLabelFile.writelines(f"a\t{posLocus}\t{locus} ({dictNames[locus]})\tblue\n")
                    


else:
    for locus in deseqDict:
        posLocus = round((int(gffDict[locus][0])+int(gffDict[locus][1]))/2)
        if deseqDict[locus][0]>l2fcTh:
            outputPosFile.writelines(f"a\t{gffDict[locus][0]}\t{gffDict[locus][1]}\t{deseqDict[locus][0]}\tblack\n")
            outputRcdFile.writelines(f"{tmpList[4]}\t{tmpList[5]}\t{tmpList[6]}\tblack\n")
            if locus not in dictNames or dictNames[locus]=="":
                outputLabelFile.writelines(f"a\t{posLocus}\t{locus}\tblack\n")
            else:
                outputLabelFile.writelines(f"a\t{posLocus}\t{locus} ({dictNames[locus]})\tblack\n")


                    
outputRcdFile.close()
outputPosFile.close()
outputLabelFile.close()

###################################################################################################








############################## writing linkfile with classification ##############################


outputClassPosFile = open(f"{outdir}/{fileName}_{ncRnaName}_links_classDeseq.txt", 'w')  # fichier pour tracer liens des gènes avec log2fc positif
outputClassRcdFile = open(f"{outdir}/{fileName}_{ncRnaName}_links_rcd_classDeseq.txt", 'w')
outputClassLabelFile = open(f"{outdir}/{fileName}_{ncRnaName}_label_classDeseq.txt", 'w')




#write header / first line
outputClassPosFile.writelines(f"chr\tsrart\tend\tvalue\tcolor\n")
outputClassRcdFile.writelines(f"chr\tsrart\tend\tcolor\n")

outputClassLabelFile.writelines(f"name\tposition\tlabel\tcolor\n")

tmpListClass = rcdLabelDict[ncRnaName.lower()]
outputClassLabelFile.writelines(f"{tmpList[0]}\t{tmpList[1]}\t{tmpList[2]}\tblack\n")



listExistingColor = []
listMulti = []

if l2fcTh == False:
    for locus in deseqDict:
        posLocus = round((int(gffDict[locus][0])+int(gffDict[locus][1]))/2)
        if len(deseqDict[locus])==3:
            colorCode="black"
        elif len(deseqDict[locus])==4:
            colorCode=colorDict[deseqDict[locus][-1]]
            if deseqDict[locus][-1] not in listExistingColor:
                listExistingColor.append(deseqDict[locus][-1])
        elif len(deseqDict[locus])>4:
            colorCode="#FF0000"
            multi=True
            listMulti.append(locus)
        outputClassPosFile.writelines(f"a\t{gffDict[locus][0]}\t{gffDict[locus][1]}\t{deseqDict[locus][0]}\t{colorCode}\n")
        outputClassRcdFile.writelines(f"{tmpList[4]}\t{tmpList[5]}\t{tmpList[6]}\tblack\n")
        if locus not in dictNames or dictNames[locus]=="":
            outputClassLabelFile.writelines(f"a\t{posLocus}\t{locus}\t{colorCode}\n")
        else:
            outputClassLabelFile.writelines(f"a\t{posLocus}\t{locus} ({dictNames[locus]})\t{colorCode}\n")
        




else:
    for locus in deseqDict:
        posLocus = round((int(gffDict[locus][0])+int(gffDict[locus][1]))/2)
        if deseqDict[locus][0]>l2fcTh:
            if len(deseqDict[locus])==3:
                colorCode="black"
            elif len(deseqDict[locus])==4:
                colorCode=colorDict[deseqDict[locus][-1]]
                if deseqDict[locus][-1] not in listExistingColor:
                    listExistingColor.append(deseqDict[locus][-1])
            elif len(deseqDict[locus])>4:
                colorCode="#FF0000"
                multi=True
                listMulti.append(locus)
            outputClassPosFile.writelines(f"a\t{gffDict[locus][0]}\t{gffDict[locus][1]}\t{deseqDict[locus][0]}\t{colorCode}\n")
            outputClassRcdFile.writelines(f"{tmpList[4]}\t{tmpList[5]}\t{tmpList[6]}\tblack\n")
            if locus not in dictNames or dictNames[locus]=="":
                outputClassLabelFile.writelines(f"a\t{posLocus}\t{locus}\t{colorCode}\n")
            else:
                outputClassLabelFile.writelines(f"a\t{posLocus}\t{locus} ({dictNames[locus]})\t{colorCode}\n")


                    
outputClassRcdFile.close()
outputClassPosFile.close()
outputClassLabelFile.close()



classificationColorFile = open(f"{outdir}/{fileName}_{ncRnaName}_classification.txt", 'w')
classificationColorFile.writelines('className,color\n')
if 'multi' in globals():
    classificationColorFile.writelines('multiple_classes,#FF0000\n')
for className in listExistingColor:
    classificationColorFile.writelines(f'{className},{colorDict[className]}\n')
classificationColorFile.writelines('unclassified,black\n')

classificationColorFile.close()


if len(listMulti)>0:
    multiFile = open(f"{outdir}/{fileName}_{ncRnaName}_listMulti.csv", 'w')
    multiFile.writelines('id,TSS\n')
    for multi in listMulti:
        onlytxt=deseqDict[multi][3]
        for v in deseqDict[multi][4:]:
            onlytxt=onlytxt+" "+v
        multiFile.writelines(f'{multi},{onlytxt}\n')

    multiFile.close()

###################################################################################################


