"""
Usage:
python gffFromInta.py ~/Desktop/22_04_transferDir/translateList.txt small_dataset rcd21 plusPositions.txt minusPositions.txt

"""

import sys

# input = intaRNA detailed output
translateFile=sys.argv[1]

inputIntaData=sys.argv[2]

# float, filter for interaction energy
suffix=sys.argv[3]

plusPosFile=sys.argv[4]
minusPosFile=sys.argv[5]

outputDict={}

translateList = open(translateFile, 'r')
for line in translateList:
    splittedLine=line.split("\t")
    outputDict[splittedLine[1].replace(">","").replace("\n","").replace(" ","")]=[int(splittedLine[0].split("-")[0].replace(">","")),int(splittedLine[0].split("-")[1])]

'''
outputdict = {
    geneID: [start, end]
}
'''

    

plusPos = open(plusPosFile, 'r')
for linePos in plusPos:
    for tl in outputDict:
        if outputDict[tl][0] == int(linePos.split(":")[1].split("-")[0]) and outputDict[tl][1] == int(linePos.split(":")[1].split("-")[1]) and outputDict[tl][1] and len(outputDict[tl])==2 :
            outputDict[tl].append("+")



minusPos = open(minusPosFile, 'r')
for lineMin in minusPos:
    for tl in outputDict:
        if outputDict[tl][0] == int(lineMin.split(":")[1].split("-")[0]) and outputDict[tl][1] == int(lineMin.split(":")[1].split("-")[1]) and outputDict[tl][1] and len(outputDict[tl])==2 :
            outputDict[tl].append("-")


'''
outputdict = {
    geneID: [start, end, strand(+/-)]
}
'''

'''
for i in outputDict:
    print(i,outputDict[i])
'''




linesIntaOut = iter(open(inputIntaData, 'r'))



linesSaveList=[] # where lines of input file are stored
outputLinesDict=[]

for line in linesIntaOut:

    linesSaveList.append(line)
    
    if line.startswith("interaction energy") and len(outputDict[linesSaveList[1].replace("\n","")])==3:

        #print(linesSaveList)
        seqName=linesSaveList[1].replace("\n","")
        
        # append interaction energy
        outputDict[seqName].append(float(line.split(" = ")[1].split(" ")[0]))
        
        posLine=linesSaveList[2].split(" ")
        check=False
        for i in posLine:
            
            if i!="":
                if outputDict[seqName][2]=="+":
                    if check==False:
                        outputDict[seqName][1]=outputDict[seqName][0]
                        outputDict[seqName][0]+=int(i)-1
                        check=True
                    elif check==True:
                        outputDict[seqName][1]+=int(i)-1
                elif outputDict[seqName][2]=="-":
                    if check==False:
                        outputDict[seqName][0]=outputDict[seqName][1]
                        outputDict[seqName][1]-=int(i)-1
                        check=True
                    elif check==True:
                        outputDict[seqName][0]-=int(i)-1


    elif line.startswith("seed Pu2"):
        linesSaveList=[]
    elif line.startswith("no favorable interaction for"):
        linesSaveList=[]


'''
outputDict = {

gene_name : [start position of interactiton, end position of interactiton, strand, interaction energy]

}

'''

# output Gff file; can add filter for interaction energy here:


outGffFile = open(f"{suffix}/linkPosition_{suffix}.gff", 'w')
for item in outputDict:
    #print(item,len(outputDict[item]))
    if len(outputDict[item])==4:
        outGffFile.writelines(f"NC_009089\t{suffix}\ttargeted region\t{outputDict[item][0]}\t{outputDict[item][1]}\t{outputDict[item][3]}\t{outputDict[item][2]}\t.\ttarget_gene_name={item}\n")

