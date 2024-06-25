"""
Usage:
python parseIntaOut.py small_dataset -11 rcd21

Output:
    - passFilter_Complete_rcd21.txt
    - passFilter_list_rcd21.tsv
    - passFilter_OnlyLink_rcd21.txt

"""

import sys

# input = intaRNA detailed output
inputIntaData=sys.argv[1]
# float, filter for interaction energy
ef=sys.argv[2]
suffix=sys.argv[3]

# check if right format of float
if "," in ef:
    ef.replace(",",".")
ef=float(ef)


lines = iter(open(inputIntaData, 'r'))


completeOutFile = open(f"{suffix}/passFilter_Complete_{suffix}.txt", 'w')
linkOutFile = open(f"{suffix}/passFilter_OnlyLink_{suffix}.txt", 'w')
listTarget = open(f"{suffix}/passFilter_list_{suffix}.tsv", 'w')
listTarget.writelines("Gene name\tInteraction energy (kcal/mol)\n")

linesSaveList=[] # where lines of input file are stored

for line in lines:

    linesSaveList.append(line)

    if line.startswith("interaction energy"):
        if float(line.split(" = ")[1].split(" ")[0])<=ef:
            checker=True
            energy=linesSaveList[-1].split(' = ')[1].split(" ")[0]
            geneName=linesSaveList[1].replace("\n","")
            listTarget.writelines(f"{geneName}\t{energy}\n")
            for i in range(len(linesSaveList)-4):
                linkOutFile.writelines(linesSaveList[i])
            linkOutFile.writelines("______________________________________________________________________________\n")
        else:
            checker=False
    
    elif line.startswith("seed Pu2"):
        if checker==True:
            for i in linesSaveList:
                completeOutFile.writelines(i)
            linesSaveList=[]
        elif checker==False:
            linesSaveList=[]
            
    elif line.startswith("no favorable interaction for"):
        linesSaveList=[]

