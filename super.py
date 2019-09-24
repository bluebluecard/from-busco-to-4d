import sys
import numpy as np

fileSpe = sys.argv[1]

filePath = sys.argv[2]

with open(fileSpe,'r') as f:
    speSeqList = []
    speIndexDict = {}
    count = 0
    for line in f:
        tmp = [line.strip(),'']
        speSeqList.append(tmp)
        speIndexDict[line.strip()] = count
        count += 1

with open(filePath,'r') as fh:
    fileContent = fh.read()

data = fileContent.split('>')
for i in data[1:]:
    info = i.split()
    IndexList = list(speIndexDict.values())
    for j in range(1,len(info),2):
        #spe = '_'.join(info[j].split('_')[1:])
        spe = info[j].split('_')[0]
        speIndex = speIndexDict[spe]
        speSeqList[speIndex][1] += info[j+1]
        IndexList.remove(speIndex)
    if IndexList:
        seqLen = len(info[2])
        for m in IndexList:
            speSeqList[m][1] += '-'*seqLen
    count += 1
a = [list(i[1]) for i  in speSeqList]
a = np.array(a)
print(a)
