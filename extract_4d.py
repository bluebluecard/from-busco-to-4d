import sys
import numpy as np

filePath = sys.argv[1]
codonPath = 'codon_4d.txt'

codonDict = {}
with open(codonPath,'r') as f:
    for line in f:
        info = line.strip().split()
        codonDict[info[0]] = info[1]

seqList = []
with open(filePath,'r') as f:
    for line in f:
        info = line.strip().split()
        seqList.append(list(info[1]))

seqArray = np.array(seqList)
seqLength = seqArray.shape[1]

firstBaseIndex = [i for i in range(0,seqLength,3)]
secondBaseIndex = [i+1 for i in range(0,seqLength,3)]
thirdBaseIndex = [i+2  for i in range(0,seqLength,3)]

firstBaseArray = seqArray[:,firstBaseIndex]
secondBaseArray = seqArray[:,secondBaseIndex]
thirdBaseArray = seqArray[:,thirdBaseIndex]

a = []
b = []

for i in range(firstBaseArray.shape[1]):
    tmp = firstBaseArray[:,i]
    if (len(set(tmp)) == 2 and '-' in tmp) or len(set(tmp)) == 1:
        a.append(i)
    tmp = secondBaseArray[:,i]
    if (len(set(tmp)) == 2 and '-' in tmp) or len(set(tmp)) == 1:
        b.append(i)

c = list(set(a) & set(b))
d = np.char.add(np.char.add(firstBaseArray[:,c],secondBaseArray[:,c]),thirdBaseArray[:,c])
temp = []

for i in range(d.shape[1]):
    tmp = d[:,i]
    for j in tmp:
        if j in codonDict:
            temp.append(i)
            break
print(thirdBaseArray[:,c][:,temp])
