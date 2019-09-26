import sys
import numpy as np

filePath = sys.argv[1]
codonDict = {'CTT':'L','CTC':'L','CTA':'L','CTG':'L','GTT':'V','GTC':'V','GTG':'V','GTA':'V', \
             'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCG':'P','CCA':'P', \
             'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCG':'A','GCA':'A', \
             'CGT':'R','CGC':'R','CGA':'R','CGG':'R','GGT':'G','GGC':'G','GGG':'G','GGA':'G'}
seqList = []
speList = []
with open(filePath,'r') as f:
    header = f.readline()
    for line in f:
        info = line.strip().split()
        seqList.append(list(info[1]))
        speList.append(info[0])

seqArray = np.array(seqList)
seqLength = seqArray.shape[1]
if seqLength%3 != 0:
    print("The length of sequences is not a multiple of 3!")
    sys.exit()

stackBool = np.stack((seqArray == 'A', seqArray == 'G', seqArray == 'C', seqArray == 'T',seqArray == '-'))    
boolAGCT =  np.sum(np.any(stackBool, axis = 1), axis = 0)
cleanPos = np.where(boolAGCT == 1)[0]

preCodonPos = np.array([i for i in range(0,seqLength,3) if (i in cleanPos) and (i+1 in cleanPos)])
codonSelect = np.char.add(np.char.add(seqArray[:,preCodonPos],seqArray[:,preCodonPos+1]),seqArray[:,preCodonPos+2])

temp = []
for i in range(codonSelect.shape[1]):
    tmp = codonSelect[:,i]
    for j in tmp:
        if j in codonDict:
            temp.append(i)
            break

a = seqArray[:,preCodonPos[temp]+2]
print('\t'+str(len(speList))+'\t'+str(len(temp)))
for i in range(0,len(speList)):
    print(speList[i]+'\t'+''.join(list(a[i])))
