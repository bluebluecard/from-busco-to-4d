import sys
import os

def splitLine(seqStr,line):

    newSeq = []
    for i in range(0,len(seqStr),line):
        newSeq.append(seqStr[i:i+line])
    return newSeq

fileList = sys.argv[1]

buscoNumber = sys.argv[2]
### species \t genome/gene \t busco full table \t busco single copy dir/ gene pep cds dir

buscoList = {}
with open(fileList,'r') as f1:
    for line in f1:
        info = line.strip().split('\t')
        proteinBusco = {}
        with open(info[2],'r') as f2:
            for buscoLine in f2:
                data = buscoLine.strip().split('\t')
                if "Complete" in data:
                    if data[0] in buscoList:
                        buscoList[data[0]].append(data[0]+'_'+info[0])
                    else:
                        buscoList[data[0]] = []
                        buscoList[data[0]].append(data[0]+'_'+info[0])

                    proteinBusco[data[2]] = data[0]+'_'+info[0] 

        for path, dirList, fileList in os.walk(info[3]):
            for fileName in fileList:
                seq = {}
                filePath = os.path.join(path,fileName)
                if info[1] == "gene" and fileName.endswith(('pep','cds')):

                    data = []
                    with open(filePath,'r') as f:
                        for line in f:
                            if '>' in line:
                                seqName = line.strip().lstrip('>').split()[0]
                                seq[seqName] = ''
                            else:
                                seq[seqName] += line.strip()

                    for i in proteinBusco:                            
                        data.append('>'+proteinBusco[i]+'\n')
                        for mSeq in splitLine(seq[i],60):
                            data.append(mSeq+'\n')

                    if fileName.endswith('pep'):
                        with open('all.pep','a') as fo:
                            fo.writelines(data)
                    elif fileName.endswith('cds'):
                        with open('all.cds','a') as fo:
                            fo.writelines(data)

                elif info[1] == "genome" and (fileName.endswith(('faa','fna'))):

                    data = []
                    with open(filePath,'r') as f:
                        for line in f:
                            if '>' in line:
                                line = line.split(':')[0]+'_'+info[0]+'\n'
                                data.append(line)        
                            else:
                                for mSeq in splitLine(line,60):
                                    data.append(mSeq+'\n')
                    if fileName.endswith('faa'):
                        with open('all.pep','a') as fo:
                            fo.writelines(data)
                    elif fileName.endswith('fna'):
                        with open('all.cds','a') as fo:
                            fo.writelines(data)
                        
count = 1
for buscoName in buscoList:
    
    if len(buscoList[buscoName]) >= int(buscoNumber):
        print(str(count)+'\t'+'\t'.join(buscoList[buscoName]))
        count += 1
