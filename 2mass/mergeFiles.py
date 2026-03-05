import glob
import sys
import os

dirPath= sys.argv[1]
print (dirPath)
overlapPath=dirPath+'/overlaps/'
txtList = glob.glob(dirPath+'/*.txt')
#print (txtList)
for chunkFile in txtList:
    print (chunkFile)
    
    overlapFile=overlapPath+os.path.basename(chunkFile).replace('.txt','_overlap.txt')
    print(overlapFile)

    if os.path.exists(overlapFile):
        f1 = open(chunkFile, 'a+')
        f2 = open(overlapFile, 'r')
        f1.write(f2.read().replace("0\n","1\n"))
        f1.close()
        f2.close()
    else :
        print ("overlap DOES NOT EXITS")
