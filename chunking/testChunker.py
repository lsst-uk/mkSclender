import numpy

import chunkerClass
import pandas as pd

chunker=chunkerClass.chunker(36,20/3600)

chunker.show()


neighbours=chunker.getNeighbours(0)
print(neighbours)
npa=chunker.getTestPoints()
print(npa)
chunker.plot()
DF = pd.DataFrame(npa)
DF.to_csv("data.csv",index=False,header=None)
neighbours=chunker.getNeighbours(1041)
print(neighbours)
print(chunker.getChunkBounds(1041),chunker.alphaStripe[chunker.getStripe(1041)])
print(chunker.getChunkBounds(969),chunker.alphaStripe[chunker.getStripe(969)])

neghbourCounts=[]
for s in range(0, chunker.numStripes):
    for c in range(0, chunker.numChunksPerStripe[s]):
        n=chunker.getNeighbours(chunker.getChunkId(s, c))
        neghbourCounts.append((chunker.getChunkId(s, c),len(n)))
        if len(n)==7:
            print(chunker.getChunkId(s, c),len(n))
neighbourArray=numpy.asarray(neghbourCounts,dtype='i,i')
#print(neighbourArray)
NF = pd.DataFrame(neighbourArray)
NF.to_csv("neighbour.csv",index=False,header=None)
print("END")


