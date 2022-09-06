import numpy

import chunkerClass
import pandas as pd

# initialize with number of stripes and overlap in arcsecs
# stripes cover 180 degrees in latitude, so 36 stripes is equivalent to 5 degree chunking

chunker=chunkerClass.chunker(36,20/3600)

chunker.show()

# find neighbours for a given chunkid

neighbours=chunker.getNeighbours(0)
print(neighbours)
neighbours=chunker.getNeighbours(1041)
print(neighbours)

# get an array of RA/Decs of corners and centres of each chunk.

npa=chunker.getTestPoints()
print(npa)
DF = pd.DataFrame(npa)
DF.to_csv("data.csv",index=False,header=None)

# show a plot of chunks

chunker.plot()

# get a count of neighbours for every chunk

neghbourCounts=[]
for s in range(0, chunker.numStripes):
    for c in range(0, chunker.numChunksPerStripe[s]):
        n=chunker.getNeighbours(chunker.getChunkId(s, c))
        neghbourCounts.append((chunker.getChunkId(s, c),len(n)))
        if len(n)==7:
            #show ones with 7 neighbours
            print(chunker.getChunkId(s, c),len(n))
neighbourArray=numpy.asarray(neghbourCounts,dtype='i,i')
#print(neighbourArray)
NF = pd.DataFrame(neighbourArray)
NF.to_csv("neighbour.csv",index=False,header=None)
print("END")


