import math

import matplotlib
import numpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon


class chunker():
    RAD_PER_DEG = math.pi / 180.0
    DEG_PER_RAD = 180.0 / math.pi
    EPSILON_DEG = 0.001 / 3600  # < 1 mas
    numStripes = 36  # 36
    numSubStripesPerStripe = 2
    overlap = 20.0 / 3600.0
    numChunksPerStripe = []
    numSubChunksPerChunk = []
    subChunkWidth = []
    alpha = []
    alphaStripe = []


    def show(self):
        print("numStripes", self.numStripes)
        print("numSubStripesPerStripe", self.numSubStripesPerStripe)
        print("overlap", self.overlap)
        print("numChunksPerStripe", self.numChunksPerStripe)
    def clampLon(self,lon):
        if lon > 360.0 - self.EPSILON_DEG:
            return 360.0
        return lon

    def clampLat(self,lat):
        if lat < -90.0:
            return -90.0
        elif lat > 90.0:
            return 90.0
        return lat

    def segments(self,latMin, latMax, width):
        lat = max(abs(latMin), abs(latMax))
        if lat > 90 - 1 / 3600.0:
            return 1
        if width >= 180:
            return 1
        elif width < 1 / 3600.0:
            width = 1 / 3600.0
        lat = lat * self.RAD_PER_DEG
        cw = math.cos(width * self.RAD_PER_DEG)
        sl = math.sin(lat)
        cl = math.cos(lat)
        x = cw - sl * sl
        u = cl * cl
        y = math.sqrt(abs(u * u - x * x))
        return math.floor(360.0 / abs(self.DEG_PER_RAD * math.atan2(y, x)))

    def segmentWidth(self, latMin, latMax, numSegments):
        lat = max(abs(latMin), abs(latMax)) * self.RAD_PER_DEG
        cw = math.cos(self.RAD_PER_DEG * (360.0 / numSegments))
        sl = math.sin(lat)
        cl = math.cos(lat)
        return math.acos(cw * cl * cl + sl * sl) * self.DEG_PER_RAD

    def getStripe(self, chunkId):
        return int(chunkId / (2 * self.numStripes))

    def getChunk(self, chunkId, stripe):
        return chunkId - stripe * 2 * self.numStripes

    def getChunkId(self, stripe, chunk):
        return stripe * 2 * self.numStripes + chunk

    # int32_t _getChunkId(int32_t stripe, int32_t chunk) const { return stripe * 2 * _numStripes + chunk; }

    def getChunkBounds(self,chunkId):
        stripe = self.getStripe(chunkId)
        chunk = self.getChunk(chunkId, stripe)
        width = 360.0 / self.numChunksPerStripe[stripe]
        lonMin = chunk * width
        lonMax = self.clampLon((chunk + 1) * width)
        latMin = self.clampLat(stripe * self.numSubStripesPerStripe * self.subStripeHeight - 90.0)
        latMax = self.clampLat((stripe + 1) * self.numSubStripesPerStripe * self.subStripeHeight - 90.0)
        return [lonMin, lonMax, latMin, latMax]

    def maxAlpha(self,r, centerLat):
        #print("r", r, centerLat)
        if r < 0.0 or r > 90.0:
            raise RuntimeError("Radius must lie in range [0, 90] deg.")
        if r == 0.0:
            return 0.0
        lat = self.clampLat(centerLat)
        if abs(lat) + r > 90.0 - 1 / 3600.0:
            return 180.0
        r = r * self.RAD_PER_DEG
        lat = lat * self.RAD_PER_DEG
        y = math.sin(r)
        x = math.sqrt(abs(math.cos(lat - r) * math.cos(lat + r)))
        return self.DEG_PER_RAD * abs(math.atan(y / x))

    def getNeighbours(self,chunkId):
        stripe = self.getStripe(chunkId)
        chunk = self.getChunk(chunkId, stripe)
        numChunks = self.numChunksPerStripe[stripe]
        neighbourIds = []
        if numChunks > 1:
            leftChunk = chunk - 1
            if leftChunk < 0:
                leftChunk = numChunks - 1
            neighbourIds.append(self.getChunkId(stripe, leftChunk))
            rightChunk = chunk + 1
            if rightChunk >= numChunks:
                rightChunk = 0
            neighbourIds.append(self.getChunkId(stripe, rightChunk))
        adjacentStripes = []
        if stripe != 0:
            adjacentStripes.append(stripe - 1)
        if stripe != self.numStripes - 1:
            adjacentStripes.append(stripe + 1)
        baseBB = self.getChunkBounds(chunkId)
        baseMinLon = baseBB[0]-self.alphaStripe[stripe]
        baseMaxLon = baseBB[1]+self.alphaStripe[stripe]
        for aStripe in adjacentStripes:
            numAChunks = self.numChunksPerStripe[aStripe]
            for ac in range(0, numAChunks):
                #  print ("AC", ac, numAChunks)
                acId = self.getChunkId(aStripe, ac)
                bb = self.getChunkBounds(acId)
                # account for overlap
                bb[0] = bb[0]-self.alphaStripe[aStripe]
                bb[1] = bb[1] + self.alphaStripe[aStripe]

                #  print("minLon",bb[0],"maxLon",bb[1])
                if bb[0] <= baseMinLon <= bb[1]:
                    neighbourIds.append(acId)
                elif bb[0] <= baseMaxLon <= bb[1]:
                    neighbourIds.append(acId)
                elif baseMinLon <= bb[0] and baseMaxLon >= bb[0]:
                    neighbourIds.append(acId)
            if chunk == 0:
                neighbourIds.append(self.getChunkId(aStripe, numAChunks - 1))
            if chunk == numChunks - 1:
                neighbourIds.append(self.getChunkId(aStripe, 0))
        #print("adjacent", adjacentStripes)
        return [*set(neighbourIds)]

    def __init__(self, numStripes, overlap, numSubStripesPerStripe=1):
        self.numStripes = numStripes
        self.overlap = overlap
        self.numSubStripesPerStripe=numSubStripesPerStripe

        if numStripes < 1 or self.numSubStripesPerStripe < 1:
            raise RuntimeError("The number of stripes and sub-stripes per stripe must be positive")
        if overlap < 0.0 or overlap > 10.0:
            raise RuntimeError("The overlap radius must be in range 0-10 deg")
        self.numSubStripes = numStripes * self.numSubStripesPerStripe
        self.stripeHeight = 180.0 / numStripes
        self.subStripeHeight = 180.0 / self.numSubStripes
        if self.subStripeHeight < overlap:
            raise RuntimeError("The overlap radius is greater than the sub-strip height")
        self.numChunksPerStripe = [0] * numStripes
        self.numSubChunksPerChunk = [0] * self.numSubStripes
        self.subChunkWidth = [0.0] * self.numSubStripes
        self.alpha = [0.0] * self.numSubStripes
        self.alphaStripe = [0] * numStripes
        maxSubChunksPerChunk = 0
        for i in range(0, numStripes):
            nc = self.segments(i * self.stripeHeight - 90.0, (i + 1) * self.stripeHeight - 90.0, self.stripeHeight)
            self.numChunksPerStripe[i] = nc
            for j in range(0, self.numSubStripesPerStripe):
                ss = i * self.numSubStripesPerStripe + j
                latMin = ss * self.subStripeHeight - 90.0
                latMax = (ss + 1) * self.subStripeHeight - 90.0
                nsc = self.segments(latMin, latMax, self.subStripeHeight) / nc
                maxSubChunksPerChunk = max(maxSubChunksPerChunk, nsc)
                self.numSubChunksPerChunk[ss] = nsc
                scw = 360.0 / (nsc * nc)
                self.subChunkWidth[ss] = scw
                a = self.maxAlpha(overlap, max(abs(latMin), abs(latMax)))
                if a > scw:
                    raise RuntimeError("The overlap radius is greater than the sub-chunk width.")
                self.alpha[ss] = a
                self.alphaStripe[i] = max(a, self.alphaStripe[i])

    def plot(self,minX=0.0, maxX=360.0, minY=-90.0, maxY=90.0 ):
        patches = []
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 9)
        plt.xlim([minX, maxX])
        plt.ylim([minY, maxY])

        for s in range(0, self.numStripes):
            for c in range(0, self.numChunksPerStripe[s]):
                bb = self.getChunkBounds(self.getChunkId(s, c))
                corners = [[bb[0], bb[2]], [bb[0], bb[3]], [bb[1], bb[3]], [bb[1], bb[2]]]
                polygon = Polygon(corners, True, linewidth=2, edgecolor='orange', facecolor='white')
                plt.text((bb[0] + bb[1]) / 2, (bb[2] + bb[3]) / 2, str(self.getChunkId(s, c)), fontsize=6,
                         horizontalalignment='center', verticalalignment='center', clip_on=True)
                if (c % 1 == 0):
                    patches.append(polygon)

        p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
        p.set_edgecolor('black')
        p.set_facecolor('white')
        p.set_linewidth(2)
        ax.add_collection(p)

        plt.show()

    def getTestPoints(self):
        points= []
        for s in range(0, self.numStripes):
            for c in range(0, self.numChunksPerStripe[s]):
                bb = self.getChunkBounds(self.getChunkId(s, c))
                cId=int(self.getChunkId(s, c))
                points.append((bb[0],bb[2],cId,-9999))
                points.append((bb[0],bb[3],cId,-9999))
                points.append((bb[1],bb[3],cId,-9999))
                points.append((bb[1], bb[2],cId,-9999))
                points.append(((bb[0]+bb[1])/2,(bb[2]+bb[3])/2,cId,0))
        pointsArray=numpy.asarray(points,dtype='f,f,i,i')
        return pointsArray
