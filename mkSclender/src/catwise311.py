import pandas as pd
import math
import logging
import os
import datetime
import healpy as hp
import numpy as np

logger = logging.getLogger('myLog')
formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                              "%Y-%m-%d %H:%M:%S")
logger.setLevel(logging.INFO )
fH     = logging.FileHandler(os.path.join(
            'logs',
            "catwiseLog_"+datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S')+'.log'))
fH.setFormatter(formatter)
logger.addHandler(fH)
logger.info("STARTINGRUN=")
# ratings_df = pd.read_csv('gaia_10000.csv', header=None)

#bcp "select cntr,source_name,ra, dec, sigra, sigdec, sigradec,ab_flags,w1mpro, w1sigmpro, w1sat,w1snr, w2mpro, w2sigmpro, w2sat, w2snr from catwise_2020   option (maxdop 1)"  queryout G:\share\mar\catwiseSkinnyish.csv -c  -S ramses19 -d WISE -U sa -P ibis1243! -t , -o G:\share\mar\catwise.log 


colNames=['cntr','designation','ra','dec','sigra','sigdec','sigradec','ab_flags','w1mpro','w1sigmpro','w1sat','w1snr','w2mpro','w2sigmpro','w2sat','w2snr']

colTypes={'cntr':int,'designation':object,'ab_flags':object}

#for csv_chunk in pd.read_csv('gaia_10000.csv', header=None,chunksize=1300,names=colNames):
#dtype = {'moon_lev':object,'size':float,
#                          'width':object,'color':object,
#                          'qty':int,'code': object})

chunkNum=0
chunkSize=1000000
skipRows=chunkNum*chunkSize
rowsWritten=0

nside=256
npix=hp.nside2npix(nside)
totCounts=np.zeros(npix,dtype=np.int32)
goodCounts=np.zeros(npix,dtype=np.int32)
badCounts=np.zeros(npix,dtype=np.int32)


possFlags=[2,1,0]
minRA=999
maxRA=-999
minDec=999
maxDec=-999
minMag=[999,999]
maxMag=[-999,-999]
minSig=999
maxSig=-999
reasons=np.zeros(shape=(10,2), dtype=np.int32)

print skipRows
logger.info("READINGCHUNK="+str(chunkNum))
for csv_chunk in pd.read_csv('/mnt/ramses25/mar/catwiseSkinnyish.csv', header=None,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):
#wiseSkinnyishGP.csv
#for csv_chunk in pd.read_csv('/mnt/ramses19/mar/wiseSkinnyish.csv', header=None,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):
#for csv_chunk in pd.read_csv('wise_10000.csv', header=None,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):
#for csv_chunk in pd.read_csv('unitWise.csv', header=None,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):

    logger.info("STARTINGCHUNK="+str(chunkNum))
    goodRows=[]
    badRows=[]
    chunkNum += 1
    print chunkNum
    radecs = csv_chunk.loc[:,('ra','dec')].values

#    print radecs[:,0]
#    print radecs[:,1]
    hpi=hp.ang2pix(nside, radecs[:,0],radecs[:,1],lonlat=True, nest=True)
    counts=np.bincount(hpi,minlength=npix)
    print counts
    print counts.shape, counts.dtype, hp.nside2npix(256)
    print counts[529707]
    totCounts=totCounts+counts
    print totCounts
    rn=-1
#    print(csv_chunk.shape)
    print (csv_chunk.dtypes)
 
#    for index, row in csv_chunk.iterrows():
#    for index, row in enumerate(csv_chunk.itertuples(),1):
    for row in csv_chunk.itertuples(index=False):
#        print row,row.source_id,"\n"
        rn += 1

        include = True
        best_flag=0
        sig_avg=-999
        flag=[2,2,2,2]
        mags=[row.w1mpro,row.w2mpro]
        sats=[row.w1sat,row.w2sat]
        magErrs=[row.w1sigmpro,row.w2sigmpro]
        snrs=[row.w1snr,row.w2snr]
        conds=[0,0]

        try:
#            print row.sigra,row.sigdec,row.sigradec

#            print row.sigra,row.sigdec,sig_maj,sig_min,flag[0],row.cc_flags
#            bitseq='0101'
#            binnum=int(bitseq,2)
#            print binnum


#            print row.cc_flags,row.ph_qual,mags
            
            for i in range(0, 2):
                if (row.sigra <=0 or row.sigdec <=0 or math.isnan(row.sigra) or math.isnan(row.sigdec)) :
                    flag[i]=-999
                    include=False
                    conds[i]=1
                    reasons[1,i] +=1
                elif (row.ab_flags[i] != "0") :
                    flag[i]=-1
                    conds[i]=2
                    reasons[2,i] +=1
                elif (math.isnan(mags[i]) or math.isnan(magErrs[i]) or magErrs[i] <=0):
                    flag[i]=-1
                    conds[i]=3
                    reasons[3,i] +=1
                elif (row.sigra > 2.0 or row.sigdec > 2.0):
                    flag[i]=-1
                    conds[i]=4
                    reasons[4,i] +=1
                elif (sats[i] > 0.1):
                    flag[i]=0
                    conds[i]=5
                    reasons[5,i] +=1
                elif (snrs[i]<5 or magErrs[i] >=0.5):
                    flag[i]=1
                    conds[i]=6
                    reasons[6,i] +=1
                       


                if (flag[i]==-999 or flag[i]==-1):
                    mags[i]=float("NaN")
                    
            covar = row.sigradec * abs(row.sigradec)
            sig_maj = math.sqrt(0.5 * (row.sigra**2 + row.sigdec**2 + math.sqrt((row.sigra**2 - row.sigdec**2)**2 + 4 * covar**2)))
            sig_min = math.sqrt(0.5 * (row.sigra**2 + row.sigdec**2 - math.sqrt((row.sigra**2 - row.sigdec**2)**2 + 4 * covar**2)))
            sig_avg = (sig_min+sig_maj)/2.0            
            
                
 
        except Exception as e:
            print e
            logger.error(e)
            logger.error(row)
            logger.error(flag)
            include=False
            for j in range(0, 2):
                flag[j]=-999
            
            
        bestFlag=-999
        bestIndex=-999
        for possFlag in possFlags:
            for j in range(0, 2):
                if (flag[j] == possFlag):
                    bestFlag=possFlag
                    bestIndex=j
                    break
            if (bestFlag >=0):
                break
                
#        print "flags", flag,conds,bestIndex,bestFlag
        
        if (bestFlag< 0) :
            include=False
        if (include):
            goodRows.append([row.cntr,row.designation,row.ra,row.dec,sig_avg,mags[0],mags[1],bestIndex,hpi[rn],bestFlag,flag[0],flag[1],row.w1mpro,row.w2mpro])
            rowsWritten += 1
            goodCounts[hpi[rn]] += 1
            if (row.ra < minRA):
                minRA=row.ra
            if (row.dec < minDec):
                minDec=row.dec
            if (sig_avg < minSig):
                minSig=sig_avg
            if minMag[0]>mags[0]:
                minMag[0]=mags[0]
            if minMag[1]>mags[1]:
                minMag[1]=mags[1]
            if (row.ra > maxRA):
                maxRA=row.ra
            if (row.dec > maxDec):
                maxDec=row.dec
            if (sig_avg > maxSig):
                maxSig=sig_avg
            if maxMag[0]<mags[0]:
                maxMag[0]=mags[0]
            if maxMag[1]<mags[1]:
                maxMag[1]=mags[1]
        else :
            badRows.append([row.cntr,row.designation,row.ra,row.dec,sig_avg,mags[0],mags[1],bestIndex,hpi[rn],bestFlag,flag[0],flag[1],row.w1mpro,row.w2mpro])
            badCounts[hpi[rn]] += 1

    goodDF=pd.DataFrame(goodRows,columns=['cntr','designation','ra','dec','sig_avg','w1mag','w2mag','bestIndex','hpi','bestFlag', \
                                          'flag0','flag1','o1mag','o2mag'])
    goodDF=goodDF.round({'ra': 9, 'dec': 8, 'sig_avg': 8, 'w1mag': 3, 'w2mag': 3, 'o1mag': 3, 'o2mag': 3})
    badDF=pd.DataFrame(badRows,columns=['cntr','designation','ra','dec','sig_avg','w1mag','w2mag','bestIndex','hpi','bestFlag', \
                                          'flag0','flag1','o1mag','o2mag'])
    badDF=badDF.round({'ra': 9, 'dec': 8, 'sig_avg': 8, 'w1mag': 3, 'w2mag': 3, 'o1mag': 3, 'o2mag': 3})
    logger.info("WRITINGCHUNK="+str(chunkNum))
    goodDF.to_csv('/unas/mar/catwiseGood.csv', mode='a', header=False, index=False)
    badDF.to_csv('/unas/mar/catwiseBad.csv', mode='a', header=False, index=False)
    logger.info("WRITTENCHUNK="+str(chunkNum))
#    for csv_row in csv_chunk.rows:
#        print row.ra,row.dec
#        print row

logger.info("MIN="+" "+str(minRA)+" "+str(minDec)+" "+str(minSig)+" "+str(minMag[0])+" "+str(minMag[1]))
logger.info("MAX="+" "+str(maxRA)+" "+str(maxDec)+" "+str(maxSig)+" "+str(maxMag[0])+" "+str(maxMag[1]))
print reasons
print "saving counts"

np.savetxt("/unas/mar/catwiseTotCounts.csv", totCounts, delimiter=",")
np.savetxt("/unas/mar/catwiseGoodCounts.csv", goodCounts, delimiter=",")
np.savetxt("/unas/mar/catwiseBadCounts.csv", badCounts, delimiter=",")

