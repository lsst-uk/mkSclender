import pandas as pd
import math
import logging
import os
import datetime
import healpy as hp
import numpy as np
import sys



from astropy.coordinates import SkyCoord
import astropy.units as u


inCSV='/disk82/panstarrs/2MASS_PSC_all_SNR.csv'

#inCSV='/mnt/ramses27/mar/2MASS_PSC_all_SNR.csv'
outGood='pscGoodSNRNan.csv'
outBad='pscBadSNRNan.csv'

numBands=3
filters=['j','h','k']
searchOrder=[2,0,1] #kjh
rowCount=0
possFlags=[2,1,0]


def flagging(row,rn):
    flags=[2,2,2]
    reasons=[0,0,0]
    useMags=[]
    useMagErrs=[]
    useMagFlags=[]


    useMags=[row.j_m,row.h_m,row.k_m]
    useMagErrs=[row.j_cmsig,row.h_cmsig,row.k_cmsig]

    for b in range(0, numBands):
        # NonDetection-1 [jhk]_m is NaN OR [jhk]_cmsig is NaN
        if  pd.isna(useMagErrs[b]) or  pd.isna(useMags[b]) :
            flags[b]=-1
            reasons[b]=1
            useMags[b]=np.nan
        # bl_flg == 0 OR rd_flg == 0 OR rd_flg == 9
        elif row.bl_flg[b]=='0' or row.rd_flg[b]=='0' or row.rd_flg[b]=='9' :
            flags[b]=-1
            reasons[b]=2
            useMags[b]=np.nan
        # Saturation 0 rd_flg == 3
        elif row.rd_flg[b]=='3' :
            flags[b]=0
            reasons[b]=3
        # Lowquality 1 rd_flg == 4 OR rd_flg == 6 OR bl_flg > 3 OR cc_flg != “0” OR ph_qual == “E” OR [jhk]_cmsig > 0.5

        elif row.rd_flg[b]=='4' or row.rd_flg[b]=='6' or int(row.bl_flg[b]) > 3 or row.cc_flg[b]!='0' or row.ph_qual[b]=='E':
            flags[b]=1
            reasons[b]=4
        elif useMagErrs[b] > 0.5 :
            flags[b]=1
            reasons[b]=5
      
    maxFlag=max(flags)
    foundBest=False
    bestBand=-1     
    if maxFlag >=0 :
        for possFlag in possFlags:
            for sb in searchOrder:
                if flags[sb]==possFlag :
                    foundBest=True
                    bestBand=sb
                    break
            if foundBest :
                break

    return bestBand,maxFlag,useMags,flags,reasons,useMagErrs





logger = logging.getLogger('myLog')
formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                              "%Y-%m-%d %H:%M:%S")
logger.setLevel(logging.INFO )
fH     = logging.FileHandler(os.path.join(
            'logs',
            "pscLog_"+datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S')+'.log'))
fH.setFormatter(formatter)
logger.addHandler(fH)
logger.info("STARTINGRUN=")


colNames=['pts_key','designation','ra','decl','err_maj','err_min','j_m','j_cmsig','h_m','h_cmsig','k_m','k_cmsig','bl_flg','rd_flg','cc_flg','ph_qual',\
          'jdate','j_snr','h_snr','k_snr']


colTypes={'pts_key':int,'designation':object,'ra':float,'decl':float,'err_maj':float,'err_min':float,\
          'j_m':float,'j_cmsig':float,'h_m':float,'h_cmsig':float,'k_m':float,'k_cmsig':float,\
          'bl_flg':object,'rd_flg':object,'cc_flg':object,'ph_qual':object,\
          'jdate':float,'j_snr':float,'h_snr':float,'k_snr':float}
          

#declare @cmd varchar(8000)
#set @cmd = 'bcp "select '+
#'pts_key,''2MASS J''+designation as designation,ra,decl,err_maj,err_min,j_m,j_cmsig,h_m,h_cmsig,k_m,k_cmsig,bl_flg,rd_flg,cc_flg,'+
#'ph_qual,jdate,j_snr,h_snr,k_snr '+
#'from twomass_psc where ext_key is null " ' +
#' queryout G:\share\mar\2MASS_PSC_all_SNR.csv -c '+
#' -S ramses27 -d TWOMASS2 -U sa -P XXXXX -t , -o G:\share\mar\2MASS_PSC.log '
#select @cmd
#exec xp_cmdshell @cmd


#sys.exit()



chunkNum=0
chunkSize=1000000
skipRows=chunkNum*chunkSize
rowsWritten=0
nside=256
npix=hp.nside2npix(nside)
totCounts=np.zeros(npix,dtype=np.int32)
goodCounts=np.zeros(npix,dtype=np.int32)
badCounts=np.zeros(npix,dtype=np.int32)
bestIndexCounts=np.zeros(numBands,dtype=np.int32)


minRA=999
maxRA=-999
minDec=999
maxDec=-999
minMag=[999,999,999,999]
maxMag=[-999,-999,-999,-999]
minSig=999
maxSig=-999
reasons=np.zeros(shape=(10,3), dtype=np.int32)




logger.info("inCSV="+inCSV)
logger.info("outBad="+outGood)
logger.info("outBad="+outBad)
#print skipRows
logger.info("READINGCHUNK="+str(chunkNum))
#print colNames
#print colTypes


#/mnt/ramses27/mar/2MASS_PSC_all.csv



for csv_chunk in pd.read_csv(inCSV, na_values=['', 'nan', 'NaN', 'inf', '-inf'],header=None,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):
    boolean_columns = csv_chunk.select_dtypes(include='boolean').columns
    #csv_chunk[boolean_columns] =  csv_chunk[boolean_columns].applymap(lambda x: True if pd.api.types.isna(x) else x)
    csv_chunk[boolean_columns] = csv_chunk[boolean_columns].fillna(True)


    logger.info("STARTINGCHUNK="+str(chunkNum))
    goodRows=[]
    badRows=[]
    chunkNum += 1
    print (chunkNum)

    ra_values=csv_chunk.loc[:,('ra')].values
    dec_values=csv_chunk.loc[:,('decl')].values        

    rn=-1

    for row in csv_chunk.itertuples(index=False):
        rn +=1
#        print(row)
        bestBand,maxFlag,mags,flags,flagReasons,magErrs=flagging(row,rn)
#        print(bestBand,maxFlag,mags,flags,flagReasons,magErrs)
        

        if rn%100000==0 :
            logger.info("ROWNUM="+str(rn))

        include=False
        if bestBand >=0 :
            include=True
        #colNames=['pts_key','designation','ra','decl','err_maj','err_min','j_m','j_cmsig','h_m','h_cmsig','k_m','k_cmsig','bl_flg','rd_flg','cc_flg','ph_qual','jdate']

        sig_avg = math.sqrt(row.err_maj * row.err_min)
        snrs=[row.j_snr,row.h_snr,row.k_snr]
        snrs = [np.nan if math.isnan(um) else snr for um, snr in zip(mags, snrs)]

        if include:
            goodRows.append([row.designation,row.ra,row.decl,sig_avg,mags[0],mags[1],mags[2],bestBand,row.jdate,
                             snrs[0],snrs[1],snrs[2],
                             magErrs[0],magErrs[1],magErrs[2],row.pts_key,maxFlag,flags[0],flags[1],flags[2]])
        else :
            badRows.append([row.designation,row.ra,row.decl,sig_avg,mags[0],mags[1],mags[2],bestBand,row.jdate,
                             snrs[0],snrs[1],snrs[2],
                             magErrs[0],magErrs[1],magErrs[2],row.pts_key,maxFlag,flags[0],flags[1],flags[2]])


    outCols=['designation','ra','dec','sig_avg','mag0','mag1','mag2','bestIndex','jdate','j_snr','h_snr','k_snr','magErr0','magErr1','magErr2','pts_key','bestFlag','pf0','pf1','pf2']                            
    goodDF=pd.DataFrame(goodRows,columns=outCols)
    badDF=pd.DataFrame(badRows,columns=outCols)
    goodDF=goodDF.round({'ra': 9, 'dec': 8, 'sig_avg': 9, 'mag0': 4, 'mag1': 4, 'mag2': 4, 'magErr0': 4, 'magErr1': 4, 'magErr2': 4, 'j_snr': 6, 'h_snr': 6, 'k_snr': 6})
    badDF=badDF.round({'ra': 9, 'dec': 8, 'sig_avg': 9, 'mag0': 4, 'mag1': 4, 'mag2': 4, 'magErr0': 4, 'magErr1': 4, 'magErr2': 4, 'j_snr': 6, 'h_snr': 6, 'k_snr': 6})
    logger.info("WRITINGCHUNK="+str(chunkNum))
    goodDF.to_csv(outGood, mode='a', header=False, index=False)
    badDF.to_csv(outBad, mode='a', header=False, index=False)
    logger.info("WRITTENCHUNK="+str(chunkNum))
