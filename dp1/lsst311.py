import pandas as pd
import math
import logging
import os
import datetime
import numpy as np
import sys
from pprint import pformat
from collections import Counter

from astropy.coordinates import SkyCoord
import astropy.units as u



latLimit=20.0
numBands=6
reasonCounts = [Counter() for _ in range(numBands)]
filters=['u','g','r','i','z','y']
#irzgyu
searchOrder=[3,2,4,1,5,0]
bestBandCounts = {f: 0 for f in filters}
numInPlane=0
rowCount=0
defaultEpoch=-9999.0
nanMags=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

def convert_to_galactic(ra_values, dec_values):
    """
    Convert RA and Dec coordinates to Galactic coordinates.
    
    Parameters:
        ra_values (array-like): Array-like object containing RA values in degrees.
        dec_values (array-like): Array-like object containing Dec values in degrees.
    
    Returns:
        array-like: Array-like object containing Galactic longitude and latitude in degrees.
    """

    ra_values = np.array(ra_values)
    dec_values = np.array(dec_values)
    # Create SkyCoord object with input RA and Dec values
    
    coords = SkyCoord(ra=ra_values*u.deg, dec=dec_values*u.deg, frame='icrs')
    
    # Convert to Galactic coordinates
    galactic_coords = coords.transform_to('galactic')
    
    # Return Galactic longitude and latitude
    return galactic_coords.l.deg, galactic_coords.b.deg


def ra_dec_to_galactic(ra, dec):
    """
    Convert Right Ascension (RA) and Declination (Dec) coordinates to 
    Galactic Longitude (l) and Latitude (b).
    
    Parameters:
        ra (float): Right Ascension in degrees.
        dec (float): Declination in degrees.
        
    Returns:
        tuple: Galactic Longitude (l) and Latitude (b) in degrees.
    """
    # Create a SkyCoord object with the RA and Dec coordinates
    sky_coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    # Convert to Galactic coordinates (l, b)
    galactic_coord = sky_coord.galactic

    # Extract Galactic longitude (l) and latitude (b)
    gLong = galactic_coord.l.degree
    gLat = galactic_coord.b.degree

    return gLong, gLat


def flagging(row,useMags,useMagErrs,useMagFlags,rn,useSNRs):
    photFlags=[2,2,2,2,2,2]
    photReasons=[0,0,0,0,0,0]
    #print('isPlane',isPlane)


    for b in range(0, numBands):
        # *_magerr 1 == NaN
        if pd.isna(useMagErrs[b]) or  pd.isna(useMags[b]) :
            photFlags[b]=-1
            photReasons[b]=1
            useMags[b]=np.nan
        #*_magerr <= 0
        elif useMagErrs[b] <= 0.0 :
            photFlags[b]=-1
            photReasons[b]=2
            useMags[b]=np.nan
        elif pd.isna(useSNRs[b]) :
            photFlags[b]=-1
            photReasons[b]=3
            
        # R *_flag is True OR merge_footprint_* is False OR *_pixelflags_edge is True OR *_pixelflags_offimage is True

        #elif useMagFlags[b] or not getattr(row,'merge_footprint_'+filters[b]) or getattr(row,filters[b]+'_pixelflags_edge')  or getattr(row,filters[b]+'_pixelflags_offimage') :
        elif useMagFlags[b] or getattr(row,filters[b]+'_pixelFlags_edge')  or getattr(row,filters[b]+'_pixelFlags_offimage') :
            photFlags[b]=-1
            photReasons[b]=4
            useMags[b]=np.nan
        # *_pixelflags_saturated is True
        elif getattr(row,filters[b]+'_pixelFlags_saturated') :
            photFlags[b]=0
            photReasons[b]=5
        elif useMagErrs[b] >= 0.5 :
            photFlags[b]=1
            photReasons[b]=6
        # *_pixelflags_bad is True OR *_pixelflags_clipped is True OR *_pixelflags_cr is True OR *_pixelflags_sensor_edge is True OR *_pixelflags_suspect is True
        elif getattr(row,filters[b]+'_pixelFlags_bad')  or  getattr(row,filters[b]+'_pixelFlags_clipped') or  getattr(row,filters[b]+'_pixelFlags_cr') or  \
             getattr(row,filters[b]+'_pixelFlags_sensor_edge') or getattr(row,filters[b]+'_pixelFlags_suspect') :
            photFlags[b]=1
            photReasons[b]=7

    maxPhotFlag=max(photFlags)
    foundBestPhot=False
    bestPhotBand=-1     

    for possFlag in possFlags:
        for sb in searchOrder:
            if photFlags[sb]==possFlag :
                foundBestPhot=True
                bestPhotBand=sb
                break
            if foundBestPhot :
                break
    return bestPhotBand,maxPhotFlag,useMags,photFlags,photReasons,useMagErrs,useSNRs

#    print (rn,'mags,magerrs',useMags,useMagErrs,useMagFlags,photFlags,photReasons)
#    print (rn,'BESTAST',bestAstBand,'BESTPHOT',bestPhotBand,photFlags,photReasons)
    

def astFlagging(row):
    global rowCount
    rowCount +=1
    astFlags=[2,2,2,2,2]
    astReasons=[0,0,0,0,0]
    for b in range(0, numBands):
#        print('flaggggTYPE',type(getattr(row,filters[b]+'_psfflux_flag')))
#        print('footprint',getattr(row,'merge_footprint_'+filters[b]))
#        if getattr(row,'merge_footprint_'+filters[b]):
#            print("TRRRRRRRUE")
#        else:
#            print("FAAALSE")
#        print('measurement',getattr(row,'merge_measurement_'+filters[b]))
#        if not getattr(row,'merge_measurement_'+filters[b]):
#            print("measurement FALSE")
#        else:
#            print("measurement TRUE")
# *_sdsscentroid_raerr == NaN OR *_sdsscentroid_decerr == NaN OR
        if pd.isna(getattr(row,filters[b]+'_sdsscentroid_raerr')) or pd.isna(getattr(row,filters[b]+'_sdsscentroid_decerr')):            
            astFlags[b]=-1
            astReasons[b]=1
            #print('flagging',b)
# *_sdsscentroid_raerr <= 0 OR *_sdsscentroid_decerr <= 0 OR
        elif getattr(row,filters[b]+'_sdsscentroid_raerr') <=0 or getattr(row,filters[b]+'_sdsscentroid_decerr') <=0:
            astFlags[b]=-1
            astReasons[b]=2
# *_sdsscentroid_flag is True OR merge_footprint_* is False OR
        elif getattr(row,filters[b]+'_sdsscentroid_flag') or not getattr(row,'merge_footprint_'+filters[b]):
            astFlags[b]=-1
            astReasons[b]=3
            #print(b,'HHHHHH',rowCount)
# *_sdsscentroid_raerr >= 1 OR *_sdsscentroid_decerr >= 1 OR
        elif getattr(row,filters[b]+'_sdsscentroid_raerr') >=1 or getattr(row,filters[b]+'_sdsscentroid_decerr') >=1:
            astFlags[b]=-1
            astReasons[b]=4
            #print(b,'GGGGG',rowCount, getattr(row,filters[b]+'_sdsscentroid_raerr'), getattr(row,filters[b]+'_sdsscentroid_decerr'))
# *_pixelflags_edge is True OR *_pixelflags_offimage is True
        elif  getattr(row,filters[b]+'_pixelflags_edge')  or getattr(row,filters[b]+'_pixelflags_offimage') :
            astFlags[b]=-1
            astReasons[b]=5
            #print(b,'JJJJJJ',rowCount,getattr(row,filters[b]+'_pixelflags_edge'), getattr(row,filters[b]+'_pixelflags_offimage'))
# *_pixelflags_saturated is True
        elif getattr(row,filters[b]+'_pixelflags_saturated') :
            astFlags[b]=0
            astReasons[b]=6
            #print(b,'KKKKK',rowCount,getattr(row,filters[b]+'_pixelflags_saturated'))
# *_sdsscentroid_raerr >= 0.5 OR *_sdsscentroid_decerr >= 0.5 OR
        elif getattr(row,filters[b]+'_sdsscentroid_raerr') >=0.5 or getattr(row,filters[b]+'_sdsscentroid_decerr') >=0.5:
            astFlags[b]=1
            astReasons[b]=6
# *_pixelflags_bad is True OR *_pixelflags_clipped is True OR
        elif getattr(row,filters[b]+'_pixelflags_bad')  or getattr(row,filters[b]+'_pixelflags_clipped') :
            astFlags[b]=1
            astReasons[b]=7
            #print(b,'LLLLLL',rowCount,getattr(row,filters[b]+'_pixelflags_bad'),getattr(row,filters[b]+'_pixelflags_clipped'))
# *_pixelflags_cr is True OR *_pixelflags_sensor_edge is True OR *_pixelflags_suspect is True
        elif getattr(row,filters[b]+'_pixelflags_cr')  or getattr(row,filters[b]+'_pixelflags_sensor_edge') or getattr(row,filters[b]+'_pixelflags_suspect') :
            astFlags[b]=1
            astReasons[b]=8
            #print(b,'MMMMM',rowCount,getattr(row,filters[b]+'_pixelflags_cr'),getattr(row,filters[b]+'_pixelflags_sensor_edge'),getattr(row,filters[b]+'_pixelflags_suspect'))
        
    #print (rowCount,'astflags',astFlags,astReasons)
    maxAstFlag=max(astFlags)
    foundBest=False
    bestAstBand=-1
    if maxAstFlag >=0 :
        for sb in searchOrder:
            if  getattr(row,'merge_peak_'+filters[sb]) and astFlags[sb]==maxAstFlag :
                foundBest=True
                bestAstBand=sb
                break
            if not foundBest :
                for possFlag in possFlags:
                    for sb in searchOrder:
                        if astFlags[sb]==possFlag :
                            foundBest=True
                            bestAstBand=sb
                            break
                    if foundBest :
                        break
    useRA=row.ra
    useDec=-row.dec
    sig_avg=np.nan

    

    #print(rowCount,'maxAstFlag',maxAstFlag,foundBest,bestAstBand,astFlags)
    if bestAstBand > -1:
        useRA=getattr(row,filters[bestAstBand]+'_sdsscentroid_ra')
        useDec=getattr(row,filters[bestAstBand]+'_sdsscentroid_dec')
        sig_avg=math.sqrt(getattr(row,filters[bestAstBand]+'_sdsscentroid_raerr')*getattr(row,filters[bestAstBand]+'_sdsscentroid_decerr'))
        #print('RA,Dec',filters[bestAstBand])
        
    #print('RA,Dec',useRA,useDec,sig_avg)

    
    return useRA,useDec,sig_avg,bestAstBand,astFlags,astReasons

#astFlags=astFlagging()
#print (astFlags)




logger = logging.getLogger('myLog')
formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                              "%Y-%m-%d %H:%M:%S")
logger.setLevel(logging.INFO )
fH     = logging.FileHandler(os.path.join(
            'logs',
            "lsstLog_"+datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S')+'.log'))
fH.setFormatter(formatter)
logger.addHandler(fH)
logger.info("STARTINGRUN=")

#objectId,coord_ra,coord_dec,coord_raErr,coord_decErr,coord_ra_dec_Cov,
#u_cModelMag,u_cModelMagErr,u_cModelFlux,u_cModelFluxErr,u_cModel_flag,
#u_psfMag,u_psfMagErr,u_psfFlux,u_psfFluxErr,u_psfFlux_flag,u_pixelFlags_edge,u_pixelFlags_offimage,
#u_pixelFlags_saturated,u_pixelFlags_bad,u_pixelFlags_clipped,u_pixelFlags_cr,u_pixelFlags_sensor_edge,u_pixelFlags_suspect,u_epoch,
#g_cModelMag,g_cModelMagErr,g_cModelFlux,g_cModelFluxErr,g_cModel_flag,
#g_psfMag,g_psfMagErr,g_psfFlux,g_psfFluxErr,g_psfFlux_flag,g_pixelFlags_edge,g_pixelFlags_offimage,
#g_pixelFlags_saturated,g_pixelFlags_bad,g_pixelFlags_clipped,g_pixelFlags_cr,g_pixelFlags_sensor_edge,g_pixelFlags_suspect,g_epoch,
#r_cModelMag,r_cModelMagErr,r_cModelFlux,r_cModelFluxErr,r_cModel_flag,
#r_psfMag,r_psfMagErr,r_psfFlux,r_psfFluxErr,r_psfFlux_flag,r_pixelFlags_edge,r_pixelFlags_offimage,
#r_pixelFlags_saturated,r_pixelFlags_bad,r_pixelFlags_clipped,r_pixelFlags_cr,r_pixelFlags_sensor_edge,r_pixelFlags_suspect,r_epoch,
#r_cModelMag,i_cModelMagErr,i_cModelFlux,i_cModelFluxErr,i_cModel_flag,
#i_psfMag,i_psfMagErr,i_psfFlux,i_psfFluxErr,i_psfFlux_flag,i_pixelFlags_edge,i_pixelFlags_offimage,
#i_pixelFlags_saturated,i_pixelFlags_bad,i_pixelFlags_clipped,i_pixelFlags_cr,i_pixelFlags_sensor_edge,i_pixelFlags_suspect,i_epoch,
#z_cModelMag,z_cModelMagErr,z_cModelFlux,z_cModelFluxErr,z_cModel_flag,
#z_psfMag,z_psfMagErr,z_psfFlux,z_psfFluxErr,z_psfFlux_flag,z_pixelFlags_edge,z_pixelFlags_offimage,
#z_pixelFlags_saturated,z_pixelFlags_bad,z_pixelFlags_clipped,z_pixelFlags_cr,z_pixelFlags_sensor_edge,z_pixelFlags_suspect,z_epoch,
#y_cModelMag,y_cModelMagErr,y_cModelFlux,y_cModelFluxErr,y_cModel_flag,
#y_psfMag,y_psfMagErr,y_psfFlux,y_psfFluxErr,y_psfFlux_flag,y_pixelFlags_edge,y_pixelFlags_offimage,
#y_pixelFlags_saturated,y_pixelFlags_bad,y_pixelFlags_clipped,y_pixelFlags_cr,y_pixelFlags_sensor_edge,y_pixelFlags_suspect,y_epoch

colNames=['objectId', 'coord_ra', 'coord_dec', 'coord_raErr', 'coord_decErr', 'coord_ra_dec_Cov',
 'u_cModelMag','u_cModelMagErr','u_cModelFlux', 'u_cModelFluxErr', 'u_cModel_flag',
 'u_psfMag', 'u_psfMagErr', 'u_psfFlux', 'u_psfFluxErr', 'u_psfFlux_flag', 'u_pixelFlags_edge', 'u_pixelFlags_offimage',
 'u_pixelFlags_saturated', 'u_pixelFlags_bad', 'u_pixelFlags_clipped', 'u_pixelFlags_cr', 'u_pixelFlags_sensor_edge', 'u_pixelFlags_suspect', 'u_epoch',
 'g_cModelMag','g_cModelMagErr','g_cModelFlux', 'g_cModelFluxErr', 'g_cModel_flag',
 'g_psfMag', 'g_psfMagErr', 'g_psfFlux', 'g_psfFluxErr', 'g_psfFlux_flag', 'g_pixelFlags_edge', 'g_pixelFlags_offimage',
 'g_pixelFlags_saturated', 'g_pixelFlags_bad', 'g_pixelFlags_clipped', 'g_pixelFlags_cr', 'g_pixelFlags_sensor_edge', 'g_pixelFlags_suspect', 'g_epoch',
 'r_cModelMag','r_cModelMagErr','r_cModelFlux', 'r_cModelFluxErr', 'r_cModel_flag',
 'r_psfMag', 'r_psfMagErr', 'r_psfFlux', 'r_psfFluxErr', 'r_psfFlux_flag', 'r_pixelFlags_edge', 'r_pixelFlags_offimage',
 'r_pixelFlags_saturated', 'r_pixelFlags_bad', 'r_pixelFlags_clipped', 'r_pixelFlags_cr', 'r_pixelFlags_sensor_edge', 'r_pixelFlags_suspect', 'r_epoch',
 'i_cModelMag','i_cModelMagErr','i_cModelFlux', 'i_cModelFluxErr', 'i_cModel_flag',
 'i_psfMag', 'i_psfMagErr', 'i_psfFlux', 'i_psfFluxErr', 'i_psfFlux_flag', 'i_pixelFlags_edge', 'i_pixelFlags_offimage',
 'i_pixelFlags_saturated', 'i_pixelFlags_bad', 'i_pixelFlags_clipped', 'i_pixelFlags_cr', 'i_pixelFlags_sensor_edge', 'i_pixelFlags_suspect', 'i_epoch',
 'z_cModelMag','z_cModelMagErr','z_cModelFlux', 'z_cModelFluxErr', 'z_cModel_flag',
 'z_psfMag', 'z_psfMagErr', 'z_psfFlux', 'z_psfFluxErr', 'z_psfFlux_flag', 'z_pixelFlags_edge', 'z_pixelFlags_offimage',
 'z_pixelFlags_saturated', 'z_pixelFlags_bad', 'z_pixelFlags_clipped', 'z_pixelFlags_cr', 'z_pixelFlags_sensor_edge', 'z_pixelFlags_suspect', 'z_epoch',
          'y_cModelMag','y_cModelMagErr','y_cModelFlux', 'y_cModelFluxErr', 'y_cModel_flag',         
 'y_psfMag', 'y_psfMagErr', 'y_psfFlux', 'y_psfFluxErr', 'y_psfFlux_flag', 'y_pixelFlags_edge', 'y_pixelFlags_offimage',
 'y_pixelFlags_saturated', 'y_pixelFlags_bad', 'y_pixelFlags_clipped', 'y_pixelFlags_cr', 'y_pixelFlags_sensor_edge', 'y_pixelFlags_suspect', 'y_epoch']




colTypes = {}

for col in colNames:

    if col == "objectId":
        colTypes[col] = "int64"

    elif col.startswith("coord_"):
        colTypes[col] = "float64"

    elif col.endswith("Flux") or col.endswith("FluxErr"):
        colTypes[col] = "float64"

    elif col.endswith("Mag") or col.endswith("MagErr"):
        colTypes[col] = "float64"

    elif col.endswith("_flag") or "pixelFlags_" in col:
        colTypes[col] = "boolean"

    elif col.endswith("_epoch"):
        colTypes[col] = "float64"   # <-- corrected

    else:
        colTypes[col] = "float64"   # safe fallback for DP1 numeric columns

# merge_measurement and merge_footprint have no nulls so can use bool


print("colTypes = " + pformat(colTypes))

#sys.exit()



chunkNum=0
chunkSize=1000000
skipRows=chunkNum*chunkSize
rowsWritten=0
nside=256

bestIndexCounts=np.zeros(numBands,dtype=np.int32)

possFlags=[2,1,0]
minRA=999
maxRA=-999
minDec=999
maxDec=-999
minMag=[999,999,999,999]
maxMag=[-999,-999,-999,-999]
minSig=999
maxSig=-999
reasons=np.zeros(shape=(10,3), dtype=np.int32)



#print skipRows
logger.info("READINGCHUNK="+str(chunkNum))
#print colNames
#print colTypes

aC = 0.0059898
bC = 8.817481E-12
mC = 7.618399
zpes=[0.0028,0.0028,0.0028]
const=-2.5/math.log(10.0)


for csv_chunk in pd.read_csv('dp1_combined_mags.csv', na_values=['', 'nan', 'NaN', 'inf', '-inf'],header=0,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):
    boolean_columns = csv_chunk.select_dtypes(include='boolean').columns
   # print(boolean_columns)
    #csv_chunk[boolean_columns] =  csv_chunk[boolean_columns].applymap(lambda x: True if pd.api.types.isna(x) else x)
    csv_chunk[boolean_columns] = csv_chunk[boolean_columns].fillna(True)


    logger.info("STARTINGCHUNK="+str(chunkNum))
    goodRows=[]
    badRows=[]
    chunkNum += 1
    print (chunkNum)
#    print(csv_chunk.shape)
    #print (csv_chunk.dtypes)
 
#    csv_chunk['g_psfflux_flag'] = csv_chunk['g_psfflux_flag'].replace({pd.NA:True}).astype('boolean')
#    csv_chunk['g_psfflux_flag'] = csv_chunk['g_psfflux_flag'].astype('boolean')
#    mapping_dict = {'True': True, 'False': False, '': pd.NA}
#    csv_chunk['g_psfflux_flag'] = csv_chunk['g_psfflux_flag'].map(mapping_dict)
    #print (csv_chunk.index)
#    vals=csv_chunk.loc[:,('g_psfflux_flag')]
    ra_values=csv_chunk.loc[:,('coord_ra')].values
    #print ('TYPE',type(ra_values))
    dec_values=csv_chunk.loc[:,('coord_dec')].values
    logger.info("STARTINGLB="+str(chunkNum))
    l_values,b_values=convert_to_galactic(ra_values, dec_values)
    logger.info("ENDLB="+str(chunkNum))

    
    #print(l_values,b_values)
    rn=-1
#    sdsscentroid_ra_idx = csv_chunk.columns.get_loc('g_sdsscentroid_ra')
    for row in csv_chunk.itertuples(index=False):
   #     print(row)

#        print('ROW',row.g_psfflux_flag,row.merge_footprint_g)

#        if row.g_psfflux_flag==row.merge_footprint_g:
#            print('ROW MATCHES')

        rn +=1
        isPlane=False
        if abs(b_values[rn]) <= latLimit:
            isPlane=True
            numInPlane=numInPlane+1
        useMags=[]
        useMagErrs=[]
        useMagFlags=[]

        if isPlane:
            useMags = [getattr(row, f"{f}_psfMag") for f in filters]
            useMagErrs = [getattr(row, f"{f}_psfMagErr") for f in filters]
            useFluxes = [getattr(row, f"{f}_psfFlux") for f in filters]
            useFluxErrs=[getattr(row, f"{f}_psfFluxErr") for f in filters]
            useMagFlags=[getattr(row, f"{f}_psfFlux_flag") for f in filters]
        else :
            useMags = [getattr(row, f"{f}_cModelMag") for f in filters]
            useMagErrs = [getattr(row, f"{f}_cModelMagErr") for f in filters]            
            useFluxes = [getattr(row, f"{f}_cModelFlux") for f in filters]
            useFluxErrs=[getattr(row, f"{f}_cModelFluxErr") for f in filters]
            useMagFlags=[getattr(row, f"{f}_cModel_flag") for f in filters]

    #    print(useFluxes)
    #    print(useFluxErrs)
    #    print(useMagFlags)

        useMagsFF = [
            (flux * u.nJy).to(u.ABmag).value
            if flux is not None and np.isfinite(flux) and flux > 0
            else np.nan
            for flux in useFluxes
        ]
        print(useMags)
        print(useMagsFF)


        
   #     print(useMags)


        useSNRs = [
            flux / err
            if (flux is not None and err is not None
                and np.isfinite(flux) and np.isfinite(err)
                and err > 0)
            else np.nan
            for flux, err in zip(useFluxes, useFluxErrs)
        ]
    #    print(useSNRs)

 #       useMagErrs = [
 #           (2.5 / np.log(10)) * (err / flux)
 #           if (flux is not None and err is not None
 #               and np.isfinite(flux) and np.isfinite(err)
 #               and flux > 0 and err >= 0)
 #           else np.nan
 #           for flux, err in zip(useFluxes, useFluxErrs)
 #       ]

        useEpochs = [getattr(row, f"{f}_epoch") for f in filters]
     #   print(useMagErrs)

        #covar = row.sigradec * abs(row.sigradec)
        #  'coord_ra', 'coord_dec', 'coord_raErr', 'coord_decErr', 'coord_ra_dec_Cov',
        sig_maj = math.sqrt(0.5 * (row.coord_raErr**2 + row.coord_decErr**2 + math.sqrt((row.coord_raErr**2 - row.coord_decErr**2)**2 + 4 * row.coord_ra_dec_Cov**2)))
        sig_min = math.sqrt(0.5 * (row.coord_raErr**2 + row.coord_decErr**2 - math.sqrt((row.coord_raErr**2 - row.coord_decErr**2)**2 + 4 * row.coord_ra_dec_Cov**2)))
        sig_avg = 3600.0*math.sqrt(sig_maj*sig_min)
      #  print(sig_avg)

        useMags_copy = useMags[:]
        if rn%100000==0 :
            logger.info("ROWNUM="+str(rn)+" "+str(numInPlane))

        bestPhotBand,maxPhotFlag,mags,photFlags,photReasons,magErrs,snrs=flagging(row,useMags,useMagErrs,useMagFlags,rn,useSNRs)
        for i, reason in enumerate(photReasons):
            if reason is not None:
                reasonCounts[i][reason] += 1
                
 #       print(bestPhotBand,maxPhotFlag,mags,photFlags,photReasons,magErrs)
 

#        sys.exit()
        #print('CSV,',ra_values[rn],',',dec_values[rn],',',l_values[rn],',',b_values[rn])


        
        include=False
        bestEpoch=defaultEpoch
        if bestPhotBand >=0 :
            include=True
            bestEpoch=useEpochs[bestPhotBand]
            wband = filters[bestPhotBand]    # convert index → band name ('u','g',...)
            bestBandCounts[wband] += 1
        if include:
            goodRows.append([row.objectId,row.coord_ra,row.coord_dec,sig_avg,mags[0],mags[1],mags[2],mags[3],mags[4],mags[5],bestPhotBand,bestEpoch,
                             snrs[0],snrs[1],snrs[2],snrs[3],snrs[4],snrs[5],
                             magErrs[0],magErrs[1],magErrs[2],magErrs[3],magErrs[4],magErrs[5],maxPhotFlag,photFlags[0],photFlags[1],photFlags[2],photFlags[3],photFlags[4],photFlags[5]])
        else :
            # reset mags
            badRows.append([row.objectId,row.coord_ra,row.coord_dec,sig_avg,useMags_copy[0],useMags_copy[1],useMags_copy[2],useMags_copy[3],useMags_copy[4],useMags_copy[5],bestPhotBand,bestEpoch,
                             snrs[0],snrs[1],snrs[2],snrs[3],snrs[4],snrs[5],
                             magErrs[0],magErrs[1],magErrs[2],magErrs[3],magErrs[4],magErrs[5],maxPhotFlag,photFlags[0],photFlags[1],photFlags[2],photFlags[3],photFlags[4],photFlags[5]])
            
#        print("RA ",row[sdsscentroid_ra_idx])

#        if pd.isna(row[sdsscentroid_ra_idx]):
#            print("NOOOOOOOOOOOO")
        
#        print(row[18]," row ",rn,type(row[18]))
        
#    for val in vals:
#        print ("valueee",val,type(val))
#        if isinstance(val, np.bool_):
#            print("here")
#            if val ==np.bool_(True):
#                print("is a boolean with value True.")
#            if val ==np.bool_(False):
#                print("is a boolean with value Falseeeee.")


    outCols=['objectId','ra','dec','sig_avg','mag0','mag1','mag2','mag3','mag4','mag5','bestIndex','epoch',
             'snr0','snr1','snr2','snr3','snr4','snr5',
             'magErr0','magErr1','magErr2','magErr3','magErr4','magErr5',\
             'bestFlag','pf0','pf1','pf2','pf3','pf4','pf5']                            
    goodDF=pd.DataFrame(goodRows,columns=outCols)
    badDF=pd.DataFrame(badRows,columns=outCols)
    goodDF=goodDF.round({'ra': 9, 'dec': 8, 'sig_avg': 8, 'mag0': 4, 'mag1': 4, 'mag2': 4, 'mag3': 4,'mag4': 4, 'mag5': 4, \
                         'snr0': 4,'snr1': 4, 'snr2': 4,'snr3': 4, 'snr4': 4, 'snr5': 4,
                         'magErr0': 4, 'magErr1': 4, 'magErr2': 4, 'magErr3': 4,'magErr4': 4, 'magErr4': 5})
    badDF=badDF.round({'ra': 9, 'dec': 8, 'sig_avg': 8, 'mag0': 4, 'mag1': 4, 'mag2': 4, 'mag3': 4,'mag4': 4, 'mag5': 4, \
                         'snr0': 4,'snr1': 4, 'snr2': 4,'snr3': 4, 'snr4': 4, 'snr5': 4,
                         'magErr0': 4, 'magErr1': 4, 'magErr2': 4, 'magErr3': 4,'magErr4': 4, 'magErr4': 5})


    logger.info("MINMAX b "+str(min(b_values))+",  "+str(max(b_values)))    
    logger.info("WRITINGCHUNK="+str(chunkNum))
    goodDF.to_csv('dp1Good_10000.csv', mode='a', header=False, index=False)
    badDF.to_csv('dp1Bad_10000.csv', mode='a', header=False, index=False)
    logger.info("WRITTENCHUNK="+str(chunkNum))

for i, flt in enumerate(filters):
    print(f"\nFilter: {flt}")
    for reason, count in sorted(reasonCounts[i].items()):
        print(f"  Reason {reason}: {count}")
    logger.info(f"\nFilter: {flt}")
    for reason, count in sorted(reasonCounts[i].items()):
        logger.info(f"  Reason {reason}: {count}")

print(bestBandCounts)
logger.info(f"\nBest band counts: {bestBandCounts}")
