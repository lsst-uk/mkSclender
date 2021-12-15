import pandas as pd
import math
import logging
import os
import datetime
import healpy as hp
import numpy as np
import sys
from gaiaFunc import  correct_flux_excess_factor
from gaiaFunc import  correct_gband

def ranrm(theta):
    '''
    Normalise an angle expressed in radians into the range 0 to 2pi.
    '''
    w = math.fmod(theta, 2.0 * math.pi)
    if w < 0.0: w = w + 2.0 * math.pi
    return w
   
def tangent_plane_to_spherical(xi, eta, raz, dcz):
    '''
    Takes the local plane coordinates (xi, eta) and de-projects back onto the celestial sphere using the given
    tangent point origin (raz, dcz). Algorithm cribbed from Pat Wallaces' dtp2s.c in SLALIB.
   
    Returns spherical coordinates RA, Dec.
   
    All arguments and results in RADIANS.
    '''
    sdecz = math.sin (dcz)
    cdecz = math.cos (dcz)
    denom = cdecz - eta * sdecz
    ra = ranrm ( math.atan2 ( xi, denom ) + raz )
    dc = math.atan2 ( sdecz + (eta * cdecz), math.sqrt( (xi * xi) + (denom * denom) ))
    return ra, dc

deg2rad = math.pi / 180.0
mas2rad = deg2rad / (3600.0 * 1000.0)

def correct(ra, dc, pmra, pmdc, tdelta):
    '''
    Proper motion corrects the coordinates [degrees] by the given proper motion components [mas/yr] and
    time [years].
   
    Returns proper motion corrected coordinates [degrees].
    '''
   
    # Offset from a (0,0) origin for local plane coordinates:
    xi = tdelta * pmra * mas2rad
    eta = tdelta * pmdc * mas2rad
   
    # project back onto celestial sphere given the tangent point (answer in rad)
    rar, dcr =  tangent_plane_to_spherical(xi, eta, ra * deg2rad, dc * deg2rad)
   
    # return coords in original units (degrees)
    return rar / deg2rad, dcr / deg2rad


#def correct_gband (bp_rp , astrometric_params_solved , phot_g_mean_mag , phot_g_mean_flux ):    
def correct_gbandBck(bp_rp , astrometric_params_solved , phot_g_mean_mag):    
    """
    Correct the G-band fluxes and magnitudes for the input list of Gaia EDR3 data.
    Parameters
    ----------
    bp_rp: float , array_like
    The (BP-RP) colour listed in the Gaia EDR3 archive.
    astrometric_params_solved: int, array_like
    The astrometric solution type listed in the Gaia EDR3 archive.
    phot_g_mean_mag: float , array_like
    The G-band magnitude as listed in the Gaia EDR3 archive.
    phot_g_mean_flux: float , array_like
    The G-band flux as listed in the Gaia EDR3 archive.
    Returns
    -------
    The corrected G-band magnitudes and fluxes. The corrections are only applied to
    sources with a 6-parameter astrometric solution fainter than G=13, for which a
    (BP-RP) colour is available.
    Example
    gmag_corr , gflux_corr = correct_gband(bp_rp , astrometric_params_solved ,
    phot_g_mean_mag , phot_g_mean_flux)
    """
    if np.isscalar(bp_rp) or np.isscalar(astrometric_params_solved) or \
                    np.isscalar(phot_g_mean_mag): # or np.isscalar(phot_g_mean_flux):
        bp_rp = np.float64(bp_rp)
        astrometric_params_solved = np.int64(astrometric_params_solved)
        phot_g_mean_mag = np.float64(phot_g_mean_mag)
#        phot_g_mean_flux = np.float64(phot_g_mean_flux)

    if not (bp_rp.shape == astrometric_params_solved.shape \
                        == phot_g_mean_mag.shape):# == phot_g_mean_flux.shape):
        raise ValueError('Function parameters must be of the same shape')

    # added  np.isnan(phot_g_mean_mag)
    print phot_g_mean_mag,bp_rp,astrometric_params_solved
    do_not_correct = np.isnan(phot_g_mean_mag) | np.isnan(bp_rp) | (phot_g_mean_mag<=13) | \
                                       (astrometric_params_solved != 95)
    print do_not_correct
    bright_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>13) & \
                                   (phot_g_mean_mag<=16)
    faint_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>16)
    bp_rp_c = np.clip(bp_rp, 0.25, 3.0)

    correction_factor = np.ones_like(phot_g_mean_mag)
    correction_factor[faint_correct] = 1.00525 - 0.02323*bp_rp_c[faint_correct] + \
        0.01740*np.power(bp_rp_c[faint_correct],2) - \
        0.00253*np.power(bp_rp_c[faint_correct],3)
    correction_factor[bright_correct] = 1.00876 - 0.02540*bp_rp_c[bright_correct] + \
        0.01747*np.power(bp_rp_c[bright_correct],2) - \
        0.00277*np.power(bp_rp_c[bright_correct],3)
    
    gmag_corrected = phot_g_mean_mag - 2.5*np.log10(correction_factor)
#    gflux_corrected = phot_g_mean_flux * correction_factor

    return gmag_corrected #, gflux_corrected





logger = logging.getLogger('myLog')
formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                              "%Y-%m-%d %H:%M:%S")
logger.setLevel(logging.INFO )
fH     = logging.FileHandler(os.path.join(
            'logs',
            "gaiaLog_"+datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S')+'.log'))
fH.setFormatter(formatter)
logger.addHandler(fH)
logger.info("STARTINGRUN=")
# ratings_df = pd.read_csv('gaia_10000.csv', header=None)

#declare @cmd varchar(8000) set @cmd = 'bcp "select '+
#'s.source_id,designation,ra,dec,ra_error,dec_error,ra_dec_corr,pmra,pmdec,phot_g_mean_mag,astrometric_chi2_al,astrometric_n_good_obs_al,'+
#'astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_matched_transits,'+
#'phot_rp_mean_mag,phot_bp_mean_mag,ruwe,phot_bp_rp_excess_factor,bp_rp,phot_bp_rp_excess_factor,'+
#'phot_bp_mean_flux_over_error,phot_g_mean_flux_over_error,phot_rp_mean_flux_over_error,'+ 
#'phot_bp_n_obs,phot_g_n_obs,phot_rp_n_obs,'+
#'phot_bp_n_blended_transits,phot_rp_n_blended_transits,astrometric_params_solved,b '+
#'from gaia_source as s " ' +
#' queryout G:\share\mar\gaiaSkinnyishEDR3.csv -c '+
#' -S ramses20 -d GAIAEDR3 -U sa -P ibis1243! -t , -o G:\share\mar\gaiaEDR3.log ' select @cmd exec xp_cmdshell @cmd



colNames=['source_id','designation','ra','dec','ra_error','dec_error','ra_dec_corr','pmra','pmdec',\
          'phot_g_mean_mag','astrometric_chi2_al','astrometric_n_good_obs_al', \
          'astrometric_excess_noise','astrometric_excess_noise_sig','astrometric_matched_transits', \
          'phot_rp_mean_mag','phot_bp_mean_mag','ruwe','phot_bp_rp_excess_factor','bp_rp','phot_bp_rp_excess_factor2',\
          'phot_bp_mean_flux_over_error','phot_g_mean_flux_over_error','phot_rp_mean_flux_over_error', \
          'phot_bp_n_obs','phot_g_n_obs','phot_rp_n_obs', \
          'phot_bp_n_blended_transits','phot_rp_n_blended_transits','astrometric_params_solved','b']

# astrometric_params_solved , phot_g_mean_mag , phot_g_mean_flux
colTypes={'source_id':int,'designation':object}
#colNames=['cntr','designation','ra','dec','sigra','sigdec','sigradec','cc_flags','ph_qual','det_bit','var_flg','nb','moon_lev','w1mpro','w1sigmpro','w1sat','w1rchi2','w2mpro','w2sigmpro','w2sat','w2rchi2','w3mpro','w3sigmpro','w3sat','w3rchi2','w4mpro','w4sigmpro','w4sat','w4rchi2']

#for csv_chunk in pd.read_csv('gaia_10000.csv', header=None,chunksize=1300,names=colNames):
#dtype = {'moon_lev':object,'size':float,
#                          'width':object,'color':object,
#                          'qty':int,'code': object})

numBands=3
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



print skipRows
logger.info("READINGCHUNK="+str(chunkNum))
print colNames
print colTypes

aC = 0.0059898
bC = 8.817481E-12
mC = 7.618399
zpes=[0.0028,0.0028,0.0028]
const=-2.5/math.log(10.0)


for csv_chunk in pd.read_csv('/disk82/panstarrs/WP3.11/gaiaSkinnyishEDR3.csv', header=None,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):
#gaiaSkinnyishEDR3.csv gaiaSkinnyish5000.csv

    #for csv_chunk in pd.read_csv('gaia_10000.csv', header=None,chunksize=chunkSize,names=colNames,dtype=colTypes,skiprows=skipRows):

    logger.info("STARTINGCHUNK="+str(chunkNum))
    goodRows=[]
    badRows=[]
    chunkNum += 1
    print chunkNum
#    print(csv_chunk.shape)
    print (csv_chunk.dtypes)


    print csv_chunk.index

#    csv_chunk.loc[0+(chunkNum-1)*chunkSize]

#    print  csv_chunk.loc[0,'phot_bp_rp_excess_factor'],csv_chunk.loc[0,'phot_bp_rp_excess_factor2']

#    for index, row in csv_chunk.iterrows():
#    for index, row in enumerate(csv_chunk.itertuples(),1):
    radecs = csv_chunk.loc[:,('ra','dec')].values
    bp_rp=csv_chunk.loc[:,('bp_rp')].values
    astrometric_params_solved=csv_chunk.loc[:,('astrometric_params_solved')].values
    phot_g_mean_mag=csv_chunk.loc[:,('phot_g_mean_mag')].values
    b=csv_chunk.loc[:,('b')].values
    print "here 100"
    logger.info("STARTINGPEF=")
    phot_bp_rp_excess_factor=csv_chunk.loc[:,('phot_bp_rp_excess_factor')].values
    corr_phot_bp_rp_excess_factor=correct_flux_excess_factor(bp_rp,phot_bp_rp_excess_factor)
    corr_phot_g_mean_mag=correct_gband(bp_rp,astrometric_params_solved,phot_g_mean_mag)
    print "here 125"
    print phot_g_mean_mag
    print corr_phot_g_mean_mag
    sigPEF=np.full_like(phot_bp_rp_excess_factor, np.nan)
    sigPEF=aC+bC*np.power(phot_g_mean_mag,mC)
    flagPEF = np.zeros_like(phot_bp_rp_excess_factor,dtype=np.bool_)
    print "here 150"
    gCopy=np.copy(phot_g_mean_mag)
    gCopy[np.isnan(gCopy)]=-999.99999
    sig3Test=np.absolute(corr_phot_bp_rp_excess_factor)-3*sigPEF
    sig3Test[np.isnan(sig3Test)]=-999.99999   
    subsetter = (np.absolute(b)<20) & (~np.isnan(sigPEF)) & (gCopy > 4) & (~np.isnan(corr_phot_bp_rp_excess_factor)) & (sig3Test >=0 )
    flagPEF[subsetter]=1
    print "here 200"
    numFlagPEF=np.count_nonzero(flagPEF)
    logger.info("ENDPEF="+str(numFlagPEF)+" out of "+str(flagPEF.shape))

#    print radecs[:,0]
#    print radecs[:,1]
#    a=np.array([16,float('NaN'),13])
#    b=np.array([95,95,95])
#    c=np.array([0.1,0.1,0.1])
#    phot_g_mean_mags_corr=correct_gband(c, b ,a)
#    print phot_g_mean_mags_corr

#    print bp_rp,"bp_rp"
#    print astrometric_params_solved,"astrometric_params_solved"
#    print phot_g_mean_mag,"phot_g_mean_mag"
#    print phot_bp_rp_excess_factor,"phot_bp_rp_excess_factor"
#    print corr_phot_bp_rp_excess_factor,"corr_phot_bp_rp_excess_factor"
#    print sigPEF,"sigPEF"
#    print sig3Test,"sig3Test"
#    print flagPEF,"flagPEF"

#    phot_g_mean_mags_corr=correct_gband(bp_rps , astrometric_params_solveds , phot_g_mean_mags)

#    print phot_g_mean_mags_corr

    hpi=hp.ang2pix(nside, radecs[:,0],radecs[:,1],lonlat=True, nest=True)
    counts=np.bincount(hpi,minlength=npix)
    print counts
    print counts.shape, counts.dtype, hp.nside2npix(256)
    print counts[529707]
    totCounts=totCounts+counts
    print totCounts
    rn=-1


    
    for row in csv_chunk.itertuples(index=False):
#        print row,row.source_id,"\n"
        rn += 1
        include = True
        best_flag=0
        sig_avg=-999
        ra=-999
        dec=-999
        flag=[2,2,2]
        conds=[0,0,0]
        mags=[row.phot_bp_mean_mag,corr_phot_g_mean_mag[rn],row.phot_rp_mean_mag]
        pbts=[row.phot_bp_n_blended_transits,-999,row.phot_rp_n_blended_transits]
        pnos=[row.phot_bp_n_obs,-999,row.phot_rp_n_obs]
        pmfoes=[row.phot_bp_mean_flux_over_error,row.phot_g_mean_flux_over_error,row.phot_rp_mean_flux_over_error]
#        print "PMFOES",pmfoes

        try:
            sig_maj = math.sqrt(0.5*(row.ra_error**2 + row.dec_error**2 + math.sqrt((row.ra_error**2 - row.dec_error**2)**2 \
                                                                                          + 4*(row.ra_dec_corr*row.ra_error*row.dec_error)**2)))
            sig_min = math.sqrt(0.5*(row.ra_error**2 + row.dec_error**2 - math.sqrt((row.ra_error**2 - row.dec_error**2)**2 \
                                                                                          + 4*(row.ra_dec_corr*row.ra_error*row.dec_error)**2)))
            sig_avg = (sig_min+sig_maj)/2.0
            

            ra=row.ra
            dec=row.dec

            if (row.ra < 360.0 and row.ra >= 0 and math.isnan(row.pmra) is False and \
                row.dec <= 90.0 and row.dec >= -90.0 and math.isnan(row.pmdec) is False):
                tra,tdec=correct(row.ra,row.dec,row.pmra,row.pmdec,-0.5)
            else:
                tra=row.ra
                tdec=row.dec

#            print row.ra, row.dec, row.pmra,row.pmdec, tra,tdec
#            print "pbts ",pbts
#            print "pnos ",pnos

            for i in range(0, numBands):
#                print i
#                print row.designation[i]
#                print  pbts[i]/pnos[i] # divide bt zero
                if (row.ruwe >= 5.0) :
                    #flag[i]=-999
                    flag=[-999,-999,-999]
                    conds=[1,1,1]
#                    print reasons[1,:]
                    reasons[1,:] +=1
                    include=False
#                    print "reason",reasons[1,:]
                    break
                elif (math.isnan(mags[i])):
                    flag[i]=-1
                    conds[i]=2
#                    print  "reason","mags",i,flag[i]
                    reasons[2,i] +=1
                elif (i==0 and mags[i] > 20.3):
                    flag[i]=-1
                    conds[i]=3
                    reasons[3,i] +=1
#                    print  "reason","mags gt 20.3",i,flag[i]
                elif (i!=1 and not math.isnan(pbts[i]) and not math.isnan(pnos[i]) and pbts[i]/pnos[i] >= 0.1):
                    flag[i]=0
                    conds[i]=4
                    reasons[4,i] +=1
#                    print  "reason pbts/pnos",i,flag[i]
                elif (i!=1 and flagPEF[rn]):
                    flag[i]=0
                    conds[i]=5
                    reasons[5,i] +=1
#                    print  "reason flagPEF",i,flag[i],"sig3test",sig3Test[rn],rn
                elif (not math.isnan(pmfoes[i]) and math.sqrt( (const/pmfoes[i])**2+zpes[i]**2 ) > 0.5  ):
                    flag[i]=1
                    conds[i]=6
                    reasons[6,i] +=1
#                    print  "reason ZPES",math.sqrt( (const/pmfoes[i])**2+zpes[i]**2) ,rn,i,flag[i]
                elif  ( (row.astrometric_excess_noise > 1 and row.astrometric_excess_noise_sig <= 2) or \
                    row.astrometric_n_good_obs_al < 60 or row.astrometric_matched_transits <= 10  ):
                    flag[i]=1
                    conds[i]=7
                    reasons[7,i] +=1
#                    print  "reason AST",i,flag[i],rn

                if (flag[i]==-999 or flag[i]==-1):
                    mags[i]=float("NaN")
                elif (i==0 and not math.isnan(mags[0]) and not math.isnan(mags[1]) and 2.0 < mags[1] < 3.94):
                    mags[0]=mags[0] - 0.9921 - 0.02598 * mags[1] + 0.1833 * mags[1]**2 - 0.02862 * mags[1]**3
                    #BP - 0.9921 - 0.02598 * G + 0.1833 * G2 - 0.02862 * G3  (2 < G < 3.94
                elif (i==1 and not math.isnan(mags[1]) and 2.0 < mags[1] < 8.0):
                    mags[1]=mags[1] - 0.09892 + 0.059 * mags[1] - 0.009775 * mags[1]**2 + 0.0004934 * mags[1]**3
                    #G - 0.09892 + 0.059 * G - 0.009775 * G2 + 0.0004934 * G3 (2 < G < 8)
                elif (i==2 and not math.isnan(mags[2]) and  2.0 < mags[2] < 3.45):
                    mags[2]=mags[2] - 14.94 + 14.41 * mags[2] - 4.657 * mags[2]**2 + 0.503 * mags[2]**3
                    #RP - 14.94 + 14.41 * RP - 4.657 * RP2 + 0.503 * RP3 (2 < RP < 3.45)
            
#            print "ZPES", math.sqrt( (const/pmfoes[i])**2+zpes[i]**2) ,rn,i

#            print "flags ",flag
#            print "checkExcess ",rn,flagPEF[rn],corr_phot_bp_rp_excess_factor[rn],sigPEF[rn]

#            if (rn >14 ):
#                sys.exit()
            
#            if (row.ra < 360.0 and row.ra > 0 and math.isnan(row.pmra) is False) :
#                ra = row.ra - 5 * row.pmra * 1e-3 / 3600 / math.cos(math.radians(row.dec))
#                if (ra > 360.0) : ra=ra-360.0
#                if (ra < 0.0) : ra=ra+360.0
#            if (row.dec <= 90.0 and row.dec >= -90.0 and math.isnan(row.pmdec) is False) :
#                dec = row.dec - 5 * row.pmdec * 1e-3 / 3600
#            flag=2 

                
#            if  ( (row.astrometric_excess_noise > 1 and row.astrometric_excess_noise_sig <= 2) \
#                   or row.astrometric_n_good_obs_al < 60 or row.astrometric_matched_observations <= 8  ):
#                flag=1

#            if ( math.sqrt(row.astrometric_chi2_al /( row.astrometric_n_good_obs_al - 5)) > \
#                 1.2 * max(1,math.exp(-0.2*( row.phot_g_mean_mag - 19.5)))    ) :
#                flag=-999

                
 
        except Exception as e:
            print e
            logger.error("ERROR: "+str(e)+" "+str(row.source_id))
            include=False
            for j in range(0, numBands):
                flag[j]=-999

        bestFlag=-999
        bestIndex=-999
        for possFlag in possFlags:
            for j in range(2, -1, -1):
                if (flag[j] == possFlag):
                    bestFlag=possFlag
                    bestIndex=j
                    bestIndexCounts[j] += 1
                    break
            if (bestFlag >=0):
                break
            
        if (bestFlag< 0) :
            include=False
            
        if (include):
            goodRows.append([row.source_id,row.designation,tra,tdec,sig_avg,mags[0],mags[1],mags[2],bestIndex,hpi[rn],bestFlag,flag[0],flag[1],flag[2],row.ra,row.dec,row.phot_bp_mean_mag,row.phot_g_mean_mag,row.phot_rp_mean_mag])
            rowsWritten += 1
            goodCounts[hpi[rn]] += 1
        else :
            badRows.append([row.source_id,row.designation,tra,tdec,sig_avg,mags[0],mags[1],mags[2],bestIndex,hpi[rn],bestFlag,flag[0],flag[1],flag[2],row.ra,row.dec,row.phot_bp_mean_mag,row.phot_g_mean_mag,row.phot_rp_mean_mag])
            badCounts[hpi[rn]] += 1
    #goodDF=pd.DataFrame(goodRows)
    outCols=['source_id','designation','tra','tdec','sig_avg','mag0','mag1','mag2','bestIndex','hpi','bestFlag','flag0','flag1','flag2','ra','dec','omag0','omag1','omag2']
    goodDF=pd.DataFrame(goodRows,columns=outCols)
    goodDF=goodDF.round({'tra': 9, 'tdec': 8, 'sig_avg': 8, 'mag0': 4, 'mag1': 4, 'mag2': 4, 'ra': 9, 'dec': 8, 'omag0': 4, 'omag1': 4, 'omag2': 4})
    badDF=pd.DataFrame(badRows,columns=outCols)
    badDF=badDF.round({'tra': 9, 'tdec': 8, 'sig_avg': 8, 'mag0': 4, 'mag1': 4, 'mag2': 4, 'ra': 9, 'dec': 8, 'omag0': 4, 'omag1': 4, 'omag2': 4})
    logger.info("WRITINGCHUNK="+str(chunkNum))
    goodDF.to_csv('/unas/mar/gaiaGood.csv', mode='a', header=False, index=False)
    badDF.to_csv('/unas/mar/gaiaBad.csv', mode='a', header=False, index=False)
    logger.info("WRITTENCHUNK="+str(chunkNum))
    
#    for csv_row in csv_chunk.rows:
#        print row.ra,row.dec
#        print row
print "reasons"
print reasons
print "saving counts"
print bestIndexCounts
np.savetxt("/unas/mar/gaiaTotCounts.csv", totCounts, delimiter=",")
np.savetxt("/unas/mar/gaiaGoodCounts.csv", goodCounts, delimiter=",")
np.savetxt("/unas/mar/gaiaBadCounts.csv", badCounts, delimiter=",")

