# -*- coding: utf-8 -*-

"""
================================================================================
Pytheas: an open-source software solution for local shear-wave splitting studies
================================================================================

Pytheas is a tool that aims to introduce a new mentality in shear-wave splitting analysis 
from local recordings, incorporating manual, semi- and fully- automatic methods under one 
Graphical User Interface. Integrating databases and offering compatibility with popular data 
and metadata file formats, Pytheas is designed with the simplification of analysis in mind, 
while, at the same time, enhanching the effectiveness of processing and quality control of results.

Pytheas is released under the GNU GPLv3 license.

Authors: Spingos I. & Kaviris G. (c) 2019
Special thanks to Millas C. for testing the software and providing valuable feedback from the 
very early stages of this endeavor!

For any issues, comments or suggestions please contact us at ispingos@geol.uoa.gr or through GitHub 
at https://www.github.com/ispingos/pytheas-splitting

"""
####################################################################
#                                                                  #
# This module includes Eigenvalues related functions.              #
#                                                                  #
####################################################################

## imports
import sys, os, logging, time, itertools
import numpy as np
from matplotlib import pyplot as plt
from obspy.signal.util import next_pow_2
from scipy.signal import correlate as xc
## function to get inverse of Cumulative Distribution Function (CDF)
from scipy.stats import f as fdistr
## pytheas
from lib import tools

## definitions ##
def SilverAndChan(stream1,bazi,pickwin,ain=False,maxDelay=0.250,method='EV'):
    """
    Apply the Silver and Chan (1991) method for determing shear-wave splitting.
    Use either 'EV' for searching for the minimum second eigenvalue or 
    'ME' to determine the optimal solution from the minimum energy in the corrected
    transverse.

    :type stream1: :class: `~obspy.core.stream.Stream`
    :param stream1: the stream object containg the waveforms
    :type bazi: float
    :param bazi: the station's backazimuth
    :type pickwin: tuple-like
    :param pickwin: a list of the signal window start and end (relative to
        the S-arrival in s).
    :type ain: float or bool, optional
    :param ain: the incidence angle at the station. If False or 0, the waveforms
        will be rotated to ZRT instead of LQT. Defaults to False.
    :type maxDelay: float, optional
    :param maxDelay: the maximum time-delay for the tests (in s). Defaults to 0.250.
    :type method: str, optional
    :param method: the method to use, either Eigenvalues ('EV') or Minimum Energy ('ME').
        Defaults to 'EV'.
    :returns: the following parameters; 
        - the estimated phi
        - the estimated time-delay
        - the array of either the second eigenvalue (for 'EV') or the minimum energy ('ME')
          for each phi/time-delay pair
        - the array containing all the phi values that trials were performed for
        - the array containing all the time-delay values that trials were performed for
        - the polarization direction of the corrected shear-wave at the NE axial system
        - the uncertainty of phi
        - the uncertainty of time-delay
        - the 95% confidence interval for either the second eigenvalue or the minimum energy
        - the correlation coefficient between the corrected horizontal channels in the FS system
        - the correlation coefficient between the corrected horizontal channels in the NE system
        - the normalized variance of the corrected transverse component
        - the number of 95% interval contours (TODO: does not work properly for now!)

    ..note:
        Silver, P.G., Chan, W.W., 1991. 
            Shear Wave Splitting and Sub continental Mantle Deformation.
            J. Geophys. Res. 96, 429–454. doi: 10.1029/91JB00899    

    """
    assert pickwin[1]-pickwin[0] > maxDelay, "Maximum time-delay larger than picked window!"
    stream=stream1.copy()
    logging.info("-------------------------------------------")
    logging.info("Initiating Silver and Chan (1991) method...")
    logging.info("Maximum time lag (s): "+str(maxDelay))
    ## get basic characteristics of data
    sps=stream[0].stats.sampling_rate
    delta=1./sps
    # convert maximum time-delay from seconds to samples
    maxDelay*=sps
    ## make the phi and td trial arrays
    step=1 # step for both phi and td
    phiTrials=np.arange(-90,90+step,step,dtype=int)
    tdTrials=np.arange(0,maxDelay+step,step,dtype=int)
    logging.debug("phi trials range is: [%i,%i]" % (min(phiTrials),max(phiTrials)))
    logging.debug("td trials range is : [%i,%i]" % (min(tdTrials),max(tdTrials)))
    ## make the arrays where the l2 and energy results will be stored. these are
    ## phi x td arrays
    phiNum=phiTrials.size
    tdNum=tdTrials.size
    Earray=np.zeros((phiNum,tdNum))
    L2array=np.zeros((phiNum,tdNum))
    L1array=np.zeros((phiNum,tdNum))
    V1array=np.zeros((phiNum,tdNum,2))
    ## rotate data to either ZRT or LQT
    ZNE=np.asarray((
         stream.select(component="Z")[0].data, 
         stream.select(component="N")[0].data,
         stream.select(component="E")[0].data
                  ))
    if ain:
        M=tools.Rmatrix3D(bazi,ain)
        ZEN=ZNE.copy()
        ZEN[1]=ZNE[2].copy() # swap N<->E
        ZEN[2]=ZNE[1].copy()
        LQT=np.dot(M,ZEN)
        #LQT=np.asarray([x.data for x in stream.rotate("LQT->ZNE",bazi,ain)]) # obspy gives different L component
        QT=LQT[1:]
    else:
        M=tools.Rmatrix2D(bazi)
        QT=np.dot(M,ZNE[1:])
    ## get selection window in samples
    selWindow=np.asarray((int(pickwin[0]*sps),int(pickwin[1]*sps)))
    ## start iterating
    logging.info("Start iterations for phi/td couples...")
    # do the tests #
    for i,j in itertools.product(range(phiNum),range(tdNum)):
        phi=phiTrials[i]; td=tdTrials[j]
        delWindow=selWindow+td
        # clockwise rotation matrix, according to IRIS services
        R=tools.Rmatrix2D(phi)
        # convert QT to FS
        FS1=np.dot(R,QT)
        # get the windows after applying time-delay
        FS2=np.asarray((
            FS1[0][selWindow[0]:selWindow[1]], # this is the fast
            FS1[1][delWindow[0]:delWindow[1]]  # this is the slow (after time-delay)
                      ))
        # rotate back to QT
        QT2=np.dot(R.T,FS2)
        # get minmum energy in transverse
        Earray[i,j]=np.sum(QT2[1]**2)
        # get the covariance matrix and the corresponding eigenvalues
        C=np.cov(QT2[0],QT2[1])
        w,v=np.linalg.eig(C) # w=eigenvalues, v=eigenvectors. v[:,i] corresponds to w[i]
        L2array[i,j]=w.min()
        L1array[i,j]=w.max()
        V1array[i,j]=v[:,np.where(w==w.max())[0][0]]
    ## get indices for both ME and EV
    Larray=L2array/L1array
    Larray=L2array
    if method == "EV":
        EVidx=np.where(Larray==Larray.min())
        cArray=Larray
        if len(EVidx) > 2:
            raise IndexError("More than 1 solutions specified!")
        V1=V1array[EVidx][0]
        phiSC=phiTrials[EVidx[0]][0]
        tdSC=tdTrials[EVidx[1]][0]
    elif method == "ME":
        MEidx=np.where(Earray==Earray.min())
        cArray=Earray
        if len(MEidx) > 2:
            raise IndexError("More than 1 solutions specified!")
        V1=V1array[MEidx][0]            
        phiSC=phiTrials[MEidx[0]][0]
        tdSC=tdTrials[MEidx[1]][0]
    logging.info("Initial SC solution is : %i,%i" % (phiSC,tdSC))
    ## phi must be converted back to ZNE
    phiCorr=bazi
    ## calculate errors ##
    ## correct waveforms
    R=tools.Rmatrix2D(phiSC)
    delWindow=selWindow+int(tdSC)
    # convert QT to FS
    FSe1=np.dot(R,QT)
    # get the windows after applying time-delay
    FSe2=np.asarray((
        FSe1[0][selWindow[0]:selWindow[1]], # this is the fast
        FSe1[1][delWindow[0]:delWindow[1]]  # this is the slow (after time-delay)
                  ))
    # rotate back to QT
    QTe=np.dot(R.T,FSe2)
    # get variance of signal in transverse
    #Tvar=tools.variance(QTe[1])
    Tvar=np.var(QTe[1]/abs(QTe[1]).max())
    logging.info("Normalized variance of corrected transverse signal: %.3E" % Tvar)
    # get number of degrees of freedom from transverse
    N=NDF(QTe[1])
    K=2 # N of parameters. here parameters are phi and td
    if N <= K:
       logging.warning("Number of Degrees of Freedom (%i) equal or less than number of parameters (2)" % N)
       if method == "EV":
         c95=cArray[EVidx[0],EVidx[1]]
       elif method == "ME":
         c95=cArray[MEidx[0],MEidx[1]]
    else:
       c95=confRegion(cArray.min(),N,K=K)
    cArrayNorm=cArray/c95
    # get the final error bounds
    try:
       phiSS,tdSS,nContours=getBounds(tdSC,phiSC,tdTrials,phiTrials,cArrayNorm)
       if phiSS == 0.:
         phiSS=np.nan
       if tdSS == 0.:
         tdSS=np.nan
    except:
       logging.exception("Could not perform error estimation!")
       phiSS=np.nan; tdSS=np.nan
       nContours=0
    ##########################################################       
    ## correct phi and td, finalize
    phiSC+=phiCorr
    phiSC%=180 # make sure it is in the 0-180 range
    ## get corrected in NE and correlate
    R=tools.Rmatrix2D(phiSC)
    delWindow=selWindow+int(tdSC)
    FSi=np.dot(R,ZNE[1:])
    FSc=np.asarray((
        FSi[0][selWindow[0]:selWindow[1]], # this is the fast
        FSi[1][delWindow[0]:delWindow[1]]  # this is the slow (after time-delay)
                  ))   
    CC_FS=tools.xcStatic(FSc[0],FSc[1])
    NEc=np.dot(R.T,FSc)
    CC_NE=tools.xcStatic(NEc[0],NEc[1])
    logging.info("Correlation coefficient for corrected FS is: %.2f" % CC_FS)
    logging.info("Correlation coefficient for corrected NE is: %.2f" % CC_NE)
    # get covariance of corrected NE
    C=np.cov(NEc[0],NEc[1])
    w,v=np.linalg.eig(C) # w=eigenvalues, v=eigenvectors. v[:,i] corresponds to w[i]
    V1=v[:,np.where(w==w.max())[0][0]]    
    # get polarization direction of the corrected S-wave
    spol=np.rad2deg(np.arctan2(V1[1],V1[0])) % 360
    logging.info("Corrected S-wave polarization: %.1f" % spol)
    ## td must be converted to seconds
    tdSC*=delta
    tdSS*=delta
    logging.info("Final SC solution is: %.1f+-%.1f,%.3f+-%.3f" % (phiSC,phiSS,tdSC,tdSS))
    logging.info("-------------------------------------------")
    return np.float(phiSC),np.float(tdSC),cArray,phiTrials,tdTrials,np.float(spol),\
           np.float(phiSS),np.float(tdSS),np.float(c95), \
           np.float(CC_FS),np.float(CC_NE),np.float(Tvar),nContours

def confRegion(cMin,N,K=2,a=0.05):
    """
    Get the 95% confidence region for the critical value (i.e.
    the second eigenvalue or the minimum energy).

    :type cMin: float
    :param cMin: the minimum critical value
    :type N: int
    :param N: the degrees of freedom
    :type K: int, optional
    :param K: the number of parameters (i.e. the phi and
        time-delay). Defaults to 2.
    :type a: float, optional
    :param a: the confidence level (1-a). Defaults to 0.05.
    :returns: the 100*(1-a) confidence interval

    """
    return cMin*(1+((K/(N-K))*fdistr.isf(a,K,N-K)))

def NDF(data):
    """
    Compute the degrees of freedom of the signal. Corrections
    according to Walsh et al. (2013) are used.

    :type data: array-like
    :param data: the data of the signal.
    :returns: the degrees of freedom.

    ..: note:
        Walsh, E., Arnold, R., Savage, M.K., 2013. 
            Silver and Chan revisited. 
            J. Geophys. Res. Solid Earth 118, 5500–5515. doi: 10.1002/jgrb.50386

    """
    # get number of points
    npts=data.size
    # get the magnitude of the spectrum
    nfft=next_pow_2(npts)
    A=abs(np.fft.fft(data,nfft))
    # get the moments
    F2=(A**2).sum()-((A[0]**2)-(A[-1]**2))/2
    #F4=(A**4).sum()-((A[0]**4)-(A[-1]**4))/2 # Splitlab
    #F4=(A**4).sum()-((A[0]**4)-(A[-1]**4))/3 # Silver and Chan (1991)
    F4=(4/3)*(A**4).sum()-((A[0]**4)-(A[-1]**4))/3 # Walsh et al. (2013) (28)
    #F4=(A**4).sum()-(A[0]**4+A[-1]**4)/2 # Walsh et al. (2013) (33)
    # get number of degrees of freedom
    N=int(round((2*(F2**2)/F4)-1))
    #N=int(round((2*(F2**2)/F4)-1)+0.5) # splitlab
    if N > npts:
        raise ValueError("Degrees of freedom of signal (%i) greater than npts (%i)"%(N,npts))
    return N

def getBounds(td,phi,tdTest,phiTest,cArrayNorm):
    """
    Get possible contours that correspond to the 95% C.I. l2/min energy value
    of the Silver and Chan (1991) test and then get the required error bounds.
    
    :type td: int
    :param td: the time-delay (in samples).
    :type phi: float
    :param phi: the Sfast polarization direction.
    :type tdTest: array-like
    :param tdTest: the array including time-delay values used in
        the grid search.
    :type phiTest: array-like
    :param phiTest: the array inclding phi values used in the grid
        search.
    :type cArrayNorm: the array of the critical value normalized to its
        95% C.I.
    :returns: (uncertainty of phi, uncertainty of time-delay, 
               number of 95% C.I. contours)

    """
    # minor corrections. td must ALWAYS be given in samples!!!
    if td < 0:
        raise ValueError("Time-delays must always be given in samples!")
    # get contours and td/phi lines
    cs=plt.contour(tdTest,phiTest,cArrayNorm,levels=[1.])
    plt.close()    # get contour lines coords
    selContour=[]
    for i,_ in enumerate(cs.collections):
        for j,pt in enumerate(cs.collections[i].get_paths()):
            if pt.contains_point((td,phi)):
                c95Contour=pt
                selContour=pt.interpolated(100).vertices
                break
    ## search for point-contour combinations
    if len(selContour) == 0:
        raise ValueError("Could not find contour that contains point")
    ## grab horizontal distances - td errors
    tdErrors=[]
    # left side
    errorLeft=selContour[selContour[:,0]<=td]
    tdErrors.append(errorLeft[np.argmin(abs(errorLeft[:,1]-phi))][0])
    # right side
    errorRight=selContour[selContour[:,0]>td]
    tdErrors.append(errorRight[np.argmin(abs(errorRight[:,1]-phi))][0])
    # grab vertical distances - phi errors
    phiErrors=[]
    # upper side
    errorUpper=selContour[selContour[:,1]>=phi]
    phiErrors.append(errorUpper[np.argmin(abs(errorUpper[:,0]-td))][1])
    # bottom side
    errorBottom=selContour[selContour[:,1]<=phi]
    phiErrors.append(errorBottom[np.argmin(abs(errorBottom[:,0]-td))][1])
    # return final errors
    if len(phiErrors) == 1:
        phiErr=abs(phiErrors[0]-phi)*0.5
    else:
        phiErr=abs(max(phiErrors)-min(phiErrors))*0.25
    if len(tdErrors) == 1:
        tdErr=abs(tdErrors[0]-td)*0.5
    else:
        tdErr=abs(max(tdErrors)-min(tdErrors))*0.25
    return phiErr,tdErr,len(cs.collections)