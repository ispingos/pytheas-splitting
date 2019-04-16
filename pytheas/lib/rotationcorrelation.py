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
# This module includes Cluster Analysis related functions/classes. #
#                                                                  #
####################################################################

## imports
import logging
import itertools
import numpy as np
from matplotlib import pyplot as plt
import scipy.stats
## pytheas 
from lib import tools
from lib import eigenvalue as EV

def crosscorrelation(stream1,pickwin,maxDelay=0.250):
    """
    Rotation - correlation method, as seen in Bowman and Ando (1989).

    :type stream1: :class: `~obspy.core.stream.Stream`
    :param stream1: the stream object containg the waveforms
    :type pickwin: tuple-like
    :param pickwin: a list of the signal window start and end (relative to
        the S-arrival in s).
    :type maxDelay: float, optional
    :param maxDelay: the maximum time-delay for the tests (in s). Defaults to 0.250.
    :returns: the following parameters; 
        - the estimated phi
        - the estimated time-delay
        - the array of the correlation coefficient for each phi/time-delay pair
        - the array containing all the phi values that trials were performed for
        - the array containing all the time-delay values that trials were performed for
        - the polarization direction of the corrected shear-wave at the NE axial system
        - the uncertainty of phi
        - the uncertainty of time-delay
        - the 95% confidence interval for either the correlation coefficient
        - the correlation coefficient between the corrected horizontal channels in the FS system
        - the correlation coefficient between the corrected horizontal channels in the NE system
        - the normalized variance of the corrected transverse component (always set to 'numpy.nan', not used)
        - the number of 95% interval contours (TODO: does not work properly for now!)    

    """
    assert pickwin[1]-pickwin[0] > maxDelay, "Maximum time-delay larger than picked window!"
    logging.debug('Initiating cross-correlation...')
    logging.debug('Maximum time lag is %s'%maxDelay)
    # grab traces
    stream=stream1.copy()
    ## get basic characteristics of data
    sps=stream[0].stats.sampling_rate
    delta=1./sps
    # convert maximum time-delay from seconds to samples
    maxDelay*=sps
    maxDelay=int(np.ceil(maxDelay))
    ## make the phi and td trial arrays
    step=1 # step for both phi and td
    phiTrials=np.arange(0,180+step,step,dtype=int)
    tdTrials=np.arange(0,maxDelay+step,step,dtype=int)
    logging.debug("phi trials range is: [%i,%i]" % (min(phiTrials),max(phiTrials)))
    logging.debug("td trials range is : [%i,%i]" % (min(tdTrials),max(tdTrials)))
    ## make the arrays where the CC
    # phi x td arrays
    phiNum=phiTrials.size
    tdNum=tdTrials.size
    cArray=np.zeros((phiNum,tdNum))
    ## get dta arrays
    ZNE=np.asarray((
         stream.select(component="Z")[0].data, 
         stream.select(component="N")[0].data,
         stream.select(component="E")[0].data
                  ))
    NE=ZNE[1:]
    nTrace=stream.select(component="N")[0]
    eTrace=stream.select(component="E")[0]
    ## get selection window in samples
    selWindow=np.asarray((int(pickwin[0]*sps),int(pickwin[1]*sps)))
    ## start iterating
    logging.info("Start iterations for phi/td couples...")
    for i,j in itertools.product(range(phiNum),range(tdNum)):
        phi=phiTrials[i]; td=tdTrials[j]
        delWindow=selWindow+td
        # clockwise rotation matrix, according to IRIS services
        R=tools.Rmatrix2D(phi)
        # convert QT to FS
        FS1=np.dot(R,NE)
        # get the windows after applying time-delay
        FS2=np.asarray((
            FS1[0][selWindow[0]:selWindow[1]], # this is the fast
            FS1[1][delWindow[0]:delWindow[1]]  # this is the slow (after time-delay)
                      ))
        # get CC
        cArray[i,j]=tools.xcStatic(FS2[0],FS2[1])
    cArray=abs(cArray) # we want anticorrelated too
    CCmax=abs(cArray).max()
    RCidx=np.where(abs(cArray)==CCmax)
    if len(RCidx) > 2:
        raise ValueError("More than one optimal td-phi pairs found")    
    phiRC=phiTrials[RCidx[0]][0]
    tdRC=tdTrials[RCidx[1]][0]            
    logging.info('Maximum CC is %.3f' % abs(cArray).max())
    ##########################################################
    ## calculate errors ##
    ##########################################################
    c95=lowCI(CCmax,0.05,cArray.size-1)
    cArrayNorm=cArray/c95
    try:
       phiRS,tdRS,nContours=EV.getBounds(tdRC,phiRC,tdTrials,phiTrials,cArrayNorm)
       if phiRS == 0.:
         phiRS=np.nan
       if tdRS == 0.:
         tdRS=np.nan
    except:
       logging.exception("Could not perform error estimation!")
       nContours=0
       phiRS=np.nan; tdRS=np.nan    ## correct phi and td, finalize
    phiRC%=180 # make sure it is in the 0-180 range
    ## get corrected in NE and correlate
    R=tools.Rmatrix2D(phiRC)
    delWindow=selWindow+int(tdRC)
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
    ## get isotropic polarization
    # get covariance of corrected NE
    C=np.cov(NEc[0],NEc[1])
    w,v=np.linalg.eig(C) # w=eigenvalues, v=eigenvectors. v[:,i] corresponds to w[i]
    V1=v[:,np.where(w==w.max())[0][0]]    
    spol=np.rad2deg(np.arctan2(V1[1],V1[0])) % 360
    ## compatibility hack
    Tvar=np.nan
    ## td must be converted to seconds
    tdRC*=delta
    tdRS*=delta
    logging.info("Final RC solution is: %.1f+-%.1f,%.3f+-%.3f" % (phiRC,phiRS,tdRC,tdRS))     
    return np.float(phiRC),np.float(tdRC),cArray,phiTrials,tdTrials,\
           np.float(spol),np.float(phiRS),np.float(tdRS),\
           np.float(c95),np.float(CC_FS),np.float(CC_NE),np.float(Tvar),nContours

def r2z(r):
    """
    Convert from R-space to CC-space

    :type r: float
    :param r: the critical value in the R space
    :returns: the critical value in the CC-space

    """
    return np.log((1+r)/(1-r))/2.0

def z2r(z):
    """
    Convert from CC-space to R-space

    :type z: float
    :param z: the critical value in the CC space
    :returns: the critical value in the R-space    

    """    
    e=np.exp(2*z)
    return((e-1)/(e+1))

def lowCI(r,a,n):
    """
    Calculate the lower Confidence Interval bound for
    the given value.

    :type r: float
    :param r: the critical value in the R-space
    :type a: float
    :param a: the confidence level (1-a).
    :type n: int
    :param n: the number of samples
    :returns: the lower bound of the 100*(1-a) confidence
        interval.

    """
    z=r2z(r)
    se=1.0/np.sqrt(n-3)
    z_crit=scipy.stats.norm.ppf(1-a)
    # get the upper and lower CI
    lo=z-z_crit*se
    hi=z+z_crit*se
    return z2r(lo) # returns only the lower CI