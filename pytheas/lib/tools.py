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
# This module includes a wide variety of complementary tools.      #
#                                                                  #
####################################################################

## imports
import os, sys, logging
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from obspy import UTCDateTime
from obspy.signal.trigger import ar_pick
from obspy.signal.cross_correlation import correlate
import scipy.signal

## archiving
def initSplittingDict(layer):
    """
    Initialize and return an empty splitting 
    dictionary based on the layer.

    :type layer: str
    :param layer: the layer id for which to init the dict.
        Accepts 'event','station' or 'method'.
    :returns: the splitting dict

    """
    if layer == "event":
        return {"origin":None,"latitude":np.nan,"longitude":np.nan,
                "depth":np.nan,"magnitude":np.nan}
    elif layer == "station":
        return {"station":None,"epicentral":np.nan,"azimuth":np.nan,
                "s_p":np.nan,"incidence":np.nan,"s_obs":np.nan,
                "s_theo":np.nan,"s_auto":np.nan,"orientation":np.nan}
    elif layer == "method":
        return {"station":None,"method":None,"phi":np.nan,"td":np.nan,
                "pol":np.nan,"CC_NE":np.nan,"var_T":np.nan,
                "err_phi":np.nan,"err_td":np.nan,"N_contours":np.nan,
                "grade":np.nan,"filter_min":np.nan,"filter_max":np.nan,
                "window_min":np.nan,"window_max":np.nan,"SNR":np.nan,
                "s_freq":np.nan,"comment":"",
                "C_array":np.empty((1,1)),"C_max":"DOUBLE","phi_test":np.empty((1,1)),"td_test":np.empty((1,1)),
                "initial_clusters":np.empty((1,1)),"initial_clusters_errors":np.empty((1,1)),"calinski_score":np.empty((1,1)),
                "clusters1":np.empty((1,1)),"clusters2":np.empty((1,1))}
    else:
        raise ValueError("Argument must be one of 'event','station' or 'method'")


def initPytheasDict(layer):
    """
    Initiate an evDict or staDict type dictionary

    :type layer: str
    :param layer: the layer code. Either 'event' or 'station'.
    :returns: a dict containing the required fields for either
        station or event information
    """
    if layer == "event":
        return {'YEAR':None,'Mo':None,'Da':None,'HR':None,
                'MN':None,'SEC':None, 'ORIGIN':None,
                'LAT':None,'LON':None,'DEPTH':None,
                'MAG':None,'evID':None}
    elif layer == "station":
        return {'STA':None,'NET':None,'DIST':None,'BAZ':None,'AN':None,
                'SECP':None,'TOBSP':None,'SECS':None,'TOBSS':None,}

## time-delay related functions
def lengthcheck(Wx,Wy):
    """
    Checks the start time and end time of two traces and clips them
    to the same start and end times.

    :type Wx: :class: `~obspy.core.trace.Trace`
    :param Wx: a trace object hosting the first trace.
    :type Wy: :class: `~obspy.core.trace.Trace`
    :param Wy: a trace object hosting the second trace.
    :returns: (first trace's corrected data array, second trace's corrected data array)

    """
    codex=str(Wx).split("|")[0].strip()
    codey=str(Wy).split("|")[0].strip()
    logging.info("Checking length for waveforms x:"+codex+" and y:"+codey)
    if Wx.stats.starttime != Wy.stats.starttime:
        sTime=min((Wx.stats.starttime,Wy.stats.starttime))
        Wx.trim(starttime=sTime,pad=True,fill_value=0,nearest_sample=False)
        Wy.trim(starttime=sTime,pad=True,fill_value=0,nearest_sample=False)
        logging.debug("Corrected starttime to "+str(sTime))
    if Wx.stats.endtime != Wy.stats.endtime:
        eTime=max((Wx.stats.endtime,Wy.stats.endtime))
        Wx.trim(endtime=eTime,pad=True,fill_value=0,nearest_sample=False)
        Wy.trim(endtime=eTime,pad=True,fill_value=0,nearest_sample=False)
        logging.debug("Corrected endtime to "+str(eTime))
    # if one sample different, the waveforms are considered synced
    sizeDiff=Wx.data.size-Wy.data.size
    if sizeDiff == 1:
        Wx.data=Wx.data[:-2]
    elif sizeDiff == -1:
        Wy.data=Wy.data[:-2]
    elif sizeDiff != 0:
        raise IndexError("Cannot sync waveforms with sizes (%i,%i)" % (Wx.data.size,Wy.data.size))
    return Wx.data,Wy.data

def timedelay(td,stream):
    """
    Adds the time delay to the slow (transverse) component.
    Pretty much adds 0s at the start of the slow component, according
    to the time delay and sampling rate.

    if np array is provided instead of stream, td must be in samples
    and stream must be an array containing ONLY the transverse.
    
    :type td: float
    :param td: the time-delay
    :type stream: :class: `~obspy.core.stream.Stream` or array-like
    :param stream: the data where the time-delay will be added to. Either 
        a stream or an array must be provided. In the first case, the stream
        must contain a 'T' component, where the time-delay will be added to. If
        an array is provide, it must contain the signal of the waveform to which
        the time-delay will be added.
    :returns: the data array of the trace with the added time-delay.

    """
    if not isinstance(stream,np.ndarray):
        trace=stream.select(component="T")[0].copy()
        sps=trace.stats.sampling_rate
        npts=trace.stats.npts
        nshift=abs(int(td*sps))
        data1=trace.data
    else:
        data1=stream
        npts=data1.size
        nshift=abs(td)
    if td > 0:
        data2=np.insert(data1,0,np.zeros(nshift))
        data2=data2[:npts]
    elif td < 0:
        data2=np.insert(data1,data1.size,np.zeros(nshift))
        data2=data2[nshift:]
    else:
        data2=data1
    return data2

## PM method related
def polarigram(stream,rotated,filterRange):
    """
    Makes the polarigram by using two horizontal channels.
    For short signals, no need to select additional window. If required,
    will reconsider.

    :type stream: :class: `~obspy.core.stream.Stream`
    :param stream: the stream containing the waveforms
    :type rotated: bool
    :param rotated: If True, the waveforms are considered to be
        in the FS system (i.e. 'R' and 'T' components will be requested).
    :type filterRange: tuple-like
    :param filterRange: the minimum and maximum bounds of the applied filter.
    :returns: (list containing the start and end point coordinates for each vector,
               dict containing the amplitudes and angles of the vectors)

    """
    stream=stream.copy()
    if rotated:
        Wy=stream.select(channel="??R")[0]
        logging.debug("POLAR RADIAL:     "+str(Wy))
        Wx=stream.select(channel="??T")[0]
        logging.debug("POLAR TRANSVERSE: "+str(Wx))
    else:
        Wy=stream.select(channel="??N")[0]
        logging.debug("POLAR NORTH: "+str(Wy))
        Wx=stream.select(channel="??E")[0]
        logging.debug("POLAR EAST:  "+str(Wx))
    # Filter?
    if filterRange != (np.nan,np.nan):
        Wy.filter(type="bandpass",freqmin=filterRange[0],freqmax=filterRange[1],zerophase=True,corners=4)
        Wx.filter(type="bandpass",freqmin=filterRange[0],freqmax=filterRange[1],zerophase=True,corners=4)
    if Wy.stats.sampling_rate != Wx.stats.sampling_rate:
        raise ValueError("Sampling rates of horizontal channels must match!")
    ## Check length of timeseries
    Ax,Ay=lengthcheck(Wx,Wy)
    Wx.data=Ax
    Wy.data=Ay
    if Wy.stats.npts > Wx.stats.npts:
        T=np.arange(0,Wy.stats.npts/Wy.stats.sampling_rate,Wy.stats.delta)
    else:
        T=np.arange(0,Wx.stats.npts/Wx.stats.sampling_rate,Wx.stats.delta)
    # Array of vector amplitudes
    A=np.sqrt((Ax**2)+(Ay**2))
    # Array of vector angles
    D=np.rad2deg(np.arctan2(Ay,Ax))
    # Make the polarigram lines
    pairList=[]
    magn=max(A)*10 # Needs more testing?
    n=0
    for t,a,d in zip(T,A,D):
        x1=t
        x2=t+(a/magn)*np.cos(np.deg2rad(d))
        y1=0
        y2=(a/magn)*np.sin(np.deg2rad(d))
        pair=((x1,y1),(x2,y2))
        pairList.append(pair)
    polarDict={"TIME":T,"AMP":A,"ANGLE":D,"Ax":Ax,"Ay":Ay}             
    return pairList, polarDict

def hodogram(A,D,norm=False):
    """
    Creates the pairs used in hodogram plotting. Input is the amplitude and
    angles arrays. Position of elements in arrays must be in correct chronological order.

    :type A: array-like
    :param A: array containing the Amplitudes of the vectors acquired from `~pytheas.tools.polarigram`.
    :type D: array-like
    :param D: array containing the Angles of the vectors acquired from `~pytheas.tools.polarigram`.    
    :returns: list containing the start and end point coordinates for each vector

    """
    pairList=[]
    magn=max(A) 
    x1=0
    y1=0
    for a,d in zip(A,D):
        x2=x1+(a/magn)*np.cos(np.deg2rad(d))
        y2=y1+(a/magn)*np.sin(np.deg2rad(d))
        pair=((x1,y1),(x2,y2))
        pairList.append(pair)
        x1=x2
        y1=y2
    pairList=np.asarray(pairList)
    # normalize?
    if norm:
        pairList=pairList/abs(pairList).max()
    return pairList

## auto grading routines
def autoGrading(null,values,bounds,grade_dict,
                sel_dict={"phi":True,"td":True,"SNR":False,"FS":True,"NE":True},weights=None):
    """
    Grade a measurement according to provided criteria

    :type null: tuple-like
    :param null: Contains the values of phi, pol and polarization offset for estimating whether 
        the measurement is null or not. This is decided by finding whether phi and pol are either 
        (sub)parallel or (sub)perpendicular.
    :type values: tuple-like
    :param values: Contains the values to use in the grading, i.e. the error of phi, error of
        td, SNR, correlation coefficient for the FS system and correlation coeficient for the NE
        system.
    :type bounds: tuple-like
    :param bounds: Contains the user-set boundaries for each value, as denoted in `values`.
    :type grade_dict: dict
    :param grade_dict: The (closed) upper bound limits for each accepted grade (``'A'``,``'B'``,
        ``'C'``,``'D'`` and ``'E'``).
    :type sel_dict: dict, optional
    :param sel_dict: Define which of the criteria (``'phi'``,``'td'``,``'SNR'``,``'FS'`` and 
        ``'NE'``) to use in the grading algorithm. Defaults to ``'True'`` for all criteria.
    :type weights: None or tuple-like, optional
    :param weights: Specify a weight to use for each criterion (same order as defined in `values`
        and `bounds`). If ``None``, all weights will be equal. Defaults to ``None``
    :returns: (score as float, grade as str)

    """
    if weights is None:
        weights=np.ones(len(values))
    if len(weights) != len(values):
        weights=np.ones(len(values))
        logging.warning("Wrong number of weights specified, switching to equal weights.")
    #
    logging.debug("Grading phi/pol: %.1f/%.1f" % (null[0],null[1]))
    ## null check
    diff=abs(null[0] - (null[1] % 180))
    clause=np.logical_or(
           np.logical_and(null[2] <= diff, diff <= 90-null[2]),
           np.logical_and(90+null[2] <= diff, diff <= 180-null[2])
                        )
    if not clause:
        logging.debug("Measurement is null")
        return np.nan, "N"
    ## scoring. abs is used mainly for getting the absolute value of CC factors
    scores=[]
    if sel_dict["phi"]:
        pScore=weights[0]*(values[0]/bounds[0])
        scores.append(pScore)
    if sel_dict["td"]:
        tScore=weights[1]*(values[1]/bounds[1])
        scores.append(tScore)
    if sel_dict["SNR"]:
        sScore=weights[2]*(bounds[2]/values[2])
        scores.append(sScore)
    if sel_dict["FS"]:
        fScore = weights[3]*(1.0-abs(values[3]))/(1.0-bounds[3])
        scores.append(fScore)
    if sel_dict["NE"]:
        nScore=weights[4]*(1.0-abs(values[4]))/(1.0-bounds[4])
        scores.append(nScore)
    # convert to array
    scores=np.asarray(scores)
    # 
    score=scores.mean()
    logging.info("total score: %.4f" % score)
    # final checks
    if np.any(scores>=1):
        return score,'E'
    if not sel_dict["SNR"]:
        if values[2] < bounds[2]:
        	return score,'E'
    # grade
    if score <= grade_dict["A"]:
        wgt="A"
    elif grade_dict["A"] < score <= grade_dict["B"]:
        wgt="B"
    elif grade_dict["B"] < score <= grade_dict["C"]:
        wgt="C"
    elif grade_dict["C"] < score <= grade_dict["D"]:
        wgt="D"
    else:
        wgt="E"
    logging.debug("Final score/grade: %.3f/%s" % (score,wgt))
    # final correction
    if np.isinf(score): score=np.nan
    return score,wgt 

def variance(signal):
    """
    Calculate the variance of the signal

    :type signal: array-like
    :param signal: the signal for which to calculate the variance
    :returns: the variance

    """
    signal/=abs(signal).max() # first normalize
    return (np.mean(abs(signal-np.mean(signal)))**2)/(signal.size-1)

def xcStatic(s1,s2):
    """
    Simple correlation of two signals for null time-shift. Uses
    `~obspy.signal.cross_correlation.correlate`

    :type s1: array-like
    :param s1: the first signal in the correlation
    :type s2: array-like
    :param s2: the second signal in the correlation
    :returns: the correlation coefficient 
    """
    return correlate(s1,s2,shift=0,demean=True,normalize='naive')[0]

## other functions ##
def errorsDelPezzo(D,td,std,sD=0.2):
    """
    Calculate the error of the normalized time-delay
    according to Del Pezzo et al., 2004). This is not
    currently used.

    :type D: float
    :param D: the epicentral distance in km
    :type td: float
    :param td: the time-delay in ms
    :type std: float
    :param std: the error of the time-delay in ms
    :type sD: float, optional
    :param sD: the error of the epicentral distance in km.
        Defaults to ``0.2``.
    :returns: error of normalized td

    ..note:
        Del Pezzo, E., Bianco, F., Petrosino, S., Saccorotti, G., 2004.
            Changes in coda decay rate and shear-wave splitting parameters 
            associated with seismic swarms at Mt. Vesuvius, Italy. 
            Bull Seismol Soc Am 94(2):439–452
    """
    return np.sqrt((D**-2)*(std**2)+(D**-4)*(td**2)*(sD**2))

def arPicker(stream,pkCNF):
    """
    Calculate the arrivals of local P and S phases with the
    AR-AIC implementation (Akazawa, 2004) of Obspy in
    `~obspy.signal.trigger.ar_pick`
   
    :type stream: s
    :param stream: The stream containing the waveforms (3 channels)
        where the arrivals will be determined.
    :type pkCNF: :class:`~pytheas.lib.parsers.parsePickerCnf`
    :param pkCNF: Configuration class containing the required settings
        for the arrival determination.
    :returns: (p arrival relative to the start time in s, 
        s arrival relative to the start time in s)

    ..note:
        Akazawa, T., 2004.
            A technique for automatic detection of onset time of P-and S-Phases 
            in strong motion records, 13th World Conference on Earthquake Engineering.

    """
    # some params
    df=stream[0].stats.sampling_rate
    z,n,e=stream
    # load from cnf file
    p,s=ar_pick( # p arrival is not used
                z.data,n.data,e.data,df,
                pkCNF.arFMIN,pkCNF.arFMAX,pkCNF.arLTAP,pkCNF.arSTAP,
                pkCNF.arLTAS,pkCNF.arSTAS,
                pkCNF.arMP,pkCNF.arMS,pkCNF.arLP,pkCNF.arLS,
                s_pick=True
                )
    logging.info("AR-AIC arrivals (p,s) in s: "+str((p,s)))
    return p,s

def getTheorArrivals(event,station,model,phases=[]):
    """
    Calculate theoretical arrivals of phases from the TauP (Crotwell et al., 2004)
    implementation of Obspy (`~obspy.taup`).

    :type event: tuple-like
    :param event: Contains the origin time, latitude, longitude and depth of the
        event.
    :type station: tuple-like
    :param station: Contains the latitude and longitude of the stations
    :type vmod: str
    :param vmod: Specify the velocity model, either as one of the accepted ones
        in `~obspy.taup` (e.g. ``'IASP91'``) or as path to a compatible file. If an
        error occurs while attempting to read the file, the ``'IASP91'`` model will
        be used.
    :type phases: tuple-like
    :param phases: Specify which phases to get output for. If empty, all phases
        will be used. Defaults to ``[]``.
    :returns: A dictionary which includes the phase hint, angle of incidence and 
        arrival time for each arrival found.

    ..note:
        Crotwell, H.P., Owens, T.J., Ritsema, J., 1999. 
        The TauP Toolkit: Flexible Seismic Travel-time and Ray-path Utilities. 
        Seismol. Res. Lett. 70, 154–160. doi: 10.1785/gssrl.70.2.154

    """
    # get arrivals
    arrivals=model.get_travel_times_geo(
                                        source_depth_in_km=event[3],
                                        source_latitude_in_deg=event[1],
                                        source_longitude_in_deg=event[2],
                                        receiver_latitude_in_deg=station[0],
                                        receiver_longitude_in_deg=station[1]
                                        #phase_list=phases
                                       )    
    # convert to a more convenient structure
    arDict={}; i=0
    for arr in arrivals:
        i+=1
        arDict.update(
            {
              i : {
                    "phase":arr.name, "ain":arr.incident_angle,
                    "time":event[0]+arr.time

              }
            }
                   )
    logging.info("Found %i phases" % i)
    return arDict

def spectrum(trace1,trace2,sps):
    """
    Calculate the spectrum for the two horizontal traces using
    `~scipy.signal.periodogram`.

    :type trace1: :class:`~obspy.core.trace.Trace`
    :param trace1: Trace of the first horizontal channel
    :type trace2: :class:`~obspy.core.trace.Trace`
    :param trace2: Trace of the second horizontal channel
    :type sps: float or int
    :param sps: The sampling rate of the signal
    :returns: Array containing the spectra of the signals in the
        two horizontal components.

    """
    specs=[]
    for data in (trace1,trace2):
        nfft=int(2**(np.ceil(np.log2(len(data))))) # ensure nfft is power of 2
        freq,spec=scipy.signal.periodogram(data,sps,window="hann",nfft=nfft,
                                           detrend=False,return_onesided=True,
                                           scaling="spectrum")
        specs.append((freq,spec))
    return specs

def centerColl(ax,coll,axPolar=False):
    """
    Get max x and y of collection and scale axes likewise

    :type ax: :class: `~matplotlib.axes.Axes`
    :param ax: the Axes instance where centering will be applied.
    :type coll: :class: `~matplotlib.collections.LineCollection` 
    :param coll: the Collection instance for which centering will be
        applied.
    :type axPolar: bool, optional
    :param axPolar: Indicate that the Collection refers to one from
        a polarigram. Defalts to False.

    """
    segs=coll.properties()['segments']
    # first for x1
    #x=[ e[0][0] for e in segs]
    x=[]; y=[]
    for i,e in enumerate(segs):
        try:
            x.append(e[0][0])
            x.append(e[1][0])
            y.append(e[0][1])
            y.append(e[1][1])
        except IndexError:
            logging.warning("IndexError in segment "+str(i))            
    # are there no values in x,y?
    if not x:
        x=[0]
    if not y:
        y=[0]
    # find min/max and adjust axes lims
    factor=1.2 # set limits to 120% of min/max
    if not axPolar:
        xmax=max(abs(min(x)),abs(max(x)))*factor
        xmin=-xmax
        ax.set_xlim(xmin,xmax)
    ymax=max(abs(min(y)),abs(max(y)))*factor
    ymin=-ymax
    ax.set_ylim(ymin,ymax)

def qcFigMAN(streamIni,phi,td,pol,sarr,baz,ain,filterRange,titleInfo,output=False,extend=0.2):
    """
    Make the Quality Control figure for the MAN method.

    :type stream: :class: `~obspy.core.stream.Stream`
    :param stream: the waveform stream object
    :type phi: float
    :param phi: the Sfast polarization direction
    :type td: float
    :param td: the time-delay (in ms)
    :type pol: float
    :param pol: the correcetd S-wave's polarization direction
    :type sarr: float
    :param sarr: the S-arrival relative to the start time of the stream (in s).
    :type baz: float
    :param baz: the backazimuth of the station.
    :type ain: float
    :param ain: the incidence angle at the station
    :type filterRange: tuple-like
    :param filterRange: the upper and lower filter bounds (in Hz).
    :titleInfo: tuple-like
    :param titleInfo: list containing information related to the title of the figure, i.e. 
        the origin time, the incidence angle, the epicentral distance, the magnitude and
        the measurement's grade.
    :type output: str or bool, optional
    :param output: If not False, the figure will be saved at the specified filepath. Defaults
        to False.
    :type extend: float, optional
    :param extend: one-sided extension of the signal window relative to the S-arrival.
        The final window will be extend*2. Defaults to 0.2.

    """
    baz%=360
    stream=streamIni.copy()
    delta=stream[0].stats.delta
    station=stream[0].stats.station
    network=stream[0].stats.network
    td=abs(td)
    spick=(UTCDateTime(sarr)-stream[0].stats.starttime)
    winstart=int(np.floor((spick-extend)/delta))
    winend=int(np.ceil((spick+extend)/delta))
    spick=extend
    # will make a 4x3 grid figure
    fig=plt.figure("Manual Method",figsize=(7,6))
    # first column is originl NE
    ax11=fig.add_subplot(4,3,1) # unfiltered NE
    ax12=fig.add_subplot(4,3,4,sharex=ax11) # filtered NE
    ax13=fig.add_subplot(4,3,7,sharex=ax11) # polarigram NE filtered
    ax14=fig.add_subplot(4,3,10) # hodogram NE filtered
    # second column is rotated FS, td is indicated by lines
    ax21=fig.add_subplot(4,3,2,sharex=ax11) # unfiltered FS
    ax22=fig.add_subplot(4,3,5,sharex=ax11) # filtered FS
    ax23=fig.add_subplot(4,3,8,sharex=ax11) # polarigram FS filtered
    ax24=fig.add_subplot(4,3,11) # hodogram FS filtered
    # third column is corrected NE
    ax31=fig.add_subplot(4,3,3,sharex=ax11) # unfiltered NE
    ax32=fig.add_subplot(4,3,6,sharex=ax11) # filtered NE
    ax33=fig.add_subplot(4,3,9,sharex=ax11) # polarigram NE filtered
    ax34=fig.add_subplot(4,3,12) # hodogram NE filtered
    ## start processing ##
    # unfiltered Initial
    uHH=np.asarray((
         stream.select(component="N")[0].data[winstart:winend],
         stream.select(component="E")[0].data[winstart:winend]
                  ))
    ax11.plot(np.arange(uHH[0].size)*delta,uHH[0],'blue',label="North")
    ax11.plot(np.arange(uHH[1].size)*delta,uHH[1],'red',label="East")
    ax11.set_ylabel("Raw")
    ax11.set_title("Initial")
    # filtered initial
    if not all((np.isnan(filterRange[0]),np.isnan(filterRange[1]))):
        stream.filter(type="bandpass",freqmin=filterRange[0],freqmax=filterRange[1],corners=4,zerophase=True)
    fHH=np.asarray((
         stream.select(component="N")[0].data[winstart:winend],
         stream.select(component="E")[0].data[winstart:winend]
                  ))
    ax12.plot(np.arange(fHH[0].size)*delta,fHH[0],'blue',label="North")
    ax12.plot(np.arange(fHH[1].size)*delta,fHH[1],'red',label="East")
    ax12.set_ylabel("Filtered")
    # polarigram initial
    stream.select(component="N")[0].data=fHH[0]
    stream.select(component="E")[0].data=fHH[1]
    pairList,polarDict=polarigram(stream,False,(np.nan,np.nan))
    colourColl=['black' for x in pairList]
    origPolar=LineCollection(pairList,colors=colourColl)
    # hodogram initial
    origPairs=hodogram(polarDict["AMP"],polarDict["ANGLE"]) 
    # unfiltered rotated
    M=Rmatrix2D((180+phi)%180) # rotation matrix for phi
    uHH=np.dot(M,uHH)
    ax21.plot(np.arange(uHH[0].size)*delta,uHH[0],'blue',label="Fast")
    ax21.plot(np.arange(uHH[1].size)*delta,uHH[1],'red',label="Slow")
    ax21.set_title("Rotated")
    # filtered rotated
    fHH=np.dot(M,fHH)
    ax22.plot(np.arange(fHH[0].size)*delta,fHH[0],'blue',label="Fast")
    ax22.plot(np.arange(fHH[1].size)*delta,fHH[1],'red',label="Slow")
    # draw td
    for ax in (ax21,ax22):
        ax.axvline(x=spick,color="blue",linestyle="dashed",linewidth=1.0)
        ax.axvline(x=spick+td,color="red",linestyle="dashed",linewidth=1.0)
        y=0.6*ax.get_ylim()[1]
        ax.arrow(spick+td,y,-td,0,color="red",length_includes_head=True,head_length=td*0.2,head_width=y*0.15)
    # polarigram rotated
    stream.select(component="N")[0].data=fHH[0]
    stream.select(component="E")[0].data=fHH[1]
    pairList,polarDict=polarigram(stream,False,(np.nan,np.nan))
    rotPolar=LineCollection(pairList,colors=colourColl)
    # hodogram rotated
    rotPairs=hodogram(polarDict["AMP"],polarDict["ANGLE"]) 
    # unfiltered corrected
    uHH[1]=timedelay(int(np.ceil(-td/delta)),uHH[1])
    uHH=np.dot(M.T,uHH)
    ax31.plot(np.arange(uHH[0].size)*delta,uHH[0],'blue',label="North")
    ax31.plot(np.arange(uHH[1].size)*delta,uHH[1],'red',label="East")
    ax31.set_title("Corrected")
    # filtered corrected
    fHH[1]=timedelay(int(np.ceil(-td/delta)),fHH[1])
    fHH=np.dot(M.T,fHH)
    ax32.plot(np.arange(fHH[0].size)*delta,fHH[0],'blue',label="North")
    ax32.plot(np.arange(fHH[1].size)*delta,fHH[1],'red',label="East")
    # polarigram corrected
    stream.select(component="N")[0].data=fHH[0]
    stream.select(component="E")[0].data=fHH[1]
    pairList,polarDict=polarigram(stream,False,(np.nan,np.nan))
    corPolar=LineCollection(pairList,colors=colourColl)
    # hodogram corrected
    corPairs=hodogram(polarDict["AMP"],polarDict["ANGLE"]) 
    ## polarigram loop
    for axPO,polCol in zip((ax13,ax23,ax33),(origPolar,rotPolar,corPolar)):
        axPO.add_collection(polCol)
        centerColl(axPO,polCol,axPolar=True)
        axPO.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
        axPO.set_aspect(1,adjustable="datalim")
        axPO.axvline(x=spick,linewidth=1.0,linestyle="dashed",color='blue')
        axPO.set_ylabel("North")
        axPO.set_xlabel("East")
    # small correction
    ax23.set_ylabel("Fast")
    ax23.set_xlabel("Slow")
    ## hodogram loop
    # normalize and add to collections
    normval=np.max((abs(origPairs),abs(rotPairs),abs(corPairs)))
    origColl=LineCollection(origPairs/normval,label="Initial")
    rotColl=LineCollection(rotPairs/normval,label="Rotated")
    corColl=LineCollection(corPairs/normval,label="Corrected")
    ## make the hodogram for the original NE
    for axHO,col in zip((ax14,ax24,ax34),(origColl,rotColl,corColl)):
        axHO.plot(0,0,marker="*",color="red")
        axHO.add_collection(col)
        axHO.set_aspect(1,adjustable="datalim")
        axHO.set_ylim(-1.3,1.3) 
        axHO.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
        axHO.set_xlabel("East")
        axHO.set_ylabel("North")
    # small correction
    ax24.set_xlabel("Slow")
    ax24.set_ylabel("Fast")
    ##
    ax11.set_xlim(0,delta*uHH[0].size)
    # plot picks to all
    for axWf in (ax11,ax12,ax21,ax22,ax31,ax32):
        axWf.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
        if axWf not in (ax21,ax22): # picks are already plotted for corrected stage waveforms
            axWf.axvline(x=spick,linewidth=1.0,linestyle='dashed',color='blue')
    # draw legends and adjust ticks
    ax11.legend(loc="upper left",handlelength=0.6,frameon=True,fontsize=8)
    ax21.legend(loc="upper left",handlelength=0.6,frameon=True,fontsize=8)
    ax31.legend(loc="upper left",handlelength=0.6,frameon=True,fontsize=8)
    ax12.tick_params(direction="in",top=False,bottom=True,left=False,right=False,labelleft=False,labelbottom=True)
    ax22.tick_params(direction="in",top=False,bottom=True,left=False,right=False,labelleft=False,labelbottom=True)
    ax32.tick_params(direction="in",top=False,bottom=True,left=False,right=False,labelleft=False,labelbottom=True)
    # adjust ticklabels and plot xlabel
    ax22.set_xlabel('rel. time from S-arrival (s)')
    for ax in (ax12,ax22,ax32):
        ax.set_xticks([0,extend,extend*2])
        ax.set_xticklabels([-extend,0,extend],rotation=-45)
    # plot title
    origtime,ain,dist,mag,grade=titleInfo
    fig.suptitle(
                    r"$\bf%s;%s.%s$" % (str(origtime[:-4]),network,station)+"\n"+
                    r"$\bfbaz$: "+"{:.1f}".format(baz)+
                    r"N$^\circ$E, $\bfain$: "+"{:.1f}".format(ain)+
                    r"$^\circ$, $\bfepi$: "+"{:.1f}".format(dist)+" km, "+
                    r"$\bfmag:$ "+"{:.1f}".format(mag)+"\n"+
                    r"$\bfφ:$ "+"{:.1f}".format(phi)+
                    r" N$^\circ$E, $\bft_d$: "+
                    "{:.1f} ms".format(abs(td)*1000.)+r", $\bfp:$ "+"{:.1f}".format(pol)+
                    r"N$^\circ$E, "+r"$\bfgrade:$ "+grade+
                    r", $\bfflt:$ %.1f | %.1f Hz" % (filterRange[0],filterRange[1]),
                    fontsize=8

                )
    for ax in fig.get_axes():
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)
    fig.tight_layout(rect=(0.1,0.1,0.9,0.9),w_pad=0.0,h_pad=0.0)
    # save to file?
    if output:
        fig.savefig(output,dpi=720,quality=95)

def qcFigSWS(stream,phi,td,pol,errors,baz,pickwin,Cmatrix,phi_test,td_test,meth,
                     cEmax,filterRange,titleInfo,output=False):
    """
    Make the Quality Control figure for one of the EV, ME and RC methods.

    :type stream: :class: `~obspy.core.stream.Stream`
    :param stream: the waveform stream object
    :type phi: float
    :param phi: the Sfast polarization direction
    :type td: float
    :param td: the time-delay (in s)
    :type pol: float
    :param pol: the correcetd S-wave's polarization direction
    :type errors: tuple-like
    :param errors: a list of the uncertainties for phi and the time-delay.
    :type baz: float
    :param baz: the backazimuth of the station.
    :type pickwin: tuple-like
    :param pickwin: a list of the signal window start and end (relative to
        the S-arrival in s).
    :type Cmatrix: array-like
    :param Cmatrix: the array of the second eigenvalue, transverse energy or correlation 
        coefficient for each phi/time-delay pair.
    :type phi_test: array-like
    :param phi_test: the array containing all the phi values that trials were performed for
    :type td_test: array-like
    :param td_test: the array containing all the time-delays values that trials were performed for ( in s)
    :type meth: str
    :param meth: the method used in the analysis. One of 'EV', 'ME' or 'RC'.
    :type cEmax: float
    :param cEmax: the critical value for selection in Cmatrix, i.e. the maximum correlation coefficient,
        the minimum second eigenvalue or the minimum energy.
    :type filterRange: tuple-like
    :param filterRange: the upper and lower filter bounds (in Hz).
    :titleInfo: tuple-like
    :param titleInfo: list containing information related to the title of the figure, i.e. 
        the origin time, the incidence angle, the epicentral distance, the magnitude and
        the measurement's grade.
    :type output: str or bool, optional
    :param output: If not False, the figure will be saved at the specified filepath. Defaults
        to False.

    """
    ##
    baz%=360 # ensure backazimuth is in the 0-359 range
    delta=stream[0].stats.delta # period rate
    station=stream[0].stats.station
    network=stream[0].stats.network
    td=abs(td)
    ##
    # convert samples to ms
    td_test=(td_test*delta)*(10**3) # convert to ms for better display
    ## extend pick window for better display
    extend=50*delta # seconds
    pickwin=list(pickwin) # convert to list for assignments
    pickwin[0]-=extend
    pickwin[1]+=extend
    ## prepare the wf data
    scStream=stream.copy()
    if filterRange != (np.nan,np.nan):
        scStream.filter(type="bandpass",freqmin=filterRange[0],freqmax=filterRange[1],
                        zerophase=True,corners=4)
    ## correct for time delay
    Wz,Wx,Wy=scStream
    # First two
    Ax,Ay=lengthcheck(Wx,Wy)
    Wx.data=Ax 
    Wy.data=Ay
    # Rest
    Az,Ay=lengthcheck(Wz,Wy)
    Wz.data=Az
    # first NE original
    tNE=(scStream.select(component="N")[0].times(),scStream.select(component="E")[0].times())
    NE=(scStream.select(component="N")[0].data,scStream.select(component="E")[0].data)
    # now RT (FS). Quick but dirty
    scStream.rotate("NE->RT",phi+180) # R->F, T->S
    data=timedelay(-td,scStream)
    scStream.select(component="T")[0].data=data
    ## correct for time delay
    Wz,Wx,Wy=scStream
    # First two
    Ax,Ay=lengthcheck(Wx,Wy)
    Wx.data=Ax 
    Wy.data=Ay
    # Rest
    Az,Ay=lengthcheck(Wz,Wy)
    Wz.data=Az
    tFS=(scStream.select(component="R")[0].times(),scStream.select(component="T")[0].times())
    FS=(scStream.select(component="R")[0].data,scStream.select(component="T")[0].data)
    # now NE corrected
    scStream.rotate("RT->NE",phi+180)
    tcNE=(scStream.select(component="N")[0].times(),scStream.select(component="E")[0].times())
    (cN,cE)=(scStream.select(component="N")[0].data,scStream.select(component="E")[0].data)
    ## make the fig
    if meth in ["EV","ME"]:
        fig=plt.figure("Silver and Chan (1991)",figsize=(7,6))
    elif meth == "RC":
        fig=plt.figure("Rotation - Correlation",figsize=(7,6))
    #fig.set_size_inches((20.0,11.5),forward=True)
    ## setup the figure
    ax1=fig.add_subplot(5,2,1) # for NE original
    ax2=fig.add_subplot(5,2,3,sharex=ax1) # for FS corrected
    ax3=fig.add_subplot(5,2,5,sharex=ax1) # for NE corrected
    ax4=fig.add_subplot(5,2,7,autoscale_on=False) # hodogram for original NE 
    ax5=fig.add_subplot(5,2,9,autoscale_on=False) # hodogram for corrected NE
    ##
    ax1.set_ylabel("NE\n(orig)")
    ax2.set_ylabel("FS\n(corr)")
    for ax in (ax1,ax2,ax3):
        ax.set_yticklabels([""])
        ax.set_xlim(pickwin)
        ax.axvspan(pickwin[0]+extend,pickwin[1]-extend,alpha=0.3,color='green')
    ## make NE
    ax1.plot(tNE[0],NE[0],'blue',label="North")
    ax1.plot(tNE[1],NE[1],'red',label="East")
    ## make FS
    ax2.plot(tFS[0],FS[0],'blue',label="Fast")
    ax2.plot(tFS[1],FS[1],'red',label="Slow")
    ## make NE corrected
    labsWf=("North","East")
    ax3.set_ylabel("NE\n(corr)")
    ax3.plot(tcNE[0],cN,'blue',label=labsWf[0])
    ax3.plot(tcNE[1],cE,'red',label=labsWf[1])
    for axWf in (ax1,ax2,ax3):
        axWf.legend(loc="lower left",handlelength=0.6,frameon=True,fontsize=8)
        axWf.tick_params(left=False)
    # small correction
    ax1.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
    ax2.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
    ax2.set_xlabel("relative time (in s)")
    ## get the polarization vectors
    rotated=False
    ## hodogram for original NE
    scStream.select(component="N")[0].data=NE[0]
    scStream.select(component="E")[0].data=NE[1]
    _,polarDict=polarigram(scStream,rotated,filterRange)
    origPairs=hodogram(polarDict["AMP"][int((pickwin[0]+extend)/delta):int((pickwin[1]-extend)/delta)],
                       polarDict["ANGLE"][int((pickwin[0]+extend)/delta):int((pickwin[1]-extend)/delta)])    
    ## hodogram for corrected NE
    scStream.select(component="N")[0].data=cN
    scStream.select(component="E")[0].data=cE
    _,polarDict=polarigram(scStream,rotated,filterRange)
    corPairs=hodogram(polarDict["AMP"][int((pickwin[0]+extend)/delta):int((pickwin[1]-extend)/delta)],
                       polarDict["ANGLE"][int((pickwin[0]+extend)/delta):int((pickwin[1]-extend)/delta)]) 
    # normalize and add to collections
    normval=np.max((abs(origPairs),abs(corPairs)))
    origColl=LineCollection(origPairs/normval,label="Original")
    corColl=LineCollection(corPairs/normval,linestyles="dashed",label="Corrected")
    ## make the hodogram for the original NE
    ax4.plot(0,0,marker="*",color="red")
    ax4.add_collection(origColl)
    ax4.set_aspect(1,adjustable="datalim")
    ax4.set_ylim(-1.3,1.3) 
    ## make the hodogram for the corrected
    ax5.plot(0,0,marker="*",color="red")
    ax5.add_collection(corColl)
    ax5.set_aspect(1,adjustable="datalim")
    ax5.set_ylim(-1.3,1.3) 
    ## final touch
    for axH in (ax4,ax5):
        #ax4.tick_params(labelleft=False,labelbottom=False,direction="in")
        axH.tick_params(top=False,bottom=False,left=False,right=False,labelleft=False,labelbottom=False)
        axH.set_xlabel("East")
        axH.set_ylabel("North")
        axH.legend(handlelength=0.6,frameon=True,fontsize=8)
    ## title
    origtime,ain,dist,mag,grade=titleInfo
    fig.suptitle(
                    r"$\bf%s;%s.%s\;%s$" % (str(origtime[:-4]),network,station,meth)+"\n"+
                    r"$\bfbaz$: "+"{:.1f}".format(baz)+
                    r"N$^\circ$E, $\bfain$: "+"{:.1f}".format(ain)+
                    r"$^\circ$, $\bfepi$: "+"{:.1f}".format(dist)+" km "+
                    r"$\bfmag:$ "+"{:.1f}".format(mag)+"\n"+
                    r"$\bfφ:$ "+"{:.1f}".format(phi)+r" $\pm$ "+"{:.1f}".format(errors[0])+
                    r" N$^\circ$E, $\bft_d$: "+
                    "{:.2f}".format(abs(td)*1000.)+r" $\pm$ "+"{:.2f}".format(errors[1]*1000.)+" ms "+
                    r", $\bfp:$ "+"{:.1f}".format(pol)+
                    r"N$^\circ$E, "+r"$\bfgrade:$ "+grade+
                    r", $\bfflt:$ %.1f | %.1f Hz" % (filterRange[0],filterRange[1]),
                    fontsize=8

                )
    ## make the contour plot
    ax=fig.add_subplot(1,2,2)
    # sort values cause reasons
    phi_idx=np.argsort(phi_test)
    phi_test=phi_test[phi_idx]
    #td_test=td_test[phi_idx]
    Cmatrix=Cmatrix[phi_idx]    
    # get values from the matrices
    if "RC" in meth:
        idxBest=np.where(Cmatrix==Cmatrix.max())
    else:
        idxBest=np.where(Cmatrix==Cmatrix.min())
    phiCR=phi_test[idxBest[0]][0]
    tdCR=td_test[idxBest[1]][0]
    # make the contour
    csRange=np.arange(1.0,2.0,0.1)
    cbRange=np.arange(1.0,2.0,0.1)
    # make the contours
    cs=ax.contour(td_test,phi_test,Cmatrix/cEmax,colors='green',linewidth=1.0)
    #ax.clabel(cs)
    # make the cross
    ax.axvline(tdCR,linestyle='--',color='blue',linewidth=1.0) 
    ax.axhline(phiCR,linestyle='--',color='blue',linewidth=1.0)
    ax.plot(tdCR,phiCR,'b+')
    # make the 95% confidence interval contour
    ax.contour(td_test,phi_test,Cmatrix/cEmax,linewidths=2,colors="black",
               levels=[1.])
    # set labels and layout
    ax.set_xlabel(r"$t_d$ (ms)")
    ax.set_ylabel(r"φ ($\circ$)",labelpad=0.)
    for ax in fig.get_axes():
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(8)
    fig.tight_layout(rect=(0.1,0.1,0.9,0.9),w_pad=0.0,h_pad=0.0)
    # save to file?
    if output:
        fig.savefig(output,dpi=720,quality=95)

def qcClusterDiagram(initial,initial_err,calinski,clusters1,clusters2,station,network,baz,
                    titleInfo,meth,filterRange,output=False):
    """
    Make the Quality Control figure for CA.

    :type initial: array-like
    :param initial: a 2D array of the initial data used in the clustering (phi/td).
    :type iniital_err: array-like
    :param intial_err: a 2D array of the initial data errors used in the clustering (phi/td).
    :type calinski: array-like
    :param calinski: an array containing the Calinski-Harabasz score for each number of clusters.
    :type clusters1: tuple-like
    :param clusters2: data and labels of clusters
    :type clusters2: tuple-like
    :param clusters2: data and labels of selected cluster
    :type station: str
    :param station: the station name
    :type network: str
    :param network: the network name
    :type baz: str
    :param baz: the station's backazimuth
    :titleInfo: tuple-like
    :param titleInfo: list containing information related to the title of the figure, i.e. 
        the origin time, the incidence angle, the epicentral distance, the magnitude and
        the measurement's grade.
    :type meth: str
    :param meth: the method used for measuring splitting
    :type filterRange: tuple-like
    :param filterRange: the upper and lower filter bounds (in Hz).
    :type output: str or bool, optional
    :param output: If not False, the figure will be saved at the specified filepath. Defaults
        to False.

    """
    # fonts
    fontTicksLabel=7
    fontAxesLabel=7
    fontAxisTitle=8
    # unpack
    phis,tds,lClear=initial
    sphis,stds=initial_err
    M,C,CHM=calinski
    O,L,Mopt=clusters1
    Dopt,td,phi=clusters2
    # make the initial figure
    fig=plt.figure("Cluster Analysis",figsize=(7,6))
    # initial data figure
    axIni=fig.add_subplot(321)
    axIni.scatter(tds*1000,phis)
    axIni.set_ylabel(r"φ (N$^\circ$E)",fontsize=fontAxesLabel)
    axIni.set_xlabel(r"$t_d$ (ms)",fontsize=fontAxesLabel)
    axIni.set_ylim(0,180)
    axIni.set_title("Initial Data",fontsize=fontAxisTitle)
    # plot the C-H criterion function
    axCH=fig.add_subplot(322)
    axCH.plot(M,C, '--bo')
    axCH.axvline(CHM)
    axCH.axhline(C.max())
    axCH.set_title("Calinski-Harabasz Criterion / M=%i" % CHM,fontsize=fontAxisTitle)
    axCH.set_xlabel("Number of Clusters",fontsize=fontAxesLabel)
    axCH.set_xticks(M[0::3])
    axCH.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=True)
    offsetFactor=axCH.yaxis.get_major_formatter().get_offset()
    axCH.yaxis.offsetText.set_visible(False)
    if not offsetFactor: offsetFactor='1'
    axCH.set_ylabel("C-H Score (x %s)" % offsetFactor,fontsize=fontAxesLabel)
    ## plot the initial clusters
    axC1=fig.add_subplot(323,sharex=axIni,sharey=axIni)
    axC1.scatter(O[:,0]*1000,O[:,1],s=60,facecolor="None",edgecolor="k",label="Rejected")
    for i in lClear:
        data=O[L[Mopt]==i]
        axC1.scatter(data[:,0]*1000,data[:,1],s=100,marker="s",edgecolor="k",label="Cluster %i"%(i+1))
    axC1.set_xlabel(r"$t_d$ (ms)",fontsize=fontAxesLabel)
    axC1.set_ylabel(r"φ (N$^\circ$E)",fontsize=fontAxesLabel)
    axC1.set_title("Final Clusters",fontsize=fontAxisTitle)
    plt.legend(handlelength=0.6,frameon=True,fontsize=8,loc='center left', bbox_to_anchor=(1, 0.5))
    ## show the final measurement
    axC2=fig.add_subplot(324,sharex=axIni,sharey=axIni)
    axC2.scatter(O[:,0]*1000,O[:,1],s=60,facecolor="None",edgecolor="k",label="Rejected Clusters")
    axC2.scatter(Dopt[:,0]*1000,Dopt[:,1],marker="o",facecolor="red",edgecolor="k",label="Selected Cluster")
    axC2.axvline(td*1000); axC2.axhline(phi)
    #axC2.set_ylabel(r"φ (N$^\circ$E)",fontsize=fontAxesLabel)
    axC2.set_xlabel(r"$t_d$ (ms)",fontsize=fontAxesLabel)
    plt.legend(handlelength=0.6,frameon=True,fontsize=8)
    axC2.set_title("Final Result",fontsize=fontAxisTitle)
    # set limits
    axIni.set_ylim(min(phis)*0.9,max(phis)*1.1)
    axIni.set_xlim(min(tds*1000)*0.9,max(tds*1000)*1.1)
    ## Teanby et al. (2004)-like plots
    windows=np.arange(1,len(phis)+1)
    # first for the phi
    axTB1=fig.add_subplot(325)
    #axTB1.scatter(windows,phis)
    axTB1.errorbar(x=windows,y=phis,yerr=sphis,fmt="bo",capsize=3)
    axTB1.set_xlabel("Window #",fontsize=fontAxesLabel)
    axTB1.set_ylabel(r"φ (N$^\circ$E)",fontsize=fontAxesLabel)
    # now for the dt
    axTB2=fig.add_subplot(326)
    #axTB1.scatter(windows,phis)
    axTB2.errorbar(x=windows,y=tds*1000,yerr=stds*1000,fmt="bo",capsize=3)
    axTB2.set_xlabel("Window #",fontsize=fontAxesLabel)
    axTB2.set_ylabel(r"$t_d$ (ms)",fontsize=fontAxesLabel)    
    # finalize
    fig.subplots_adjust(wspace=0.35,hspace=0.25)
    origtime,ain,dist,mag,grade=titleInfo
    fig.suptitle(
                    #r"$\bf%s$" % str(origtime[:-4])+r" $\bf%s.%s\;%s$" % (network,station,meth),
                    "%s %s.%s %s (Windows: %s)" % (str(origtime[:-4]),network,station,meth,len(windows)),
                    fontweight="bold",fontsize=8
                )
    # change font size for tick labels
    for tax in fig.get_axes():
        tax.tick_params(labelsize=fontTicksLabel,direction="in")
    fig.tight_layout(rect=[0,0.03,1.00,0.97])
    # save to file?
    if output:
        fig.savefig(output,dpi=720,quality=95)

def Rmatrix2D(baz):
    """
    Returns the clockwise rotation matrix (according to IRIS services) for
    the provided azimuth (measured from the North clockwise).

    :type baz: float
    :param baz: the backazimuth used for the rotation (in degrees).
    :returns: the 2D rotation matrix

    """
    baz=np.deg2rad(baz)
    return np.asarray([
            [np.cos(baz), np.sin(baz)],
            [-np.sin(baz),np.cos(baz)]
                    ])

def Rmatrix3D(baz,ain):
    """
    Returns the clockwise rotation matrix for the provided
    azimuth and angle of incidence. Rotation according to
    Plešinger et al. (1986).0

    :type baz: float
    :param baz: the backazimuth used for the rotation (in degrees).
    :type ain: float
    :param ain: the incidence angle at the station (in degrees)>
    :returns: the 3D rotation matrix
    
    ..note:
        Plesinger, A., Hellweg, M., Seidl, D., 1986. 
            Interactive high-resolution polarization analysis of broad-band seismograms. 
            J. Geophys. 59, 129–139.
    
    """
    baz=np.deg2rad(baz)
    ain=np.deg2rad(ain)
    return np.asarray([
            [np.cos(baz),-np.sin(ain)*np.sin(baz),-np.sin(ain)*np.cos(baz)],
            [np.sin(ain), np.cos(ain)*np.sin(baz), np.cos(ain)*np.cos(baz)],
            [0,          -np.cos(baz),             np.sin(baz)]
                     ])

def isreal(string):
    """
    Validates that the string is a real number. We're using
    `~numpy.float` to perform the check to include 'nan' strings.

    :type string: str
    :param string: The string for which the check will be performed
    :returns: True if the string is real, otherwise false.

    """
    try:
        ck = np.float(string)
        return True
    except ValueError:
        return False

def isdate(string):
    """
    Validates that the string is a date. We're using `~obspy.core.utcdatetime.UTCDateTime` 
    to perform the check.

    :type string: str
    :param string: The string for which the check will be performed
    :returns: True if the string is datetime, otherwise false.

    """
    try:
        ck = UTCDateTime(string)
        return True
    except ValueError:
        return False

def isnone(value):
    """
    Validates that the value is None, `~numpy.nan` or `~numpy.inf`.

    :type value: any type
    :param value: The string for which the check will be performed
    :returns: True if the string is datetime, otherwise false.

    """
    if value is None:
        return True
    elif np.isnan(value):
        return True
    elif np.isinf(value):
        return True
    else:
        return False      