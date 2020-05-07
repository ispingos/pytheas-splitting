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

Authors: Spingos I. & Kaviris G. (c) 2019-2020
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
from PyQt5 import QtCore
import os, logging, time
import itertools
import collections
import numpy as np
from obspy import read
from matplotlib import pyplot as plt
## scikit-learn imports 
from sklearn import metrics
from sklearn.cluster import AgglomerativeClustering as AC
# pytheas
from lib import eigenvalue as SC
from lib import rotationcorrelation as RC

class clustering(QtCore.QThread):
    """
    The main clustering class. This is the implementation of
    the method proposed by Teanby et al. (2004), with slight
    alterations.

    Uses :class: `~PyQt5.QtCore.QThread`

    ..note:
        Teanby, N., Kendall, J.-M., van der Baan, M., 2004.
            Automation of shear-wave splitting measurements using cluster analysis. 
            Bull. Seismol. Soc. Am. 94, 453–463. doi: 10.1785/0120030123

    """

    # setup custom signals
    iterDone=QtCore.pyqtSignal(str)
    iterFail=QtCore.pyqtSignal(Exception)

    def __init__(self,stream,method,sPick,bazi,ain,tbCNF,parent=None):
        """
        Initialize

        :type stream: :class: ~obspy.core.stream.Stream`
        :param stream: the waveform data in a Stream object.
        :type method: str
        :param method: analysis method to use with CA (EV, ME or RC).
        :type sPick: float
        :param sPick: the S-arrival relative from the stream's start time (in s).
        :type bazi: float
        :param bazi: the backazimuth of the station.
        :type ain: float
        :param ain: the incidence angle at the station.
        :type tbCNF: :class:`~pytheas.lib.parsers.parsePickerCnf`
        :param tbCNF: Configuration class containing parameters related to clustering.

        """
        QtCore.QThread.__init__(self,parent)
        # store arguments to self
        self.stream=stream
        self.method=method
        self.sPick=sPick
        self.bazi=bazi
        self.ain=ain
        self.caCNF=tbCNF
        self.excFlag=False
        self._isRunning=True

    def __del__(self):
        """Change the __del__ method."""
        self.wait()

    def windowing(self,TbegA,TendA,DTbeg,DTend,Nbeg,Nend):
        """
        make the windows to be used for measurements. According
        to Teanby et al. (2004):

        (1) Window should be representative of the S-wave
        (2) Window should be long enough to include several S periods
        (3) Window shouldn't be very long to not include secondary phases
        (4) Window should be at least one period long
        (5) Window should start slightly before the S arrival

        :type TbegA: tuple-like
        :param TbegA: list containing the start and end points for generating
            the beginning of the windows (relative to the S-arrival in s).
        :type TendA: tuple-like
        :param TbegA: list containing the start and end points for generating
            the ending of the windows (relative to the S-arrival in s).
        :type DTbeg: float
        :param DTbeg: the step for window beginnings.
        :type DTend: float
        :param DTend: the step for window endings.
        :type Nbeg: float
        :param Nbeg: the number of window beginnings.
        :type Nend: float
        :param Nend: the number of window endings.
        :returns: tuple with window endings-beginnings combinations,
            in ascending length order.

        """
        # unpack #
        Tbeg0,Tbeg1=TbegA
        Tend0,Tend1=TendA  
        # get ready #
        windows=[]
        # run # 
        for i,j in itertools.product(range(Nbeg),range(Nend)):
            i+=1; j+=1
            Tbeg=Tbeg1-(i-1)*DTbeg
            Tend=Tend0+(j-1)*DTend
            if (Tbeg < Tbeg0) or (Tend > Tend1):
                continue
            windows.append((Tbeg,Tend))
        windows=np.asarray(windows)
        # sort windows based on length
        dur=windows[:,1]-windows[:,0]
        windows=windows[np.argsort(dur)]
        return windows

    def standardZ4(self,X):
        """
        Perform standarization of variables per the Z4 method
        described in Milligan and Cooper (1988). If all the values
        are the same, the returned array is comprised by ones.
        
        :type X: array-like
        :param X: the 1D array to standarize.
        :returns: the standarized array
        
        Milligan, G.W., Cooper, M.C., 1988. 
            A study of standardization of variables in cluster analysis. 
            J. Classif. 5, 181–204. doi: 10.1007/BF01897163
        
        """
        # what if only one value in array? if not addressed, this
        # will result to inf and subsequent exception in clustering!
        val=X[1] # any value from the array will do
        # test
        if np.all(X==val):
            return X/X # this is equal to np.ones(X.shape)
        else:
            return X/(np.nanmax(X)-np.nanmin(X))

    def getClustersSK(self,X,method="single"):
        """
        Get the model and labels for all possible clusters from
        the built-in sklearn function for agglomerative clustering.
        No connectivity matrix is used, since data are (apparently)
        unstructured.
        k=[1,N-1], where N the number of observations

        :type X: array-like
        :param X: 2D array containing the x,y coordinates of the points
            to be clustered.
        :type method: str, optional
        :param method: method for linking clusters. Defaults to 'single'. 

        """
        Mmax=len(X)-1
        M=np.arange(1,Mmax+1)
        L={}
        for k in M:
            model=AC(n_clusters=k,linkage=method,affinity="euclidean").fit(X)
            L.update({k:model.labels_})
        return L

    def getClusterDataSK(self,X,L):
        """
        Get the data for each k and each cluster

        :type X: array-like
        :param X: a 2D array containing the x,y 
            coordinates of the data.
        :type L: array-like
        :param L: a 1D array containing the labels
            of the data.
        :returns: a dict containing the data points for
            each cluster at each number of clusters.

        """
        ZData={k:{} for k in L.keys()}
        for k in sorted(L):
            kLabels=L[k]
            kDict={}
            for i in np.arange(k):
                data=X[kLabels==i]
                kDict.update({i:data})
            ZData[k]=kDict
        return ZData

    def getClusterCentresSK(self,ZData): 
        """
        Get the centre of each cluster for each k
        
        :type ZData: dict
        :param ZData: a dict containing all the data points
            for each cluster for each number of clusters k
        :returns: a dict which includes the centre of each
            cluster for each k

        """
        C={k:{} for k in ZData}
        for k in sorted(ZData):
            kDict=ZData[k]
            tempDict={}
            for i in sorted(kDict):
                data=kDict[i]
                dt=np.asarray([x[0] for x in data],dtype="float").mean()
                phi=np.asarray([x[1] for x in data],dtype="float").mean()
                tempDict.update({i:(dt,phi)})
            C[k]=tempDict
        return C
      
    def critCalHara(self,X,L,Mmax):
        """
        Calculate the Calinski and Harabasz (1974) score from
        `~sklearn.metrics.calinski_harabaz_score`.

        :type X: array-like
        :param X: the x,y coordinates of the data.
        :type L: array-like
        :param L: the data labels.
        :type Mmax: int
        :param Mmax: the maximum accepted number of clusters.
        :returns: (the number of clusters, the CH scores)

        ..note:
            Caliński, T., Harabasz, J., 1974. 
                A dendrite method for cluster analysis. 
                Commun. Stat. 3, 1–27. doi: 10.1080/03610927408827101

        """     
        M=np.arange(2,Mmax+1)
        C=np.zeros(M.shape)
        for k in sorted(M): # skip k=1
            labels=L[k]
            c=metrics.calinski_harabaz_score(X,labels)
            C[k-2]=c
        return M,C
        
    def critDudaHartSK(self,X,Z,k,ccrit=3.2):   
        """
        Compute Ccritical for all clusters from the
        agglomerative clustering, according to 
        Duda and Hart (1973). Use a sklearn dataset
        
        :type X: array-like
        :param X: the x,y coordinates of the data.
        :type Z: dict
        :param Z: a dict containing all the data points for 
            each cluster for each number of clusters k
        :type k: int
        :param k: the number of clusters
        :type ccrit: float, optional
        :param ccrit: the critical value for which the 
            null hypothesis will be checked. Defaults to 3.2.
        :returns: the selected number of clusters and the respective
            score

        ..note:
            Duda, R.O., Hart, P.E., 1973. 
                Pattern Classification and Scene Analysis,
                1st ed. Wiley, New York.        

        """
        # prepare
        p=2 # number of parameters (i.e. phi and dt)
        N=len(X)-1
        M=np.flip(np.arange(1,N+1),0)
        C=np.zeros(N)
        # let's roll
        logging.debug("Duda and Hart trial for k = %i" % k)
        combos=itertools.combinations(range(k),2)
        for k1,k2 in combos:
            # now convert to array
            d1=Z[k][k1]; d2=Z[k][k2]
            da=np.concatenate((d1,d2))
            dt1=d1[:,0]; phi1=d1[:,1]
            dt2=d2[:,0]; phi2=d2[:,1]
            dta=da[:,0]; phia=da[:,1]
            # calculate sigma values for the DH criterion
            s22 = np.sum(((dt1-dt1.mean())**2)+((phi1-phi1.mean())**2))+\
                  np.sum(((dt2-dt2.mean())**2)+((phi2-phi2.mean())**2))
            s21 = np.sum(((dta-dta.mean())**2)+((phia-phia.mean())**2))
            # calculate the DH criterion
            c=(1-(s22/s21)-(2/(np.pi*p)))*np.sqrt((N*p)/(2*(1-(8/(np.pi*np.pi*p)))))
            if c <= ccrit:
                logging.debug("Duda and Hart score for k = %i: %.4f"%(k,c))
                return k,c
        return None,None

    def clusterVar(self,C,V):
        """
        Calculate the within-cluster and mean data variances
        for the given cluster.
        
        :type C: array-like
        :param C: the x,y coordinates of the cluster's data
        :type V: array-like
        :param V: the errors of the x and y parameters.
        :returns: (within-cluster variance, 
                   mean data variance, 
                   overall variance)

        """
        dt=C[:,0]; phi=C[:,1]; Nj=C.size
        sdt=V[:,0]; sphi=V[:,1]
        # within-cluster
        s2cj = np.nansum(((dt-np.nanmean(dt))**2)+(phi-np.nanmean(phi))**2)/Nj
        # mean data
        s2dj = 1 / (np.nansum(1/(sdt**2))) + 1 / (np.nansum(1/(sphi**2)))
        # overall variance
        s2oj = np.nanmax((s2cj,s2dj))
        return s2cj,s2dj,s2oj

    def stop(self):
        """Set flag to stop iterating over windows."""
        self._isRunning=False
        self.wait()

    def run(self):
        """Main function to perform the cluster analysis."""
        try:
            ## run ##
            logging.info("starting clustering...")
            # define the windows # 
            logging.info("..: define time windows for analysis")
            TbegA=(self.caCNF.Tbeg0,self.caCNF.Tbeg1); TendA=(self.caCNF.Tend0,self.caCNF.Tend1)
            windows=self.windowing(TbegA,TendA,self.caCNF.DTbeg,self.caCNF.DTend,self.caCNF.Nbeg,self.caCNF.Nend)
            windows+=self.sPick # in seconds
            logging.info("..: number of windows: %i"%len(windows))
            # perform measurements
            measDict={}; n=0
            tic=time.time()
            nWindows=len(windows); tempResList=[]
            for nw,pickwin in enumerate(windows):
                if not self._isRunning:
                    raise IOError("User cancelled the iterations.")
                logging.info("===========================================")
                logging.info("(%i/%i) Working on window: %s" % (nw+1,nWindows,pickwin))
                if pickwin.sum() == 0: # if no values were assigned to array
                    continue        
                try:
                    if self.method == "EV":
                        tempRes=SC.SilverAndChan(self.stream,self.bazi,pickwin,maxDelay=self.caCNF.maxTd/1000.,method='EV',ain=self.ain)
                    elif self.method == "ME":
                        tempRes=SC.SilverAndChan(self.stream,self.bazi,pickwin,maxDelay=self.caCNF.maxTd/1000.,method='ME',ain=self.ain)
                    elif self.method  == "RC":
                        tempRes=RC.crosscorrelation(self.stream,pickwin,maxDelay=self.caCNF.maxTd/1000.)
                    phi,dt,cArray,pTest,dTest,\
                    pol,sphi,sdt,cEmax,\
                    CC_FS,CC_NE,Tvar,nContours=tempRes
                    # if no errors were estimated it is a bad measurement (i.e.
                    # open 95% CL contour or point not in contour. Should be discarded
                    # regardless
                    if np.isnan(sphi) or np.isnan(sdt):
                        continue
                    n+=1
                    measDict.update({n:{"phi":phi,"dt":dt,
                                        "sphi":sphi,"sdt":sdt,
                                        "win":pickwin,"tempRes":tempRes}})
                except KeyboardInterrupt:
                    break
                except:
                    logging.exception("Error while processing window")
                self.iterDone.emit(str(round(100*nw/nWindows)))
            toc=time.time()
            logging.info("Time to apply window analysis: %.3f"%(toc-tic))
            # pack values in np arrays #
            phis=np.asarray([measDict[x]['phi'] for x in sorted(measDict)],dtype="float")
            sphis=np.asarray([measDict[x]['sphi'] for x in sorted(measDict)],dtype="float")
            dts=np.asarray([measDict[x]['dt'] for x in sorted(measDict)],dtype="float")
            sdts=np.asarray([measDict[x]['sdt'] for x in sorted(measDict)],dtype="float")
            idxes=np.asarray([x for x in sorted(measDict)],dtype="int")
            # standarize per the ranges #
            nphis=self.standardZ4(phis)
            nsphis=self.standardZ4(sphis)
            ndts=self.standardZ4(dts)
            nsdts=self.standardZ4(sdts)
            # store original data
            O=np.array((dts,phis),dtype="float").T # original data
            U=np.array((sdts,sphis),dtype="float").T # original errors
            # combine and prepare for clustering
            X=np.array((ndts,nphis),dtype="float").T # normalized data
            S=np.array((nsdts,nsphis),dtype="float").T # normalized errors
            # get labels
            L=self.getClustersSK(X,method=self.caCNF.linkage)
            # get data for every cluster and k
            ZData=self.getClusterDataSK(X,L)
            ## apply the two criteria to obtain optimal number of clusters
            # first Duda and Hart
            Mmax=np.inf
            for k in sorted(ZData):
                MD,CD=self.critDudaHartSK(X,ZData,k,ccrit=self.caCNF.ccrit)    
                if MD != None:
                    Mmax=k
                    break    
            # is Mmax from Duda and Hart greater than accepted k ?                    
            Mmax=self.caCNF.kmax if Mmax > self.caCNF.kmax else Mmax
            # now Calinski-Harabasz
            MO,CO=self.critCalHara(X,L,Mmax)
            # keep k/score pairs for k less than user specified
            MH=MO[MO<=self.caCNF.kmax]
            CH=CO[MO<=self.caCNF.kmax]
            # use both criteria to select number of clusters
            idx=CH.argsort()[::-1][:] # get C-H score descending indices
            Mopt=MH[idx][0]
            M=MH; C=CH; CHM=Mopt 
            ## apply the two criteria to obtain optimal cluster
            # filter out spurious clusters with number of points selection
            lList=np.arange(Mopt); lIdx=L[Mopt] 
            lCount=collections.Counter(lIdx)
            delList=[]
            for l in sorted(lCount):
                if lCount[l] < self.caCNF.Ncmin:
                    delList.append(l)
            if len(lList) != len(delList):
                lClear=np.setdiff1d(lList,delList)
                logging.info("Removed %i clusters due to the minimum N = %i points limit"%(len(delList),self.caCNF.Ncmin))
            else:
                raise ValueError("Did not apply minimum N criterion, none of the clusters would pass")
            # calculate the within-cluster and mean data variances
            varDict={x:None for x in lClear}
            for l in lClear:
                D=X[L[Mopt]==l]
                _,_,varDict[l]=self.clusterVar(D,S)
            Copt=min(varDict,key=varDict.get)
            # get original data for optimal cluster
            Dopt=O[L[Mopt]==Copt] # optimal data
            Uopt=U[L[Mopt]==Copt] # optimal uncertainties
            Iopt=idxes[L[Mopt]==Copt] # optimal data indices
            # combine normalized errors for optimal cluster
            Sopt=S[L[Mopt]==Copt]
            Scomb=Sopt[:,1]+Sopt[:,0]
            optIdx=np.where(Scomb==np.nanmin(Scomb))[0]
            # get the values
            if optIdx.size == 1:
                phi=Dopt[optIdx[0],1]
                dt=Dopt[optIdx[0],0]
                sphi=Uopt[optIdx[0],1]
                sdt=Uopt[optIdx[0],0]
                ix=int(Iopt[optIdx[0]])
                optWindow=measDict[ix]["win"]
                tempResFinal=measDict[ix]["tempRes"]
            else:
                raise ValueError("More than 1 optimal measurements!")
            logging.info("Final results: phi: %.3f +/- %.3f, dt: %.4f +/- %.4f"%(phi,sphi,dt,sdt))
            logging.info("Number of points: %i"%Dopt.size)  
            toctoc=time.time()                
            logging.debug("Elapsed time for whole process (m): %.3f / windows: %i"%((toctoc-tic)/60,len(windows)))
            # prepare for the QC figure
            self.initial=phis,dts,lClear
            self.initial_err=sphis,sdts
            self.calinski=M,C,CHM
            self.clusters1=O,L,Mopt
            self.clusters2=Dopt,dt,phi
            #
            self.phi=phi; self.dt=dt; self.sphi=sphi; self.sdt=sdt
            self.optWindow=optWindow; self.tempRes=tempResFinal
        except Exception as exc:
            self.excFlag=True
            self.iterFail.emit(exc)