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
# This module includes related to I/O of Preferences.              #
#                                                                  #
####################################################################

## imports 
import logging
from configparser import ConfigParser

class parseGeneralCnf():
    """
    class that contains all required settings
    from the config file.

    """

    def __init__(self,filename):
        """
        Reads a configuration file and returns the settings.
        
        :type filename: str
        :param filename: the path to the file        

        """
         # initiate the configparser object and read the file
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        # parse
        self.cleanLogs=cnf.getfloat('GENERAL','cleanLogs')
        self.matching=cnf.getfloat('WAVEFORMS','matching')
        self.chanPref=cnf.get('WAVEFORMS','prefOrder')
        self.orientFlag=cnf.getboolean('WAVEFORMS','orientation')
        self.trimStart=cnf.getfloat('WAVEFORMS','trimStart')
        self.trimEnd=cnf.getfloat('WAVEFORMS','trimEnd')
        self.trimFlag=cnf.getboolean('WAVEFORMS','trim')
        self.snrStart=cnf.getfloat('SNR','snrStart')
        self.snrEnd=cnf.getfloat('SNR','snrEnd')
        self.maxTd=cnf.getfloat('SPLITTING','maxTd')


    def write(self,filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file

        """
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        cnf.set('GENERAL','cleanLogs',str(self.cleanLogs))
        cnf.set('WAVEFORMS','matching',str(self.matching))
        cnf.set('WAVEFORMS','trimStart',str(self.trimStart))        
        cnf.set('WAVEFORMS','trimEnd',str(self.trimEnd))      
        cnf.set('WAVEFORMS','trim',str(self.trimFlag))          
        cnf.set('WAVEFORMS','prefOrder',str(self.chanPref))
        cnf.set('WAVEFORMS','orientation',str(self.orientFlag))
        cnf.set('SNR','snrStart',str(self.snrStart))
        cnf.set('SNR','snrEnd',str(self.snrEnd))
        cnf.set('SPLITTING','maxTd',str(self.maxTd))
        with open(filename,"w") as fid: cnf.write(fid)
        logging.debug("Saved General settings to %s" % filename)        

class parsePickerCnf():
    """
    Class that contains all required settings
    from the config file

    """

    def __init__(self,filename):
        """
        Reads an ini file and returns the settings

        :type filename: str
        :param filename: the path to the file

        """
        # initiate the configparser object and read the file
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        # parse the defaults first
        self.arFMIN=cnf.getfloat('ARPICKER','f1')
        self.arFMAX=cnf.getfloat('ARPICKER','f2')
        self.arLTAP=cnf.getfloat('ARPICKER','lta_p')
        self.arSTAP=cnf.getfloat('ARPICKER','sta_p')
        self.arLTAS=cnf.getfloat('ARPICKER','lta_s')
        self.arSTAS=cnf.getfloat('ARPICKER','sta_s')
        self.arMP=cnf.getint('ARPICKER','m_p')
        self.arMS=cnf.getint('ARPICKER','m_s')
        self.arLP=cnf.getfloat('ARPICKER','l_p')
        self.arLS=cnf.getfloat('ARPICKER','l_s')

    def write(self,filename):
        """
        Write values to file

        :type filename: str
        :param filename: the path to the file

        """
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        cnf.set('ARPICKER','f1',str(self.arFMIN))        
        cnf.set('ARPICKER','f2',str(self.arFMAX))
        cnf.set('ARPICKER','lta_p',str(self.arLTAP))
        cnf.set('ARPICKER','sta_p',str(self.arSTAP))
        cnf.set('ARPICKER','lta_s',str(self.arLTAS))
        cnf.set('ARPICKER','sta_s',str(self.arSTAS))
        cnf.set('ARPICKER','m_p',str(self.arMP))
        cnf.set('ARPICKER','m_s',str(self.arMS))
        cnf.set('ARPICKER','l_p',str(self.arLP))
        cnf.set('ARPICKER','l_p',str(self.arLS))
        with open(filename,"w") as fid: cnf.write(fid)
        logging.debug("Saved AR-AIC settings to %s" % filename)

class parseClusteringCnf():
    """
    Class that contains all required settings
    from the config file of Cluster Analysis.

    """

    def __init__(self,filename):
        """
        Reads an ini file and returns the settings

        :type filename: str
        :param filename: the path to the file

        """
        # initiate the configparser object and read the file
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        # parse the windows settings
        self.Tbeg1=cnf.getfloat('WINDOWS','Tbeg1')
        self.DTbeg=cnf.getfloat('WINDOWS','DTbeg')
        self.DTend=cnf.getfloat('WINDOWS','DTend')
        self.Tend0=cnf.getfloat('WINDOWS','Tend0')
        self.tptsMax=cnf.getfloat('WINDOWS','tptsMax')
        self.specWindow=cnf.getfloat('WINDOWS','specWindow')
        self.multPeriod=cnf.getfloat('WINDOWS','multPeriod')
        self.minPeriod=cnf.getfloat('WINDOWS','minPeriod')
        self.maxPeriod=cnf.getfloat('WINDOWS','maxPeriod')
        # parse the clustering settings
        self.ccrit=cnf.getfloat('CLUSTERING','ccrit')
        self.kmax=cnf.getint('CLUSTERING','kmax')
        self.Ncmin=cnf.getint('CLUSTERING','Ncmin')
        self.linkage=cnf.get('CLUSTERING','linkage')

    def write(self,filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file

        """
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        cnf.set('WINDOWS','Tbeg1',str(self.Tbeg1))        
        cnf.set('WINDOWS','DTbeg',str(self.DTbeg))
        cnf.set('WINDOWS','DTend',str(self.DTend))
        cnf.set('WINDOWS','Tend0',str(self.Tend0))
        cnf.set('WINDOWS','tptsMax',str(self.tptsMax))
        cnf.set('WINDOWS','specWindow',str(self.specWindow))
        cnf.set('WINDOWS','multPeriod',str(self.multPeriod))
        cnf.set('WINDOWS','minPeriod',str(self.minPeriod))
        cnf.set('WINDOWS','maxPeriod',str(self.maxPeriod))
        cnf.set('CLUSTERING','ccrit',str(self.ccrit))
        cnf.set('CLUSTERING','kmax',str(self.kmax))
        cnf.set('CLUSTERING','Ncmin',str(self.Ncmin))
        cnf.set('CLUSTERING','linkage',str(self.linkage))
        with open(filename,"w") as fid: cnf.write(fid)
        logging.debug("Saved CLUSTERING settings to %s" % filename)

class parseTaupCnf():
    """
    Class that contains all required settings
    from the config file.

    """

    def __init__(self,filename):
        """
        Reads an ini file and returns the settings

        :type filename: str
        :param filename: the path to the file

        """
        # initiate the configparser object and read the file
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        # parse the defaults first
        self.model=cnf.get('PATHS','model').strip()
        self.stations=cnf.get('PATHS','stations').strip()     
        # options
        self.ainFlag=cnf.getboolean('OPTIONS','ain')

    def write(self,filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file

        """
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        cnf.set('PATHS','model',str(self.model))        
        cnf.set('PATHS','stations',str(self.stations))
        cnf.set('OPTIONS','ain',str(self.ainFlag))
        with open(filename,"w") as fid: cnf.write(fid)
        logging.debug("Saved CLUSTERING settings to %s" % filename)

class parseGradeCnf():
    """
    Class that contains all required settings
    from the config file.
    
    """
    def __init__(self,filename):
        """
        Reads an ini file and returns the settings
        
        :type filename: str
        :param filename: the path to the file
        
        """
        # initiate the configparser object and read the file
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        # angle
        self.polOff=cnf.getfloat('ANGLE','diff')
        # error bounds
        self.snr_bound=cnf.getfloat('BOUNDS','SNR')
        self.error_bounds=[cnf.getfloat('BOUNDS','phi'),cnf.getfloat('BOUNDS','td')]
        self.CC_FS_bound=cnf.getfloat('BOUNDS','CC_FS')
        self.CC_NE_bound=cnf.getfloat('BOUNDS','CC_NE')
        # grading
        A=cnf.getfloat('GRADING','A')
        B=cnf.getfloat('GRADING','B')
        C=cnf.getfloat('GRADING','C')
        D=cnf.getfloat('GRADING','D')
        self.gradeDict={"A":A,"B":B,"C":C,"D":D}

    def write(self,filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file
                
        """
        cnf=ConfigParser(inline_comment_prefixes="!")
        with open(filename,"r") as fid: cnf.read_file(fid)
        cnf.set('ANGLE','diff',str(self.polOff))        
        cnf.set('BOUNDS','SNR',str(self.snr_bound))
        cnf.set('BOUNDS','phi',str(self.error_bounds[0]))
        cnf.set('BOUNDS','td',str(self.error_bounds[1]))
        cnf.set('BOUNDS','CC_FS',str(self.CC_FS_bound))
        cnf.set('BOUNDS','CC_NE',str(self.CC_NE_bound))
        cnf.set('GRADING','A',str(self.gradeDict["A"]))
        cnf.set('GRADING','B',str(self.gradeDict["B"]))
        cnf.set('GRADING','C',str(self.gradeDict["C"]))
        cnf.set('GRADING','D',str(self.gradeDict["D"]))
        with open(filename,"w") as fid: cnf.write(fid)
        logging.debug("Saved CLUSTERING settings to %s" % filename)