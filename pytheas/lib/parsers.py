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

Authors: Spingos I. & Kaviris G. (c) 2019-2021
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
import os
import logging
import numpy as np
from obspy import read_events, read_inventory
from configparser import *

class parseGeneralCnf():
    """
    class that contains all required settings
    from the config file.

    """

    def __init__(self, filename):
        """
        Reads a configuration file and returns the settings.
        
        :type filename: str
        :param filename: the path to the file        

        """
        self.default_sections = ['GENERAL', 'WAVEFORMS', 'SNR', 'SPLITTING']
         # initiate the configparser object and read the file
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)
        # parse
        try:
            self.cleanLogs = cnf.getfloat('GENERAL', 'cleanLogs')
        except NoOptionError:
            self.cleanLogs = 1.
        try:
            self.max_log_size = cnf.getfloat('GENERAL', 'max_log_size')
        except NoOptionError:
            self.max_log_size = 50.
        try:
            self.matching = cnf.getfloat('WAVEFORMS', 'matching')
        except NoOptionError:
            self.matching = 5.0
        try:
            self.chanPref = cnf.get('WAVEFORMS', 'prefOrder')
        except NoOptionError:
            self.chanPref = 'HH,EN,HN,EH,BH,BN'
        try:
            self.orientFlag = cnf.getboolean('WAVEFORMS', 'orientation')
        except NoOptionError:
            self.orientFlag = True
        try:
            self.trimStart = cnf.getfloat('WAVEFORMS', 'trimStart')
        except NoOptionError:
            self.trimStart = -5.
        try:
            self.trimEnd = cnf.getfloat('WAVEFORMS', 'trimEnd')
        except NoOptionError:
            self.trimEnd = 5.
        try:
            self.trimFlag = cnf.getboolean('WAVEFORMS', 'trim')
        except NoOptionError:
            self.trimFlag = False
        try:
            self.snrStart = cnf.getfloat('SNR', 'snrStart')
        except NoOptionError:
            self.snrStart = -0.2
        try:
            self.snrEnd = cnf.getfloat('SNR', 'snrEnd')
        except NoOptionError:
            self.snrEnd = 0.2
        try:
            self.maxTd = cnf.getfloat('SPLITTING', 'maxTd')
        except NoOptionError:
            self.maxTd = 250.


    def write(self, filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file

        """
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)
        cnf.set('GENERAL', 'cleanLogs', str(self.cleanLogs))
        cnf.set('GENERAL','max_log_size', str(self.max_log_size))
        cnf.set('WAVEFORMS', 'matching', str(self.matching))
        cnf.set('WAVEFORMS', 'trimStart', str(self.trimStart))        
        cnf.set('WAVEFORMS', 'trimEnd', str(self.trimEnd))      
        cnf.set('WAVEFORMS', 'trim', str(self.trimFlag))          
        cnf.set('WAVEFORMS', 'prefOrder', str(self.chanPref))
        cnf.set('WAVEFORMS', 'orientation', str(self.orientFlag))
        cnf.set('SNR', 'snrStart', str(self.snrStart))
        cnf.set('SNR', 'snrEnd', str(self.snrEnd))
        cnf.set('SPLITTING', 'maxTd', str(self.maxTd))
        with open(filename, "w") as fid:
            cnf.write(fid)
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
        self.default_sections = ['ARPICKER']
        # initiate the configparser object and read the file
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)
        # parse the defaults first
        try:
            self.arFMIN = cnf.getfloat('ARPICKER', 'f1')
        except NoOptionError:
            self.arFMIN = 1.
        try:
            self.arFMAX = cnf.getfloat('ARPICKER', 'f2')
        except NoOptionError:
            self.arFMAX = 10.
        try:
            self.arLTAP = cnf.getfloat('ARPICKER','lta_p')
        except NoOptionError:
            self.arLTAP = 1.
        try:
            self.arSTAP = cnf.getfloat('ARPICKER','sta_p')
        except NoOptionError:
            self.arSTAP = 0.5
        try:
            self.arLTAS = cnf.getfloat('ARPICKER','lta_s')
        except NoOptionError:
            self.arLTAS = 2.0
        try:
            self.arSTAS = cnf.getfloat('ARPICKER','sta_s')
        except NoOptionError:
            self.arSTAS = 1.0
        try:
            self.arMP = cnf.getint('ARPICKER','m_p')
        except NoOptionError:
            self.arMP = 2
        try:
            self.arMS = cnf.getint('ARPICKER','m_s')
        except NoOptionError:
            self.arMS = 8
        try:
            self.arLP = cnf.getfloat('ARPICKER','l_p')
        except NoOptionError:
            self.arLP = 0.2
        try:
            self.arLS = cnf.getfloat('ARPICKER','l_s')
        except NoOptionError:
            self.arLS = 0.2

    def write(self,filename):
        """
        Write values to file

        :type filename: str
        :param filename: the path to the file

        """
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)        
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
        with open(filename,"w") as fid:
            cnf.write(fid)
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
        self.default_sections = ['WINDOWS', 'CLUSTERING']
        # initiate the configparser object and read the file
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)        # parse the windows settings
        try:
            self.Tbeg1 = cnf.getfloat('WINDOWS','Tbeg1')
        except NoOptionError:
            self.Tbeg1 = -0.1
        try:
            self.DTbeg = cnf.getfloat('WINDOWS','DTbeg')
        except NoOptionError:
            self.DTbeg = 0.2
        try:
            self.DTend = cnf.getfloat('WINDOWS','DTend')
        except NoOptionError:
            self.DTend = 0.02
        try:
            self.Tend0 = cnf.getfloat('WINDOWS','Tend0')
        except NoOptionError:
            self.Tend0 = 0.1
        try:
            self.tptsMax = cnf.getfloat('WINDOWS','tptsMax')
        except NoOptionError:
            self.tptsMax = 1.5
        try:
            self.specWindow = cnf.getfloat('WINDOWS','specWindow')
        except NoOptionError:
            self.specWindow = 0.5
        try:
            self.multPeriod = cnf.getfloat('WINDOWS','multPeriod')
        except NoOptionError:
            self.multPeriod = 3.
        try:
            self.minPeriod = cnf.getfloat('WINDOWS','minPeriod')
        except NoOptionError:
            self.minPeriod = 0.1
        try:
            self.maxPeriod = cnf.getfloat('WINDOWS','maxPeriod')
        except NoOptionError:
            self.maxPeriod = 2.0
        # parse the clustering settings
        try:
            self.ccrit = cnf.getfloat('CLUSTERING','ccrit')
        except NoOptionError:
            self.ccrit = 3.2
        try:
            self.kmax = cnf.getint('CLUSTERING','kmax')
        except NoOptionError:
            self.kmax = 10
        try:
            self.Ncmin = cnf.getint('CLUSTERING','Ncmin')
        except NoOptionError:
            self.Ncmin = 5
        try:
            self.linkage = cnf.get('CLUSTERING','linkage')
        except NoOptionError:
            self.linkage = 'ward'

    def write(self,filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file

        """
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)        
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
        with open(filename,"w") as fid:
            cnf.write(fid)
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
        self.default_sections = ['PATHS', 'OPTIONS']
        # initiate the configparser object and read the file
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)        # parse the defaults first
        try:
            self.model = cnf.get('PATHS','model').strip()
        except NoOptionError:
            self.model = 'example%srigo_crustal.nd' % os.sep
        try:
            self.stations = cnf.get('PATHS','stations').strip()
        except NoOptionError:
            self.stations = 'example%swgoc_stations_eida.xml' % os.sep     
        # options
        try:
            self.ainFlag = cnf.getboolean('OPTIONS','ain')
        except NoOptionError:
            self.ainFlag = True

    def write(self,filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file

        """
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section)        
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
        self.default_sections = ['ANGLE', 'BOUNDS', 'GRADING']
        # initiate the configparser object and read the file
        cnf = ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section) 
        # angle
        try:
            self.polOff = cnf.getfloat('ANGLE','diff')
        except NoOptionError:
            self.polOff = 10.
        # error bounds
        try:
            self.snr_bound = cnf.getfloat('BOUNDS','SNR')
        except NoOptionError:
            self.snr_bound = 1.5
        try:
            self.error_bounds = [cnf.getfloat('BOUNDS','phi'),cnf.getfloat('BOUNDS','td')]
        except NoOptionError:
            self.error_bounds = [10., 10.]
        try:
            self.CC_FS_bound = cnf.getfloat('BOUNDS','CC_FS')
        except NoOptionError:
            self.CC_FS_bound = 0.6
        try:
            self.CC_NE_bound = cnf.getfloat('BOUNDS','CC_NE')
        except NoOptionError:
            self.CC_NE_bound = 0.6
        # grading
        try:
            A = cnf.getfloat('GRADING','A')
        except NoOptionError:
            A = 0.25
        try:
            B = cnf.getfloat('GRADING','B')
        except NoOptionError:
            B = 0.50
        try:
            C = cnf.getfloat('GRADING','C')
        except NoOptionError:
            C = 0.75
        try:
            D = cnf.getfloat('GRADING','D')
        except NoOptionError:
            D = 1.00
        self.gradeDict={"A":A,"B":B,"C":C,"D":D}

    def write(self,filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file
                
        """
        cnf=ConfigParser(inline_comment_prefixes="!")
        try:
            with open(filename, 'r') as fid:
                cnf.read_file(fid)
            # check for missing sections
            current_sections = cnf.sections()
            diff_sections = list(set(self.default_sections).difference(set(current_sections)))
            for section in diff_sections:
                cnf.add_section(section)
            # remove sections that exist in the file, but not in the defaults
            diff_sections_file = list(set(current_sections).difference(set(self.default_sections)))
            for section_r in diff_sections_file:
                cnf.remove_section(section)
        except FileNotFoundError:
            for section in self.default_sections:
                cnf.add_section(section) 
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
        with open(filename,"w") as fid:
            cnf.write(fid)
        logging.debug("Saved CLUSTERING settings to %s" % filename)


class ParseFilterCnf():
    """
    class that contains all required settings
    from the config file concerning filters.

    """

    def __init__(self, filename):
        """
        Reads a configuration file and returns the settings.
        
        :type filename: str
        :param filename: the path to the file        

        """
        # initiate the configparser object and read the file
        if not os.path.exists(filename):
        	self.make_default_file(filename)
        with open(filename,"r") as fid:
            cnf_filter_lines = fid.readlines()
        # parse
        filter_array = np.loadtxt(filename, dtype=np.float, delimiter=',')
        self.filter_ranges = filter_array[2:]
        logging.debug('Found %i filters' % self.filter_ranges.shape[0])
        self.filter_preset_1 = filter_array[0][1:]
        logging.debug('Filter preset 1: %s' % str(self.filter_preset_1))
        self.filter_preset_2 = filter_array[1][1:]
        logging.debug('Filter preset 2: %s' % str(self.filter_preset_2))

    def write(self, filename):
        """
        Write values to file.

        :type filename: str
        :param filename: the path to the file

        """
        logging.debug('Writing #%i filters to filter file' % (self.filter_ranges.shape[0] + 2))
        filter_array = np.concatenate(
                                    [
            np.insert(self.filter_preset_1, 0, -2).reshape((1, 3)),
            np.insert(self.filter_preset_2, 0, -1).reshape((1, 3)),
            self.filter_ranges
                                    ])
        np.savetxt(filename, filter_array, fmt='%.0f,%.8f,%.8f')
        logging.debug("Saved Filter settings to %s" % filename)

    def make_default_file(self, filename):
    	"""
        Write a default filters file

        :type filename: str
        :param filename: the path to the file

    	"""
    	logging.debug('Creating filters file %s' % filename)
    	with open(filename, 'w') as fid:
    		fid.writelines(
    						['-2,1.00000000,20.00000000\n',
							 '-1,1.00000000,10.00000000\n',
							 '0,nan,nan\n',
							 '1,0.40000000,4.00000000\n',
							 '2,0.50000000,5.00000000\n',
							 '3,0.20000000,3.00000000\n',
							 '4,0.30000000,3.00000000\n',
							 '5,0.50000000,4.00000000\n',
							 '6,0.60000000,3.00000000\n',
							 '7,0.80000000,6.00000000\n',
							 '8,1.00000000,3.00000000\n',
							 '9,1.00000000,5.00000000\n',
							 '10,1.00000000,8.00000000\n',
							 '11,2.00000000,3.00000000\n',
							 '12,2.00000000,6.00000000\n',
							 '13,3.00000000,8.00000000\n',
							 '14,4.00000000,10.00000000\n',
							 '15,1.00000000,20.00000000\n']
						  )
def path_parser(paths_file, mode="read", catFile=None, dataPath=None, dbPath=None):
    """
    Check whether a file containing the previously used data directory, catalogue file and
    database file paths exists and then either read or write.
    :type paths_file: str
    :param paths_file: Path to the ASCII file containing paths (or the one where the paths
        will be stored in 'read' mode).
    :type mode: str, optional
    :param mode: One of 'read' or 'write', selects the mode with which the function will
        operate. 'read' means reading previous paths and 'write' means (over)writing the 
        new ones to the file. Defaults to 'read'.
    :type catFile: str or None, optional
    :param catFile: path to the catalogue file. If None, this path will be skipped. Defaults to 
        None.
    :type dataPath: str or None, optional
    :param dataPath: path to data directory. If None, this path will be skipped. Defaults to
        None.
    :type dbPath: str or None, optional
    :param dbPath: path to database file. If None, this path will be skipped. Defaults to
        None.
    """
    # set name for paths file
    if mode == "read": # read from existing file
        try:
            with open(paths_file,"r") as fid:
                temp=fid.readlines()
            paths=[x.split(",") for x in temp]
            for path in paths:
                if path[0] == "ERRPATH": # catalogue. 'ERR' is a remnant from old times...
                    catFile=path[1].strip()
                elif path[0] == "DATAPATH": # data directory
                    dataPath=path[1].strip()
                elif path[0] == "DBPATH": # database directory
                    dbPath=path[1].strip()
            # if the file doesn't exist, set current directory
            if not catFile or not os.path.exists(catFile):
                catFile=os.getcwd()
            if not dataPath or not os.path.exists(dataPath):
                dataPath=os.getcwd()
            if not dbPath or not os.path.exists(dbPath):
                dbPath=os.getcwd()
            return dataPath,catFile,dbPath
        except IOError:
            catFile=os.getcwd()
            dataPath=os.getcwd()
            dbPath=os.getcwd()
            return dataPath,catFile,dbPath
    elif mode == "write": # write to the file
        if catFile and dataPath and dbPath:
            with open(paths_file,"w") as fid:
                fid.write("ERRPATH"+","+catFile+"\n")
                fid.write("DATAPATH"+","+dataPath+"\n")
                fid.write("DBPATH"+","+dbPath+"\n")


def read_stations(path_station_file):
    """
    Read the stations information file. Must be either a
    StationXML or an ASCII file in the format:
    station network latitude longitude elevation
    :type path_station_file: str
    :param path_station_file: path to the station information file
    :returns: if the file is StationXML returns a tuple with 
        (a dict containing the station information, 
            the :class:`~obspy.core.stream.Inventory` object),
        otherwise returns (a dict containing the station information, None)
    """
    try:
        # simple check to validate file is of xml type, as obspy's
        # read_inventory sometimes parses non-xml files as xml
        with open(path_station_file, 'r') as fid:
            check_line = fid.readline()
        if not '?xml version' in check_line:
            raise TypeError('Input file is not StationXML')
        # read the xml
        xml = read_inventory(path_station_file)
        inventory = {}
        for net in xml:
            for sta in net:
                inventory.update(
                    {
                        sta.code:{
                        "network":net.code,
                        "latitude":sta.latitude,
                        "longitude":sta.longitude,
                        "elevation":sta.elevation
                                }
                    }
                                )
        return inventory, xml
    except TypeError:
        xml = None
        with open(path_station_file, "r") as fid:
            staLines=fid.readlines()
        return {x.split()[0]:
                    {
                    "network":x.split()[1],
                    "latitude":float(x.split()[2]),
                    "longitude":float(x.split()[3]),
                    "elevation":float(x.split()[4])
                    } for x in staLines
                }, xml
    except:
        logging.exception('Caught exception')
        raise IOError("Could not read station file %s" % path_station_file)