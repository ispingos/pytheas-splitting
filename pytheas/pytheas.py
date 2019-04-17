#!/usr/bin/python
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

LIC_STR=("""

     Pytheas // PYThon sHear - wavE Analysis Suite

  Copyright (C) 2019 Spingos I. & Kaviris G.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
   
    For additional information, questions or suggestions 
    please contact: ispingos@geol.uoa.gr

""")

global VERSION
global VERDATE
VERSION="0.1.0"
VERDATE="15/04/2019"

## first get script's path
import os
global WORKDIR
WORKDIR=os.path.dirname(os.path.realpath(__file__))

## splash screen ##
print("..: Starting module imports...")
# using tk cause it's easier/faster for this (and is in the PSL)
import tkinter as tk
splroot=tk.Tk()
splroot.overrideredirect(True)
# get screen width/height and set as global
global SCWIDTH 
global SCHEIGHT
SCWIDTH=splroot.winfo_screenwidth()
SCHEIGHT=splroot.winfo_screenheight()
# grab image file and get its dimensions
splfile=WORKDIR+os.sep+"etc/images/splash_color.png"
splimage=tk.PhotoImage(file=splfile)
iwidth=splimage.width(); iheight=splimage.height()
# scale it to 40% of screen for better display
scale=SCWIDTH*0.4/iwidth
splimage=splimage.subsample(int(round(1/scale)))
iwidth=splimage.width(); iheight=splimage.height()
# center image
x=(SCWIDTH-iwidth)/2
y=(SCHEIGHT-iheight)/2
splroot.geometry('%dx%d+%d+%d' % (iwidth, iheight, x, y))
# init canvas
canvas=tk.Canvas(splroot,height=iheight,width=iwidth)
canvas.create_image(iwidth/2,iheight/2,image=splimage)
canvas.pack()
# show splash
splroot.update()

## imports ##
print("..: Loading modules...")    
import sys, shutil
import logging, traceback 
import time, datetime
from glob import glob
import numpy as np
import scipy
from PyQt5 import QtWidgets, QtCore, QtGui
# ensure we're using Qt5 as the MPL backend
import matplotlib
matplotlib.use("Qt5Agg")
# continue with imports
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import LinearLocator, NullLocator
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_qt5agg \
import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backend_bases import MouseEvent
import configparser
from obspy import read, Stream, UTCDateTime, read_events, read_inventory
from obspy.geodetics.base import degrees2kilometers, gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.taup.taup_create import get_builtin_model_files as taup_get_builtin_model_files
from obspy.taup.velocity_model import VelocityModel
# pytheas imports
from lib import eigenvalue as SC
from lib import clustering as TB
from lib import rotationcorrelation as RC
from lib import db_handler as DB
from lib import tools
from lib import parsers
# make 'logs' directory
try:
    os.makedirs(WORKDIR+os.sep+"logs")
except: # if folder already exists
    pass

# init the log file
logfile=WORKDIR+os.sep+"logs%s%s_Pytheas.log" % (os.sep,UTCDateTime().strftime("%Y%m%d%H%M%S"))
logging.basicConfig(
                        filename=logfile,
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(message)s'
                   )

# setup logging to terminal
terminal=logging.StreamHandler()
terminal.setLevel(logging.DEBUG)
terminal.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
logging.getLogger().addHandler(terminal)

# kill splash
splroot.destroy()

## Main Class ##
class Pytheas(QtWidgets.QMainWindow):
    """
    Main class of the Pytheas software. Is responsible for building the
    GUI and handling all functionality connected to it.

    Uses a :class:`~PyQt5.QtWidgets.QMainWindow` 

    """
    def __init__(self):
        """Initiliaze main app window and setups labels, flags etc."""
        logging.info("Start Application...")
        # init the main window
        QtWidgets.QMainWindow.__init__(self)
        # check for default dirs. create them if needed
        self.checkDefaultDirs()
        # set icon and initial GUI params
        self.appIcon=WORKDIR+os.sep+"etc/images/Pytheas_ICO.ico"
        self.globalName="Pytheas"
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setWindowTitle(self.globalName+" v"+str(VERSION))
        self.setWindowIcon(QtGui.QIcon(self.appIcon))
        # set validators
        self.validateInt=QtGui.QIntValidator()
        self.validateFloat=QtGui.QDoubleValidator()
        ## set toolbar menus and actions
        self.activateMenus=[]
        # set file menu
        logging.debug("Set File Menu...")
        self.fileMenu=QtWidgets.QMenu("&File",self)
        self.fileMenu.addAction("&Open Catalogue",self.openCatWindow,QtCore.Qt.SHIFT+QtCore.Qt.Key_O)
        self.fileMenu.addAction("&Preferences",self.prefWindow,QtCore.Qt.SHIFT+QtCore.Qt.Key_P)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction("&Quit",self.onlyQuit)        
        self.menuBar().addMenu(self.fileMenu) # Add Open menu to bar
        # set database menu
        logging.debug("Set Database Menu...")
        self.dbMenu=QtWidgets.QMenu("&Database",self)
        self.dbMenu.addAction("Save As...",self.saveAsDB)
        self.dbMenu.addAction("Load New",self.loadNewDB)
        self.dbMenu.addAction("Save Solution to Database",self.saveSolution,QtCore.Qt.CTRL+QtCore.Qt.Key_S)
        self.dbMenu.addAction("Load Solution from Database",self.getSolution)
        self.dbMenu.addAction("Export...",self.dbWindowQuery)
        self.dbMenu.addAction("Export Active",self.exportActive)
        self.dbMenu.addAction("Export Figure",self.exportFigure)
        self.menuBar().addMenu(self.dbMenu) # add db menu bar
        self.dbMenu.setEnabled(False)
        self.activateMenus.append(self.dbMenu)
        # set view menu
        logging.debug("Set View Menu...")
        self.viewMenu=QtWidgets.QMenu("&View",self)
        self.viewMenu.addAction("&Grid",self.changeGrid)
        self.actionGrid=self.viewMenu.actions()[0] # Since only one item so far. Maybe add an if check?
        self.actionGrid.setCheckable(True)
        self.actionGrid.setChecked(False)
        self.viewMenu.addAction("&Center S",self.centerS)
        self.viewMenu.addAction("&Set X Axis Limits",self.applyXAxisLimits)
        self.menuBar().addMenu(self.viewMenu) # add view menu bar
        self.viewMenu.setEnabled(False)
        self.activateMenus.append(self.viewMenu)        
        # set visual inspection menu 
        logging.debug("Set Visual Inspection Menu...")
        self.visualMenu=QtWidgets.QMenu("&Manual",self)
        self.visualMenu.addAction("&Initial Stage",self.goToInitial,QtCore.Qt.CTRL+QtCore.Qt.Key_1)
        self.visualMenu.addAction("&Rotated Stage",self.goToRotated,QtCore.Qt.CTRL+QtCore.Qt.Key_2)
        self.visualMenu.addAction("&Corrected Stage",self.goToCorrected,QtCore.Qt.CTRL+QtCore.Qt.Key_3)
        self.visualMenu.addSeparator()
        self.visualMenu.addAction("&Set phi",self.getPhi)
        self.visualMenu.addAction("&Set time-delay",self.applyTimedelay,QtCore.Qt.CTRL+QtCore.Qt.Key_D)
        self.visualMenu.addAction("&Set grade",self.getGrade,QtCore.Qt.CTRL+QtCore.Qt.Key_X)
        self.visualMenu.addAction("&Set comment",self.getComment,QtCore.Qt.CTRL+QtCore.Qt.Key_G)
        self.menuBar().addMenu(self.visualMenu) # Add Visual menu to bar
        self.visualMenu.setEnabled(False)
        self.activateMenus.append(self.visualMenu)           
        # set other methods menu
        self.splitMenu=QtWidgets.QMenu("&Splitting",self)
        self.splitMenu.addAction("&Set Shear-wave Window",self.getMaxAin)
        self.splitMenu.addSeparator()
        self.actRT=QtWidgets.QAction("&ZRT", checkable=True, checked=True)
        self.actQT=QtWidgets.QAction("&LQT", checkable=True, checked=False)
        self.rotGrp=QtWidgets.QActionGroup(self, exclusive=True)
        self.rotGrp.addAction(self.actRT)
        self.rotGrp.addAction(self.actQT)
        self.splitMenu.addAction(self.actRT)
        self.splitMenu.addAction(self.actQT)
        self.splitMenu.addSeparator()
        self.splitMenu.addAction("&Rotation-Correlation",self.applyRC,QtCore.Qt.CTRL+QtCore.Qt.Key_K)
        self.splitMenu.addAction("&Eigenvalue",lambda: self.applyEV('EV'),QtCore.Qt.CTRL+QtCore.Qt.Key_L)
        self.splitMenu.addAction("&Minimum Energy",lambda: self.applyEV('ME'),QtCore.Qt.CTRL+QtCore.Qt.Key_E)
        self.splitMenu.addSeparator()        
        self.splitMenu.addAction("&Cluster Analysis (EV)",lambda: self.applyCA('EV',self.stream,self.filtered,self.freqmin,self.freqmax,
                                                                      self.sPick,self.baz,self.ain,True,False
                                                                    ),QtCore.Qt.CTRL+QtCore.Qt.Key_T)
        self.splitMenu.addAction("&Cluster Analysis (ME)",lambda: self.applyCA('ME',self.stream,self.filtered,self.freqmin,self.freqmax,
                                                                      self.sPick,self.baz,self.ain,True,False
                                                                    ),QtCore.Qt.CTRL+QtCore.Qt.Key_U)
        self.splitMenu.addAction("&Cluster Analysis (RC)",lambda: self.applyCA('RC',self.stream,self.filtered,self.freqmin,self.freqmax,
                                                                      self.sPick,self.baz,self.ain,True,False
                                                                    ),QtCore.Qt.CTRL+QtCore.Qt.Key_Y)
        self.splitMenu.addSeparator()        
        self.splitMenu.addAction("&Catalogue Cluster Analysis",self.caMultWindow,QtCore.Qt.SHIFT+QtCore.Qt.Key_A)
        self.menuBar().addMenu(self.splitMenu) # Add Splitting menu to bar
        self.splitMenu.setEnabled(False)
        self.activateMenus.append(self.splitMenu)          
        # set navigation menu
        logging.debug("Set Navigation Menu...")
        self.navMenu=QtWidgets.QMenu("&Navigate",self)
        self.navMenu.addAction("&Station List",self.stationsListWindow,QtCore.Qt.SHIFT+QtCore.Qt.Key_S)
        self.navMenu.addAction("&Next Station",self.nextStation,QtCore.Qt.CTRL+QtCore.Qt.Key_W)
        self.navMenu.addAction("&Previous Station",self.prevStation,QtCore.Qt.CTRL+QtCore.Qt.Key_Q)
        self.navMenu.addSeparator()
        self.navMenu.addAction("&Event List",self.eventsListWindow,QtCore.Qt.SHIFT+QtCore.Qt.Key_E)
        self.navMenu.addAction("&Next Event",self.nextEvent,QtCore.Qt.SHIFT+QtCore.Qt.Key_N)
        self.navMenu.addAction("&Previous Event",self.prevEvent,QtCore.Qt.SHIFT+QtCore.Qt.Key_B)
        self.menuBar().addMenu(self.navMenu)
        self.navMenu.setEnabled(False)
        self.activateMenus.append(self.navMenu)            
        # set tools menu
        logging.debug("Set Tools Menu...")
        self.toolsMenu=QtWidgets.QMenu("&Tools",self)
        self.toggleAR=QtWidgets.QAction("&AR-AIC", checkable=True, checked=False)
        self.toggleAR.triggered.connect(self.changeAR)
        self.toolsMenu.addAction(self.toggleAR)
        self.toolsMenu.addSeparator()
        self.toolsMenu.addAction("&Bandpass Filter",self.applyBandFilter,QtCore.Qt.CTRL+QtCore.Qt.Key_B)
        self.toolsMenu.addAction("&Bandpass 1-20 Hz",self.applyFixed120Filter,QtCore.Qt.CTRL+QtCore.Qt.Key_F)
        self.toolsMenu.addAction("&Bandpass 1-10 Hz",self.applyFixed110Filter,QtCore.Qt.CTRL+QtCore.Qt.Key_H)
        self.toolsMenu.addAction("&Recommend Filter",self.applyRecFilter,QtCore.Qt.CTRL+QtCore.Qt.Key_A)
        self.toolsMenu.addAction("&Remove filter",self.removeFilter,QtCore.Qt.CTRL+QtCore.Qt.Key_R)
        self.toolsMenu.addAction("&Show Horizontal Spectra",self.plotSpectrum)
        self.menuBar().addMenu(self.toolsMenu)
        self.toolsMenu.setEnabled(False)
        self.activateMenus.append(self.toolsMenu)     
        # disable actions to avoid errors before loading any data
        for menu in self.activateMenus:
            for action in menu.actions():
                action.setEnabled(False)    
        # set the main widget and layout
        logging.debug("Set Main Widget and Layout...")
        self.mainWidg=QtWidgets.QWidget(self)
        self.vLayout=QtWidgets.QVBoxLayout(self.mainWidg)
        self.vLayout.setAlignment(QtCore.Qt.AlignCenter)
        self.mainWidg.setStyleSheet("background-color:white;")
        # Qt magic
        self.dumWidg=QtWidgets.QWidget()
        # set grid layout for the information widget
        self.gLayout=QtWidgets.QGridLayout()
        self.gLayout.setAlignment(QtCore.Qt.AlignCenter)
        # set fonts for the information widget
        self.lblFont=QtGui.QFont("Calibri",11,QtGui.QFont.Bold)
        self.inpFont=QtGui.QFont("Calibri",11,QtGui.QFont.Normal)
        # add the labels where information will be shown 
        self.lblOrigin=QtWidgets.QLabel("Pytheas v.%s" % VERSION,font=QtGui.QFont("Calibri",14,QtGui.QFont.Bold))
        self.lblOrigin.setAlignment(QtCore.Qt.AlignCenter)
        self.lblStation=QtWidgets.QLabel("V. Date: %s" % VERDATE,font=QtGui.QFont("Calibri",14,QtGui.QFont.Bold))
        self.lblStation.setAlignment(QtCore.Qt.AlignCenter)
        self.lblAIN=QtWidgets.QLabel("incidence ("+u"\u00b0"+"):",font=self.lblFont)
        self.inpAIN=QtWidgets.QLabel("45",font=self.inpFont)
        self.lblBAZ=QtWidgets.QLabel("backazimuth (N"+u"\u00b0"+"E):",font=self.lblFont)
        self.inpBAZ=QtWidgets.QLabel("270",font=self.inpFont)
        self.lblEPI=QtWidgets.QLabel("epicentral (km):",font=self.lblFont)
        self.inpEPI=QtWidgets.QLabel("30",font=self.inpFont)
        self.lblMAG=QtWidgets.QLabel("magnitude:",font=self.lblFont)
        self.lblMAG.setFixedSize(64,20)
        self.inpMAG=QtWidgets.QLabel("1.0",font=self.inpFont)
        self.inpMAG.setFixedSize(32,20)
        self.lblSNR=QtWidgets.QLabel("SNR:",font=self.lblFont)
        self.inpSNR=QtWidgets.QLabel("5.9",font=self.inpFont)
        self.lblSYS=QtWidgets.QLabel("System",font=self.lblFont)
        self.inpSYS=QtWidgets.QLabel("ZNE",font=self.inpFont)
        self.lblV2H=QtWidgets.QLabel("V<H",font=self.lblFont)
        self.inpV2H=QtWidgets.QCheckBox()
        self.inpV2H.setChecked(False)
        self.inpV2H.toggled.connect(self.disableV2HToggle)
        self.lblPHI=QtWidgets.QLabel("φ (N"+u"\u00b0"+"E):",font=self.lblFont)
        self.inpPHI=QtWidgets.QLabel("130",font=self.inpFont)
        self.lblTD=QtWidgets.QLabel("t<sub>d</sub> (ms): ",font=self.lblFont)
        self.inpTD=QtWidgets.QLabel("50",font=self.inpFont)
        self.lblPOL=QtWidgets.QLabel("pol (N"+u"\u00b0"+"E):",font=self.lblFont)
        self.inpPOL=QtWidgets.QLabel("30",font=self.inpFont)
        self.lblGRD=QtWidgets.QLabel("grade: ",font=self.lblFont)
        self.lblGRD.setFixedSize(64,20)
        self.inpGRD=QtWidgets.QLabel("X",font=self.inpFont)
        self.inpGRD.setFixedSize(32,20)
        self.lblFLT=QtWidgets.QLabel("Filter (Hz):",font=self.lblFont)
        self.inpFLT=QtWidgets.QLabel("1.0|20.0",font=self.inpFont)
        self.lblORN=QtWidgets.QLabel("Orient",font=self.lblFont)
        self.inpORN=QtWidgets.QLabel("0.0",font=self.inpFont)
        self.lblMET=QtWidgets.QLabel("Method:",font=self.lblFont)
        self.inpMET=QtWidgets.QLabel("MAN",font=self.inpFont)
        # add labels to grid
        # first int is row, second int is column
        self.gLayout.addWidget(self.lblOrigin,1,5)
        self.gLayout.addWidget(self.lblStation,1,7)
        self.gLayout.addWidget(self.lblAIN,2,1)
        self.gLayout.addWidget(self.inpAIN,2,2)
        self.gLayout.addWidget(self.lblBAZ,2,3)
        self.gLayout.addWidget(self.inpBAZ,2,4)
        self.gLayout.addWidget(self.lblEPI,2,5)
        self.gLayout.addWidget(self.inpEPI,2,6)
        self.gLayout.addWidget(self.lblMAG,2,7)
        self.gLayout.addWidget(self.inpMAG,2,8)
        self.gLayout.addWidget(self.lblSNR,2,9)
        self.gLayout.addWidget(self.inpSNR,2,10)
        self.gLayout.addWidget(self.lblSYS,2,11)
        self.gLayout.addWidget(self.inpSYS,2,12)
        self.gLayout.addWidget(self.lblV2H,2,13)
        self.gLayout.addWidget(self.inpV2H,2,14)
        self.gLayout.addWidget(self.lblPHI,3,1)
        self.gLayout.addWidget(self.inpPHI,3,2)
        self.gLayout.addWidget(self.lblTD,3,3)
        self.gLayout.addWidget(self.inpTD,3,4)
        self.gLayout.addWidget(self.lblPOL,3,5)
        self.gLayout.addWidget(self.inpPOL,3,6)
        self.gLayout.addWidget(self.lblGRD,3,7)
        self.gLayout.addWidget(self.inpGRD,3,8)
        self.gLayout.addWidget(self.lblORN,3,9)
        self.gLayout.addWidget(self.inpORN,3,10)
        self.gLayout.addWidget(self.lblFLT,3,11)
        self.gLayout.addWidget(self.inpFLT,3,12)
        self.gLayout.addWidget(self.lblMET,3,13)
        self.gLayout.addWidget(self.inpMET,3,14)
        # add the grid layout to the dummy widget and add it to
        # the main layout
        self.dumWidg.setLayout(self.gLayout)
        self.vLayout.addWidget(self.dumWidg)
        # make the placeholder MPL plot and add it to the main widget
        logging.debug("Add Matplotlib plot...")
        self.figDict=self.makeFig()
        self.activeFig=self.figDict["FIGURE"]
        self.axZ=self.figDict["TRACE"]["Z"]
        self.axN=self.figDict["TRACE"]["N"]
        self.axE=self.figDict["TRACE"]["E"]
        self.axPolar=self.figDict["POLAR"]
        self.axPart1=self.figDict["HOD1"]
        self.axPart2=self.figDict["HOD2"]
        self.axPart3=self.figDict["HOD3"]
        self.colPal=self.figDict["PALETTE"]
        self.addFig(self.activeFig)
        # read configuration files
        self.generalCNF=parsers.parseGeneralCnf(WORKDIR+os.sep+"etc%soptions%sgeneral.cnf"%(os.sep,os.sep))
        self.pkCNF=parsers.parsePickerCnf(WORKDIR+os.sep+"etc%soptions%spicker.cnf"%(os.sep,os.sep))
        self.caCNF=parsers.parseClusteringCnf(WORKDIR+os.sep+"etc%soptions%sclustering.cnf"%(os.sep,os.sep))
        self.tpCNF=parsers.parseTaupCnf(WORKDIR+os.sep+"etc%soptions%staup.cnf"%(os.sep,os.sep))
        self.gradeCNF=parsers.parseGradeCnf(WORKDIR+os.sep+"etc%soptions%sgrading.cnf"%(os.sep,os.sep))
        # small hack
        self.caCNF.maxTd=self.generalCNF.maxTd
        # clean logs
        self.cleanLogs()
        # read inventory
        try:
            self.inventory,self.xmlInventory=self.readStations(self.tpCNF.stations)
        except: # no station file! show a warning to the user...
            logging.exception("Could not find station file %s" % self.tpCNF.stations)
            winTitle="StationXML Warning"
            genText="Could not read the StationXML!"
            infText="Could not read %s " % self.tpCNF.stations
            self.warnMsgBox(genText,infText,winTitle)
            self.inventory=None
        # set flags and various initial parameters
        self.flagReset()
        self.maxAin=np.inf
        # all done!
        logging.debug("Set focus, central widget and navigation toolbar...")
        self.setCentralWidget(self.mainWidg)

    ## SCIENCE FUNCTIONS ##
    def calcV2H(self,pick,start=-0.5,end=0.5):
        """
        Calculate the vertical to horizontal ratio of amplitudes for both
        channels.

        :type pick: float
        :param pick: the arrival of the S-wave relative to the trace's start time (in s)
        :type start: float, optional
        :param start: start of the window for which the V/H ratio will be calculated,
            relative to the S-arrival (in s). Defaults to -0.5.
        :type end: float, optional
        :param end: end of the window for which the V/H ratio will be calculated,
            relative to the S-arrival (in s). Defaults to 0.5.
        :returns: True if V<H and False if V>H

        """
        ## adjust windows to S pick
        sps=self.stream[0].stats.sampling_rate
        winstart=pick+start
        idxstart=int(np.floor(sps*winstart))
        winend=pick+end
        idxend=int(np.ceil(sps*winend))
        if (idxstart < 0) or (idxend > self.stream[0].stats.npts):
            logging.warning("Not enough samples for V2H based on pick")
            return False
        ## get arrays
        H1="N"; H2="E"
        if self.stageRotated: H1="R"; H2="T"
        Z=self.stream.select(component="Z")[0].data[idxstart:idxend]
        N=self.stream.select(component=H1)[0].data[idxstart:idxend]
        E=self.stream.select(component=H2)[0].data[idxstart:idxend]
        ## return the ratios
        return any((np.abs(Z).max()<np.abs(N).max(),np.abs(Z).max()<np.abs(E).max()))

    def calcSNR(self,pick,winNoise=None,winSignal=None):
        """
        Calculates the Signal-to-Noise Ratio (SNR) around the S-arrival (if picked)
        or the middle of the selected signal window. The equation used is:

        SNR = ( RMSsignal / RMSnoise)**2, where RMS = sqrt(sum(A**2))

        The SNR is acquired for each horizontal component and their mean is defined
        as the final value.
        
        :type pick: float
        :param pick: the arrival of the S-wave relative to the trace's start time (in s)
        :type winNoise: float or None, optional
        :param winNoise: start of the noise window (its end is the S-arrival/middle of window).
            If None, its value is acquired from the Preferences. Defaults to None.
        :type winSignal: float or None, optional
        :param winSignal: end of the signal window (its start is the S-arrival/middle of window).
            If None, its value is acquired from the Preferences. Defaults to None.
        :returns: the average SNR of the horizontal waveforms
    
        """
        # get windows bounds from Preferences?
        if winNoise is None: winNoise=self.generalCNF.snrStart
        if winSignal is None: winSignal=self.generalCNF.snrEnd
        # grab the two horizontal waveforms
        if self.stageRotated:
            M1=self.stream.select(component="R")[0].data
            M2=self.stream.select(component="T")[0].data
        else:
            M1=self.stream.select(component="N")[0].data
            M2=self.stream.select(component="E")[0].data
        # get indices of the signal and noise windows
        idxNoise=int(np.floor((pick+winNoise)*self.stream[0].stats.sampling_rate))
        idxSignl=int(np.ceil((pick+winSignal)*self.stream[0].stats.sampling_rate))
        # however, if the bounds are beyond the length of the waveform...
        if (idxNoise < 0) or (idxSignl > self.stream[0].stats.npts):
            logging.warning("Not enough samples for SNR based on pick")
            return 0
        idxPick=int(pick*self.stream[0].stats.sampling_rate)
        # calculate SNR for N/F channels
        noise=M1[idxNoise:idxSignl]
        signl=M1[idxPick:idxSignl]        
        RMSnoise=np.sqrt(np.mean(noise**2))
        RMSsignl=np.sqrt(np.mean(signl**2))
        SNR1=(RMSsignl/RMSnoise)**2
        # calculate SNR for E/S channels
        noise=M2[idxNoise:idxSignl]
        signl=M2[idxPick:idxSignl]
        RMSnoise=np.sqrt(np.mean(noise**2))
        RMSsignl=np.sqrt(np.mean(signl**2))
        SNR2=(RMSsignl/RMSnoise)**2
        # get average SNR
        SNR=(SNR1+SNR2)/2.0
        if tools.isnone(SNR): # final check, just to be sure
            SNR=0.0
        return SNR
  
    def tdCheck(self,td):
        """
        Checks that the input time delay value is proper according to the
        sampling rate and corrects it if it isn't.

        :type td: float
        :param td: the input time-delay (in s)
        :returns: the corrected time-delay

        """
        logging.debug("input td: "+str(td))
        # get values straight from the trace object
        df=self.stream[0].stats.delta
        sps=self.stream[0].stats.sampling_rate
        # convert dt to samples
        td*=sps 
        # if no modulo then is fine
        if td % sps:
            td=np.ceil(td)
            td*=df # convert back to s
            logging.debug("input td was corrected to "+str(td))
        else:
            td*=df # convert back to s
            logging.debug("input td is correct!")
        return td

    def rotation(self,phi,stream,method):
        """
        Rotate the traces according to the measured Sfast polarization.
        
        :type phi: float
        :param phi: the Sfast polarization direction
        :type stream: :class:`~obspy.core.stream.Stream`
        :param stream: the stream containing the waveforms to be rotated
        :type method: str
        :param method: determines which rotation to perform. Must be either
            "NE->RT" or "RT->NE".
        :returns: the Stream object with the rotated waveforms.

        """
        # first check start times. Differences of 1 sample or less
        # are allowed.
        delta=stream[0].stats.delta
        try:
            startNorth=stream.select(component="N")[0].stats.starttime
            startEast=stream.select(component="E")[0].stats.starttime
        except IndexError:
            startNorth=stream.select(component="R")[0].stats.starttime
            startEast=stream.select(component="T")[0].stats.starttime
        if abs(startNorth - startEast) > delta:
            raise ValueError("Components are not synced!")
        # make the rotation matrix
        M=tools.Rmatrix2D(phi)
        #  according to requested rotation
        if method.upper() == "NE->RT":
            ne=np.asarray((
                 stream.select(component="N")[0].data,
                 stream.select(component="E")[0].data
                          ))
            rt=np.dot(M,ne)
            # assign new values to stream
            stream.select(component="N")[0].data=rt[0]
            stream.select(component="N")[0].stats.channel=\
                      stream.select(component="N")[0].stats.channel[:2]+"R"
            stream.select(component="E")[0].data=rt[1]
            stream.select(component="E")[0].stats.channel=\
                      stream.select(component="E")[0].stats.channel[:2]+"T"
        elif method.upper() ==  "RT->NE":
            rt=np.asarray((
                 stream.select(component="R")[0].data,
                 stream.select(component="T")[0].data
                          ))
            ne=np.dot(M.T,rt)
            # assign new values to stream
            stream.select(component="R")[0].data=ne[0]
            stream.select(component="R")[0].stats.channel=\
                      stream.select(component="R")[0].stats.channel[:2]+"N"
            stream.select(component="T")[0].data=ne[1]
            stream.select(component="T")[0].stats.channel=\
                      stream.select(component="T")[0].stats.channel[:2]+"E"
        return stream
        

    def timedelay(self,td,stream,warning=True):
        """
        Wrapper for adding the time-delay to the slow (transverse) component. Changes
        occur in-place, so a new stream is not returned.

        :type td: float
        :param td: the time-delay
        :type stream: :class:`~obspy.core.stream.Stream`
        :param stream: stream where the time-delay will be added to
        :type warning: bool, optional
        :param warning: selects whether to log a warning if the time-delay
            is attempted to be added in the initial stage. Defaults to True.

        """
        if warning:
            if not self.stageRotated and not self.stageCorrected:
                logging.warning("Traces aren't rotated!")
                return
        try:
            return tools.timedelay(td,stream)
        except IndexError:
            logging.warning("Traces have not been rotated!")

    ## ARCHIVING FUNCTIONS ##
    def checkDefaultDirs(self):
        """Make the default directories."""
        if not os.path.isdir(WORKDIR+os.sep+"etc%sindex"%os.sep):
            os.makedirs(WORKDIR+os.sep+"etc%sindex"%os.sep)
            logging.debug("Created etc%sindex."%os.sep)
        if not os.path.isdir(WORKDIR+os.sep+"etc%soptions"%os.sep):
            os.makedirs(WORKDIR+os.sep+"etc%soptions"%os.sep)
            logging.debug("Created etc%soptions."%os.sep)
        if not os.path.isdir(WORKDIR+os.sep+"etc%sevcor"%os.sep):
            os.makedirs(WORKDIR+os.sep+"etc%sevcor"%os.sep)
            logging.debug("Created etc%sevcor"%os.sep)

    def getTree(self,path):
        """
        Get tree of directories that contain waveforms. It is mandatory for the
        final branch to be named after the event code, e.g. %Y-%m-%d-%H-%M-%S

        :type path: str
        :param path: master directory that contains separate event folders.
        :returns: a list of the event directories found

        """
        fullEvents={}; n=0
        for root,dirs,files in os.walk(path):
            if files and not dirs:
                n+=1
                logging.debug("Adding path to event %i"%n)
                fullEvents[os.path.split(root)[-1]]=root
        return fullEvents

    def getActiveStations(self):
        """
        Gets station names that exist both in the waveform directory of the 
        active event and have arrivals in the catalogue.
        
        :returns: a list of the station names

        """
        evFolder=self.fullEvents[self.activeEvent]
        staList=sorted(
                        self.statsDict[self.activeEvent],
                        key=lambda x:self.statsDict[self.activeEvent][x]["AN"]
                      )
        activeStations=[]
        for station in staList:
            path=evFolder+os.sep+"*."+station+".*"
            try:
                if glob(path)[0]:
                    activeStations.append(station)
            except:
                pass
        return activeStations

    def getActiveEvents(self):
        """
        Gets events that exist both in the data directory and the
        catalogue

        :returns: a list of the active events

        """
        activeEvents=[]
        for event in sorted(self.evsDict):
            try:
                path=self.fullEvents[event]
                if glob(path)[0]:
                    activeEvents.append(event)
            except:
                continue
        return activeEvents

    def readStations(self,staFile):
        """
        Read the stations information file. Must be either a
        StationXML or an ASCII file in the format:
        station network latitude longitude elevation

        :type staFile: str
        :param staFile: path to the station information file
        :returns: if the file is StationXML returns a tuple with 
            (a dict containing the station information, 
             the :class:`~obspy.core.stream.Inventory` object),
            otherwise returns (a dict containing the station information, None)

        """
        try:
            xml=read_inventory(staFile)
            inventory={}
            for net in xml:
                for sta in net:
                    inventory.update({sta.code:{"network":net.code,"latitude":sta.latitude,
                                      "longitude":sta.longitude,"elevation":sta.elevation}})
            return inventory, xml
        except TypeError:
            xml=None
            with open(staFile,"r") as fid: staLines=fid.readlines()
            return {x.split()[0]:
                     {
                       "network":x.split()[1],"latitude":float(x.split()[2]),
                        "longitude":float(x.split()[3]),"elevation":float(x.split()[4])
                     } for x in staLines
                   }, xml
        except:
            raise IOError("Could not read station file %s" % staFile)


    def readEventCat(self,catFile):
        """
        Read a simple catalogue file in the format of:
        (origin time must be in the first 23 columns)
        year mo da hr mn second latit longi depth magnitude
                              ^
                              23
        :type catFile: str
        :param catFile: path to the catalogue file
        :returns: (dict with event information, dict with arrival information)
        
        """    
        # initialize dictionaries
        evsDict={}; statsDict={}
        frmt="%Y-%m-%d-%H-%M-%S"
        # read the file
        with open(catFile,"r") as fid: catLines=fid.readlines()
        for evid,line in enumerate(catLines):
          try:
              evDict=tools.initPytheasDict("event")
              origin=UTCDateTime(line[:23])
              dateKey=origin.strftime(frmt)
              evlat=float(line[23:].strip().split()[0])
              evlon=float(line[23:].strip().split()[1])
              evdep=float(line[23:].strip().split()[2])
              mag=float(line[23:].strip().split()[3])
              evDict["YEAR"]=origin.year
              evDict["Mo"]=origin.month
              evDict["Da"]=origin.day
              evDict["HR"]=origin.hour
              evDict["MN"]=origin.minute
              evDict["SEC"]=origin.second+origin.microsecond/10**6
              evDict["ORIGIN"]=origin
              evDict["LAT"]=evlat; evDict["LON"]=evlon; evDict["DEPTH"]=evdep
              evDict["MAG"]=mag; evDict["evID"]=str(evid)
              # start adding stations
              statsDict.update({dateKey:{}})
              for station in sorted(self.inventory):
                try:
                    staDict=tools.initPytheasDict("station")                
                    # get distance and azimuth
                    dist,az,baz=gps2dist_azimuth(evlat,evlon,
                                                  self.inventory[station]["latitude"],
                                                  self.inventory[station]["longitude"])
                    # calculate ain for direct linear raypath #
                    ain=np.rad2deg(np.arctan((dist/1000.)/evdep))
                    # store values to dictionary
                    staDict['STA']=station
                    staDict['NET']=self.inventory[station]['network']
                    staDict['DIST']=dist/1000.
                    staDict['BAZ']=baz; staDict['AN']=ain
                    statsDict[dateKey].update({station:staDict})
                except:
                    logging.exception("Could not process %s for %s" % (station,dateKey))
                    continue
              # store values
              evsDict.update({dateKey:evDict})
          except:
              logging.exception("Could not process %s" % dateKey)
              continue
        return evsDict,statsDict


    def readQML(self,qmlFile):
        """
        Read a QuakeML file.

        :type catFile: str
        :param catFile: path to the QuakeML file
        :returns: (dict with event information, dict with arrival information)

        """
        # Set event ID before running any loops
        evID=0
        # Read the actual QuakeML file.
        logging.info("Attempting to read %s"%qmlFile)
        try:
          events=read_events(qmlFile)
        except:
          logging.exception("Could not read %s"%qmlFile)
          return
        # Create a dictionary for the event and arrival information.
        evsDict={}; statsDict={}
        for evid,event in enumerate(events):
            try:
                # first, update the evDict object
                prefOrig=event.preferred_origin()
                rid=prefOrig.resource_id
                origin=prefOrig.time
                originZero=UTCDateTime(year=origin.year,month=origin.month,day=origin.day,
                                       hour=origin.hour,minute=origin.minute,second=0)
                lat=prefOrig.latitude
                lon=prefOrig.longitude
                try:
                    dep=prefOrig.depth/1000. # depth should be in km
                except TypeError: # is depth None?
                    dep=1.0
                # set pref mag before searching for corresponding r-id
                mag=event.preferred_magnitude()
                if mag is None:
                    magVal=np.nan
                else:
                    for magRef in event.magnitudes:
                        if magRef.origin_id == rid:
                            mag=magRef
                    magVal=mag.mag
                # parse to evDict object
                evDict={'YEAR':origin.year,'Mo':origin.month,'Da':origin.day,'HR':origin.hour,
                 'MN':origin.minute,'SEC':origin.second+origin.microsecond/10**6, 'ORIGIN':origin,
                 'LAT':lat,'LON':lon,'DEPTH':dep,
                 'MAG':magVal,'evID':str(evid)}
                frmt="%Y-%m-%d-%H-%M-%S"
                dateKey=origin.strftime(frmt)
                # now, need to make the staDict object
                picksRef=event.picks; staDict={}; picks={}
                for arr in prefOrig.arrivals:
                    pid=arr.pick_id
                    phase=arr.phase
                    # get the epicentral distance
                    try: # distance should be written in the file
                        dist=degrees2kilometers(arr.distance)
                    except:
                        dist=None
                    # get the backazimuth
                    try:
                        baz=arr.backazimuth
                    except AttributeError: # if baz is not there, maybe az?
                        try:
                            baz=(arr.azimuth+180) % 360
                        except: # no az either? 
                            baz=None
                    # get the angle of incidence
                    try: # if ain is not written
                        ain=arr.incidence_angle
                    except AttributeError: # maybe the takeoff is
                        try:
                            # calcualte the incidence angle from the takeoff
                            # as seen in Kapetanidis (2017).
                            ih=arr.takeoff_angle # takeoff angle
                            vTop=np.float(self.vmodel.evaluate_below(0,'s')) # velocity at station
                            vSrc=np.float(self.vmodel.evaluate_below(dep,'s')) # velocity at source
                            ain=np.rad2deg(np.arcsin((vTop/vSrc)*np.sin(np.deg2rad(ih))))
                        except:
                            ain=np.nan
                    ## TODO: this should be added as a user-switch
                    ##if dist > degrees2kilometers(1): # accept ONLY local phases. should this be a user option?
                    ##    continue
                    
                    # check for phase and select
                    if phase not in ["p","P","s","S"]: # "P"/"S" accepted only cause of erroneous registration
                        continue                       # of direct upward local phases in the QuakeML
                    for tmp in picksRef: # get the station info for the arrival
                        if tmp.resource_id == pid: # match arrival to pick
                            station=tmp.waveform_id.station_code
                            network=tmp.waveform_id.network_code
                            component=tmp.waveform_id.channel_code
                            ptime=tmp.time
                            if not network:   # can this happen?
                                network="XX"  # should skip in that case?
                            if not station:
                                station="XXX"
                                continue
                            if not component:
                                component="XXX"
                            if any((tools.isnone(ain), tools.isnone(baz), tools.isnone(dist))):
                                # get distance and azimuth
                                if station not in self.inventory:
                                    continue
                                gps=gps2dist_azimuth(
                                                      lat,lon,
                                                      self.inventory[station]["latitude"],
                                                      self.inventory[station]["longitude"]
                                                             )
                                # replace values where needed
                                if tools.isnone(dist):
                                    dist=gps[0]/1000.
                                if tools.isnone(baz):
                                    baz=gps[2]
                                if tools.isnone(ain):
                                    # calculate ain for direct linear raypath
                                    ain=np.rad2deg(np.arctan((dist)/dep))
                                    if tools.isnone(ain): ain = 90.0
                    # get phase information
                    if phase in ['p','P']:
                        temp={'STA':station,'NET':network,'COM':component,
                             'DIST':dist,'BAZ':baz,'AN':ain,
                             'SECP':ptime-originZero,'TOBSP':ptime-origin}
                        if station in staDict:
                            staDict[station].update(temp)                       
                        else:
                          staDict.update({station:temp}) 
                    elif phase in ['s','S']:
                        temp={'STA':station,'NET':network,'COM':component, # re-adding these in case no
                             'DIST':dist,'BAZ':baz,'AN':ain,               # P-arrival was found in the 
                             'SECS':ptime-originZero,'TOBSS':ptime-origin} # QuakeML
                        if station in staDict:
                            staDict[station].update(temp)                       
                        else:
                          staDict.update({station:temp})   
                if not staDict:
                  logging.warning("No picks found for %s" % dateKey)
                  for station in sorted(self.inventory):
                    try:
                        tempDict=tools.initPytheasDict("station")                
                        # get distance and azimuth
                        dist,az,baz=gps2dist_azimuth(lat,lon,
                                                      self.inventory[station]["latitude"],
                                                      self.inventory[station]["longitude"])
                        # calculate ain for direct linear raypath #
                        if dep == 0.0: dep = 1.0 # hack to permit ain calc
                        ain=np.rad2deg(np.arctan((dist/1000.)/dep))
                        if tools.isnone(ain): ain = 0.0
                        # store values to dictionary
                        tempDict['STA']=station
                        tempDict['NET']=self.inventory[station]['network']
                        tempDict['DIST']=dist/1000.
                        tempDict['BAZ']=baz
                        tempDict['AN']=ain
                        staDict.update({station:tempDict})
                    except:
                        logging.exception("Could not process %s for %s" % (station,dateKey))
                        continue
                evsDict.update({dateKey:evDict})                            
                evsDict[dateKey]['NSTA']=len(staDict)
                statsDict.update({dateKey:staDict})
            except:
                logging.exception("Could not parse %s from catalogue, skipping..." % evid)
        return evsDict,statsDict

    def eventCorrect(self,evPath,evsDict,statsDict):
       """
       Corrects the keys of events and stations dictionaries to be
       in accordance with event ids acquired from the folder search.

       :type evPath: tuple-like
       :param evPath: tuple-like of event codes found in the data directory
       :type evsDict: dict
       :param evsDict: dict containing all the event information acquired from the
            catalogue
        :type statsDict: dict
        :param statsDict: dict containing all the arrival information acquired from the
            catalogue
        :returns: (event information dict, arrival information dict) both with corrected
            event code keys (i.e. the origin times from the catalogue were replaced with
            the matched ones from the data directory)

       """
       tOff=self.generalCNF.matching
       keyTrash=[]
       logging.debug("Performing keys correction...")
       # get previous corrections
       try: # if it exists
           with open(WORKDIR+os.sep+"etc%sevcor%sevcor.dat"%(os.sep,os.sep)) as fid:
               corList=fid.readlines()
               corDict={x.split()[0].strip():x.split()[1].strip() for x in corList}
       except IOError: # if file doesn't exist
           with open(WORKDIR+os.sep+"etc%sevcor%sevcor.dat"%(os.sep,os.sep),"w") as fid:
               fid.write("# INITIAL CORRECTED \n")
           corDict={}
       # get keys. catalogue keys are rounded to 0 microseconds
       catalogueKeys=np.asarray([round(UTCDateTime(x).timestamp) for x in sorted(evsDict)])
       directoryKeys=np.asarray([UTCDateTime(y).timestamp for y in sorted(evPath)])
       keyFrmt="%Y-%m-%d-%H-%M-%S"
       # we want to match the CATALOGUE keys to the DIRECTORY keys
       for eTms in catalogueKeys:
           diff=directoryKeys-eTms
           if abs(diff).min() == 0:
               continue
           # get indices
           didx=np.abs(diff)<=tOff
           if np.count_nonzero(didx) == 1: # only 1 solution found, good!
               pTms=directoryKeys[didx]
           elif np.count_nonzero(didx) > 1:
               didx2=np.argmin(diff[didx])
               pTms=directoryKeys[didx][didx2]
           else: # no match found
               continue
           pKey=UTCDateTime(pTms).strftime(keyFrmt)              
           eKey=UTCDateTime(eTms).strftime(keyFrmt)
           evsDict[pKey]=evsDict[eKey]
           statsDict[pKey]=statsDict[eKey]
           keyTrash.append(eKey)            
           with open(WORKDIR+os.sep+"etc%sevcor%sevcor.dat"%(os.sep,os.sep),"a") as fid:
               fid.write(str(eKey)+" "+str(pKey)+"\n")
           logging.debug("Corrected key: "+eKey+" -> "+pKey)
       # Clear dictionaries of unused keys
       for tKey in keyTrash:
           if tKey in evsDict:
               _=evsDict.pop(tKey)
           if tKey in statsDict:
               _=statsDict.pop(tKey)
       return evsDict,statsDict

    def pathParser(self,mode="read",catFile=None,dataPath=None,dbPath=None):
        """
        Check whether a file containing the previously used data directory, catalogue file and
        database file paths exists and then either read or write.

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
        exPathFile=WORKDIR+os.sep+"etc%soptions%spaths.exp"%(os.sep,os.sep)
        if mode == "read": # read from existing file
            try:
                with open(exPathFile,"r") as fid:
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
                with open(exPathFile,"w") as fid:
                    fid.write("ERRPATH"+","+catFile+"\n")
                    fid.write("DATAPATH"+","+dataPath+"\n")
                    fid.write("DBPATH"+","+dbPath+"\n")

    def updateSplittingDict(self,splitting=True):
        """
        Update the splitting dictionary. Values are held in ``self``.

        :type splitting: bool, optional
        :param splitting: enables updating splitting values (i.e. the method layer)
            in the dictionary. Defaults to True.

        """
        logging.debug("Updating splitting dict with combination: %s/%s/%s" % (self.activeEvent,self.station,self.method))
        # for each layer (event,station,method) initialize a default dictionary, if the layer does not exist
        # in the dictionary
        if self.activeEvent not in self.splittingDict: 
            self.splittingDict[self.activeEvent]=tools.initSplittingDict("event")
        if self.station not in self.splittingDict[self.activeEvent]: 
            self.splittingDict[self.activeEvent][self.station]=tools.initSplittingDict("station")
        if self.method not in self.splittingDict[self.activeEvent][self.station]:
            self.splittingDict[self.activeEvent][self.station][self.method]=tools.initSplittingDict("method")
        # update the event layer
        self.splittingDict[self.activeEvent]["origin"]=self.evsDict[self.activeEvent]["ORIGIN"].datetime
        self.splittingDict[self.activeEvent]["latitude"]=float(self.evsDict[self.activeEvent]["LAT"])
        self.splittingDict[self.activeEvent]["longitude"]=float(self.evsDict[self.activeEvent]["LON"])
        self.splittingDict[self.activeEvent]["depth"]=float(self.evsDict[self.activeEvent]["DEPTH"])
        self.splittingDict[self.activeEvent]["magnitude"]=float(self.evsDict[self.activeEvent]["MAG"])
        # update the station layer
        self.splittingDict[self.activeEvent][self.station]["station"]=self.station
        self.splittingDict[self.activeEvent][self.station]["network"]=self.statsDict[self.activeEvent][self.station]["NET"]
        self.splittingDict[self.activeEvent][self.station]["epicentral"]=float(self.statsDict[self.activeEvent][self.station]["DIST"])
        self.splittingDict[self.activeEvent][self.station]["azimuth"]=float(self.statsDict[self.activeEvent][self.station]["BAZ"])
        self.splittingDict[self.activeEvent][self.station]["incidence"]=float(self.statsDict[self.activeEvent][self.station]["AN"])
        self.splittingDict[self.activeEvent][self.station]["s_p"]=self.tpts
        self.splittingDict[self.activeEvent][self.station]["orientation"]=self.orCorrection
        # update arrivals and picks
        if not self.sPick:
            try:
                self.sArr=self.spickDict[self.activeEvent][self.station]
            except:
                pass
        else:
            self.sArr=self.stream[0].stats.starttime+self.sPick  
        if self.sArr not in [self.stream[0].stats.starttime,0,np.nan]:      
            self.splittingDict[self.activeEvent][self.station]["s_obs"]=self.sArr.datetime
        if self.sArrTheor not in [self.stream[0].stats.starttime,0,np.nan]:      
            self.splittingDict[self.activeEvent][self.station]["s_theo"]=self.sArrTheor.datetime
        if self.sArrAuto not in [self.stream[0].stats.starttime,0,np.nan]:      
            self.splittingDict[self.activeEvent][self.station]["s_auto"]=self.sArrAuto.datetime
        # update the method layer
        if splitting:
            self.splittingDict[self.activeEvent][self.station][self.method]["station"]=self.station
            self.splittingDict[self.activeEvent][self.station][self.method]["method"]=self.method
            self.splittingDict[self.activeEvent][self.station][self.method]["network"]=self.statsDict[self.activeEvent][self.station]["NET"]
            self.splittingDict[self.activeEvent][self.station][self.method]["SNR"]=self.SNR
            self.splittingDict[self.activeEvent][self.station][self.method]["s_freq"]=self.dPer
            self.splittingDict[self.activeEvent][self.station][self.method]["phi"]=self.phi%180
            self.splittingDict[self.activeEvent][self.station][self.method]["td"]=abs(self.td)*(10**3)
            self.splittingDict[self.activeEvent][self.station][self.method]["pol"]=self.pol%180
            self.splittingDict[self.activeEvent][self.station][self.method]["CC_FS"]=self.CC_FS
            self.splittingDict[self.activeEvent][self.station][self.method]["CC_NE"]=self.CC_NE
            self.splittingDict[self.activeEvent][self.station][self.method]["var_T"]=self.Tvar
            self.splittingDict[self.activeEvent][self.station][self.method]["err_phi"]=self.errors[0]
            self.splittingDict[self.activeEvent][self.station][self.method]["err_td"]=self.errors[1]*(10**3)
            self.splittingDict[self.activeEvent][self.station][self.method]["N_contours"]=self.nContours
            self.splittingDict[self.activeEvent][self.station][self.method]["grade_score"]=self.grade_score
            self.splittingDict[self.activeEvent][self.station][self.method]["grade"]=self.grade
            self.splittingDict[self.activeEvent][self.station][self.method]["filter_min"]=self.freqmin
            self.splittingDict[self.activeEvent][self.station][self.method]["filter_max"]=self.freqmax
            if self.method[1:] in ["EV","RC","ME"]:
                try:
                    self.splittingDict[self.activeEvent][self.station][self.method]["window_min"]=(self.stream[0].stats.starttime+self.minPick).datetime
                    self.splittingDict[self.activeEvent][self.station][self.method]["window_max"]=(self.stream[0].stats.starttime+self.maxPick).datetime
                except:
                    pass
            if self.comment != "":
                self.splittingDict[self.activeEvent][self.station][self.method]["comment"]=self.comment
            self.splittingDict[self.activeEvent][self.station][self.method]["phi_test"]=np.asarray(self.pTest)
            self.splittingDict[self.activeEvent][self.station][self.method]["td_test"]=np.asarray(self.dTest)
            self.splittingDict[self.activeEvent][self.station][self.method]["C_max"]=self.cEmax
            self.splittingDict[self.activeEvent][self.station][self.method]["C_array"]=np.asarray(self.cArray)
            self.splittingDict[self.activeEvent][self.station][self.method]["initial_clusters"]=np.asarray(self.initialClusters)
            self.splittingDict[self.activeEvent][self.station][self.method]["initial_clusters_errors"]=np.asarray(self.initialClustersErrors)
            self.splittingDict[self.activeEvent][self.station][self.method]["calinski_score"]=np.asarray(self.calinski)
            self.splittingDict[self.activeEvent][self.station][self.method]["clusters1"]=np.asarray(self.clusters1)
            self.splittingDict[self.activeEvent][self.station][self.method]["clusters2"]=np.asarray(self.clusters2)


    def updateFromSplittingDict(self,method="MAN"):
        """
        Store values from splitting dictionary to the appropriate objects for the
        given method.

        :type method: str, optional
        :param method: one of the available splitting methods. Defaults to 'MAN'
        
        """
        # checks
        if self.activeEvent not in self.splittingDict:
            return
        else:
            if self.station not in self.splittingDict[self.activeEvent]:
                return
            else:
                if self.method not in self.splittingDict[self.activeEvent][self.station]:
                    return
        # station layer
        self.tpts=self.splittingDict[self.activeEvent][self.station]["s_p"]
        # method layer
        metDict=self.splittingDict[self.activeEvent][self.station][method]
        # get observed arrival
        try:
            sArr=UTCDateTime(self.splittingDict[self.activeEvent][self.station]["s_obs"])
        except:
            sArr=0.
        if sArr not in [self.stream[0].stats.starttime,0,np.nan]:
            self.sArr=UTCDateTime(sArr)
            self.spickDict[self.activeEvent][self.station]=self.sArr
            self.sPick=self.sArr-self.stream[0].stats.starttime
        else:
            try:
                self.sArr=self.spickDict[self.activeEvent][self.station]
                self.sPick=self.sArr-self.stream[0].stats.starttime
            except:
                pass
        # get theoretical arrival
        if self.sArrTheor in [self.stream[0].stats.starttime,0,np.nan]:
            sArrTheor=self.splittingDict[self.activeEvent][self.station]["s_theo"]
            if sArrTheor not in [self.stream[0].stats.starttime,0,np.nan]:
                self.sArrTheor=UTCDateTime(sArrTheor)
                self.sPickTheor=self.sArrTheor-self.stream[0].stats.starttime
        # get automatic arrrival
        if self.sArrAuto in [self.stream[0].stats.starttime,0,np.nan]:
            sArrAuto=self.splittingDict[self.activeEvent][self.station]["s_auto"]
            if sArrAuto not in [self.stream[0].stats.starttime,0,np.nan]:
                self.sArrAuto=UTCDateTime(sArrAuto)
                self.sPickAuto=self.sArrAuto-self.stream[0].stats.starttime
        # get rest
        self.SNR=metDict["SNR"]
        self.dPer=metDict["s_freq"]
        self.phi=metDict["phi"]
        self.td=-metDict["td"]/(10**3)
        self.pol=metDict["pol"]
        self.CC_FS=metDict["CC_FS"]
        self.CC_NE=metDict["CC_NE"]
        self.Tvar=metDict["var_T"]
        self.errors=(metDict["err_phi"],metDict["err_td"]/(10**3))
        self.nContours=metDict["N_contours"]
        self.grade_score=metDict["grade_score"]
        self.grade=metDict["grade"]
        self.freqmin=metDict["filter_min"]
        self.freqmax=metDict["filter_max"]
        if not all((tools.isnone(self.freqmin),tools.isnone(self.freqmax))):
            self.filtered=True
        try:
            self.minPick=UTCDateTime(metDict["window_min"])-self.stream[0].stats.starttime
            self.maxPick=UTCDateTime(metDict["window_max"])-self.stream[0].stats.starttime
        except:
            self.minPick=np.nan; self.maxPick=np.nan
        self.comment=metDict["comment"]
        self.pTest=self.splittingDict[self.activeEvent][self.station][self.method]["phi_test"]
        self.dTest=self.splittingDict[self.activeEvent][self.station][self.method]["td_test"]
        self.cEmax=self.splittingDict[self.activeEvent][self.station][self.method]["C_max"]
        self.cArray=self.splittingDict[self.activeEvent][self.station][self.method]["C_array"]
        self.initialClusters=self.splittingDict[self.activeEvent][self.station][self.method]["initial_clusters"]
        self.initialClustersErrors=self.splittingDict[self.activeEvent][self.station][self.method]["initial_clusters_errors"]
        self.calinski=self.splittingDict[self.activeEvent][self.station][self.method]["calinski_score"]
        self.clusters1=self.splittingDict[self.activeEvent][self.station][self.method]["clusters1"]
        self.clusters2=self.splittingDict[self.activeEvent][self.station][self.method]["clusters2"]

    def writeSplittingDict(self,splittingDict,showInfo=True):
        """
        Write the splitting information to an ASCII format.

        :type splittingDict: dict
        :param splittingDict: the dict containing the splitting information that will be written
            into the file.
        :type showInfo: bool, optional
        :param showInfo: enables showing an info dialog to the user when the writing process is
            completed. Defaults to True.

        """
        # set the database's path as default
        splname=os.path.splitext(self.dbPath)[0]
        # ask the user for a new name
        splname,_=QtWidgets.QFileDialog.getSaveFileName(self,"Select the destination file:",splname,
            "Comma-separated file (*.csv);; Pytheas splitting file (*.spl);; All files (*)")        
        # initiate an empty list and gradually fill it with the required info
        logging.debug("Writing results to "+splname)
        writeList=[]
        # add header to writelist
        if splname.lower().endswith("spl"):
            hdr=str(
                        "YY mm dd HH MM SS.SSS".rjust(23)+"Latitude".rjust(9)+"Longitude".rjust(10)+
                        "Depth".rjust(6)+
                        "Mag".rjust(4)+"Delta".rjust(8)+"Baz".rjust(6)+"Inc".rjust(6)+
                        "phi".rjust(7)+"td".rjust(7)+"pol".rjust(7)+"Stat".rjust(5)+
                        "Meth".rjust(5)+"SNR".rjust(5)+"err_phi".rjust(10)+"err_td".rjust(10)+
                        "CC_FS".rjust(8)+"CC_NE".rjust(8)+"Score".rjust(7)+"Orient".rjust(6)+
                        3*" "+"G"+" # "+"Comment"
                    )
        elif splname.lower().endswith("csv"):
            hdr=",".join(("Origin Time","Latitude","Longitude",
                           "Depth","Magnitude","Epicentral (km)","Backazimuth","Incidence",
                           "phi","td (ms)","pol","Station","Method",
                           "SNR","phi error","td error (ms)","CC FS","CC NE",
                           "Score","Grade","Orientation","Comment"))
        writeList.append(hdr+"\n")
        # event details
        frmt="%Y %m %d %H %M %S.%f" # format for origin information
        for activeEvent in sorted(splittingDict):
            evLayer=splittingDict[activeEvent]
            date=UTCDateTime(evLayer["origin"]).strftime(frmt)[:-3]
            latitude="{:.4f}".format(evLayer["latitude"])
            longitude="{:.4f}".format(evLayer["longitude"])
            depth=evLayer["depth"]
            magnitude="{:.1f}".format(evLayer["magnitude"])
            # station details
            for station in sorted(evLayer):
                stLayer=evLayer[station]
                if not isinstance(stLayer,dict): continue # if object is not a dictionary, it's not a station layer
                depth="{:.1f}".format(float(depth)) # reset it for next station
                logging.debug("Writing station "+station)
                network=stLayer["network"]
                distance="{:.2f}".format(stLayer["epicentral"])
                backazimuth="{:.1f}".format(stLayer["azimuth"])
                ain="{:.1f}".format(stLayer["incidence"])
                orientation="{:.1f}".format(stLayer["orientation"])
                for method in sorted(stLayer):
                    meLayer=stLayer[method]
                    if not isinstance(meLayer,dict): continue # same as above
                    grade=meLayer["grade"]
                    if grade in ("F","X"): continue # skip failed attempts
                    phi="{:.1f}".format(meLayer["phi"])
                    td="{:.1f}".format(meLayer["td"])
                    pol="{:.1f}".format(meLayer["pol"])
                    SNR="{:.1f}".format(meLayer["SNR"])
                    cc_fs="{:.2f}".format(meLayer["CC_FS"])
                    cc_ne="{:.2f}".format(meLayer["CC_NE"])
                    sphi="{:.2f}".format(meLayer["err_phi"])
                    sdt="{:.3f}".format(meLayer["err_td"])
                    score="{:.3f}".format(meLayer["grade_score"])
                    comment=meLayer["comment"]
                    # convert to str
                    if splname.lower().endswith("spl"):
                        line=str(
                                    date.rjust(23)+latitude.rjust(9)+longitude.rjust(10)+depth.rjust(6)+
                                    magnitude.rjust(4)+distance.rjust(8)+backazimuth.rjust(6)+ain.rjust(6)+
                                    phi.rjust(7)+td.rjust(7)+pol.rjust(7)+station.rjust(5)+
                                    method.rjust(5)+SNR.rjust(5)+sphi.rjust(10)+sdt.rjust(10)+
                                    cc_fs.rjust(8)+cc_ne.rjust(8)+score.rjust(7)+orientation.rjust(6)+
                                    3*" "+grade+" # "+comment
                                )
                    elif splname.lower().endswith("csv"):
                        line=",".join((str(UTCDateTime(evLayer["origin"]))[:-3],latitude,longitude,depth,magnitude,distance,backazimuth,ain,phi,
                                        td,pol,station,method,SNR,sphi,sdt,cc_fs,cc_ne,score,grade,orientation,comment))
                    # append to list for writing in event file
                    writeList.append(line+"\n")
        # write to event file
        with open(splname,"w") as fid: fid.writelines(writeList)
        # let the user know
        winTitle="Database Export"
        genText="Successfully exported results from database!"
        infText="Saved results to %s" % splname
        self.infoMsgBox(genText,infText,winTitle)

    ## FIGURE FUNCTIONS ##
    def makeFig(self):
        """Makes the initial Figure instance"""
        # TODO: add different themes for the user to selet in the Preferences.
        #       Who wouldn't like a dark/night theme for working?
        colPal={"fcl":"white","traceZ":"red","traceN":"blue","traceE":"black",
                "polarcl":"blue","partcl":"magenta"}
        # init the figure
        fig=Figure()
        fig.set_facecolor(colPal["fcl"])
        ## add subplots. these share the x axis
        # waveforms
        axZ=fig.add_subplot(5,1,1) # vertical trace 
        axN=fig.add_subplot(5,1,2,sharex=axZ,sharey=axZ) # horizontal 1
        axE=fig.add_subplot(5,1,3,sharex=axZ,sharey=axZ) # horizontal 2
        axTrace={"Z":axZ,"N":axN,"E":axE}
        # particle motion diagrams
        axPolar=fig.add_subplot(5,1,4,sharex=axZ) # polarigram, sharing x with waveforms
        axPart1=fig.add_subplot(5,3,13,autoscale_on=False) # left hodogram
        axPart2=fig.add_subplot(5,3,14,sharey=axPart1,autoscale_on=False) # mid hodogram
        axPart3=fig.add_subplot(5,3,15,sharey=axPart1,autoscale_on=False) # right hodogram
        # adjust the figure's margins to better fit in the widget
        fig.subplots_adjust(left=0.02,right=0.98,top=0.99,bottom=0.02)
        # make some dummy values for the initial plot. Straight lines @ y=0
        x=range(0,80); y=np.zeros(len(x))
        ## plot dummies and set initial labels
        # For the Trace subplots
        fullName={"Z":"VERTICAL","N":"NORTH","E":"EAST"}
        for k in ["Z","N","E"]:
            ax=axTrace[k]
            ax.plot(x,y,colPal["trace"+k])
            ax.set_ylabel(fullName[k])
        # For the Polarigram subplot
        axPolar.plot(x,y,colPal["polarcl"],label="POLAR")
        axPolar.set_ylabel("POLAR")
        axPolar.set_xlabel("East")
        # For the hodogram subplots
        axPart1.plot(x,y,colPal["partcl"],label="HOD1")
        axPart2.plot(x,y,colPal["partcl"],label="HOD2")
        axPart3.plot(x,y,colPal["partcl"],label="HOD3")
        # set all ticks inside and disable labels at the left
        for ax in fig.get_axes():
            ax.tick_params(labelleft=False,left=False,which="both",direction="in")
            if ax not in (axPolar,axE): # disable bottom labels for vertical and horizontal 1 plots
                ax.tick_params(labelbottom=False)
            if self.actionGrid.isChecked(): # show a grid?
                ax.grid()
        # set labels and adjust ticks for particle motion diagrams
        for ax in (axPart1,axPart2,axPart3,axPolar):
            ax.set_ylabel("North")
            ax.set_xlabel("East")
            ax.tick_params(top=False,bottom=False,left=False,right=False)
        logging.debug("Figure Created.")
        # add dummy picks
        for ax in (axZ,axN,axE,axPolar):
            # 1 -> obs, 2 -> auto, 3 -> theor
            ax.axvline(x=0,linewidth=2.0,color='red',label='observed')
            ax.axvline(x=0,color='darkviolet',linewidth=1.5,linestyle='-',label='autopicker') # picker line
            ax.axvline(x=0,color='darkviolet',linewidth=1.5,linestyle='--',label='theoretical') # theoretical line
            ax.set_xlim(0,60)
        # adjust legend of picks. only shown at the vertical trace!
        axZ.legend(fancybox=True,loc='center left',shadow=True,framealpha=1.0)
        # store figure objects for easier retrieval later
        figDict={
                    "FIGURE":fig, "TRACE":axTrace, "POLAR":axPolar,
                    "HOD1":axPart1,"HOD2":axPart2,"HOD3":axPart3,
                    "PALETTE":colPal
                }
        self.flagIniPlot=True # this means that the figure is the initial placeholder
        return figDict

    def changeGrid(self):
        """
        Changes the grid, depending of whether the grid action
        is checked or not.
        
        """
        if self.actionGrid.isChecked(): # show grid
            for ax in self.activeFig.get_axes():
                ax.grid(True)
        elif not self.actionGrid.isChecked(): # don't show grid
            for ax in self.activeFig.get_axes():
                ax.grid(False)
        self.canvas.draw()

    def centerS(self):
        """Centers view on S pick"""
        xmin=self.sPick-1. # arrival +/- 1 second
        xmax=self.sPick+1.
        self.axPolar.set_xlim((xmin,xmax)) # setting only to the polarigram, means 
        self.canvas.draw()                 # it'll change for all plots. gotta love sharex!

    def addFig(self,fig):
        """
        Adds an existing figure to the main widget canvas.
        
        :type fig: :class: `~matplotlib.figure.Figure`
        :param fig: the figure object to be added to the main widget

        """
        # set the figure to the canvas
        self.canvas=FigureCanvas(fig)
        # make it resizable
        self.canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
        self.canvas.updateGeometry()
        # add the new canvas to the widget
        self.vLayout.addWidget(self.canvas)
        # connect mouse actions
        self.mplConnect(fig)
        # final draw
        self.canvas.draw()
        logging.info("Added Figure.")

    def updateFig(self,center=True,keepAxisLimits=False):
        """
        Updates the current figure. This includes reprocessing the initial figure
        and making the particle motion diagrams.

        :type center: bool, optional
        :param center: center the figure. Defaults to True.
        :type keepAxisLimits: bool, optional
        :param keepAxisLimits: do not change the limits of the x axis.
            Defaults to False.

        """
        # backup stream or load the initial one
        if self.streamIni:
            self.stream=self.streamIni.copy()
        else:
            self.streamIni=self.stream.copy()
        # should we add time-delay? Also calculate the correlation coefficient
        # for corrected FS or NE while at it!
        if self.stageRotated:
            self.stream=self.rotation(self.phi,self.stream,"NE->RT")
            self.streamRotated=self.stream.copy()
            if self.td != 0 and (self.stageRotated or self.stageCorrected):
                self.applyTimedelay(arrowInput=True)
            # if there's a pick and the method is MAN, calculate the correlation
            # coefficient for the corrected F/S
            if self.sPick not in [0,np.nan] and self.method == "MAN":
                ccIdx1=int(np.floor((self.sPick-0.2)*self.stream[0].stats.sampling_rate))
                ccIdx2=int(np.ceil((self.sPick+0.2)*self.stream[0].stats.sampling_rate))
                self.CC_FS=tools.xcStatic(self.stream.select(component="R")[0][ccIdx1:ccIdx2],\
                                          self.stream.select(component="T")[0][ccIdx1:ccIdx2])
                logging.info("Correlation in the FS system is %.3f" % self.CC_FS)
        elif self.stageCorrected:
            # this works as follows: starting from R/T -> add time-delay ->
            # calculate CC @ F/S -> rotate to N/E -> calculate CC @ N/E
            if self.td != 0 and (self.stageRotated or self.stageCorrected):
                self.applyTimedelay(arrowInput=True)
            # if there's a pick and the method is MAN, calculate the correlation
            # coefficient for the corrected F/S and then N/E.
            if self.sPick not in [0,np.nan] and self.method == "MAN":
                ccIdx1=int(np.floor((self.sPick-0.2)*self.stream[0].stats.sampling_rate))
                ccIdx2=int(np.ceil((self.sPick+0.2)*self.stream[0].stats.sampling_rate))
                self.CC_FS=tools.xcStatic(self.stream.select(component="R")[0][ccIdx1:ccIdx2],\
                                          self.stream.select(component="T")[0][ccIdx1:ccIdx2])
            logging.info("Correlation in the FS system is %.3f" % self.CC_FS)
            # rotate to N/E
            self.stream=self.rotation(self.phi,self.stream,"RT->NE")
            # calculate the CC @ N/E
            if self.sPick and self.method == "MAN":
                ccIdx1=int(np.floor((self.sPick-0.2)*self.stream[0].stats.sampling_rate))
                ccIdx2=int(np.ceil((self.sPick+0.2)*self.stream[0].stats.sampling_rate))
                self.CC_NE=tools.xcStatic(self.stream.select(component="N")[0][ccIdx1:ccIdx2],\
                                          self.stream.select(component="E")[0][ccIdx1:ccIdx2])
            logging.info("Correlation in the NE system is %.3f" % self.CC_NE)
            self.streamCorrected=self.stream.copy()
        # change axis limits?
        if keepAxisLimits:
            bottom,top=self.axN.get_ybound()
            left,right=self.axN.get_xbound()
        # store trace for each channel
        if self.stageRotated:
            vert=self.stream.select(component="Z")[0]
            nort=self.stream.select(component="R")[0]
            east=self.stream.select(component="T")[0]
        else:
            vert=self.stream.select(component="Z")[0]
            nort=self.stream.select(component="N")[0]
            east=self.stream.select(component="E")[0]           
        # Filter?
        if self.filtered:
            vert.filter(type="bandpass",freqmin=self.freqmin,freqmax=self.freqmax,zerophase=True,corners=4)
            nort.filter(type="bandpass",freqmin=self.freqmin,freqmax=self.freqmax,zerophase=True,corners=4)
            east.filter(type="bandpass",freqmin=self.freqmin,freqmax=self.freqmax,zerophase=True,corners=4)
            self.V2H=np.asarray(self.calcV2H(self.sPick))
        # get relative times
        tV=vert.times(type="relative")
        tN=nort.times(type="relative")
        tE=east.times(type="relative")
        ## update data in traces subplots
        # vertical
        self.axZ.lines[0].set_xdata(tV)
        self.axZ.lines[0].set_ydata(vert.data)
        self.axZ.set_ylabel(vert.stats.channel)
        # recenter plots?
        if center:
            self.centerAxes(self.axZ)
        # horizontal 1
        self.axN.lines[0].set_xdata(tN)
        self.axN.lines[0].set_ydata(nort.data)
        # set suitable label (i.e. if R set F)
        if not self.stageRotated:
            self.axN.set_ylabel(nort.stats.channel)
        else:
            self.axN.set_ylabel(nort.stats.channel[:2]+"F")
        self.axN.set_xlabel("Relative Time (s)")
        if center:
            self.centerAxes(self.axN)
        # horizontal 2
        self.axE.lines[0].set_xdata(tE)
        self.axE.lines[0].set_ydata(east.data)
        # set suitable label (i.e. if T set S)
        if not self.stageRotated:
            self.axE.set_ylabel(east.stats.channel)
        else:
            self.axE.set_ylabel(east.stats.channel[:2]+"S")       
        if center:
            self.centerAxes(self.axE)
        # make the polarigram
        self.updatePolarigram()
        # make the hodogram plot
        self.updateHodogram()
        # draw picks
        logging.debug("Drawing picks...")
        # first draw (if there's one) the selected signal window
        for vspan in self.vspanList:
            try:
                vspan.remove()
            except ValueError:
                continue
        # now draw the picks proper
        for ax in (self.axZ,self.axN,self.axE):
            ax.lines[1].set_xdata(self.sPick)
            ax.lines[2].set_xdata(self.sPickAuto)
            ax.lines[3].set_xdata(self.sPickTheor)
            self.vspanList.append(ax.axvspan(self.minPick,self.maxPick,color='green',alpha=0.3))
        for i,arr in enumerate((self.sPick,self.sPickAuto,self.sPickTheor)):
            self.axPolar.lines[i].set_xdata(arr)
        self.vspanList.append(self.axPolar.axvspan(self.minPick,self.maxPick,color='green',alpha=0.3))
        # calculate SNR
        if self.method == "MAN":
            self.SNR=self.calcSNR(self.sPick)
        else:
            try:
                if self.sPick:
                    self.SNR=self.calcSNR(self.sPick)
                else:
                    raise ValueError("No S pick set.")
            except:
                self.SNR=self.calcSNR(round(int((self.minPick+self.maxPick)/2)))
            finally:
                self.SNR=np.nan    
        # update the title
        self.updateTitle()
        # fix the limits, if need
        if keepAxisLimits: # Only in axN because both y n x -axis are shared
            self.axN.set_ylim(bottom,top)
            self.axN.set_xlim(left,right)
        else:
            left=0
            right=max(tN)
            self.axN.set_xlim(left,right)
        # set background color to be sure
        self.activeFig.set_facecolor('white')
        for ax in self.activeFig.get_axes():
            ax.set_facecolor('white')
        ## TODO: maybe have it as a user option at some point?
        ##       set background color according to incidence value
        #bkgClr='navajowhite' # navajowhite, ivory or blanchedalmond?
        #warnClr='lightgreen'
        #dangClr='orange'
        #bkgClr=warnClr=dangClr="white"
        #if 45 < float(self.ain) < 60:
        #    self.activeFig.set_facecolor(warnClr)
        #    for ax in self.activeFig.get_axes():
        #        ax.set_facecolor(warnClr)
        #elif 60 <= float(self.ain):
        #    self.activeFig.set_facecolor(dangClr)
        #    for ax in self.activeFig.get_axes():
        #        ax.set_facecolor(dangClr)
        #else:
        #    for ax in self.activeFig.get_axes():
        #        ax.set_facecolor(bkgClr)           
        #    self.activeFig.set_facecolor(bkgClr)
        self.canvas.draw()

    def updateTitle(self):
        """Updates the title of the figure"""
        try: # can we get the latest values from the splitting dict?
            self.updateFromSplittingDict(method=self.method)
        except KeyError: # if the selected method doesn't exist, (re)set the
            pass         # ones of the last method?
        # set the values
        self.lblOrigin.setText(str(UTCDateTime(self.evsDict[self.activeEvent]["ORIGIN"]))[:-3])
        self.lblStation.setText(self.station)
        self.inpAIN.setText("{:.1f}".format(self.ain))
        self.inpBAZ.setText("{:.1f}".format(self.baz))
        self.inpEPI.setText("{:.2f}".format(self.dist))
        self.inpMAG.setText("{:.1f}".format(self.mag))
        self.inpSNR.setText("{:.2f}".format(self.SNR))
        self.setSystem() # updates self.inpSYS
        self.inpV2H.setChecked(self.V2H)
        self.inpPHI.setText("{:.1f}".format(self.phi))
        self.inpTD.setText("{:.2f}".format(abs(self.td)*1000))
        self.inpPOL.setText("{:.1f}".format(self.pol))
        self.inpGRD.setText(self.grade)
        self.inpORN.setText("{:.1f}".format(self.orCorrection))
        self.inpFLT.setText("%.2f|%.2f" % (self.freqmin,self.freqmax))
        self.inpMET.setText(self.method)

    def updatePolarigram(self,draw=False):
        """
        Updates the polarigram
        
        :type draw: bool, optional
        :param draw: redraw the active figure's canvas. Defaults to False.

        """
        # clear the polarigram plot
        self.axPolar.cla()
        # draw dummy picks. this is done to ensure they're always in the correct order.
        self.axPolar.axvline(x=0,linewidth=2.0,color='red',label='observed')
        self.axPolar.axvline(x=0,color='darkviolet',linewidth=1.5,linestyle='-',label='autopicker') # picker line
        self.axPolar.axvline(x=0,color='darkviolet',linewidth=1.5,linestyle='--',label='theoretical') # theoretical line
        # show the grid?
        if self.actionGrid.isChecked():
            self.axPolar.grid(True)
        # get filter and rotation
        filterRange=(self.freqmin,self.freqmax)
        rotated=True if self.stageRotated else False
        # make the polarigram!
        polarPairs,polarDict=tools.polarigram(self.stream,rotated,filterRange)
        self.polarDict=polarDict
        # recolor the collection of vectors and convert it to LineCollection
        colourColl=['black' for x in polarPairs]
        polarColl=LineCollection(polarPairs,colors=colourColl)
        self.axPolar.add_collection(polarColl)
        # center the polarigram
        tools.centerColl(self.axPolar,polarColl,axPolar=True)
        # this is extremely important to NOT distort the image and, as a result, 
        # the visual result of the polarigram. Ratio is set as 1:1, i.e. the data
        # in the x axis have the same scale as the data in the y axis.
        self.axPolar.set_aspect(1,adjustable="datalim")
        ## TODO: set labels?
        #if not self.stageRotated:
        #    self.axPolar.set_ylabel("North")
        #    self.axPolar.set_xlabel("East")
        #else:
        #    self.axPolar.set_ylabel("Fast")
        #s    self.axPolar.set_xlabel("Slow")
        # temp hack, till we find a better way to display the above labels
        self.axPolar.set_ylabel("")
        self.axPolar.set_xlabel("")            
        ## get either aux or pol
        T=polarDict["TIME"] # time dimension of vectors
        D=polarDict["ANGLE"] # angle of each vector
        # get the angle at the pick
        # get the vector CLOSEST to the pick
        pAngle=D[np.abs(T-self.sPick).argmin()] % 180
        # convert to N/E-system (or F/S) angle
        pAngle=self.mpl2ns(pAngle)        
        # get polarization angle?
        if self.stageRotated and self.method is "MAN":
            self.aux=pAngle
            self.pol=(self.aux+self.phi) % 180
            self.updateSplittingDict()
        if self.stageCorrected and self.method is "MAN":
            self.pol=pAngle
            self.updateSplittingDict()
        if draw:
            self.canvas.draw()

    def updateHodogram(self,draw=False):
        """
        Updates the hodograms
        
        :type draw: bool, optional
        :param draw: redraw the active figure's canvas. Defaults to False.

        """
        n=3 # Offset of samples for hodogram, e.g. 5 means a total of 10 samples to show
        # get the vector CLOSEST to the pick
        x=np.abs(self.stream[0].times()-self.sPick).argmin()
        hodoList=[]; maxvals=[]
        for ax,j in zip((self.axPart1,self.axPart3,self.axPart2),(-1,1,0)):
            ax.cla()
            if self.actionGrid.isChecked():
                ax.grid(True)
            ax.plot(0,0,marker="*",color="red")
            # set plots around the pick fo the S
            npts=self.stream[1].stats.npts
            x1=int(x-n+2*j*n)
            if x1 < 0:
                logging.warning("x1 < 0! Using 0 as x1")
                x1=0
            elif x1 > npts:
                logging.warning("x1 > npts! Using npts as x1")
                x1=int(npts)
            x2=int(x+n+2*j*n)+1 # +1 cause of np arrays indexing
            if x2 > npts:
                logging.warning("x2 > npts! Using npts as x2")
                x2=int(npts)
            if x2 < 0:
                logging.warning("x2 < 0! Using 0 as x2")
                x2=0
            t1="{:.2f}".format(float(x1)/self.stream[1].stats.sampling_rate)
            t2="{:.2f}".format(float(x2-1)/self.stream[1].stats.sampling_rate)
            tStr=t1+" - "+t2+" s"
            ax.text(0.01,0.8,tStr,transform=ax.transAxes,fontsize=9,zorder=100)
            try:
                partPairs=tools.hodogram(self.polarDict["AMP"][x1:x2],self.polarDict["ANGLE"][x1:x2])
                hodoList.append(partPairs)
                maxvals.append(abs(partPairs[:,:,1]).max())
            except ValueError: # if cannot get a non 0-element array
                logging.exception("0-element array for "+str(ax))
                partPairs=np.asarray([0,0])
                continue
        normval=max(maxvals)
        ## plot                
        for ax,partPairs in zip((self.axPart1,self.axPart3,self.axPart2),hodoList):
            partColl=LineCollection(partPairs/normval)
            ax.add_collection(partColl)
            ax.set_aspect(1,adjustable="datalim")
            ax.set_ylim(-1.2,1.2)
            #tools.centerColl(ax,partColl)    
            if not self.stageRotated:
                ax.set_ylabel("North")
                ax.set_xlabel("East")
            else:
                ax.set_ylabel("Fast")
                ax.set_xlabel("Slow")
        if draw:
            self.canvas.draw()

    def centerAxes(self,ax):
        """
        Centers data in axes to y=0.
        
        :type ax: :class: `~matplotlib.axes.Axes`
        :param ax: axes for which to center the plot

        """
        y=ax.lines[0].get_ydata()
        ymax=max(abs(min(y)),abs(max(y)))
        ax.set_ylim(-ymax,ymax)

    def setAxisLimits(self,ax,axis,min,max,draw=True):
        """
        Set limits of given axis
        
        :type ax: :class: `~matplotlib.axes.Axes`
        :param ax: axes for which to set the limits
        :type axis: str
        :param axis: axis for which to set the limtis. Either 'x' or 'y'
        :type min: float or int
        :param min: the bottom limit to be set
        :type max: float or int
        :param max: the upper limit to be set
        :type draw: bool, optional
        :param draw: redraw the active figure's canvas. Defaults to False.

        """
        if axis == "y":
            ax.set_ylim((min,max))
        elif axis == "x":
            ax.set_xlim((min,max))
        if draw:
            self.canvas.draw()

    def mpl2ns(self,phi):
        """
        Converts an angle from the plot system to the NS system
        
        :type phi: float or int
        :param phi: the angle in the matplotlib system to convert
        :returns: the angle in the N/S (or F/S) system

        """
        if 0 <= phi < 90: # Q1
            phi=90-phi  
        elif -90 <= phi < 0: # Q2
            phi=90+abs(phi)
        elif -180 <= phi < -90: # Q3
            phi=abs(phi)-90
        elif 90 <= phi < 180: # Q4
            phi=270-phi
        return phi

    def ns2mpl(self,phi):
        """
        Converts an angle from the plot system to the NS system
        
        :type phi: float or int
        :param phi: the angle in the N/S (or F/S) system to convert
        :returns: the angle in the matplotlib system

        """
        if 0 <= phi < 90: # Q1
            phi=90-phi
        elif 90 <= phi < 180: # Q2
            phi=90-phi
        elif 180 <= phi < 270: # Q3
            phi=-phi-90
        elif 270 <= phi < 360: # Q4
            phi=270-phi
        return phi

    def mplConnect(self,fig):
        """
        Connect the existing canvas to mouse functions
        
        :type fig: :class: `~matplotlib.figure.Figure`
        :param fig: the figure for which the connections will be made

        """
        fig.canvas.mpl_connect('button_press_event',self.mplClick)
        fig.canvas.mpl_connect('scroll_event',self.mplScrollWheel)

    def mplClick(self,event):
        """
        Connect left/right click. Middle button?
        
        :type event: `matplotlib event`
        :param event: the mouse event

        """
        if event.button == 1: # left click
            self.mplLeftClick(event)
        elif event.button == 2: # middle click
            self.mplMidClick(event)
        elif event.button == 3: # right click
            self.mplRightClick(event)

    def mplLeftClick(self,event):
        """
        Function for left mouse click
           
        :type event: `matplotlib event`
        :param event: the mouse event     

        """
        if event.inaxes in (self.axZ, self.axE, self.axN, self.axPolar):
            qMods=QtWidgets.QApplication.keyboardModifiers()
            if qMods == QtCore.Qt.ControlModifier:
                # ask the user if they want to repick
                if self.td and not any((self.stageCorrected,self.stageRotated)) \
                    and ((self.method == "MAN") or ("MAN" in self.splittingDict[self.activeEvent][self.station])):

                    ok = QtWidgets.QMessageBox.question(
                                self,"Repick","You are attempting to measure phi in an already manually processed pair.\
                                               Are you sure? (Proceeding will reset any manual measurement for this pair)",
                                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
                                                    )
                    if ok == QtWidgets.QMessageBox.No:
                        return
                    else:
                        self.method="RESET"
                # reset values first if method wasn't MNA
                if self.method != "MAN":
                    self.flagReset(general=False,splitting=True,indexing=False,titling=False,splitDict=False)
                    self.grade="X"
                    self.method="MAN"
                #
                T=self.polarDict["TIME"]
                D=self.polarDict["ANGLE"]
                # Get t closest to pick
                idx=np.abs(self.stream[0].times()-event.xdata).argmin()
                self.sPick=self.stream[0].times()[idx]
                self.sArr=self.stream[0].stats.starttime+self.sPick
                self.spickDict[self.activeEvent][self.station]=self.sArr
                logging.info("S picked @ %s" % self.spickDict[self.activeEvent][self.station])
                # get SNR
                try:
                    self.SNR=self.calcSNR(self.sPick)
                except:
                    self.SNR=np.nan
                self.V2H=self.calcV2H(self.sPick)
                self.errors=(1.,self.stream[0].stats.delta/2.)             
                # get angles
                pAngle=D[np.abs(T-event.xdata).argmin()] % 180
                # Convert to Fast-Slow according to quadrant.
                pAngle=self.mpl2ns(pAngle)
                # Draw pick lines
                for ax in (self.axZ,self.axE,self.axN):
                    ax.lines[1].set_xdata(self.sPick)
                self.axPolar.lines[0].set_xdata(self.sPick)
                # set measured angles
                if not self.stageRotated and not self.stageCorrected:
                    self.phi=pAngle
                if self.stageRotated:
                    self.aux=pAngle # Auxiliary
                    self.pol=(self.phi+self.aux) % 180
                if self.stageCorrected:
                    self.pol=pAngle
                self.updateSplittingDict()
                self.updateTitle()
                self.updateHodogram()             
                self.canvas.draw()
            elif qMods == QtCore.Qt.ShiftModifier:
                self.minPick=event.xdata
                if not tools.isnone(self.maxPick):
                    if self.vspanList: # first delete the old ones
                        for vspan in self.vspanList:
                            try:
                                vspan.remove()
                            except ValueError:
                                logging.warning("ValueError while removing vspans")
                    for ax in (self.axZ, self.axE, self.axN, self.axPolar):
                        self.vspanList.append(ax.axvspan(self.minPick,self.maxPick,color='green',alpha=0.3))
                    self.canvas.draw()

    def mplMidClick(self,event):
        """
        Function for middle mouse click. 
           
        :type event: `matplotlib event`
        :param event: the mouse event     

        """
        if event.inaxes in (self.axZ, self.axE, self.axN, self.axPolar):
            qMods=QtWidgets.QApplication.keyboardModifiers()
            if qMods == QtCore.Qt.ShiftModifier:
                if self.vspanList: # first delete the old ones
                    for vspan in self.vspanList:
                        try:
                            vspan.remove()
                        except ValueError:
                            continue
                    self.minPick=np.nan; self.maxPick=np.nan
                    self.canvas.draw()

    def mplRightClick(self,event):
        """
        Function for right mouse click
           
        :type event: `matplotlib event`
        :param event: the mouse event     

        """
        if event.inaxes in (self.axZ, self.axE, self.axN, self.axPolar):
            qMods=QtWidgets.QApplication.keyboardModifiers()
            if qMods == QtCore.Qt.ShiftModifier:
                self.maxPick=event.xdata
                if not tools.isnone(self.minPick):
                    if self.vspanList: # first delete the old ones
                        for vspan in self.vspanList:
                            try:
                                vspan.remove()
                            except ValueError:
                                pass
                    for ax in (self.axZ, self.axE, self.axN, self.axPolar):
                        self.vspanList.append(ax.axvspan(self.minPick,self.maxPick,color='green',alpha=0.3))
                    self.canvas.draw()
            elif qMods == QtCore.Qt.ControlModifier:
                # get t closest to pick
                idx=np.abs(self.stream[0].times()-event.xdata).argmin()
                self.sPick=self.stream[0].times()[idx]
                self.sArr=self.stream[0].stats.starttime+self.sPick
                self.spickDict[self.activeEvent][self.station]=self.sArr
                logging.info("S picked @ %s" % self.spickDict[self.activeEvent][self.station])
                # get SNR
                try:
                    self.SNR=self.calcSNR(self.sPick)
                except:
                    self.SNR=np.nan
                self.V2H=self.calcV2H(self.sPick)
                # Draw pick lines
                for ax in (self.axZ,self.axE,self.axN):
                    ax.lines[1].set_xdata(self.sPick)
                self.axPolar.lines[0].set_xdata(self.sPick)
                self.updateSplittingDict(splitting=False)
                self.updateTitle()
                self.canvas.draw()

    def keyPressEvent(self,event):
        """
        Reroute Qt Key input for left and right arrows
           
        :type event: `Qt Event`
        :param event: the keyboard event     

        """
        if self.stageRotated or self.stageCorrected:
            step=self.stream[0].stats.delta
            qMods=QtWidgets.QApplication.keyboardModifiers()
            if qMods == QtCore.Qt.ShiftModifier:
                step*=5 # TODO: make the number of steps a user's option
            if event.key() == QtCore.Qt.Key_Left:
                self.td-=step
                if self.method != "MAN":
                    self.flagReset(general=False,basicSplit=False,indexing=False,titling=False,splitDict=False)
                    self.grade="X"
                    self.method="MAN"
                    self.grade="X"
                    self.method="MAN"
                self.updateSplittingDict()
                self.updateFig(center=False,keepAxisLimits=True)
            elif event.key() == QtCore.Qt.Key_Right:
                self.td+=step
                if self.method != "MAN":
                    self.grade="X"
                    self.method="MAN"
                self.updateSplittingDict()                
                self.updateFig(center=False,keepAxisLimits=True)

    def mplScrollWheel(self,event):
       """
       Controls the scroll wheel
           
       :type event: `matplotlib event`
       :param event: the wheel event     

       """
       ax=event.inaxes
       if ax not in (None,self.axPart1,self.axPart2,self.axPart3):
        qMods=QtWidgets.QApplication.keyboardModifiers()
        if qMods == QtCore.Qt.ShiftModifier:
            (left,right)=ax.get_xbound()
            zoomFactor= 2.4
            if event.button == 'up':
                left -= (event.xdata - left) / zoomFactor
                right += (right - event.xdata) / zoomFactor
            elif event.button == 'down':
                left += (event.xdata - left) / zoomFactor
                right -= (right - event.xdata) / zoomFactor
            ax.set_xlim(left,right)
        else:
            if ax == self.axPolar:
                return
            (bottom,top)=ax.get_ybound()
            zoomFactor= 1.5
            if event.button == 'up':
                top *= zoomFactor
                bottom *= zoomFactor
            elif event.button == 'down':
                top /= zoomFactor
                bottom /= zoomFactor
            ax.set_ylim(bottom,top)
        self.axPolar.set_aspect(1,adjustable="datalim")
        self.canvas.draw()

    def nextStation(self,autoAinSkip=True):
        """
        Go to next station
        
        :type autoAinSkip: bool, optional
        :param autoAinSkip: enables automatically skipping the station due to
            the shear-wave window. Defaults to True.

        """
        # set inital params
        foundNextStation=False
        # store the initial station index, in case of a reset
        iniStIdx=self.stIdx
        # start iterating for station indexing, from the current one onwards
        for stIdx,_ in enumerate(self.staList[self.stIdx:]):
            self.stIdx+=stIdx+1
            if self.stIdx >= len(self.staList): break       
            event=sorted(self.evsDict.keys())[self.evIdx]
            station1=self.staList[self.stIdx] 
            # shear-wave window?
            if (self.statsDict[event][station1]["AN"] > self.maxAin) and autoAinSkip:
                logging.warning("Skipping %s - %s cause of AIN criterion..." % (event,station1))
                continue
            # fetch the next event/station and if all goes well, break the loop
            try:
                self.fetchEvent(self.evIdx,self.stIdx,eventMode=False)
                logging.info("Found stream: "+str(self.stream))
                foundNextStation=True
                break
            except: # woops!
                logging.exception("Error while opening waveforms for "+event)
                winTitle=sorted(self.evsDict.keys())[self.evIdx]
                if self.stIdx+1 < len(self.staList):
                    station2=self.staList[self.stIdx+1]
                    genText="Can't find station "+station1+". Try to open "+station2+"?"
                    reply=self.replyMsgBox(genText,winTitle)
                    if reply == QtWidgets.QMessageBox.No:
                        return
                    else:
                        self.stIdx=iniStIdx
                        continue
        if not foundNextStation:
            logging.info("Reached last station!")
            genTxt=""
            infTxt="You reached the last available station for this event."
            winTit="Last Station"
            self.warnMsgBox(genTxt,infTxt,winTit)
            self.stIdx=iniStIdx
            self.updateFromSplittingDict()
            return
        self.undoStream=self.stream.copy()
        self.streamIni=self.stream.copy()
        splitRes=False
        try:
            self.updateFromSplittingDict()
            self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
        except KeyError: # means that one of event/station/method is not in the DB
            try:
                self.sArr=self.spickDict[self.activeEvent][self.station]
                self.sPick=self.sArr-self.stream[0].stats.starttime
                self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
            except:
                pass 
        self.updateFig()

    def prevStation(self,autoAinSkip=True):
        """
        Go to previous station
        
        :type autoAinSkip: bool, optional
        :param autoAinSkip: enables automatically skipping the station due to
            the shear-wave window. Defaults to True.

        """
        foundPrevStation=False
        iniStIdx=self.stIdx
        for stIdx in range(iniStIdx-1,-1,-1):
            self.stIdx=stIdx
            if self.stIdx < 0: break       
            event=sorted(self.evsDict.keys())[self.evIdx]
            station1=self.staList[self.stIdx] 
            if (self.statsDict[event][station1]["AN"] > self.maxAin) and autoAinSkip:
                logging.warning("Skippnig %s - %s cause of AIN criterion..." % (event,station1))
                continue
            try:
                self.fetchEvent(self.evIdx,self.stIdx,eventMode=False)
                logging.info("Found stream: "+str(self.stream))
                foundPrevStation=True
                break
            except:
                logging.exception("Error while opening waveforms for "+event)
                winTitle=sorted(self.evsDict.keys())[self.evIdx]
                if self.stIdx-1 > 0:
                    station2=self.staList[self.stIdx-1]
                    genText="Can't find station "+station1+". Try to open "+station2+"?"
                    reply=self.replyMsgBox(genText,winTitle)
                    if reply == QtWidgets.QMessageBox.No:
                        return
                    else:
                        self.stIdx=iniStIdx
                        continue
        if not foundPrevStation:
            logging.info("Reached first station!")
            genTxt=""
            infTxt="You reached the first available station for this event."
            winTit="Last Station"
            self.warnMsgBox(genTxt,infTxt,winTit)
            self.stIdx=iniStIdx
            return
        self.undoStream=self.stream.copy()
        self.streamIni=self.stream.copy()
        splitRes=False
        try:
            self.updateFromSplittingDict()
            self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
        except KeyError: # means that one of event/station/method is not in the DB
            try:
                self.sArr=self.spickDict[self.activeEvent][self.station]
                if self.sArr not in [self.stream[0].stats.starttime,0,np.nan]: 
                    self.sPick=self.sArr-self.stream[0].stats.starttime
                #splitRes=True
                self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
            except KeyError:
                pass         
        self.updateFig()

    def nextEvent(self):
        """Go to next event"""
        inEvIdx=self.evIdx
        inStIdx=self.stIdx
        foundNextEvent=False
        while not foundNextEvent:
            self.evIdx+=1
            if self.evIdx >= len(self.evsDict):
                logging.info("This is the last event!")
                self.evIdx-=1
                return            
            event1=sorted(self.evsDict.keys())[self.evIdx]
            ## skip cause of ain?
            self.staList=sorted(self.statsDict[event1],key=lambda x:self.statsDict[event1][x]["AN"])
            station=self.staList[0]
            if self.statsDict[event1][station]["AN"] > self.maxAin:
                logging.warning("Skipping %s - %s cause of AIN criterion..." % (event1,station))
                continue
            if self.evIdx == len(self.evsDict)-1:
                event2=event1
            else:
                event2=sorted(self.evsDict.keys())[self.evIdx+1]
            try:
                self.stIdx=0
                self.fetchEvent(self.evIdx,self.stIdx)
                logging.info("Found stream: "+str(self.stream))
                foundNextEvent=True
            except StationException:
                self.nextStation()
                foundNextEvent=True
                return
            except:
                if event1 != event2:
                    logging.exception("Error while opening waveforms for "+event1)
                    winTitle=sorted(self.evsDict.keys())[self.evIdx]
                    genText=str("Cannot open event "+event1+". Try to open "+event2+"?")
                    reply=self.replyMsgBox(genText,winTitle)
                    if reply == QtWidgets.QMessageBox.No:
                        self.stIdx=inStIdx
                        self.evIdx=inEvIdx
                        return
                else:
                    logging.exception("Error while opening waveforms for "+event1)
                    logging.debug("this was the last event!")
        self.undoStream=self.stream.copy()
        self.streamIni=self.stream.copy()
        self.station=self.stream[0].stats.station
        splitRes=False
        try:
            self.updateFromSplittingDict()
            self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
        except KeyError: # means that one of event/station/method is not in the DB
            try:
                #splitRes=True
                self.sArr=self.spickDict[self.activeEvent][self.station]
                if self.sArr not in [self.stream[0].stats.starttime,0,np.nan]: 
                    self.sPick=self.sArr-self.stream[0].stats.starttime
                self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
            except KeyError:
                pass 
        self.updateFig()

    def prevEvent(self):
        """Go to previous event"""
        inEvIdx=self.evIdx
        inStIdx=self.stIdx
        foundPrevEvent=False
        while not foundPrevEvent:
            self.evIdx-=1
            if self.evIdx < 0:
                logging.info("This is the first event!")
                self.evIdx+=1
                return            
            event1=sorted(self.evsDict.keys())[self.evIdx]
            ## skip cause of ain?
            self.staList=sorted(self.statsDict[event1],key=lambda x:self.statsDict[event1][x]["AN"])
            station=self.staList[0]
            if self.statsDict[event1][station]["AN"] > self.maxAin:
                logging.warning("Skipping %s - %s cause of AIN criterion..." % (event1,station))
                continue            
            if self.evIdx == 0:
                event1=event2
            else:
                event2=sorted(self.evsDict.keys())[self.evIdx-1]
            try:
                self.stIdx=0
                self.fetchEvent(self.evIdx,self.stIdx)
                logging.info("Found stream: "+str(self.stream))
                foundPrevEvent=True
            except StationException:
                self.prevStation()
                foundPrevEvent=True
                return
            except:
                if event1 != event2:
                    logging.exception("Error while opening waveforms for "+event1)
                    winTitle=sorted(self.evsDict.keys())[self.evIdx]
                    genText=str("Cannot open event "+event1+". Try to open "+event2+"?")
                    reply=self.replyMsgBox(genText,winTitle)
                    if reply == QtWidgets.QMessageBox.No:
                        self.stIdx=inStIdx
                        self.evIdx=inEvIdx
                        return
                else:
                    logging.exception("Error while opening waveforms for "+event1)
                    logging.debug("this was the first event!")
        self.undoStream=self.stream.copy()
        self.streamIni=self.stream.copy()
        self.station=self.stream[0].stats.station
        splitRes=False
        try:
            self.updateFromSplittingDict()
            self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
        except KeyError: # means that one of event/station/method is not in the DB
            try:
                #splitRes=True
                self.sArr=self.spickDict[self.activeEvent][self.station]
                if self.sArr not in [self.stream[0].stats.starttime,0,np.nan]: 
                    self.sPick=self.sArr-self.stream[0].stats.starttime
                self.flagReset(general=False,indexing=False,titling=False,splitDict=False,splitting=splitRes)
            except KeyError:
                pass  
        self.updateFig()

    ## GUI FUNCTIONS ##
    def getWave(self,fileName,inventory=None,prefChannel="HH,EH,BH,HN",orientCorr=True,showWarn=True):
        """
        Ask for required waveforms. Order is Vertical -> North -> East
        
        :type fileName: str
        :param fileName: query string for waveforms
        :type inventory: :class: `~obspy.core.inventory` or None, optional
        :param inventory: inventory with stored station information or None.
            Defaults to None.
        :type prefChannel: str, optional
        :param prefChannel: string with the order of the channel codes (per SEED) that the
            program will preferentialy look for, comma-separated. Defaults to 'HH,EH,BH,HN'.
        :type orientCorr: bool, optional
        :param orientCorr: perform orientation correction if required. Defaults to True.
        :type showWarn: bool, optional
        :param showWarn: display a warning dialog if something goes wrong. Defaults to True.
        :returns: (obspy Stream object, instrument orientation correction)

        """
        # convert preferred channels object to list
        prefChannel=prefChannel.split(",")
        # clean up file name just in case
        fileName=str(fileName).replace("/",os.sep)
        logging.debug("Trying to open "+str(fileName)+"...")
        # start looking for the station/event
        try:
            if not glob(fileName):
                fileNamePath=os.path.dirname(fileName)
                if os.path.exists(fileNamePath):
                    raise StationException("Could not find station "+fileName)
                else:
                    raise EventException("Could not find event "+fileNamePath)
            stream=read(str(fileName))
            logging.info("Read "+str(stream))
            # look for the channels as provided by the user. if '??' is reached,
            # the first ones will be used
            for prefCode in prefChannel+["??"]:
                if prefCode == "??":
                    prefCode=stream[0].stats.channel[:2]
                numChan=[x.stats.channel for x in stream.select(channel=prefCode+"?")]
                # did we find 3 channels?
                if len(numChan) == 3:
                    stream=stream.select(channel=prefCode+"?")
                    logging.debug("Found channels %s?" % stream[0].stats.channel[:2])
                    break
            # if no 3 channels could be found, raise an exception
            if len(numChan) != 3:
                raise StationException("Station does not have 3 channels!")
            # correct calibration factor if null
            for tr in stream:
                if float(tr.stats.calib) == 0.0:
                    tr.stats.calib=1.0
            # get waveforms and correct orientation from inventory
            comps=sorted(list(set([tr.stats.channel[-1] for tr in stream])))
            if "3" in comps:
                inpComponents="Z23"
                if "1" in comps:
                    inpComponents="123"
            elif "2" in comps:
                inpComponents="Z12"
            elif "N" in comps:
                inpComponents="ZNE"
            streamBk=stream.copy() # save a stream backup, apparently Obspy (1.1.1) deletes the
                                   # data if an exception is raised during rotation. see 
                                   # obspy/obspy #2372
            try:
                if orientCorr:
                    if inventory: # is the station in the inventory?
                        stream.rotate("->ZNE",inventory=inventory,components=(inpComponents))     
                        # the azimuth of the 'north' (second as defined above) is the
                        # orientation correction
                        orCorrection=inventory.select(
                            network=stream[0].stats.network,
                            station=stream[0].stats.station,
                            location=stream[0].stats.location,
                            channel=stream[0].stats.channel[:2]+inpComponents[1],
                            time=stream[0].stats.starttime
                                                      )[0][0][0].azimuth
                    else:
                        raise IOError("No available inventory")
                else:
                    raise ValueError("Instrument orientation correction disabled.")
            except:
                stream=streamBk.copy()
                logging.exception("Could not perform orientation correction. Renaming channels...")
                stream.select(component=inpComponents[0])[0].stats.channel=\
                                stream.select(component=inpComponents[0])[0].stats.channel[:2]+"Z"
                stream.select(component=inpComponents[1])[0].stats.channel=\
                                stream.select(component=inpComponents[1])[0].stats.channel[:2]+"N"
                stream.select(component=inpComponents[2])[0].stats.channel=\
                                stream.select(component=inpComponents[2])[0].stats.channel[:2]+"E"
                orCorrection=0.
            # get final combo of channels
            comps=sorted(list(set([tr.stats.channel[-1] for tr in stream])))
            # preprocess. we need to remove mean/trend to be able to merge.
            logging.info("Applying detrend/demean...")
            stream.detrend("linear")
            stream.detrend("demean")
            stream.merge(method=1,fill_value=0,interpolation_samples=0)
        except Exception as exc:
            logging.exception("Could not read stream")
            if showWarn:
                genTxt="Could not read "+str(fileName)
                infTxt=str(exc)
                winTit="Stream Error"
                self.warnMsgBox(genTxt,infTxt,winTit)
            return
        return stream, orCorrection

             
    def getCat(self,dataPath,catFile,dbPath):
        """
        Ask for data path and catalogue file
        
        :type dataPath: str
        :param dataPath: path to the data directory
        :type catFile: str
        :param catFile: path to the catalogue file
        :type dbPath: str
        :param dbPath: path to the database file

        """
        # initial params
        indexed=False
        idxPath=WORKDIR+os.sep+"etc%sindex%s"%(os.sep,os.sep)
        idxFile=WORKDIR+os.sep+"etc%sindex%sindex.cat"%(os.sep,os.sep)
        # search for an existing index to the folder
        try:
            self.pDial.setLabelText("Getting event directories...")
            self.pDial.setValue(30)  
            QtCore.QCoreApplication.processEvents()
        except:
            pass
        # start searching in the folders for available events and archive them
        # to an index file for faster access the next time.
        logging.debug("Initiating folder searching...")
        try:
            with open(idxFile,"r") as fid:
                indexCat=fid.readlines()
                indexCat=[x.strip() for x in indexCat]
            for idx in indexCat[1:]: # First line is comment
                mFile,mPath=idx.split(",")
                if mPath == dataPath:
                    with open(idxPath+mFile,"r") as fid:
                        temp=fid.readlines()
                        temp=[x.strip() for x in temp]
                    self.fullEvents={
                                        x.strip().split(os.sep)[-1]:x.strip() for x in temp

                                    }
                    logging.debug("Read path information from existing index file "+mFile)
                    indexed=True      
        except IOError:
            with open(WORKDIR+os.sep+"etc/index/index.cat","w") as fid:
                fid.write("# Index file, main path\n")
            logging.debug("Created index catalogue file.")
        if not indexed:
            self.fullEvents=self.getTree(dataPath)
            logging.info("Found "+str(len(self.fullEvents))+" events in given path.")
            idxNew=UTCDateTime.now().strftime("%Y%m%d%H%M%S")+".idx"
            with open(idxPath+idxNew,"w") as fid:
                fid.writelines([self.fullEvents[x]+"\n" for x in self.fullEvents])
            logging.debug("Created new index file: "+idxPath+idxNew)
            with open(idxFile,"a") as fid:
                fid.write(idxNew+","+dataPath+"\n")
        logging.info("Found "+str(len(self.fullEvents))+" events in given path.")
        ## read the catalogue file
        try:
            self.pDial.setLabelText("Reading the catalogue file...")
            self.pDial.setValue(40)  
            QtCore.QCoreApplication.processEvents() # maybe we should replace this.
        except:                                     # can it lead to unstable behaviour?
            pass
        try:
            ## TODO: implement a better method for validating the catalogue is in QuakeML
            ##       format
            
            # read the first line and check if it contains the xml tag 
            with open(catFile,'r') as fid: xCheck=fid.readline()
            if not '<?xml' in xCheck: raise TypeError('Given catalogue is not in XML format...')
            self.evsDict,self.statsDict=self.readQML(catFile)
        except TypeError: # then it's not a QuakeML!
            try:
                self.pDial.setLabelText("Calculating angles of incidence...")
                self.pDial.setValue(50)  
                QtCore.QCoreApplication.processEvents()
            except:
                pass
            logging.info("QUAKEML file not identified! Trying to read file as catalogue...")
            self.evsDict,self.statsDict=self.readEventCat(catFile)
        # if no events/catlogues were found, let the user know something went wrong
        if len(self.evsDict) == 0:
            logging.warning("Could not find events in %s" % catFile)
            try:
                self.pDial.close()
            except:
                pass
            winTitle="No Catalogue Found"
            genText="Could not find the catalogue!"
            infText="No corresponding catalogue found in %s"%catFile
            self.warnMsgBox(genText,infText,winTitle)
            return
        logging.info("Found "+str(len(self.evsDict))+" events in given catalogue.")
        ## intialize the velocity model
        # set model object
        try:
            if self.tpCNF.model.endswith("tvel") or self.tpCNF.model.endswith("nd"):
                from obspy.taup.taup_create import build_taup_model
                modroot=os.path.split(self.tpCNF.model)[0]
                _=build_taup_model(self.tpCNF.model,output_folder=modroot)
                # get vel model name
                vmod=os.path.split(os.path.splitext(self.tpCNF.model)[0])[1]
                vmod=os.path.splitext(self.tpCNF.model)[0]+'.npz'
            # set taup model
            self.tmodel=TauPyModel(model=vmod)
        except:
            logging.exception("Could not use custom model! Using IASP91 instead...")
            vmod="IASP91"
            # set taup model
            self.tmodel=TauPyModel(model=vmod)            
        # set velocity model
        try:
            self.vmodel=VelocityModel.read_velocity_file(self.tpCNF.model)
        except:
            # if the vel model doesn't exist it might mean this
            # is one of the builtins, so search for it
            bltin=taup_get_builtin_model_files()
            vmod=[x for x in bltin if vmod.lower() in x.lower()]
            try:
                vmod=vmod[0]
            except IndexError: # means no results were returned
                vmod='IASP91'
                vmod=[x for x in bltin if vmod.lower() in os.path.split(x)[1].lower()][0]
            self.vmodel=VelocityModel.read_velocity_file(vmod)
        logging.debug('Using velocity model %s' % vmod)
        ## read the database if it exists, otherwise create it
        try: # create
            self.dbConn,self.dbCur=DB.init(self.dbPath)
            self.splittingDict={}
            self.spickDict={}
        except IOError: # read
            self.dbConn,self.dbCur=DB.open(self.dbPath)
        # correct keys
        try:
            self.pDial.setLabelText("Matchin directories to events...")
            self.pDial.setValue(60)  
            QtCore.QCoreApplication.processEvents()
        except:
            pass
        self.evsDict,self.statsDict=self.eventCorrect(
                                                sorted(self.fullEvents.keys()),
                                                self.evsDict,self.statsDict
                                                        )

    def openCat(self):
        """Function for Open Catalogue action."""
        if not self.loadCatalogue:
            return
        # get the files
        self.pDial=self.progressDialog("Opening catalogue")
        self.pDial.setFixedSize(400,100)
        logging.debug("Selected catalogue file: "+self.catFile)
        self.pDial.setLabelText("Loading...")
        try: # get the catalogue, files etc
            self.getCat(self.dataPath,self.catFile,self.dbPath)
        except FileNotFoundError: # woops!
            logging.warning("Could not find file %s" % self.catFile)
            self.pDial.close()
            winTitle="No Catalogue Found"
            genText="Could not find the catalogue!"
            infText="No corresponding catalogue found in %s" % self.catFile
            self.warnMsgBox(genText,infText,winTitle)
            return
        self.pDial.setLabelText("Fetch S arrivals...")
        self.pDial.setValue(50)  
        QtCore.QCoreApplication.processEvents()
        # get s picks from DB/QuakeML
        self.fetchSPicks()
        #
        self.pDial.setLabelText("Finding first available event-station pair...")
        self.pDial.setValue(70)
        # assign self parameters
        self.evIdx=0 # start from zero index and keep going
        logging.info("Now reading waveforms...")
        for n in range(len(self.evsDict)):
            try: # keep going til we find the first available pair
                logging.debug("Fetching event #%i" % n)
                self.fetchEvent(n)
                self.evIdx=n
                logging.info("Found stream: "+str(self.stream))
                break
            except:
                logging.exception("Error while opening waveforms for "+sorted(self.evsDict.keys())[n])
        try:
            if len(self.stream) == 0:
                raise AttributeError('Stream is empty')
            self.streamIni=self.stream.copy()
        except AttributeError: # means self.stream wasn't created, i.e. could find no pairs
            logging.exception("self.stream was not created")
            self.pDial.close()
            winTitle="No Waveforms Found"
            genText="Could not find any waveforms!"
            infText="No corresponding waveforms found in %s" % self.dataPath
            self.warnMsgBox(genText,infText,winTitle)
            return
        self.updateFig()
        # enable menus/actions
        for menu in self.activateMenus: 
            menu.setEnabled(True)
            for action in menu.actions():
                action.setEnabled(True)
        self.pDial.close()

    def exportFigure(self):
        """ Export figure(s) for the selected available method."""
        ok=True
        _,vals=DB.retrieve(self.dbCur,"method","method","%s/%s/*" % (self.activeEvent,self.station))
        accMethods=sorted(list(set([x[0] for x in vals])))
        if not accMethods:
            genTxt="Could not find methods in database"
            infTxt="There are no available solutions in the database"
            winTit="Solution Not Found"
            self.warnMsgBox(genTxt,infTxt,winTit)
            return
        try:
            idx=accMethods.index(self.method)
        except:
            idx=0
        method, ok = QtWidgets.QInputDialog.getItem(
            self,"Method to load","Select method: ",accMethods,idx,False
                                                  )
        if not ok:
            return
        # selected output file
        root=os.path.split(self.dbPath)[0]
        figFile="%s%s%s_%s_%s_%s" % (root,os.sep,
                                         self.activeEvent,
                                         self.stream[0].stats.network,
                                         self.stream[0].stats.station,
                                         method)
        figFile,_=QtWidgets.QFileDialog.getSaveFileName(self,"Select the database file:",figFile,
            "Portable Network Graphics (*.png);; JPEG compressed (*.jpg);; PostScript (*.ps);; Encapsulated PostScript (*.eps);; \
             Scalable Vector Graphics (*.svg);; Portable Document Format (*.pdf);;All files (*)")
        if len(figFile) == 0:
            return
        # get dict with solution
        metDict=self.splittingDict[self.activeEvent][self.station][method]
        # make the figure(s)
        titleInfo=str(self.origtime),self.ain,self.dist,self.mag,metDict["grade"]
        if method== "MAN":
            tools.qcFigMAN(self.streamIni.copy(),
                           metDict["phi"],metDict["td"]/(10**3),metDict["pol"],
                           self.splittingDict[self.activeEvent][self.station]["s_obs"],
                           self.baz,
                           self.ain,
                           (metDict["filter_min"],metDict["filter_max"]),
                           titleInfo,output=figFile,extend=0.2)
            plt.show()
        elif method in ["MEV","MME","MRC"]:
            tools.qcFigSWS(
                           self.streamIni.copy(),metDict["phi"],metDict["td"]/(10**3),metDict["pol"],
                           (metDict["err_phi"],metDict["err_td"]/(10**3)),
                           self.baz,
                           (
                            UTCDateTime(metDict["window_min"])-self.stream[0].stats.starttime,
                            UTCDateTime(metDict["window_max"])-self.stream[0].stats.starttime
                            ),
                           metDict["C_array"],metDict["phi_test"],
                           metDict["td_test"],method[1:],
                           metDict["C_max"],(metDict["filter_min"],metDict["filter_max"]),titleInfo,
                           output=figFile
                          )  
            plt.show()      
        elif method in ["AEV","AME","ARC"]:
            tools.qcFigSWS(
                           self.streamIni.copy(),metDict["phi"],metDict["td"]/(10**3),metDict["pol"],
                           (metDict["err_phi"],metDict["err_td"]/(10**3)),
                           self.baz,
                           (
                            UTCDateTime(metDict["window_min"])-self.stream[0].stats.starttime,
                            UTCDateTime(metDict["window_max"])-self.stream[0].stats.starttime
                            ),
                           metDict["C_array"],metDict["phi_test"],
                           metDict["td_test"],method[1:],
                           metDict["C_max"],(metDict["filter_min"],metDict["filter_max"]),titleInfo,
                           output=figFile
                          ) 
            root=os.path.split(figFile)[0]
            figFile="%s%s%s_%s_%s_%s.%s" % (root,os.sep,
                                 self.activeEvent,
                                 self.stream[0].stats.network,
                                 self.stream[0].stats.station,
                                 "CA",os.path.splitext(figFile)[1])
            figFile,_=QtWidgets.QFileDialog.getSaveFileName(self,"Select the database file:",figFile,
                "Portable Network Graphics (*.png);; JPEG compressed (*.jpg);; PostScript (*.ps);; Encapsulated PostScript (*.eps);; \
                 Scalable Vector Graphics (*.svg);; Portable Document Format (*.pdf);;All files (*)")
            if len(figFile) == 0:
                return
            tools.qcClusterDiagram(
             metDict["initial_clusters"],
             metDict["initial_clusters_errors"],
             metDict["calinski_score"],
             metDict["clusters1"],metDict["clusters2"],
             self.stream[0].stats.station,
             self.stream[0].stats.network,self.baz,
             titleInfo, method[1:],
             (metDict["filter_min"],metDict["filter_max"]),
             output=figFile)
            plt.show()

    def exportActive(self):
        """Export spl file for active event-station pair."""
        splittingDict={self.activeEvent:{
                        self.station:{
                           self.method:self.splittingDict[self.activeEvent][self.station][self.method]
                                     }
                                      }
                                        }
        for key in sorted(self.splittingDict[self.activeEvent]):
            if not isinstance(self.splittingDict[self.activeEvent][key],dict):
                splittingDict[self.activeEvent][key]=self.splittingDict[self.activeEvent][key]
        for key in sorted(self.splittingDict[self.activeEvent][self.station]):
            if not isinstance(self.splittingDict[self.activeEvent][self.station][key],dict):
                splittingDict[self.activeEvent][self.station][key]=self.splittingDict[self.activeEvent][self.station][key]
        self.writeSplittingDict(splittingDict)

    def saveAsDB(self):
        """Save existing database to different path."""
        # ask for the db location
        dbPath,_=QtWidgets.QFileDialog.getSaveFileName(self,"Select the database file:",self.dbPath,
            "Database v3 (*.db3);;Database (*.db);;SQLITE3 (*.sqlite3) ;;All files (*)")
        if len(dbPath) == 0:
            return
        # close connection to database and copy file
        DB.close(self.dbConn)
        shutil.copy2(self.dbPath,dbPath)
        # connect to new database
        self.dbPath=dbPath
        self.dbConn,self.dbCur=DB.open(self.dbPath)
        logging.info("Saved database as %s" % self.dbPath)

    def loadNewDB(self):
        """Load a new database."""
        try:
            # ask for the db location
            dbPath,_=QtWidgets.QFileDialog.getOpenFileName(self,"Select the database file:",self.dbPath,
                "Database v3 (*.db3);;Database (*.db);;SQLITE3 (*.sqlite3) ;;All files (*)")
            if len(dbPath) == 0:
                return
            self.dbPath=dbPath
            logging.info("Changes will take place in %s" % dbPath)
            #
            DB.close(self.dbConn)
            self.dbConn,self.dbCur=DB.open(self.dbPath)
            self.splittingDict=DB.load(self.dbCur,[self.activeEvent,])
            self.updateFromSplittingDict()
            self.updateFig()
            self.canvas.draw()
        except:
            tbStack=traceback.format_exc()
            genTxt="Loading New Database"
            infTxt="Loading new database exited with exception:\n %s" % tbStack
            winTit="Could not load new database!"
            self.warnMsgBox(genTxt,infTxt,winTit) 

    def fetchEvent(self,evIdx,stIdx=0,eventMode=True):
        """
        Fetch event stream from imported catalogue
        
        :type evIdx: int
        :param evIdx: the event index
        :type stIdx: int, optional
        :param stIdx: the station index. Defaults to 0.
        :type eventMode: bool, optional
        :param eventMode: Enable opening stations AND events. Defaults to True.

        """
        self.flagReset(splitDict=False,indexing=False)
        self.activeEvent=sorted(self.evsDict)[evIdx]
        logging.info("~"*25)
        logging.info("fetchEvent: "+str(self.activeEvent))
        self.staDict=self.statsDict[self.activeEvent]
        self.staList=sorted(self.staDict,key=lambda x: self.staDict[x]["AN"])
        self.evDict=self.evsDict[self.activeEvent]
        # get required waveforms
        station=self.staList[stIdx]
        evFolder=self.fullEvents[self.activeEvent]
        # iterate through ZNE
        while True:
            self.stream=Stream()
            tempStream=Stream()
            stationFile=evFolder+os.sep+"*."+station+".*"
            try:
                tempStream,orCorrection=self.getWave(stationFile,
                                                     inventory=self.xmlInventory,
                                                     prefChannel=self.generalCNF.chanPref,
                                                     orientCorr=self.generalCNF.orientFlag,
                                                     showWarn=False)
            except:
                logging.exception('Exception while trying to read %s' % stationFile)
                pass # pass cause stream will be empty and handled below
            if not tempStream:
                if not eventMode:
                    raise StationException("Could not find %s" %station)
                stIdx+=1
                if stIdx >= len(self.staList):
                    raise StationException("Could not find any stations for %s" % self.activeEvent)
                station=self.staList[stIdx]
            else:
                break
        self.stIdx=stIdx
        self.stream=tempStream
        self.station=self.stream[0].stats.station
        self.orCorrection=orCorrection
        logging.info("Orientation correction: %.0f" % self.orCorrection)
        # update splitting dict
        self.splittingDict=DB.load(self.dbCur,[self.activeEvent])
        if self.splittingDict:
            logging.debug("Found splitting entries for %s in the DB!" % self.station)
        else:
            logging.warning("No splitting entries for %s in the DB!" % self.station)
        # check if all timeseries have the same length
        Wz=self.stream.select(component="Z")[0]
        Wy=self.stream.select(component="N")[0]
        Wx=self.stream.select(component="E")[0]
        # First horizontals
        Ax,Ay=tools.lengthcheck(Wx,Wy)
        self.stream.select(component="E")[0].data=Ax 
        self.stream.select(component="N")[0].data=Ay
        # Now the vertical
        Az,Ay=tools.lengthcheck(Wz,Wy)
        self.stream.select(component="Z")[0].data=Az
        self.stream.select(component="N")[0].data=Ay
        logging.debug("Corrected stream before trim %s" % str(self.stream))
        # Horizontals second time 
        Ax,Ay=tools.lengthcheck(Wx,Wy)
        self.stream.select(component="E")[0].data=Ax 
        self.stream.select(component="N")[0].data=Ay 
        ## sync starttimes if less than a sample difference
        delta=self.stream[0].stats.delta
        if 0 < abs(Wx.stats.starttime-Wy.stats.starttime) <= delta:
            Wx.stats.starttime=max((Wy.stats.starttime,Wx.stats.starttime))
            Wy.stats.starttime=max((Wy.stats.starttime,Wx.stats.starttime))
        if 0 < abs(Wz.stats.starttime-Wy.stats.starttime) <= delta:
            Wz.stats.starttime=max((Wy.stats.starttime,Wz.stats.starttime))
            Wy.stats.starttime=max((Wy.stats.starttime,Wz.stats.starttime))            
        if 0 < abs(Wx.stats.starttime-Wy.stats.starttime) <= delta:
            Wx.stats.starttime=max((Wy.stats.starttime,Wx.stats.starttime))
            Wy.stats.starttime=max((Wy.stats.starttime,Wx.stats.starttime))            
        # initiliaze picks
        self.sPick=0; self.sArr=0
        self.sPickTheor=0; self.sArrTheor=0
        self.sPickAuto=0; self.sArrAuto=0
        # get theoretical information from TauP
        try:
          taupdat=tools.getTheorArrivals( # TODO: display secondary phases too
                                   (self.evsDict[self.activeEvent]["ORIGIN"],self.evDict["LAT"],
                                                                  self.evDict["LON"],
                                                                  self.evDict["DEPTH"]),
                                   (self.inventory[station]["latitude"],
                                    self.inventory[station]["longitude"]),
                                   self.tmodel
                                        )
          # grab measurements for the s arrival (not S)
          foundS=False
          for reg in taupdat.values():
            if reg["phase"] == "s": 
                foundS=True
                self.sArrTheor=reg["time"]
                self.sPickTheor=self.sArrTheor-self.stream[0].stats.starttime
                if self.tpCNF.ainFlag:
                    if tools.isnone(reg['ain']):
                        raise ValueError("TauP couldn't estimate a valid incidence angle.")
                    logging.debug('Estimated TauP incidence angle at %.1f' % reg['ain'])
                    self.ain=reg["ain"]
                    self.statsDict[self.activeEvent][station]["AN"]=self.ain
                else:
                    self.ain=self.statsDict[self.activeEvent][station]["AN"]
                break # can there be more than 1 "s" arrivals??
          if not foundS:
            raise ValueError('Could not find a valid s arrival with TauP!')
        except:
          logging.exception("Could not calculate theoretical arrivals for %s in %s pair!" % (station,self.activeEvent))
          self.sArrTheor=0; self.sPickTheor=0
          self.ain=self.statsDict[self.activeEvent][station]["AN"]
        logging.info("Theoretical arrival @ %s" % self.sArrTheor)
        logging.info("Theoretical pick @ %.3f (before trim)" % self.sPickTheor)
        # update other values
        self.baz=self.statsDict[self.activeEvent][station]["BAZ"]
        self.dist=self.statsDict[self.activeEvent][station]["DIST"]
        self.mag=self.evsDict[self.activeEvent]["MAG"]
        self.origtime=self.evsDict[self.activeEvent]["ORIGIN"]        
        # load values from DB
        try:
            self.updateFromSplittingDict()
            if self.sArr in [self.stream[0].stats.starttime,0,np.nan]:
                raise KeyError("Could not find arrival in splitting dictionary")
        except KeyError: # means that one of event/station/method is not in the DB
            try:
                self.sArr=self.spickDict[self.activeEvent][self.station]
                if self.sArr not in [self.stream[0].stats.starttime,0,np.nan]: 
                    self.sPick=self.sArr-self.stream[0].stats.starttime
            except KeyError:
                pass  
        logging.info("Observed arrival @ %s" % self.sArr)
        logging.info("Observed pick @ %.3f (before trim)" % self.sPick)            
        ## use AR-AIC?
        if self.toggleAR.isChecked():
            try:
                logging.debug("Applying the AR-AIC autopicker...")
                _,self.sPickAuto=tools.arPicker(self.stream,self.pkCNF)
                self.sArrAuto=self.stream[0].stats.starttime+self.sPickAuto
            except:
                logging.exception("Error while getting automatic pick.")
                self.sPickAuto=0; self.sArrAuto=0
        else:
            logging.debug("automatic picker is deactivated")
            self.sPickAuto=0; self.sArrAuto=0
        logging.info("Automatic arrival @ %s" % self.sArrAuto)
        logging.info("Automatic pick @ %.3f (before trim)" % self.sPickAuto)   
        # trim waveforms for faster processing
        if self.generalCNF.trimFlag and (self.sArr not in [self.stream[0].stats.starttime,0,np.nan]):
            self.stream.trim(
                             starttime=self.sArr+self.generalCNF.trimStart,
                             endtime=self.sArr+self.generalCNF.trimEnd,
                             pad=True,nearest_sample=False,fill_value=0
                             )  
            # check if all timeseries have the same length
            Wz=self.stream.select(component="Z")[0]
            Wy=self.stream.select(component="N")[0]
            Wx=self.stream.select(component="E")[0]
            # First horizontals
            Ax,Ay=tools.lengthcheck(Wx,Wy)
            self.stream.select(component="E")[0].data=Ax 
            self.stream.select(component="N")[0].data=Ay
            # Now the vertical
            Az,Ay=tools.lengthcheck(Wz,Wy)
            self.stream.select(component="Z")[0].data=Az
            self.stream.select(component="N")[0].data=Ay
            # Horizontals second time 
            Ax,Ay=tools.lengthcheck(Wx,Wy)
            self.stream.select(component="E")[0].data=Ax 
            self.stream.select(component="N")[0].data=Ay       
            logging.debug("Corrected stream after trim %s" % str(self.stream))
            # recalc picks
            if self.sArr not in [self.stream[0].stats.starttime,0,np.nan]: 
                self.sPick=self.sArr-self.stream[0].stats.starttime
            if self.sArrTheor not in [self.stream[0].stats.starttime,0,np.nan]:
                self.sPickTheor=self.sArrTheor-self.stream[0].stats.starttime
            if self.sArrAuto not in [self.stream[0].stats.starttime,0,np.nan]:
                self.sPickAuto=self.sArrAuto-self.stream[0].stats.starttime
            logging.info("Observed pick @ %.3f (after trim)" % self.sPick) 
            logging.info("Theoretical pick @ %.3f (after trim)" % self.sPickTheor)
            logging.info("Automatic pick @ %.3f (after trim)" % self.sPickAuto) 
        # calc V/H 
        self.V2H=self.calcV2H(self.sPick)
        # calc SNR
        if self.sPick not in [0,np.nan]:
            self.SNR=self.calcSNR(self.sPick)
        else:
            self.SNR=np.nan


    def changeAR(self):
        """Turn AR-AIC picker on and off."""
        try:
            if self.toggleAR.isChecked():
                logging.debug("Applying the AR-AIC autopicker...")
                _,s=tools.arPicker(self.stream,self.pkCNF)

            else:
                logging.debug("automatic picker is deactivated")
                s=0
            self.sPickAuto=s
            if self.sPickAuto:
                self.sArrAuto=self.stream[0].stats.starttime+self.sPickAuto
            for ax in (self.axZ,self.axN,self.axE):
                # plot s
                ax.lines[2].set_xdata(s)
            self.axPolar.lines[1].set_xdata(s)  
            self.updateSplittingDict(splitting=False)              
            self.canvas.draw()
        except AttributeError:
            return

    def fetchSPicks(self):
        """
        Store all individual S picks per event per station
        in a comprehensive dictionary.

        """
        # make the s pick dictionary and store the arrival times
        self.spickDict={}
        for ev in sorted(self.evsDict):
            self.spickDict.update({ev:{}})
            for sta in sorted(self.statsDict[ev]):
                try: # DB picks takes precedent
                    _,vals=DB.retrieve(self.dbCur,"station","s_obs","%s/%s" % (ev,sta))
                    if vals[0] in [0,np.nan,None]:
                        raise IndexError("No arrival specified through Pytheas")
                    if len(vals[0]) > 1:
                        raise IndexError("More than 1 arrival found specified in database!")
                    self.spickDict[ev].update({sta:UTCDateTime(vals[0][0])})
                except: # fetch from DB
                    try:
                        pickval=self.evsDict[ev]['ORIGIN']+self.statsDict[ev][sta]['TOBSS']
                        self.spickDict[ev].update({sta:pickval})
                    except:
                        self.spickDict[ev].update({sta:np.nan})

    def setSystem(self):
        """Get the axial system from waveforms in the stream"""
        comps=[x.stats.channel[2] for x in self.stream]
        if "N" in comps: # if North exists, then it's ZNE
            self.inpSYS.setText("ZNE")
        elif "T" in comps: # if Transverse exists, then it's ZFS
            self.inpSYS.setText("ZFS")


    def saveSolution(self):
        """Save current solution to DB"""
        DB.addValues(self.dbConn,self.dbCur,self.splittingDict)

    def getSolution(self):
        """Get solution values from DB."""
        method="MAN"; ok=True
        _,vals=DB.retrieve(self.dbCur,"method","method","%s/%s/*" % (self.activeEvent,self.station))
        accMethods=sorted(list(set([x[0] for x in vals])))
        if not accMethods:
            genTxt="Could not find methods in database"
            infTxt="There are no available solutions in the database"
            winTit="Solution Not Found"
            self.warnMsgBox(genTxt,infTxt,winTit)
            return
        try:
            idx=accMethods.index(method)
        except:
            idx=0
        method, ok = QtWidgets.QInputDialog.getItem(
            self,"Method to load","Select method: ",accMethods,idx,False
                                                  )
        if not ok:
            return
        self.method=method
        self.updateFromSplittingDict(method)
        self.updateFig(keepAxisLimits=True)
        self.updateTitle()

    def getMaxAin(self):
        """
        Get a maximum Angle of incidence for processing. This is true for both
        event-station viewing and Cluster Analysis
        """
        ain=False; ok=True
        while not ain and ok:
            ain, ok = QtWidgets.QInputDialog.getText(
                        self,"Shear-wave Window",
                        "Current Max: %.1f\nUser defined angle of incidence (degrees): " % self.maxAin
                                                    )
        if not ok:
            return
        if ain.upper() == "INF":
            ain=np.inf
        self.maxAin=float(ain)

    def getPhi(self):
        """Get phi value from user input."""
        phi, ok = QtWidgets.QInputDialog.getDouble(
            self,"S Fast Polarity","User defined phi (degrees):",
            self.phi,0,180,2
            )
        if not ok:
            return
        self.phi=float(phi) # fix for cases like int("110.0")
        self.updateSplittingDict()
        DB.addValues(self.dbConn,self.dbCur,self.splittingDict)        
        self.updateTitle()

    def getGrade(self):
        """Get weight value from user input."""
        weight="X"; ok=True
        accGrades=("A","B","C","D","E","N","X")
        idx=accGrades.index(self.grade)
        grade, ok = QtWidgets.QInputDialog.getItem(
            self,"Measurement Grade (A-E or N)","Set grade: ",accGrades,idx,False
                                                  )
        if not ok:
            return
        self.grade=grade
        if self.method == "MAN":
            self.grade_score=np.nan
        self.updateSplittingDict()
        DB.addValues(self.dbConn,self.dbCur,self.splittingDict)
        self.updateTitle()

    def getComment(self):
        """Get comment value from user input."""
        comment, ok = QtWidgets.QInputDialog.getText(
            self,"Measurement Comment","Comment:",
            QtWidgets.QLineEdit.Normal,self.comment
                                                    )
        if not ok:
            return
        self.comment=comment
        self.updateSplittingDict()
        DB.addValues(self.dbConn,self.dbCur,self.splittingDict)        
        self.updateTitle()

    def applyXAxisLimits(self):
        """Apply the change in axis limits."""
        axInp, ok = QtWidgets.QInputDialog.getText(self,"X Axis Limits","Limits (Min,Max):")
        axLimits=str(axInp).split(",")
        while len(axLimits) != 2 and ok:
            axLimits, ok = QtWidgets.QInputDialog.getText(self,"X Axis Limits","Limits (Min,Max):")
            axLimits=str(axInp).split(",")
        if not ok:
            return
        self.setAxisLimits(self.axPolar,"x",float(axLimits[0]),float(axLimits[1]))


    def applyFixed110Filter(self):
        """Apply a fixed 1-10 Hz filter."""
        self.applyBandFilter(freqmin=1,freqmax=10)

    def applyFixed120Filter(self):
        """Apply a fixed 1-20 Hz filter."""
        self.applyBandFilter(freqmin=1,freqmax=20)

    def applyRecFilter(self):
        """
        Find the dominant frequency in the window and recommend 
        filter for analysis,

        """
        # get unfiltered data
        if self.streamIni:
            stream=self.streamIni
        else:
            stream=self.stream
        if tools.isnone(self.minPick) or tools.isnone(self.maxPick):
            logging.warning("Both start and end of window must be picked!")
            return
        sps=stream[0].stats.sampling_rate
        picks=(self.minPick,self.maxPick)
        pickwin=(int(round(min(picks)*sps)),int(round(max(picks)*sps)))
        dFreq=self.getDomFreq(stream,pickwin)
        # set the filter boundaries
        self.freqmin=dFreq*0.5
        self.freqmax=dFreq*1.5
        self.filtered=True
        if self.method.startswith("M"):
            self.updateSplittingDict()
        self.updateFig(keepAxisLimits=True)

    def applyBandFilter(self,freqmin=False,freqmax=False):
        """
        Apply bandpass filter

        :type freqmin: float or bool, optional
        :param freqmin: the minimum frequency of the filter. If False, it will not be used.
            Defaults to False.
        :type freqmax: float or bool, optional
        :param freqmax: the maximum frequency of the filter. If False, it will not be used.
            Defaults to False.

        """
        if not freqmin and not freqmax:
            filtCorns, ok = QtWidgets.QInputDialog.getText(self,"Bandpass Filter","Filter Corners (Min,Max): ")
            freqCorns=str(filtCorns).split(",")
            while len(freqCorns) != 2 and ok:
                filtCorns, ok = QtWidgets.QInputDialog.getText(self,"Bandpass Filter","Filter Corners (Min,Max): ")
                freqCorns=str(filtCorns).split(",")        
            if not ok:
                return
            self.freqmin=float(freqCorns[0]); self.freqmax=float(freqCorns[1])
        else:
            self.freqmin=freqmin
            self.freqmax=freqmax
        self.filtered=True
        if self.method.startswith("M"):
            self.updateSplittingDict()
        self.updateFig(keepAxisLimits=True)

    def removeFilter(self):
        """Removes any filter applied in the traces."""
        self.filtered=False
        self.freqmin=np.nan; self.freqmax=np.nan     
        if self.method.startswith("M"):
            self.updateSplittingDict()        
        self.updateFig(keepAxisLimits=True)

    def goToInitial(self):
        """Goes to the initial stage."""
        self.stageRotated=False 
        self.stageCorrected=False
        self.updateFig(center=False,keepAxisLimits=True)

    def goToRotated(self):
        """Goes to the Rotated stage."""
        if tools.isnone(self.phi):
            winTitle="Invalid value for S<sub>fast</sub> pol direction!"
            genText="Invalid value"
            infText="You must provide a valid phi value to rotate."
            self.warnMsgBox(genText,infText,winTitle)
            return
        self.stageRotated=True
        self.stageCorrected=False
        self.updateFig(center=False,keepAxisLimits=True)

    def goToCorrected(self):
        """Goes to the Corrected stage."""
        self.stageRotated=False
        self.stageCorrected=True
        self.updateFig(center=False,keepAxisLimits=True)

    def applyTimedelay(self,arrowInput=False,td=False):
        """
        Apply the time delay correction

        :type arrowInput: bool, optional
        :param arrowInput: declares whether the arrow keys were used
            to change the time-delay. If False, a user dialog will
            ask for a value. Defaults to False.
        :type td: float or bool, optional
        :param td: the time-delay to be added. Defaults to False.
        """
        if not arrowInput:
            td, ok = QtWidgets.QInputDialog.getDouble(
                self,"Time-delay","Input time-delay (ms):",
                abs(self.td),0,10**6,3
                                                    )
            if not ok:
                    return
            self.td=-float(td)/1000.
        if self.td > 0.0: 
            self.td=0.0
            #return
        self.td=self.tdCheck(self.td)
        if self.stageCorrected:
            self.stream=self.rotation(self.phi,self.stream,"NE->RT")
        data=self.timedelay(self.td,self.stream)
        self.stream.select(component="T")[0].data=data
        # Correct lengths
        Wz=self.stream.select(component="Z")[0]
        Wy=self.stream.select(component="R")[0]
        Wx=self.stream.select(component="T")[0]
        # First horizontals
        Ax,Ay=tools.lengthcheck(Wx,Wy)
        Wx.data=Ax 
        Wy.data=Ay
        # Now the vertical
        Az,Ay=tools.lengthcheck(Wz,Wy)
        Wz.data=Az
        Wy.data=Ay
        # Horizontals second time
        Ax,Ay=tools.lengthcheck(Wx,Wy)
        Wx.data=Ax 
        Wy.data=Ay
        # get either aux or pol
        T=self.polarDict["TIME"]
        D=self.polarDict["ANGLE"]                
        # get angles
        pAngle=D[np.abs(T-self.sPick).argmin()] % 180
        # Convert to Fast-Slow according to quadrant.
        pAngle=self.mpl2ns(pAngle)        
        # get polarization angle?
        if self.stageRotated and self.method is "MAN":
            self.aux=pAngle
            self.pol=(self.aux+self.phi) % 180
        if self.stageCorrected and self.method is "MAN":
            self.pol=pAngle
        # update splitting dict
        self.updateSplittingDict()
        if not arrowInput:
            self.updateFig(center=False,keepAxisLimits=True)

    def tdCheck(self,td):
        """
        Checks whether provided time-delay is a multiple of the 
        sampling rate.

        TODO: for now this is disabled.

        :type td: float
        :param td: the time-delay
        :returns: the corrected time-delay
        """
        delta=self.stream[0].stats.delta
        #return np.ceil(td/delta)*delta
        return td

    ## Splitting Methods ##
    def applyCA(self,meth,stream,filtered,freqmin,freqmax,sPick,
                baz,ain,showProg,runMult,eventCode=None):
        """
        Function to apply Cluster Analysis on a single station-event pair.

        :type meth: str
        :param meth: the method used for analyzing splitting (EV, RC or ME)
        :type stream: :class: `~obspy.core.stream.Stream`
        :param stream: the waveforms to be used in the analysis
        :type filtered: bool
        :param filtered: declare whether the waveforms should be filtered.
        :type freqmin: float
        :param freqmin: the minimum bound of the filter
        :type freqmax: float
        :param freqmax: the maximum bound of the filter
        :type sPick: float
        :param sPick: the arrival time of the S-wave relative to the stream's
            start-time (in s).
        :type baz: float
        :param baz: the backazimuth
        :type ain: float
        :param ain: the angle of incidence
        :type showProg: bool
        :param showProg: show progress dialog
        :type eventCode: str or None, optional
        :param eventCode: the event code. If it is None, the self.activeEvent value
            will be used instead. Defaults to None.

        """
        if eventCode is None:
            eventCode=self.activeEvent
        station=stream[0].stats.station
        self.station=station
        ## assign the proper prefix to the method
        if meth == "EV":
            self.method="AEV"
        elif meth == "ME":
            self.method="AME"
        elif meth == "RC":
            self.method="ARC"
        ## CA prep
        logging.debug("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        logging.debug("Event code: %s - Station: %s" % (eventCode,station))
        # get S-P
        sArr=self.sArr
        if sPick == 0 and tools.isnone(sArr):
            try:
                sArr=UTCDateTime(self.splittingDict[eventCode][station]["s_obs"])
                sPick=sArr-self.stream[0].stats.starttime
            except KeyError:
                try:
                    sArr=self.spickDict[eventCode][station]
                    sPick=sArr-self.stream[0].stats.starttime
                except KeyError:
                    logging.warning("No manual arrival")
                    return
            except:
                logging.warning("Must add theoretical here?")
                return
        self.sArr=sArr; self.sPick=sPick
        if sPick:
            try:
                tpts=((stream[1].stats.starttime+sPick)-self.evsDict[eventCode]["ORIGIN"])-\
                     self.statsDict[eventCode][stream[0].stats.station]["TOBSP"]
                # apply a maximum S-P limit to avoid overextending the window in cases of deep earthquakes,
                # where S-P can be long for small incidences
                self.tpts=tpts
                tpts=self.caCNF.tptsMax if tpts > self.caCNF.tptsMax else tpts
            except:
                logging.exception("Error in S-P calculation, using default (%.2f s)" % self.caCNF.tptsMax)
            tpts=self.caCNF.tptsMax
            self.tpts=tpts                 
        else:
            logging.warning("No S pick! Using fixed Tp-Ts = %.2f s" % self.caCNF.tptsMax)
            tpts=self.caCNF.tptsMax 
            self.tpts=tpts
        logging.info("Ts-Tp is %.3f s" % tpts)
        # get dominant period for S
        pickwin=int(np.floor(sPick*stream[0].stats.sampling_rate)),\
                int(np.ceil((self.caCNF.specWindow+sPick)*stream[0].stats.sampling_rate))
        dPer=1./self.getDomFreq(stream,pickwin)
        self.dPer=dPer
        # limit period to accepted limits
        if dPer < self.caCNF.minPeriod:
            dPer=self.caCNF.minPeriod
        elif dPer > self.caCNF.maxPeriod:
            dPer=self.caCNF.maxPeriod
        # constrain dominant period
        logging.debug("Dominant period is %.4f s" % dPer)
        # calculate window parameters
        self.caCNF.Tbeg0=-tpts/2.
        self.caCNF.Tend0=dPer#+self.caCNF.Tbeg1
        self.caCNF.Tend1=(self.caCNF.multPeriod*dPer)#+self.caCNF.Tbeg0
        self.caCNF.Nbeg=int(np.ceil((self.caCNF.Tbeg1-self.caCNF.Tbeg0)/self.caCNF.DTbeg))
        self.caCNF.Nend=int(np.ceil((self.caCNF.Tend1-self.caCNF.Tend0)/self.caCNF.DTend))
        logging.debug("CA beg windows: [%.4f , %.4f]" % (self.caCNF.Tbeg0,self.caCNF.Tbeg1))
        logging.debug("CA end windows: [%.4f , %.4f]" % (self.caCNF.Tend0,self.caCNF.Tend1))
        logging.debug("N beg/end windows: %i x %i" % (self.caCNF.Nbeg,self.caCNF.Nend))
        ## 
        tbStream=stream.copy()
        if filtered:
            tbStream.filter(type="bandpass",corners=4,freqmin=freqmin,freqmax=freqmax,zerophase=True)
        if not runMult:
            self.progBar=self.prgBar("Processing time windows, this may take a few minutes...",minVal=0,maxVal=100)
            self.progBar.setWindowTitle("Cluster Analysis")
            self.progBar.setWindowIcon(QtGui.QIcon(self.appIcon))
            self.progBar.setValue(0)
            self.progBar.show()
        if self.actQT.isChecked():
            logging.debug("Will rotate to LQT")
        else:
            logging.debug("Will rotate to ZRT")
            ain=0.            
        self.CAthread=TB.clustering(tbStream,meth,sPick,baz,ain,self.caCNF)
        if not runMult:
            if showProg:
                self.CAthread.iterDone.connect(self.updateProgBar)
            self.CAthread.iterFail.connect(self.failCA)            
            self.CAthread.finished.connect(self.doneCA)
        self.CAthread.start()
       
    def failCA(self,exc):
        """
        Handling exception in the CA threads

        :type exc: Exception
        :param exc: the Exception that was raised during CA

        """
        try: # if bar wasnt created
            self.progBar.hide()
        except:
            pass
        try:
            raise exc
        except:
            logging.exception("Teanby process failed!")
            tbStack=traceback.format_exc()
            genTxt="Clustering process failure"
            infTxt="Clustering process exited with exception:\n %s" % tbStack
            winTit="Clustering Failure"
            self.warnMsgBox(genTxt,infTxt,winTit)            
        

    def doneCA(self):
        """Handling success in the CA threads."""
        if self.CAthread.excFlag:
            return
        try: # if progBar wasnt created
            self.progBar.hide()
        except:
            pass
        self.phi=self.CAthread.phi % 180
        self.td=self.CAthread.dt
        # grab all required values
        self.minPick,self.maxPick=self.CAthread.optWindow
        self.errors=(self.CAthread.sphi,self.CAthread.sdt)
        _,_,self.cArray,self.pTest,self.dTest,\
        self.pol,_,_,self.cEmax,\
        self.CC_FS,self.CC_NE,self.Tvar,self.nContours=self.CAthread.tempRes 
        self.grade_score,self.grade=tools.autoGrading(
                    (self.phi,self.pol,self.gradeCNF.polOff),
                    (self.errors[0],self.errors[1],self.SNR,self.CC_FS,self.CC_NE), # values
                    (self.gradeCNF.error_bounds[0],
                     self.gradeCNF.error_bounds[1],
                     self.gradeCNF.snr_bound,
                     self.gradeCNF.CC_FS_bound, self.gradeCNF.CC_NE_bound),  # bounds
                    self.gradeCNF.gradeDict
                                    )
        # update CA related stuff
        self.initialClusters=self.CAthread.initial; self.initialClustersErrors=self.CAthread.initial_err
        self.calinski=self.CAthread.calinski
        self.clusters1=self.CAthread.clusters1; self.clusters2=self.CAthread.clusters2          
        # write results to splitting dict
        self.updateSplittingDict()
        DB.addValues(self.dbConn,self.dbCur,self.splittingDict)
        # get dict with solution
        method=self.method
        metDict=self.splittingDict[self.activeEvent][self.station][method]
        # make the figure(s)
        titleInfo=str(self.origtime),self.ain,self.dist,self.mag,metDict["grade"]
        tools.qcFigSWS(
                       self.streamIni.copy(),metDict["phi"],metDict["td"]/(10**3),metDict["pol"],
                       (metDict["err_phi"],metDict["err_td"]/(10**3)),
                       self.baz,
                       (
                        UTCDateTime(metDict["window_min"])-self.stream[0].stats.starttime,
                        UTCDateTime(metDict["window_max"])-self.stream[0].stats.starttime
                        ),
                       metDict["C_array"],metDict["phi_test"],
                       metDict["td_test"],method[1:],
                       metDict["C_max"],(metDict["filter_min"],metDict["filter_max"]),titleInfo
                      ) 
        tools.qcClusterDiagram(
         metDict["initial_clusters"],
         metDict["initial_clusters_errors"],
         metDict["calinski_score"],
         metDict["clusters1"],metDict["clusters2"],
         self.stream[0].stats.station,
         self.stream[0].stats.network,self.baz,
         titleInfo, method[1:],
         (metDict["filter_min"],metDict["filter_max"])
                                )
        plt.show()
        #
        self.updateFig(center=False,keepAxisLimits=True)


    ## cluster analysis multi
    def applyCCA(self,meth,selectedStations,maxAin,filterBounds,skipPairs,skipFailed):
        """
        Function to apply Cluster Analysis on the whole catalogue.

        :type meth: str
        :param meth: analysis method to use with CA (EV, ME or RC)
        :type selectedStations: tuple-like
        :param selectedStations: list of stations to be used in the analysis
        :type maxAin: float
        :param maxAin: the incidence angle corresponding to the shear-wave window
        :type filterBounds: tuple-like
        :param filterBounds: the minimum and maximum filter boundaries
        :type skipPairs: bool
        :param skipPairs: select whether to skip pairs that already exist in the database
        :type skipFailed: bool
        :param skipFailed: select whether to skip pairs that previously led to errors during CA

        """
        logging.debug("===================================================")
        logging.debug("Starting Catalogue CA process with method %s" % meth)
        # save original values to reload after finishing here 
        self.defaults=(
                  self.activeEvent,self.station,self.method,self.SNR,self.sPick,
                  self.dPer,self.tpts,self.phi,self.td,self.pol,self.CC_FS,self.CC_NE,self.Tvar,
                  self.errors,self.nContours,self.grade_score,self.grade,self.freqmin,self.freqmax,
                  self.minPick,self.maxPick,self.comment,self.maxAin,self.stream.copy()
                )
        # set values
        self.freqmin=filterBounds[0]
        self.freqmax=filterBounds[1]
        self.maxAin=maxAin
        # init custom logging file for process
        tic=UTCDateTime()
        actEvList=sorted(self.getActiveEvents())
        if meth == "EV":
            self.method="AEV"
        elif meth == "ME":
            self.method="AME"
        elif meth == "RC":
            self.method="ARC"
        ######################################
        self.progBar=self.prgBar("Initiating Catalogue Cluster Analysis...",minVal=0,maxVal=100,runMult=True)
        self.progBar.setWindowTitle("Cluster Analysis")
        self.progBar.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.progBar.setValue(0)
        self.progBar.show()          
        ## thread
        self.CCAthread=applyCAMult(self,self.method,selectedStations,maxAin,filterBounds,skipPairs,skipFailed,actEvList)
        self.CCAthread.iterDone.connect(self.updateProgBar)
        self.CCAthread.iterFail.connect(self.failCCA)
        self.CCAthread.finished.connect(self.doneCCA)
        self.CCAthread.start()

    def failCCA(self,exc):
        """
        Handling exception in the Catalogue CA thread

        :type exc: Exception
        :param exc: the Exception that was raised during CA

        """
        try: # if bar wasnt created
            self.progBar.hide()
        except:
            pass
        try:
            raise exc
        except UserCancelException:
            logging.exception("Catalogue CA process cancelled!")
            logging.debug("===================================================")
            genTxt="Clustering process"
            infTxt="The user cancelled the operation!"
            winTit="Clustering Cancelled"
            self.warnMsgBox(genTxt,infTxt,winTit)  
        except:
            logging.exception("Catalogue CA process failed!")
            logging.debug("===================================================")
            tbStack=traceback.format_exc()
            genTxt="Clustering process failure"
            infTxt="Clustering process exited with exception:\n %s" % tbStack
            winTit="Clustering Failure"
            self.warnMsgBox(genTxt,infTxt,winTit)  

    def doneCCA(self):
        """Handling exception in the Catalogue CA thread."""
        if self.CCAthread.excFlag:
            return
        ######################################
        tic=self.CCAthread.tic
        toc=UTCDateTime()
        diff=str(datetime.timedelta(seconds=toc-tic))
        with open(self.CCAthread.logfileCA,"a") as fid:
            fid.write("# Process finished @ %s\n" % toc)   
            fid.write("# Number of pairs analyzed (done/fail): %i/%i\n" % (self.CCAthread.ispairs,self.CCAthread.ipairs))
            fid.write("# Total time elapsed (s): %s" % diff)
        # reload previous settings
        self.activeEvent,self.station,self.method,self.SNR,self.sPick,\
        self.dPer,self.tpts,self.phi,self.td,self.pol,self.CC_FS,self.CC_NE,self.Tvar,\
        self.errors,self.nContours,self.grade_score,self.grade,self.freqmin,self.freqmax,\
        self.minPick,self.maxPick,self.comment,self.maxAin,self.stream=self.defaults
        # finish off
        try:
            self.progBar.hide()
        except:
            pass
        logging.debug("Catalogue CA completed successfully!")
        logging.debug("===================================================")
        # let the user know the process is completed
        winTitle="Catalogue CA"
        genText="CA Process completed!"
        infText="Analyzed %i pairs in %s !" % (self.CCAthread.ipairs,diff.split(".")[0])
        self.infoMsgBox(genText,infText,winTitle)

    def applyEV(self,meth):
        """
        Function to apply the Eigenvalue or Minimum Energy method.

        :type meth: str
        :param meth: select the method to use (EV or ME)

        """ 
        try:
            if tools.isnone(self.minPick) or tools.isnone(self.maxPick):
                raise ValueError("Both start and end of window must be picked!")           
            picks=(self.minPick,self.maxPick)
            pickwin=(min(picks),max(picks))
            if np.nan in (self.baz,self.ain,self.minPick,self.maxPick):
                logging.warning("Not all required values for EV/ME set")
                raise ValueError("Not all required values (e.g. picks)")
            if pickwin[1]-pickwin[0] < self.generalCNF.maxTd/1000.:
                logging.exception("Signal window shorter than maximum time-delay!")
                winTitle="Signal Window Warning"
                genText="Signal window length is invalid"
                infText="Signal window for analysis should be longer than %.3f s (current: %.3f s)" \
                         % ((self.generalCNF.maxTd/1000.),(pickwin[1]-pickwin[0]))
                self.warnMsgBox(genText,infText,winTitle)
                return
            # calculate SNR
            try:
                self.SNR=self.calcSNR(self.sPick)
            except: 
                self.SNR=self.calcSNR(round(int((self.minPick+self.maxPick)/2))) 
            # perform the grid search
            if self.actQT.isChecked():
                logging.debug("Will rotate to LQT")
                ain=self.statsDict[self.activeEvent][self.stream[0].stats.station]["AN"]
            else:
                logging.debug("Will rotate to ZRT")
                ain=0.
            # reset waveforms
            if self.streamIni:
                self.stream=self.streamIni.copy()
            else:
                self.streamIni=self.stream.copy()
            self.stageRotated=False; self.stageCorrected=False
            if self.filtered:
                self.stream.filter(type="bandpass",freqmin=self.freqmin,freqmax=self.freqmax,zerophase=True,corners=4)
            tempRes=SC.SilverAndChan(self.stream,self.baz,pickwin,ain=ain,method=meth,maxDelay=self.generalCNF.maxTd/1000.)
            self.phi,self.td,self.cArray,self.pTest,self.dTest,\
            self.pol,stdPhi,stdTd,self.cEmax,\
            self.CC_FS,self.CC_NE,self.Tvar,self.nContours=tempRes
            if meth == "ME":
                self.method="MME"
            elif meth == "EV":
                self.method="MEV"
            self.errors=(stdPhi,stdTd)
            # invert sign
            self.td=-self.td
            self.phi%=180
            #
            self.grade_score,self.grade=tools.autoGrading(
                        (self.phi,self.pol,self.gradeCNF.polOff),
                        (self.errors[0],self.errors[1],self.SNR,self.CC_FS,self.CC_NE), # values
                        (self.gradeCNF.error_bounds[0],
                         self.gradeCNF.error_bounds[1],
                         self.gradeCNF.snr_bound,
                         self.gradeCNF.CC_FS_bound, self.gradeCNF.CC_NE_bound),  # bounds
                        self.gradeCNF.gradeDict
                                        )
            # write results to dict
            self.updateSplittingDict()
            DB.addValues(self.dbConn,self.dbCur,self.splittingDict)
            ## plot QC figure
            # get dict with solution
            method="M"+meth
            metDict=self.splittingDict[self.activeEvent][self.station][method]
            # make the figure(s)
            titleInfo=str(self.origtime),self.ain,self.dist,self.mag,metDict["grade"]
            tools.qcFigSWS(
                           self.streamIni.copy(),metDict["phi"],metDict["td"]/(10**3),metDict["pol"],
                           (metDict["err_phi"],metDict["err_td"]/(10**3)),
                           self.baz,
                           (
                            UTCDateTime(metDict["window_min"])-self.stream[0].stats.starttime,
                            UTCDateTime(metDict["window_max"])-self.stream[0].stats.starttime
                            ),
                           metDict["C_array"],metDict["phi_test"],
                           metDict["td_test"],method[1:],
                           metDict["C_max"],(metDict["filter_min"],metDict["filter_max"]),titleInfo
                          )  
            plt.show()            
        except:
            logging.exception("SC failure!")
            tbStack=traceback.format_exc()
            genTxt="%s Application failed" % meth
            infTxt="%s failed with exception:\n %s" % (meth,tbStack)
            winTit="%s Application" % meth
            self.warnMsgBox(genTxt,infTxt,winTit)
        self.updateFig(center=False,keepAxisLimits=True)

    def applyRC(self):
        """Function to apply the Rotation-Correlation method."""
        try:
            if tools.isnone(self.minPick) or tools.isnone(self.maxPick):
                raise ValueError("Both start and end of window must be picked!") 
            picks=(self.minPick,self.maxPick)
            pickwin=(min(picks),max(picks))
            if np.nan in (self.baz,self.ain,self.minPick,self.maxPick):
                logging.warning("Not all required values for EV/ME set")
                raise ValueError("Not all required values (e.g. picks)")
            if pickwin[1]-pickwin[0] < self.generalCNF.maxTd/1000.:
                logging.exception("Signal window shorter than maximum time-delay!")
                winTitle="Signal Window Warning"
                genText="Signal window length is invalid"
                infText="Signal window for analysis should be longer than %.3f s (current: %.3f s)" \
                         % ((self.generalCNF.maxTd/1000.),(pickwin[1]-pickwin[0]))
                self.warnMsgBox(genText,infText,winTitle)
                return
            # calculate SNR
            try:
                self.SNR=self.calcSNR(self.sPick)
            except: 
                self.SNR=self.calcSNR(round(int((self.minPick+self.maxPick)/2)))             
            # set method
            self.method="MRC"
            meth="RC"
            # reset waveforms
            if self.streamIni:
                self.stream=self.streamIni.copy()
            else:
                self.streamIni=self.stream.copy()
            self.stageRotated=False; self.stageCorrected=False
            if self.filtered:
                self.stream.filter(type="bandpass",freqmin=self.freqmin,freqmax=self.freqmax,zerophase=True,corners=4)
            # perform the cross-correlation
            ain=self.statsDict[self.activeEvent][self.stream[0].stats.station]["AN"]
            tempRes=RC.crosscorrelation(self.stream,pickwin,maxDelay=self.generalCNF.maxTd/1000.)
            self.phi,self.td,self.cArray,self.pTest,self.dTest,\
            self.pol,stdPhi,stdTd,self.cEmax,\
            self.CC_FS,self.CC_NE,self.Tvar,self.nContours=tempRes
            if meth == "ME":
                self.method="MME"
            elif meth == "EV":
                self.method="MEV"
            self.errors=(stdPhi,stdTd)
            # invert sign
            self.td=-self.td
            self.phi%=180
            #
            self.grade_score,self.grade=tools.autoGrading(
                        (self.phi,self.pol,self.gradeCNF.polOff),
                        (self.errors[0],self.errors[1],self.SNR,self.CC_FS,self.CC_NE), # values
                        (self.gradeCNF.error_bounds[0],
                         self.gradeCNF.error_bounds[1],
                         self.gradeCNF.snr_bound,
                         self.gradeCNF.CC_FS_bound, self.gradeCNF.CC_NE_bound),  # bounds
                        self.gradeCNF.gradeDict
                                        )
            # write results to dict
            self.updateSplittingDict()
            DB.addValues(self.dbConn,self.dbCur,self.splittingDict)
            ## plot QC figure
            method="M"+meth
            # get dict with solution
            metDict=self.splittingDict[self.activeEvent][self.station][method]
            # make the figure(s)
            titleInfo=str(self.origtime),self.ain,self.dist,self.mag,metDict["grade"]
            tools.qcFigSWS(
                           self.streamIni.copy(),metDict["phi"],metDict["td"]/(10**3),metDict["pol"],
                           (metDict["err_phi"],metDict["err_td"]/(10**3)),
                           self.baz,
                           (
                            UTCDateTime(metDict["window_min"])-self.stream[0].stats.starttime,
                            UTCDateTime(metDict["window_max"])-self.stream[0].stats.starttime
                            ),
                           metDict["C_array"],metDict["phi_test"],
                           metDict["td_test"],method[1:],
                           metDict["C_max"],(metDict["filter_min"],metDict["filter_max"]),titleInfo
                          )  
            plt.show() 
        except:
            logging.exception("RC failure")
            tbStack=traceback.format_exc()
            genTxt="%s Application failed" % meth
            infTxt="%s failed with exception:\n %s" % (meth,tbStack)
            winTit="%s Application" % meth
            self.warnMsgBox(genTxt,infTxt,winTit)  
        self.updateFig(center=False,keepAxisLimits=True)


    def getDomFreq(self,stream,pickwin):
        """
        Calculate the dominant frequency in given window

        :type stream: :class: `~obspy.core.stream.Stream`
        :param stream: the stream containing the waveforms
        :type pickwin: tuple-like
        :param pickwin: a list of the two picks defining the signal
            window
        :returns: the dominant frequency
        """
        sps=stream[0].stats.sampling_rate
        north=stream.select(component="N")[0].data[pickwin[0]:pickwin[1]]
        east=stream.select(component="E")[0].data[pickwin[0]:pickwin[1]]
        specs=tools.spectrum(north,east,sps)
        nSpec=specs[0]; eSpec=specs[1]
        # filter out frequencies less than 1Hz. 
        # TODO: add as user option
        minFreqIdx=(nSpec[0]>=1,eSpec[0]>=1)
        # get max frequency
        nFreq=nSpec[0][minFreqIdx[0]][np.where(nSpec[1]==nSpec[1].max())][0]
        eFreq=eSpec[0][minFreqIdx[1]][np.where(eSpec[1]==eSpec[1].max())][0]
        dFreq=(nFreq+eFreq)/2.
        return dFreq

    def plotSpectrum(self):
        """Plot the spectra used for recommending the filter."""
        if self.streamIni:
            stream=self.streamIni
        else:
            stream=self.stream
        if tools.isnone(self.minPick) or tools.isnone(self.maxPick):
            logging.warning("Both start and end of window must be picked!")
            return
        sps=stream[0].stats.sampling_rate
        picks=(self.minPick,self.maxPick)
        pickwin=(int(round(min(picks)*sps)),int(round(max(picks)*sps)))
        #
        north=stream.select(component="N")[0].data[pickwin[0]:pickwin[1]]
        east=stream.select(component="E")[0].data[pickwin[0]:pickwin[1]]
        specs=tools.spectrum(north,east,sps)
        nSpec=specs[0]; eSpec=specs[1]
        # filter out frequencies less than 1Hz. 
        # TODO: add as user option
        minFreqIdx=([nSpec[0]>=1],[eSpec[0]>=1])
        # get max frequency
        nFreq=nSpec[0][minFreqIdx[0]][np.where(nSpec[1]==nSpec[1].max())]
        eFreq=eSpec[0][minFreqIdx[1]][np.where(eSpec[1]==eSpec[1].max())]
        # get max frequency
        nFreq=nSpec[0][np.where(nSpec[1]==nSpec[1].max())]
        eFreq=eSpec[0][np.where(eSpec[1]==eSpec[1].max())]
        dFreq=(nFreq+eFreq)/2.
        freqmin=dFreq*0.5
        freqmax=dFreq*1.5
        # make the plot
        fig,axes=plt.subplots(2)
        axes[0].set_ylabel("North")
        axes[1].set_ylabel("East")
        fig.suptitle("f0: %.2f\nRecommended boundaries: %.2f - %.2f Hz" % (dFreq,freqmin,freqmax))
        for ax,spec,freq in zip(axes,(nSpec,eSpec),(nFreq,eFreq)):
            ax.plot(spec[0],spec[1])
            ax.axvline(freq,linestyle='--')
            ax.set_xlabel("Frequency (Hz)")
        # show (full screen please!)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()

    def flagReset(self,general=True,splitting=True,basicSplit=True,secSplit=True,
                       indexing=True,titling=True,splitDict=True):
        """
        Resets basic parameters and other flags

        :type general: bool, optional
        :param general: resets general values and flags. Defaults to True.
        :type splitting: bool, optional
        :param splitting: resets splitting values and flags. Defaults to True.
        :type basicSplit: bool, optional
        :param basicSplit: resets basic splitting values and flags. Defaults to True.
        :type secSplit: bool, optional
        :param secSplit: resets secondary splitting values and flags. Defaults to True.
        :type indexing: bool, optional
        :param indexing: resets indexing values and flags. Defaults to True.
        :type titling: bool, optional
        :param titling: resets values and flags related to the information widget. Defaults to True.
        :type splitDict: bool, optional
        :param splitDict: sets the splittingDict to empty. Defaults to True

        """
        if general:
            self.V2H=False
            self.filtered=False
            self.freqmin=np.nan
            self.freqmax=np.nan
            self.stageRotated=False
            self.stageCorrected=False
            self.streamIni=None
            self.streamRotated=False
            self.streamCorrected=False
            try:
                for vspan in self.vspanList:
                    try:
                        vspan.remove()
                    except ValueError:
                        continue
            except AttributeError: # if not yet created
                self.vspanList=[]
            for ax in self.activeFig.get_axes():
              if ax == self.axPolar:
                try:
                  ax.lines[0].set_xdata(0)
                except:
                  pass
              else:
                try:
                  ax.lines[1].set_xdata(0)
                except:
                  pass
        # splitting parameters
        if splitting:
            if basicSplit:
                self.phi=np.nan
                self.td=0.
                self.sArr=np.nan
                self.sPick=0.
                self.aux=np.nan
                self.pol=np.nan
                self.SNR=np.nan
            if secSplit:
                self.Tvar=np.nan
                self.nContours=0
                self.errors=(np.nan,np.nan)
                self.grade="X"
                self.grade_score=np.nan
                self.comment=""
                self.method="MAN"
                self.dPer=np.nan
                self.tpts=np.nan
                self.minPick=np.nan
                self.maxPick=np.nan
                self.CC_NE=np.nan
                self.CC_FS=np.nan
                self.pTest=np.zeros((1,1))
                self.dTest=np.zeros((1,1))
                self.cEmax=np.nan
                self.cArray=np.zeros((1,1))
                self.initialClusters=np.zeros((1,1))
                self.initialClustersErrors=np.zeros((1,1))
                self.calinski=np.zeros((1,1))
                self.clusters1=np.zeros((1,1))
                self.clusters2=np.zeros((1,1))
        if splitDict:
            self.splittingDict={}
        # indexing parameters
        if indexing:
            self.evIdx=0
            self.stIdx=0
        # title parameters
        if titling:
            self.station=np.nan
            self.orCorrection=np.nan
            self.ain=np.nan
            self.dist=np.nan
            self.baz=np.nan
            self.mag=np.nan
            self.origtime="None"

    ## General Purpose ##
    def prgBar(self,labelText,runMult=False,cancelButtonText="Cancel",minVal=0,maxVal=100):
        """
        Progress dialog wrapper function for PyQt4 progress

        :type labelText: str
        :param labelText: sets the text at the bar's label
        :type runMult: bool, optional
        :param runMult: enables CCA mode for the progress bar. Defaults to False.
        :type cancelButtonText: str, optional
        :param cancelButtonText: sets the text at the dialog's cancel button. Defaults
            to 'Cancel'.
        :type minVal: int, optional
        :param minVal: sets the minimum progress value. Defaults to 0.
        :type maxVal: int, optional
        :param maxVal: sets the maximum progress value. Defaults to 100.

        """
        cancBtn=QtWidgets.QPushButton("Cancel")
        if runMult:
            cancBtn.clicked.connect(lambda: self.CCAthread.stop())
        else:
            cancBtn.clicked.connect(lambda: self.CAthread.stop())
        prgDiag=QtWidgets.QProgressDialog(
                   labelText,cancelButtonText,
                   minVal, maxVal
                                        )
        prgDiag.setCancelButton(cancBtn)
        return prgDiag

    def updateProgBar(self,i):
        """
        Updates the progress bar
        
        :type i: list or float or int
        :param i: either a tuple of (label text, progress number) or
            a progress number.

        """
        if isinstance(i,list):
            self.progBar.setLabelText(i[0])
            self.progBar.setValue(int(round(i[1])))
        else:
            self.progBar.setValue(int(float(i))) # the int(float()) thing is for compatibility with Python 2

    def infoMsgBox(self,genText,infText,winTitle):
        """
        Information message box

        :type genText: str
        :param genText: sets the general text of the dialog window.
        :type infText: str
        :param infText: sets the informative text of the dialog window.
        :type winTitle: str
        :param winTitle: sets the title of the dialog window.

        """
        logging.debug("Generating info message box with title %s" % winTitle)
        self.iMsg=QtWidgets.QMessageBox(self)
        self.iMsg.setIcon(QtWidgets.QMessageBox.Information)
        self.iMsg.setText(genText)
        self.iMsg.setInformativeText(infText)
        self.iMsg.setWindowTitle(winTitle)
        self.iMsg.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.iMsg.show()

    def progressDialog(self,title,cancel=None):
        """
        Progress dialog window

        :type title: str
        :param title: sets the title of the dialog window

        """
        logging.debug("Generating title message box with title %s" % title)
        cancel=None
        self.pDial=QtWidgets.QProgressDialog("Please wait...",cancel,0,100,self)
        self.pDial.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.pDial.setWindowTitle(title)
        self.pDial.setWindowIcon(QtGui.QIcon(self.appIcon))
        #self.pDial.setWindowModality(QtCore.Qt.WindowModal)
        self.pDial.show()
        return self.pDial

    def warnMsgBox(self,genText,infText,winTitle):
        """
        Warning message box

        :type genText: str
        :param genText: sets the general text of the dialog window.
        :type infText: str
        :param infText: sets the informative text of the dialog window.
        :type winTitle: str
        :param winTitle: sets the title of the dialog window.

        """
        logging.debug("Generating warning message box with title %s" % winTitle)
        self.wMsg=QtWidgets.QMessageBox(self)
        self.wMsg.setIcon(QtWidgets.QMessageBox.Warning)
        self.wMsg.setText(genText)
        self.wMsg.setInformativeText(infText)
        self.wMsg.setWindowTitle(winTitle)
        self.wMsg.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.wMsg.show()

    def replyMsgBox(self,infText,winTitle):
        """
        Accept/decline message box.

        :type infText: str
        :param infText: sets the informative text of the dialog window.
        :type winTitle: str
        :param winTitle: sets the title of the dialog window.

        """
        reply=QtWidgets.QMessageBox.question(
                                        self,winTitle,infText,
                                        QtWidgets.QMessageBox.Yes,QtWidgets.QMessageBox.No
                                      )
        return reply

    def cleanAll(self):
        """ Cleans all existing files (logs,indices,etc)/"""
        ok = QtWidgets.QMessageBox.question(
                    self,"Quit","Are you sure you want reset? (The program will close)",
                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
                                        )
        if ok == QtWidgets.QMessageBox.No:
            return
        # now evcor
        evcdir=WORKDIR+os.sep+"etc%sevcor%s" % (os.sep,os.sep)
        for fl in os.listdir(evcdir):
            os.remove(evcdir+fl)
        # now indices
        inddir=WORKDIR+os.sep+"etc%sindex%s" % (os.sep,os.sep)
        for fl in os.listdir(inddir):
            os.remove(inddir+fl)
        # now logs
        self.cleanLogs()
        # now quit
        QtWidgets.qApp.quit()


    def cleanLogs(self):
        """Cleans logs created in periods greater than specified days."""
        logging.debug("Cleaning logs...")
        days=self.generalCNF.cleanLogs
        days_in_sec=3600.*24*days
        logs=os.listdir(WORKDIR+os.sep+"logs"+os.sep)
        logs=[WORKDIR+os.sep+"logs"+os.sep+x for x in logs]
        now=UTCDateTime()
        for log in logs:
            if  (now - os.path.getmtime(log)) > days_in_sec:
                logging.debug("Removing logfile %s" % log)
                os.remove(log)

    def onlyQuit(self):
        """Quit Pytheas."""
        ok = QtWidgets.QMessageBox.question(
                    self,"Quit","Are you sure you want to quit?",
                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
                                        )
        if ok == QtWidgets.QMessageBox.Yes:
            try:
                DB.close(self.conn)
            except:
                pass
            QtWidgets.qApp.quit()
        else:
            return

    def closeEvent(self,event):
        """Quit event."""
        ok = QtWidgets.QMessageBox.question(
                    self,"Quit","Are you sure you want to quit?",
                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
                                        )
        if ok == QtWidgets.QMessageBox.Yes:
            QtWidgets.qApp.quit()
        else:
            event.ignore()

    ## Open Catalogue window
    def openCatWindow(self):
        """Open catalogue window."""
        #
        # ask for events master directory
        self.exDataPath,self.exCatFile,self.exDbPath=self.pathParser(mode="read")
        self.dataPath=self.exDataPath; self.catFile=self.exCatFile; self.dbPath=self.exDbPath
        # various parameters
        self.loadCatalogue=False
        self.ocWin=QtWidgets.QDialog(self)
        self.ocWin.setModal(False)
        self.ocWin.setWindowTitle("Open Catalogue")
        self.ocWin.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.ocWin.setMinimumSize(600,300)
        self.ocWin.setsFont1=QtGui.QFont("Calibri",12,QtGui.QFont.Bold)
        self.ocWin.setsFont2=QtGui.QFont("Calibri",11,QtGui.QFont.Normal)
        self.ocWin.gBox=QtWidgets.QGridLayout(self.ocWin)
        ## add labels/edits
        # first for the data directory
        self.ocWin.gBox.addWidget(QtWidgets.QLabel("Data Directory",font=self.ocWin.setsFont1),1,1)
        self.ocWin.dataInp=QtWidgets.QLineEdit(self.exDataPath,font=self.ocWin.setsFont2)
        self.ocWin.dataBtn=QtWidgets.QPushButton("Browse...",font=self.ocWin.setsFont2)
        self.ocWin.dataBtn.setFixedSize(110,32)
        self.ocWin.dataBtn.clicked.connect(self.ocWinData)
        self.ocWin.gBox.addWidget(self.ocWin.dataInp,2,1)
        self.ocWin.gBox.addWidget(self.ocWin.dataBtn,2,2)
        # then for the catalogue
        self.ocWin.gBox.addWidget(QtWidgets.QLabel("Catalogue File",font=self.ocWin.setsFont1),3,1)
        self.ocWin.catInp=QtWidgets.QLineEdit(self.exCatFile,font=self.ocWin.setsFont2)
        self.ocWin.catBtn=QtWidgets.QPushButton("Browse...",font=self.ocWin.setsFont2)
        self.ocWin.catBtn.setFixedSize(110,32)
        self.ocWin.catBtn.clicked.connect(self.ocWinCat)
        self.ocWin.gBox.addWidget(self.ocWin.catInp,4,1)
        self.ocWin.gBox.addWidget(self.ocWin.catBtn,4,2)
        # last for the database
        self.ocWin.gBox.addWidget(QtWidgets.QLabel("Database File",font=self.ocWin.setsFont1),5,1)
        self.ocWin.dbInp=QtWidgets.QLineEdit(self.exDbPath,font=self.ocWin.setsFont2)
        self.ocWin.dbBtn=QtWidgets.QPushButton("Browse...",font=self.ocWin.setsFont2)
        self.ocWin.dbBtn.setFixedSize(110,32)
        self.ocWin.dbBtn.clicked.connect(self.ocWinDb)
        self.ocWin.gBox.addWidget(self.ocWin.dbInp,6,1)
        self.ocWin.gBox.addWidget(self.ocWin.dbBtn,6,2)
        # Final touches
        self.ocWin.diagBox=QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok |
                                                      QtWidgets.QDialogButtonBox.Cancel)
        self.ocWin.diagBox.accepted.connect(self.ocWinOk)
        self.ocWin.diagBox.rejected.connect(self.ocWinCancel)
        self.ocWin.gBox.addWidget(self.ocWin.diagBox,7,1)
        self.ocWin.show()

    def ocWinData(self):
        """Prompt for data path."""
        dataPath=str(QtWidgets.QFileDialog.getExistingDirectory(self,"Select the master data path:",self.exDataPath))
        if len(dataPath) == 0:
            return
        else:
            self.dataPath=dataPath
            self.ocWin.dataInp.setText(self.dataPath)
            logging.debug("Selected data path: %s" % self.dataPath)

    def ocWinCat(self):
        """Prompt for catalogue path. """
        catFile,_=QtWidgets.QFileDialog.getOpenFileName(self,"Select the catalogue file:",self.exCatFile,
            "QuakeML (*.xml);;Catalogue (*.cat);;Text (*.txt) ;;All files (*)")
        if len(catFile) == 0:
            return
        else:
            self.catFile=catFile
            self.ocWin.catInp.setText(self.catFile)
            logging.debug("Selected catalogue file: %s" % self.catFile)

    def ocWinDb(self):
        """Prompt for db path/"""
        if not os.path.exists(self.exDbPath) or self.exCatFile != self.catFile:
            if not os.path.isdir(self.catFile):
                self.exDbPath=os.path.splitext(self.catFile)[0]+".db3"
            else:
                self.exDbPath=os.getcwd()
        dbPath,_=QtWidgets.QFileDialog.getSaveFileName(self,"Select the database file:",self.exDbPath,
            "Database v3 (*.db3);;Database (*.db);;SQLITE3 (*.sqlite3) ;;All files (*)")
        if len(dbPath) == 0:
            return
        else:
            self.dbPath=dbPath
            self.ocWin.dbInp.setText(self.dbPath)
            logging.debug("Selected database file: %s" % self.dbPath)



    def ocWinOk(self):
        """Set the new paths and close."""
        self.dataPath=self.ocWin.dataInp.text()
        self.catFile=self.ocWin.catInp.text()
        self.dbPath=self.ocWin.dbInp.text() 
        self.ocWin.close()
        if not os.path.exists(self.dataPath):
            genTxt="Input Does Not Exist"
            infTxt="Could not locate data directory!"
            winTit="Open Catalogue"
            self.warnMsgBox(genTxt,infTxt,winTit)
            return
        if not os.path.exists(self.catFile):
            genTxt="Input Does Not Exist"
            infTxt="Could not locate catalogue!"
            winTit="Open Catalogue"
            self.warnMsgBox(genTxt,infTxt,winTit)
            return
        if os.path.isdir(self.dbPath):  
            genTxt="Input Does Not Exist"
            infTxt="Invalid database path given!"
            winTit="Open Catalogue"
            self.warnMsgBox(genTxt,infTxt,winTit)
            return   
        self.loadCatalogue=True
        try:
            # write new paths?
            if self.exDataPath != self.dataPath or self.exCatFile != self.catFile or self.exDbPath != self.dbPath:
                self.pathParser(mode="write",catFile=self.catFile,dataPath=self.dataPath,dbPath=self.dbPath)
            self.openCat()
        except:
            logging.exception("Could not open catalouge!")
            try:
                self.pDial.close()
            except:
                pass
            tbStack=traceback.format_exc()
            genTxt="Opening Catalogue Failure"
            infTxt="Opening catalogue exited with exception:\n %s" % tbStack
            winTit="Open Catalogue"
            self.warnMsgBox(genTxt,infTxt,winTit)  

    def ocWinCancel(self):
        """Close the Open Catalogue window."""
        self.ocWin.close()

    ## Fully automatic Cluster analysis Window
    def caMultWindow(self):
        """Open window for CA Mult analysis."""
        stCodes=sorted(self.inventory)
        self.cmWin=QtWidgets.QDialog(self)
        self.cmWin.setModal(True)
        self.cmWin.setWindowTitle("Automated Cluster Analysis Preferences")
        self.cmWin.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.cmWin.setsFont1=QtGui.QFont("Calibri",14,QtGui.QFont.Bold)
        self.cmWin.setsFont2=QtGui.QFont("Calibri",11,QtGui.QFont.Normal)
        self.cmWin.vBox=QtWidgets.QVBoxLayout(self.cmWin)
        ## add widgets ##
        # event selection
        self.cmStLbl=QtWidgets.QLabel("Selected stations: * ",font=self.cmWin.setsFont1)
        self.cmStList=QtWidgets.QListWidget()
        self.cmStList.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.cmStList.addItems(stCodes)
        self.cmStList.itemSelectionChanged.connect(self.cmListGet)
        self.cmWin.vBox.addWidget(self.cmStLbl)
        self.cmWin.vBox.addWidget(self.cmStList)
        ## other options
        self.cmWinDum=QtWidgets.QWidget()
        self.cmWinDum.gBox=QtWidgets.QGridLayout()
        self.cmWinDum.setLayout(self.cmWinDum.gBox)
        # maximum angle of incidence
        self.cmMaxAinLbl=QtWidgets.QLabel("Maximum angle of incidence ("+u"\u00b0"+")",font=self.cmWin.setsFont2)
        self.cmMaxAinInp=QtWidgets.QLineEdit("{:.1f}".format(self.maxAin),font=self.cmWin.setsFont2)
        self.cmMaxAinInp.setValidator(self.validateFloat)
        self.cmWinDum.gBox.addWidget(self.cmMaxAinLbl,1,1)
        self.cmWinDum.gBox.addWidget(self.cmMaxAinInp,1,2)
        # filters
        self.cmF1Lbl=QtWidgets.QLabel("Filter lower bound (Hz)",font=self.cmWin.setsFont2)
        self.cmF1Inp=QtWidgets.QLineEdit("{:.2f}".format(self.freqmin),font=self.cmWin.setsFont2)
        self.cmF2Lbl=QtWidgets.QLabel("Filter upper bound (Hz)",font=self.cmWin.setsFont2)
        self.cmF2Inp=QtWidgets.QLineEdit("{:.2f}".format(self.freqmax),font=self.cmWin.setsFont2)        
        self.cmWinDum.gBox.addWidget(self.cmF1Lbl,2,1)
        self.cmWinDum.gBox.addWidget(self.cmF1Inp,2,2)
        self.cmWinDum.gBox.addWidget(self.cmF2Inp,3,1)
        self.cmWinDum.gBox.addWidget(self.cmF2Inp,3,2)
        # selected method
        self.cmMetLbl=QtWidgets.QLabel("Method",font=self.cmWin.setsFont2)
        self.cmMetInp=QtWidgets.QComboBox()
        accMethods=sorted(("RC", "EV", "ME"))
        for itm in accMethods: self.cmMetInp.addItem(itm)
        self.cmWinDum.gBox.addWidget(self.cmMetLbl,4,1)
        self.cmWinDum.gBox.addWidget(self.cmMetInp,4,2)
        # checkbox for skipping pairs
        self.cmWinDum.skipFlag=QtWidgets.QCheckBox('Skip existing pairs',font=self.cmWin.setsFont2)        
        self.cmWinDum.skipFlag.setChecked(True)
        self.cmWinDum.gBox.addWidget(self.cmWinDum.skipFlag,5,1)
        # checkbox for skipping failed pairs
        self.cmWinDum.skipFailedFlag=QtWidgets.QCheckBox('Skip previously failed pairs',font=self.cmWin.setsFont2)        
        self.cmWinDum.skipFailedFlag.setChecked(True)
        self.cmWinDum.gBox.addWidget(self.cmWinDum.skipFailedFlag,5,2)
        #
        self.cmWin.vBox.addWidget(self.cmWinDum)
        #
        self.cmWin.diagBox=QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|
                                                      QtWidgets.QDialogButtonBox.Cancel)
        self.cmWin.diagBox.accepted.connect(self.cmWinOk)
        self.cmWin.diagBox.rejected.connect(self.cmWinCancel)
        self.cmWin.vBox.addWidget(self.cmWin.diagBox)

        self.cmWin.show()

    def cmWinOk(self):
        """Finalize CA Mult selection."""
        selectedStations=self.cmListGet(export=True)
        winTitle="Automated CA Confirmation"
        infText="The process will now begin for %i events and %i stations (potentially %i pairs). Are you sure?" \
                % (len(self.evsDict),len(selectedStations),len(self.evsDict)*len(selectedStations))
        rep=self.replyMsgBox(infText,winTitle)
        if rep == QtWidgets.QMessageBox.No:
            return
        filterBounds=(float(self.cmF1Inp.text()),float(self.cmF2Inp.text()))
        if not any(filterBounds): filterBounds=(np.nan,np.nan)
        maxAin=self.cmMaxAinInp.text()
        maxAin=np.inf if maxAin.upper() == "INF" else float(maxAin)
        method=self.cmMetInp.currentText()
        skipPairs=self.cmWinDum.skipFlag.isChecked()
        skipFailed=self.cmWinDum.skipFailedFlag.isChecked()
        self.cmWin.close()
        self.applyCCA(method,selectedStations,maxAin,filterBounds,skipPairs,skipFailed)
        

    def cmWinCancel(self):
        """Close CA Mult window."""
        self.cmWin.close()

    def cmListGet(self,export=False):
        """
        Get selected items and update labels

        :type export: bool, optional
        :param export: select whether selected items will be 
            exported as list. Defaults to False.
        """
        # update stations
        stations=self.cmStList.selectedItems()
        self.cmStLbl.setText(str("Selected stations: %i" % len(stations)).replace("0","*"))
        if export:
            if len(stations) == 0:
                return sorted([self.cmStList.item(idx).text() for idx in range(self.cmStList.count())])
            else:
                return [x.text() for x in stations]

    ## Preferences window related
    def prefWindow(self):
        """ Open up and set the preferences window."""
        self.pfWin=QtWidgets.QDialog(self)
        self.pfWin.setModal(True)
        self.pfWin.setWindowTitle("Preferences")
        self.pfWin.setWindowIcon(QtGui.QIcon(self.appIcon))
        #self.pfWin.resize(200,500)      
        ## set tabs ##
        self.pfTabs=QtWidgets.QTabWidget()
        # set default fonts
        self.pfWin.setsFont1=QtGui.QFont("Calibri",14,QtGui.QFont.Bold)
        self.pfWin.setsFont2=QtGui.QFont("Calibri",11,QtGui.QFont.Normal)
        # general
        self.gnTab=QtWidgets.QWidget()
        self.gnTab.gBox=QtWidgets.QGridLayout()
        self.gnTab.setLayout(self.gnTab.gBox)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("General",font=self.pfWin.setsFont1),1,1)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Clean Logs (days)",font=self.pfWin.setsFont2),2,1)
        self.gnTab.cleanLogs=QtWidgets.QLineEdit(str(self.generalCNF.cleanLogs),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)        
        self.gnTab.cleanLogs.setValidator(self.validateFloat)
        self.gnTab.gBox.addWidget(self.gnTab.cleanLogs,2,2)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Waveforms (appropriate values in s)",font=self.pfWin.setsFont1),3,1)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Offset for matching event folders to origin times",font=self.pfWin.setsFont2),4,1)
        self.gnTab.matching=QtWidgets.QLineEdit(str(self.generalCNF.matching),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)        
        self.gnTab.matching.setValidator(self.validateFloat)
        self.gnTab.gBox.addWidget(self.gnTab.matching,4,2)
        self.gnTab.trimFlag=QtWidgets.QCheckBox('Trim waveforms',font=self.pfWin.setsFont2)        
        self.gnTab.trimFlag.setChecked(self.generalCNF.trimFlag)
        self.gnTab.gBox.addWidget(self.gnTab.trimFlag,5,1)        
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Start time from S-arrival",font=self.pfWin.setsFont2),6,1)
        self.gnTab.trimStart=QtWidgets.QLineEdit(str(self.generalCNF.trimStart),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)        
        self.gnTab.trimStart.setValidator(self.validateFloat)
        self.gnTab.gBox.addWidget(self.gnTab.trimStart,6,2)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("End time from S-arrival",font=self.pfWin.setsFont2),7,1)
        self.gnTab.trimEnd=QtWidgets.QLineEdit(str(self.generalCNF.trimEnd),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)        
        self.gnTab.trimEnd.setValidator(self.validateFloat)
        self.gnTab.gBox.addWidget(self.gnTab.trimEnd,7,2)        
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Channel Code Preference",font=self.pfWin.setsFont2),8,1)
        self.gnTab.chanPref=QtWidgets.QLineEdit(str(self.generalCNF.chanPref),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.gnTab.gBox.addWidget(self.gnTab.chanPref,8,2)
        self.gnTab.orientFlag=QtWidgets.QCheckBox('Instrument Orientation Correction',font=self.pfWin.setsFont2)
        self.gnTab.gBox.addWidget(self.gnTab.orientFlag,9,1)
        self.gnTab.orientFlag.setChecked(self.generalCNF.orientFlag)        
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("SNR (values in s)",font=self.pfWin.setsFont1),10,1)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Start of noise window",font=self.pfWin.setsFont2),11,1)
        self.gnTab.snrStart=QtWidgets.QLineEdit(str(self.generalCNF.snrStart),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)        
        self.gnTab.snrStart.setValidator(self.validateFloat)
        self.gnTab.gBox.addWidget(self.gnTab.snrStart,11,2)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("End of signal window",font=self.pfWin.setsFont2),12,1)
        self.gnTab.snrEnd=QtWidgets.QLineEdit(str(self.generalCNF.snrEnd),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)        
        self.gnTab.snrEnd.setValidator(self.validateFloat)
        self.gnTab.gBox.addWidget(self.gnTab.snrEnd,12,2)  
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Splitting",font=self.pfWin.setsFont1),13,1)
        self.gnTab.gBox.addWidget(QtWidgets.QLabel("Maximum accepted t<sub>d</sub> (ms)",font=self.pfWin.setsFont2),14,1)
        self.gnTab.maxTd=QtWidgets.QLineEdit(str(self.generalCNF.maxTd),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)        
        self.gnTab.maxTd.setValidator(self.validateFloat)
        self.gnTab.gBox.addWidget(self.gnTab.maxTd,15,2)
        self.gnTab.resetBtn=QtWidgets.QPushButton("Reset Settings",font=self.pfWin.setsFont2)
        self.gnTab.resetBtn.clicked.connect(self.cleanAll)
        self.gnTab.resetBtn.setFixedSize(110,32)
        self.gnTab.gBox.addWidget(self.gnTab.resetBtn,16,1)
        # grading
        self.grTab=QtWidgets.QWidget()
        self.grTab.setMaximumSize(1024,350)
        self.grTab.gBox=QtWidgets.QGridLayout()
        self.grTab.setLayout(self.grTab.gBox)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("Thresholds",font=self.pfWin.setsFont1),1,1)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("Polarization Difference for Nulls (<sup>o</sup>)",font=self.pfWin.setsFont2),2,1)
        self.grTab.polOff=QtWidgets.QLineEdit(str(self.gradeCNF.polOff),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.polOff.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.polOff,2,2)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("SNR<sub>min</sub>",font=self.pfWin.setsFont2),3,1)
        self.grTab.snr_bound=QtWidgets.QLineEdit(str(self.gradeCNF.snr_bound),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.snr_bound.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.snr_bound,3,2)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("δφ<sub>max</sub>",font=self.pfWin.setsFont2),4,1)
        self.grTab.error_phi=QtWidgets.QLineEdit(str(self.gradeCNF.error_bounds[0]),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.error_phi.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.error_phi,4,2)        
        self.grTab.gBox.addWidget(QtWidgets.QLabel("δt<sub>d</sub><sub>max</sub>",font=self.pfWin.setsFont2),5,1)
        self.grTab.error_td=QtWidgets.QLineEdit(str(self.gradeCNF.error_bounds[0]),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.error_td.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.error_td,5,2)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("Minimum Correlation Coefficient (FxS)",font=self.pfWin.setsFont2),6,1) 
        self.grTab.CC_FS_bound=QtWidgets.QLineEdit(str(self.gradeCNF.CC_FS_bound),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.CC_FS_bound.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.CC_FS_bound,6,2)        
        self.grTab.gBox.addWidget(QtWidgets.QLabel("Minimum Correlation Coefficient (NxE)",font=self.pfWin.setsFont2),7,1) 
        self.grTab.CC_NE_bound=QtWidgets.QLineEdit(str(self.gradeCNF.CC_NE_bound),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.CC_NE_bound.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.CC_NE_bound,7,2)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("Grades (Specify only the UPPER closed bound for each grade)",font=self.lblFont),8,1)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("A",alignment=QtCore.Qt.AlignCenter,font=self.pfWin.setsFont2),9,1) 
        self.grTab.gradeA=QtWidgets.QLineEdit(str(self.gradeCNF.gradeDict["A"]),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.gradeA.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.gradeA,9,2)
        self.grTab.gBox.addWidget(QtWidgets.QLabel("B",alignment=QtCore.Qt.AlignCenter,font=self.pfWin.setsFont2),10,1) 
        self.grTab.gradeB=QtWidgets.QLineEdit(str(self.gradeCNF.gradeDict["B"]),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.gradeB.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.gradeB,10,2)  
        self.grTab.gBox.addWidget(QtWidgets.QLabel("C",alignment=QtCore.Qt.AlignCenter,font=self.pfWin.setsFont2),11,1)
        self.grTab.gradeC=QtWidgets.QLineEdit(str(self.gradeCNF.gradeDict["C"]),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.gradeC.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.gradeC,11,2) 
        self.grTab.gBox.addWidget(QtWidgets.QLabel("D",alignment=QtCore.Qt.AlignCenter,font=self.pfWin.setsFont2),12,1)
        self.grTab.gradeD=QtWidgets.QLineEdit(str(self.gradeCNF.gradeDict["D"]),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.grTab.gradeD.setValidator(self.validateFloat)
        self.grTab.gBox.addWidget(self.grTab.gradeD,12,2)                
        # cluster analysis
        self.caTab=QtWidgets.QWidget()
        self.caTab.setMaximumSize(800,1024)        
        self.caTab.gBox=QtWidgets.QGridLayout()
        self.caTab.setLayout(self.caTab.gBox)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Windows (all values in seconds, relative to S arrival)",font=self.pfWin.setsFont1),1,1)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Minimum window start point (T<sub>beg1</sub>)",font=self.pfWin.setsFont2),2,1)
        self.caTab.Tbeg1=QtWidgets.QLineEdit(str(self.caCNF.Tbeg1),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.Tbeg1.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.Tbeg1,2,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Step for windows start points (DT<sub>beg</sub>)",font=self.pfWin.setsFont2),3,1)
        self.caTab.DTbeg=QtWidgets.QLineEdit(str(self.caCNF.DTbeg),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.DTbeg.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.DTbeg,3,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Step for windows end points (DT<sub>end</sub>)",font=self.pfWin.setsFont2),4,1)
        self.caTab.DTend=QtWidgets.QLineEdit(str(self.caCNF.DTend),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.DTend.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.DTend,4,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Minimum window end point (T<sub>end0</sub>)",font=self.pfWin.setsFont2),5,1)
        self.caTab.Tend0=QtWidgets.QLineEdit(str(self.caCNF.Tend0),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.Tend0.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.Tend0,5,2) 
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Windows Bounds",font=self.pfWin.setsFont1),6,1)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Maximum t<sub>s</sub>-t<sub>p</sub> (s)",font=self.pfWin.setsFont2),7,1)
        self.caTab.tptsMax=QtWidgets.QLineEdit(str(self.caCNF.tptsMax),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.tptsMax.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.tptsMax,7,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Time from S-arrival for period determination (s)",font=self.pfWin.setsFont2),8,1)
        self.caTab.specWindow=QtWidgets.QLineEdit(str(self.caCNF.specWindow),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.specWindow.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.specWindow,8,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Minimum S-wave period (s)",font=self.pfWin.setsFont2),9,1)
        self.caTab.minPeriod=QtWidgets.QLineEdit(str(self.caCNF.minPeriod),alignment=QtCore.Qt.AlignCenter,font=self.pfWin.setsFont2)
        self.caTab.minPeriod.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.minPeriod,9,2)                
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Maximum S-wave period (s)",font=self.pfWin.setsFont2),10,1)
        self.caTab.maxPeriod=QtWidgets.QLineEdit(str(self.caCNF.maxPeriod),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.maxPeriod.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.maxPeriod,10,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Period factor",font=self.pfWin.setsFont2),11,1)
        self.caTab.multPeriod=QtWidgets.QLineEdit(str(self.caCNF.multPeriod),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.multPeriod.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.multPeriod,11,2)        
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Clustering",font=self.pfWin.setsFont1),12,1)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("C<sub>critical</sub>",font=self.pfWin.setsFont2),13,1)
        self.caTab.ccrit=QtWidgets.QLineEdit(str(self.caCNF.ccrit),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.ccrit.setValidator(self.validateFloat)
        self.caTab.gBox.addWidget(self.caTab.ccrit,13,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Maximum number of clusters (k<sub>max</sub>)",font=self.pfWin.setsFont2),14,1)
        self.caTab.kmax=QtWidgets.QLineEdit(str(self.caCNF.kmax),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.kmax.setValidator(self.validateInt)
        self.caTab.gBox.addWidget(self.caTab.kmax,14,2)
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Minimum number of points per cluster",font=self.pfWin.setsFont2),15,1)
        self.caTab.Ncmin=QtWidgets.QLineEdit(str(self.caCNF.Ncmin),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.caTab.Ncmin.setValidator(self.validateInt)
        self.caTab.gBox.addWidget(self.caTab.Ncmin,15,2)  
        self.caTab.gBox.addWidget(QtWidgets.QLabel("Linkage Criterion",font=self.pfWin.setsFont2),16,1)
        self.caTab.linkage=QtWidgets.QComboBox()
        accLinks=sorted(("ward", "complete", "average","single"))
        for itm in accLinks: self.caTab.linkage.addItem(itm)
        idx=accLinks.index(self.caCNF.linkage)
        self.caTab.linkage.setCurrentIndex(idx)
        self.caTab.gBox.addWidget(self.caTab.linkage,16,2)
        # TauP
        self.tauTab=QtWidgets.QWidget()
        self.tauTab.setMaximumSize(1024,200)
        self.tauTab.gBox=QtWidgets.QGridLayout()
        self.tauTab.setLayout(self.tauTab.gBox)
        self.tauTab.gBox.addWidget(QtWidgets.QLabel("Model (specify default model name if no file is provided)",font=self.pfWin.setsFont1),1,1)
        self.tauTab.model=QtWidgets.QLineEdit(self.tpCNF.model,font=self.pfWin.setsFont2)
        self.tauTab.modelBtn=QtWidgets.QPushButton("Browse...",font=self.pfWin.setsFont2)
        self.tauTab.modelBtn.clicked.connect(self.pfWinGetModel)
        self.tauTab.gBox.addWidget(self.tauTab.model,2,1)
        self.tauTab.gBox.addWidget(self.tauTab.modelBtn,2,2)
        self.tauTab.gBox.addWidget(QtWidgets.QLabel("Station Information File (requires restart to take effect)",font=self.pfWin.setsFont1),3,1)        
        self.tauTab.stations=QtWidgets.QLineEdit(self.tpCNF.stations,font=self.pfWin.setsFont2)
        self.tauTab.stationsBtn=QtWidgets.QPushButton("Browse...",font=self.pfWin.setsFont2)
        self.tauTab.stationsBtn.clicked.connect(self.pfWinGetStations)
        self.tauTab.gBox.addWidget(self.tauTab.stations,4,1)
        self.tauTab.gBox.addWidget(self.tauTab.stationsBtn,4,2)
        self.tauTab.ainFlag=QtWidgets.QCheckBox('Calculate incidence',font=self.pfWin.setsFont2)
        self.tauTab.gBox.addWidget(self.tauTab.ainFlag,5,1)
        self.tauTab.ainFlag.setChecked(self.tpCNF.ainFlag)
        # AR-AIC
        self.arTab=QtWidgets.QWidget()
        self.arTab.gBox=QtWidgets.QGridLayout()
        self.arTab.setLayout(self.arTab.gBox)
        self.arTab.gBox.addWidget(QtWidgets.QLabel("Lower bandpass window (Hz)",font=self.pfWin.setsFont2),1,1)
        self.arTab.arFMIN=QtWidgets.QLineEdit(str(self.pkCNF.arFMIN),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arFMIN.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arFMIN,1,2)
        self.arTab.gBox.addWidget(QtWidgets.QLabel("Upper bandpass window (Hz)",font=self.pfWin.setsFont2),2,1)
        self.arTab.arFMAX=QtWidgets.QLineEdit(str(self.pkCNF.arFMAX),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arFMAX.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arFMAX,2,2)   
        self.arTab.gBox.addWidget(QtWidgets.QLabel("LTA Length for P (s)",font=self.pfWin.setsFont2),3,1)
        self.arTab.arLTAP=QtWidgets.QLineEdit(str(self.pkCNF.arLTAP),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arLTAP.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arLTAP,3,2)   
        self.arTab.gBox.addWidget(QtWidgets.QLabel("STA Length for P (s)",font=self.pfWin.setsFont2),4,1)
        self.arTab.arSTAP=QtWidgets.QLineEdit(str(self.pkCNF.arLTAP),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arSTAP.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arSTAP,4,2)  
        self.arTab.gBox.addWidget(QtWidgets.QLabel("LTA Length for S (s)",font=self.pfWin.setsFont2),5,1)
        self.arTab.arLTAS=QtWidgets.QLineEdit(str(self.pkCNF.arLTAS),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arLTAS.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arLTAS,5,2)   
        self.arTab.gBox.addWidget(QtWidgets.QLabel("STA Length for S (s)",font=self.pfWin.setsFont2),6,1)
        self.arTab.arSTAS=QtWidgets.QLineEdit(str(self.pkCNF.arSTAS),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arSTAS.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arSTAS,6,2)          
        self.arTab.gBox.addWidget(QtWidgets.QLabel("Number of AR coefficients for P",font=self.pfWin.setsFont2),7,1)
        self.arTab.arMP=QtWidgets.QLineEdit(str(self.pkCNF.arMP),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arMP.setValidator(self.validateInt)
        self.arTab.gBox.addWidget(self.arTab.arMP,7,2)   
        self.arTab.gBox.addWidget(QtWidgets.QLabel("Number of AR coefficients for S",font=self.pfWin.setsFont2),8,1)
        self.arTab.arMS=QtWidgets.QLineEdit(str(self.pkCNF.arMS),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arMS.setValidator(self.validateInt)
        self.arTab.gBox.addWidget(self.arTab.arMS,8,2)  
        self.arTab.gBox.addWidget(QtWidgets.QLabel("Length of variance window for P",font=self.pfWin.setsFont2),9,1)
        self.arTab.arLP=QtWidgets.QLineEdit(str(self.pkCNF.arLP),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arLP.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arLP,9,2)   
        self.arTab.gBox.addWidget(QtWidgets.QLabel("Length of variance window for S",font=self.pfWin.setsFont2),10,1)
        self.arTab.arLS=QtWidgets.QLineEdit(str(self.pkCNF.arLS),font=self.pfWin.setsFont2,alignment=QtCore.Qt.AlignCenter)
        self.arTab.arLS.setValidator(self.validateFloat)
        self.arTab.gBox.addWidget(self.arTab.arLS,10,2)
        # add tabs
        self.pfTabs.addTab(self.gnTab,"General")
        self.pfTabs.addTab(self.grTab,"Grading")
        self.pfTabs.addTab(self.caTab,"Cluster Analysis")
        self.pfTabs.addTab(self.tauTab,"TauP")
        self.pfTabs.addTab(self.arTab,"AR-AIC")        
        # finalize
        self.pfWin.vBox=QtWidgets.QVBoxLayout(self.pfWin)
        self.pfWin.vBox.addWidget(self.pfTabs)
        self.pfWin.diagBox=QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|
                                                      QtWidgets.QDialogButtonBox.Apply|
                                                      QtWidgets.QDialogButtonBox.Cancel)
        self.pfWin.diagBox.accepted.connect(self.pfWinOk)
        self.pfWin.diagBox.clicked.connect(self.pfWinSave)
        self.pfWin.diagBox.rejected.connect(self.pfWinCancel)
        self.pfWin.vBox.addWidget(self.pfWin.diagBox)
        #
        self.pfWin.show()

    def pfWinGetModel(self):
        """Get a new location for the model file."""
        modFile,_=QtWidgets.QFileDialog.getOpenFileName(
            self,"Select the model file",
            str(self.tauTab.model.text()),
            "Named Discontinuities (*.nd);;All files (*)"
                                                        )
        if len(modFile) == 0:
            return       
        self.tauTab.model.setText(modFile)

    def pfWinGetStations(self):
        """Get a new location for the model file."""
        staFile,_=QtWidgets.QFileDialog.getOpenFileName(
            self,"Select the stations file",
            str(self.tauTab.stations.text()),
            "StationXML (*.xml);;Text(*.txt);;All files (*)"
                                                        )
        if len(staFile) == 0:
            return       
        self.tauTab.stations.setText(staFile)

    def pfWinCancel(self):
        """Cancel Preferences."""
        self.pfWin.close()

    def pfWinSave(self):
        """Save preferences to files."""
        # general
        self.generalCNF.cleanLogs=float(self.gnTab.cleanLogs.text())
        self.generalCNF.matching=float(self.gnTab.matching.text())
        self.generalCNF.chanPref=str(self.gnTab.chanPref.text())
        self.generalCNF.orientFlag=self.gnTab.orientFlag.isChecked()
        self.generalCNF.trimFlag=self.gnTab.trimFlag.isChecked()
        self.generalCNF.trimStart=float(self.gnTab.trimStart.text())
        self.generalCNF.trimEnd=float(self.gnTab.trimEnd.text())
        self.generalCNF.snrStart=float(self.gnTab.snrStart.text())
        self.generalCNF.snrEnd=float(self.gnTab.snrEnd.text())
        self.generalCNF.maxTd=float(self.gnTab.maxTd.text())
        self.generalCNF.write(WORKDIR+os.sep+"etc%soptions%sgeneral.cnf"%(os.sep,os.sep))
        # grading 
        self.gradeCNF.polOff=float(self.grTab.polOff.text())
        self.gradeCNF.snr_bound=float(self.grTab.snr_bound.text())
        self.gradeCNF.error_bounds[0]=float(self.grTab.error_phi.text())
        self.gradeCNF.error_bounds[1]=float(self.grTab.error_td.text())
        self.gradeCNF.CC_FS_bound=float(self.grTab.CC_FS_bound.text())
        self.gradeCNF.CC_NE_bound=float(self.grTab.CC_NE_bound.text())
        self.gradeCNF.gradeDict["A"]=float(self.grTab.gradeA.text())
        self.gradeCNF.gradeDict["B"]=float(self.grTab.gradeB.text())
        self.gradeCNF.gradeDict["C"]=float(self.grTab.gradeC.text())
        self.gradeCNF.gradeDict["D"]=float(self.grTab.gradeD.text())
        self.gradeCNF.write(WORKDIR+os.sep+"etc%soptions%sgrading.cnf"%(os.sep,os.sep))
        # cluster analysis
        self.caCNF.Tbeg1=float(self.caTab.Tbeg1.text())
        self.caCNF.DTbeg=float(self.caTab.DTbeg.text())
        self.caCNF.DTend=float(self.caTab.DTend.text())
        self.caCNF.Tend0=float(self.caTab.Tend0.text())
        self.caCNF.tptsMax=float(self.caTab.tptsMax.text())
        self.caCNF.specWindow=float(self.caTab.specWindow.text())
        self.caCNF.multPeriod=float(self.caTab.multPeriod.text())
        self.caCNF.minPeriod=float(self.caTab.minPeriod.text())
        self.caCNF.maxPeriod=float(self.caTab.maxPeriod.text())
        self.caCNF.ccrit=float(self.caTab.ccrit.text())
        self.caCNF.kmax=int(self.caTab.kmax.text())
        self.caCNF.Ncmin=int(self.caTab.Ncmin.text())
        self.caCNF.linkage=self.caTab.linkage.currentText()
        # small hack
        self.caCNF.maxTd=self.generalCNF.maxTd
        # write the file now
        self.caCNF.write(WORKDIR+os.sep+"etc%soptions%sclustering.cnf"%(os.sep,os.sep))
        # taup
        self.tpCNF.model=str(self.tauTab.model.text())
        self.tpCNF.stations=str(self.tauTab.stations.text())
        self.tpCNF.ainFlag=self.tauTab.ainFlag.isChecked()
        self.tpCNF.write(WORKDIR+os.sep+"etc%soptions%staup.cnf"%(os.sep,os.sep))
        # AR-AIC
        self.pkCNF.arFMIN=float(self.arTab.arFMIN.text())
        self.pkCNF.arFMAX=float(self.arTab.arFMAX.text())
        self.pkCNF.arLTAP=float(self.arTab.arLTAP.text())
        self.pkCNF.arSTAP=float(self.arTab.arSTAP.text())
        self.pkCNF.arLTAS=float(self.arTab.arLTAS.text())
        self.pkCNF.arSTAS=float(self.arTab.arSTAS.text())
        self.pkCNF.arMP=int(self.arTab.arMP.text())
        self.pkCNF.arMS=int(self.arTab.arMS.text())
        self.pkCNF.arLP=float(self.arTab.arLP.text())
        self.pkCNF.arLS=float(self.arTab.arLS.text())
        self.pkCNF.write(WORKDIR+os.sep+"etc%soptions%spicker.cnf"%(os.sep,os.sep))
        logging.info("Saved preferences")

    def pfWinOk(self):
        """Ok response."""
        self.pfWinSave()
        self.pfWin.close()

    ## GUI DB query related functions
    def dbWindowQuery(self):
        """Open a window that permits database querying."""
        # prepwork
        self.dbWin=QtWidgets.QDialog(self)
        self.dbWin.setModal(True)
        self.dbWin.setWindowTitle("Database Search")
        self.dbWin.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.dbWin.vBox=QtWidgets.QVBoxLayout(self.dbWin)
        ## add widgets ##
        # event selection
        self.evLbl=QtWidgets.QLabel("Selected events: * ")
        self.evList=QtWidgets.QListWidget()
        self.evList.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        _,vals=DB.retrieve(self.dbCur,"event","evid","*")
        origins=sorted(list(set([x[0] for x in vals])))
        self.evList.addItems(origins)
        self.evList.itemSelectionChanged.connect(self.dbListGet)
        self.dbWin.vBox.addWidget(self.evLbl)
        self.dbWin.vBox.addWidget(self.evList)
        # station selection
        self.stLbl=QtWidgets.QLabel("Selected stations: * ")
        self.stList=QtWidgets.QListWidget()
        self.stList.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        _,vals=DB.retrieve(self.dbCur,"station","station","*/*")
        stCodes=sorted(list(set([x[0] for x in vals])))
        self.stList.addItems(stCodes)
        self.stList.itemSelectionChanged.connect(self.dbListGet)
        self.dbWin.vBox.addWidget(self.stLbl)
        self.dbWin.vBox.addWidget(self.stList)
        # method selection
        self.meLbl=QtWidgets.QLabel("Selected methods: * ")
        self.meList=QtWidgets.QListWidget()
        self.meList.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        _,vals=DB.retrieve(self.dbCur,"method","method","*/*/*")
        meCodes=sorted(list(set([x[0] for x in vals])))
        self.meList.addItems(meCodes)
        self.meList.itemSelectionChanged.connect(self.dbListGet)
        self.dbWin.vBox.addWidget(self.meLbl)
        self.dbWin.vBox.addWidget(self.meList)
        # export/clear button
        self.clBtn=QtWidgets.QPushButton("Clear selection")
        self.clBtn.clicked.connect(self.dbListClear)
        self.dbWin.vBox.addWidget(self.clBtn)
        self.exBtn=QtWidgets.QPushButton("Export")
        self.exBtn.clicked.connect(self.dbWinFinish)
        self.dbWin.vBox.addWidget(self.exBtn)
        # fin
        self.dbWin.resize(300,600)
        self.dbWin.show()

    def dbListClear(self):
        """Clear lists from selected items."""
        for listwidget in (self.evList,self.stList,self.meList):
            for i in range(listwidget.count()):
                item=listwidget.item(i)
                item.setSelected(False)

    def dbListGet(self,export=False):
        """
        Get selected items and update labels

        :type export: bool, optional
        :param export: select whether selected items will be 
            exported as list. Defaults to False.

        """
        # update events
        events=self.evList.selectedItems()
        self.evLbl.setText(str("Selected events: %i" % len(events)).replace("0","*"))
        # update stations
        stations=self.stList.selectedItems()
        self.stLbl.setText(str("Selected stations: %i" % len(stations)).replace("0","*"))
        # update methods
        methods=self.meList.selectedItems()
        self.meLbl.setText(str("Selected methods: %i" % len(methods)).replace("0","*"))
        if export:
            return events,stations,methods

    def dbWinFinish(self):
        """Query in db and save in file."""
        selected=self.dbListGet(export=True)
        selected=[[x.text() for x in y] if len(y) > 0 else "*" for y in selected]
        splittingDict=DB.load(self.dbCur,evids=selected[0],stids=selected[1],meids=selected[2])
        self.dbWin.close()
        self.writeSplittingDict(splittingDict)

    ## Evens list window
    def eventsListWindow(self):
        """
        Open a window containing selectable events. In case the events
        doesn't exist? Maybe events from catalogue that do not correspond
        to an event in folder with warnings (and not selectable).
        """
        self.evWin=QtWidgets.QDialog(self)
        self.evWin.setModal(True)
        self.evWin.setWindowTitle("Events' List")
        self.evWin.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.evWin.vBox=QtWidgets.QVBoxLayout(self.evWin)
        self.evWin.setLayout(self.evWin.vBox)
        self.evWin.selectEv=QtWidgets.QTableWidget()
        self.evWin.selectEv.setColumnCount(4)
        self.evWin.selectEv.verticalHeader().setVisible(False)
        self.evWin.selectEv.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.evWin.vBox.addWidget(self.evWin.selectEv)
        allEvList=sorted(self.evsDict)
        actEvList=sorted(self.getActiveEvents())
        maxLen=len(allEvList)
        self.evWin.selectEv.setRowCount(maxLen)
        self.evWin.selectEv.setHorizontalHeaderLabels(["#","Event","Magnitude","Depth"])
        hdr=self.evWin.selectEv.horizontalHeader()
        for hi in range(self.evWin.selectEv.columnCount()): 
            hdr.setSectionResizeMode(hi, QtWidgets.QHeaderView.ResizeToContents)
        for j,i in enumerate(allEvList):
            j+=1
            minAin=min([self.statsDict[i][x]["AN"] for x in self.statsDict[i]])
            font=QtGui.QFont()
            if (i in actEvList) and (minAin <= self.maxAin):
                font.setWeight(QtGui.QFont.Bold)
            elif i not in actEvList:
                font.setStrikeOut(True)
            # add event code
            tItems=[]
            for ii,jj in enumerate((str(j),str(UTCDateTime(i)).split(".")[0],
                                             "{:.1f}".format(self.evsDict[i]["MAG"]),
                                             "{:.2f}".format(self.evsDict[i]["DEPTH"]))):
                if tools.isreal(jj):
                    itm=QTableWidgetNumItem(jj)
                elif tools.isdate(jj):
                    itm=QTableWidgetDateItem(jj)
                else:
                    itm=QtWidgets.QTableWidgetItem(jj)                
                itm.setTextAlignment(QtCore.Qt.AlignCenter)
                itm.setFont(font)
                self.evWin.selectEv.setItem(j-1,ii,itm)
        self.evWin.buttons=QtWidgets.QDialogButtonBox(
           QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
                                    )
        self.evWin.buttons.accepted.connect(self.runEvWin)
        self.evWin.buttons.rejected.connect(self.closeEvWin)
        self.evWin.vBox.addWidget(self.evWin.buttons)
        self.evScroll()
        # resize
        self.evWin.setMinimumSize(420,600)
        # final properties
        self.evWin.selectEv.setSortingEnabled(True)
        self.evWin.selectEv.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.evWin.show()

    def evScroll(self):
        """Scroll to current event."""
        try:
            lookup=str(UTCDateTime(self.activeEvent)).split(".")[0]
            item=self.evWin.selectEv.findItems(lookup,QtCore.Qt.MatchExactly)[0]
            self.evWin.selectEv.selectRow(item.row())
            self.evWin.selectEv.scrollToItem(item,QtWidgets.QAbstractItemView.PositionAtTop)
        except:
            logging.exception("Couldn't scroll to active event")

    def closeEvWin(self):
        """Close the events' list window."""
        self.evWin.close()

    def runEvWin(self):
        """Close events' list window and run."""
        tm=self.evWin.selectEv.item(self.evWin.selectEv.currentRow(),1)
        self.activeEvent=UTCDateTime(tm.text()).strftime("%Y-%m-%d-%H-%M-%S")
        self.evWin.close()
        evList=sorted(self.evsDict)
        self.evIdx=evList.index(self.activeEvent)-1
        self.nextEvent()

    def stationsListWindow(self):
        """
        Open a window containing selectable stations. In case the events
        doesn't exist? Maybe events from catalogue that do not correspond
        to an event in folder with warnings (and not selectable).
        """
        self.stWin=QtWidgets.QDialog(self)
        self.stWin.setModal(False)
        self.stWin.setWindowTitle("Station List")
        self.setWindowIcon(QtGui.QIcon(self.appIcon))
        self.stWin.vBox=QtWidgets.QVBoxLayout(self.stWin)
        self.stWin.setLayout(self.stWin.vBox)
        self.stWin.selectSt=QtWidgets.QTableWidget()
        self.stWin.selectSt.setColumnCount(4)
        self.stWin.selectSt.verticalHeader().setVisible(False)
        self.stWin.selectSt.setSelectionBehavior(QtWidgets.QTableView.SelectRows)        
        self.stWin.vBox.addWidget(self.stWin.selectSt)
        actStaList=self.getActiveStations()
        allStaList=sorted(
                        self.statsDict[self.activeEvent],
                        key=lambda x:self.statsDict[self.activeEvent][x]["AN"]
                      )
        maxLen=len(allStaList)
        self.stWin.selectSt.setRowCount(maxLen)
        self.stWin.selectSt.setHorizontalHeaderLabels(["#","Station","Incidence","Backazimuth"])
        hdr=self.stWin.selectSt.horizontalHeader()
        for hi in range(self.stWin.selectSt.columnCount()): 
            hdr.setSectionResizeMode(hi, QtWidgets.QHeaderView.ResizeToContents)
        for j,i in enumerate(allStaList):
            j+=1
            ain=self.statsDict[self.activeEvent][i]["AN"]
            font=QtGui.QFont()
            if (i in actStaList) and (ain <= self.maxAin):
                font.setWeight(QtGui.QFont.Bold)
            elif i not in actStaList:
                font.setStrikeOut(True)
            # add event code
            tItems=[]
            for ii,jj in enumerate((str(j),i,
                                      "{:.1f}".format(self.statsDict[self.activeEvent][i]["AN"]),
                                      "{:.2f}".format(self.statsDict[self.activeEvent][i]["BAZ"]))):
                if tools.isreal(jj):
                    itm=QTableWidgetNumItem(jj)
                else:
                    itm=QtWidgets.QTableWidgetItem(jj)
                itm.setTextAlignment(QtCore.Qt.AlignCenter)
                itm.setFont(font)
                if ii != 1: # index 1 is the station name
                    itm.setData(QtCore.Qt.EditRole, QtCore.QVariant(float(jj)))
                self.stWin.selectSt.setItem(j-1,ii,itm)       
        self.stWin.buttons=QtWidgets.QDialogButtonBox(
           QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
                                    )
        self.stWin.buttons.accepted.connect(self.runStWin)
        self.stWin.buttons.rejected.connect(self.closeStWin)
        self.stWin.vBox.addWidget(self.stWin.buttons)
        # resize
        self.stWin.setMinimumSize(355,500)
        # final properties
        self.stWin.selectSt.setSortingEnabled(True)
        self.stWin.selectSt.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.stScroll()
        self.stWin.show()

    def stScroll(self):
        """Scroll to active station."""
        try:
            item=self.stWin.selectSt.findItems(self.station,QtCore.Qt.MatchExactly)[0]
            self.stWin.selectSt.selectRow(item.row())
            self.stWin.selectSt.scrollToItem(item,QtWidgets.QAbstractItemView.PositionAtTop)
        except:
            logging.exception("Couldn't scroll to active event")

    def closeStWin(self):
        """Close the events' list window."""
        self.stWin.close()

    def runStWin(self):
        """Close events' list window and run."""
        tm=self.stWin.selectSt.item(self.stWin.selectSt.currentRow(),1)
        self.station=tm.text()
        self.stWin.close()
        staList=sorted(
                        self.statsDict[self.activeEvent],
                        key=lambda x:self.statsDict[self.activeEvent][x]["AN"]
                      )        
        self.stIdx=staList.index(self.station)-1
        self.nextStation(autoAinSkip=False)

    # other
    def disableV2HToggle(self):
        """Hack to indirectly disable this checkbox."""
        self.inpV2H.setChecked(self.V2H)

## class for threading Catalogue CA
class applyCAMult(QtCore.QThread):
    """
    Class for applying Catalogue CA. Admittedly, it's a bit messy, 
    might have to rewrite it at some point and try and make it 
    thread-safe (also enabling parallelism). Is it a good practice to use
    the active class that controls the GUI as an argument?

    Uses :class: `~PyQt5.QtCore.QThread`
 
    """
    # setup custom signals
    iterDone=QtCore.pyqtSignal(list)
    iterFail=QtCore.pyqtSignal(Exception)

    def __init__(self,pytheas,meth,selectedStations,maxAin,filterBounds,skipPairs,skipFailed,actEvList,parent=None):
        """
        Initilization

        :type pytheas: :class: `~pytheas.Pytheas`
        :param pytheas: the Pytheas' main class that contains all the required information and calls to
            run the Catalogue CA.
        :type meth: str
        :param meth: analysis method to use with CA (EV, ME or RC)
        :type selectedStations: tuple-like
        :param selectedStations: list of stations to be used in the analysis
        :type maxAin: float
        :param maxAin: the incidence angle corresponding to the shear-wave window
        :type filterBounds: tuple-like
        :param filterBounds: the minimum and maximum filter boundaries
        :type skipPairs: bool
        :param skipPairs: select whether to skip pairs that already exist in the database
        :type skipFailed: bool
        :param skipFailed: select whether to skip pairs that previously led to errors during CA
        :type actEvList: tuple-like
        :param actEvList: list of active events

        """
        QtCore.QThread.__init__(self,parent)
        self.pytheas=pytheas
        self.meth=meth
        self.selectedStations=selectedStations
        self.maxAin=maxAin
        self.filterBounds=filterBounds
        self.freqmin=filterBounds[0]
        self.freqmax=filterBounds[1]
        self.skipPairs=skipPairs
        self.skipFailed=skipFailed
        self.actEvList=actEvList
        self.excFlag=False
        self._isRunning=True

    def stop(self):
        """Set flag to stop iterating over windows."""
        self._isRunning=False
        try:
            self.pytheas.CAthread._isRunning=False
        except:
            pass
        self.wait()

    def run(self):
        """Run the full catalogue."""
        try:
            # get start time
            tic=UTCDateTime()
            self.tic=tic
            # setup the log for the CA
            logfileCA=WORKDIR+os.sep+"logs"+os.sep+"CA_%s.log" % tic.strftime("%Y%m%d_%H%M%S")
            self.logfileCA=logfileCA
            logging.info("Will save CA related information to %s" % logfileCA)
            with open(logfileCA,"w") as fid:
                fid.write("# Shear-wave Window: %.1f deg | Filter (Hz): %.1f - %.1f\n" % (self.maxAin,self.freqmin,self.freqmax))
                fid.write("# Process started @ %s\n" % tic)
            # pairs analyzed
            self.ipairs=0; self.ispairs=0
            maxI=len(self.actEvList)
            ## connecting to DB
            self.dbConn,self.dbCur=DB.open(self.pytheas.dbPath)
            # find all failed attempts
            if self.skipFailed:
                _,uids=DB.find(self.dbCur,"method","grade","F")
                failedPairs=["%s/%s/%s" % x for x in uids]
            # start iterating events
            for event in sorted(self.actEvList):
                if not self._isRunning:
                    raise UserCancelException("User cancelled the iterations.")
                i=sorted(self.actEvList).index(event)+1
                # get info first #
                self.pytheas.activeEvent=event
                logging.info("Iterating %i/%i events"%(i,maxI))
                self.pytheas.multStDict=stDict=self.pytheas.statsDict[event]
                self.pytheas.multEvDict=evDict=self.pytheas.evsDict[event]
                maxJ=len(stDict.keys())
                with open(logfileCA,"a") as fid:
                    fid.write("===== Starting event %s =====\n" % event)
                # start iterating stations
                for station in sorted(stDict.keys()):
                    ## checks ##
                    skip=False; state=False                    
                    # reset measurements first
                    self.pytheas.flagReset(general=False,splitting=True,indexing=False,titling=False,splitDict=False)
                    # keep going
                    if not self._isRunning:
                        raise UserCancelException("User cancelled the iterations.")
                    # was this a failed attempt before? skip?
                    if self.skipFailed:
                        uid="%s/%s/%s" % (event,station,self.meth)
                        if uid in failedPairs:
                            logging.info("%s has previously failed, skipping..." % uid)
                            continue
                    # has this combo already been analyzed? skip?
                    if self.skipPairs:
                        uid="%s/%s/%s" % (event,station,self.meth)
                        _,vals=DB.retrieve(self.dbCur,"method","method",uid)
                        if vals:
                            logging.info("%s exists in database" % uid)
                            continue
                    # station was not one of the selected ones. skip!
                    if station not in self.selectedStations:
                        logging.debug("%s was not in the selected stations, skipping..." % station)
                        continue
                    # get station's #
                    j=sorted(stDict.keys()).index(station)+1
                    # send a signal to main thread and update progress bar
                    self.iterDone.emit(["Processing event %s (%i/%i)\nProcessing station %s (%i/%i)" \
                            % (
                                event,i,maxI, 
                                station,j,maxJ
                              ),100*i/maxI])
                    self.pytheas.station=station
                    try:
                        logging.info("Iterating %i/%i stations"%(j,maxJ))
                        # calculate incidence angle?
                        if self.pytheas.tpCNF.ainFlag:
                            try:
                              taupdat=tools.getTheorArrivals( # TODO: display secpndary phases too
                                                       (evDict["ORIGIN"],evDict["LAT"],
                                                        evDict["LON"],evDict["DEPTH"]),
                                                       (self.pytheas.inventory[station]["latitude"],
                                                        self.pytheas.inventory[station]["longitude"]),
                                                        self.pytheas.tmodel
                                                            )
                              # grab measurements for the s arrival (not S)
                              for reg in taupdat.values():
                                if reg["phase"] == "s": 
                                    ain=reg["ain"]
                                    stDict[station]["AN"]=ain
                                    break # can there be more than 1 "s" arrivals??
                            except:
                                logging.debug("Could not calculate TauP ain for %s -%s" % (event,station))            
                        # is the pair outside the shear-wave window?
                        if float(stDict[station]['AN']) > self.pytheas.maxAin:
                            logging.warning("Angle of incidence greater than %.1f, skipping..." % self.pytheas.maxAin)
                            skip=True
                            state=" SKIP AIN\n"
                        # does the pair have an active S arrival?
                        if station not in self.pytheas.spickDict[event]:
                            logging.warning("No S pick for %s , skipping..."%station)     
                            skip=True
                            if not state:
                                state=" SKIP NO_S_PICKS\n"                        
                        # check if arrival is valid
                        self.pytheas.sArr=self.pytheas.spickDict[event][station]
                        if self.pytheas.sArr in [None,0,np.nan]:
                            logging.warning("No S pick for %s , skipping..."%station)     
                            skip=True
                            if not state:
                                state=" SKIP NO_S_PICKS\n"   
                        # are we skipping?                                          
                        if skip:
                            with open(logfileCA,"a") as fid:
                                fid.write(event.ljust(25)+station.rjust(6)+state)
                            # WRITE DUMMY DICT WITH F GRADE TO DB
                            continue                 
                        # get origin time
                        otime=evDict["ORIGIN"]
                        # find the event folder and fetch the waveforms
                        evFolder=self.pytheas.fullEvents[event]
                        # iterate through ZNE
                        try:
                            stationFile=evFolder+os.sep+"*."+station+".*"
                            tempStream,orCorrection=self.pytheas.getWave(stationFile,
                                                                         inventory=self.pytheas.xmlInventory,
                                                                         prefChannel=self.pytheas.generalCNF.chanPref,
                                                                         orientCorr=self.pytheas.generalCNF.orientFlag,
                                                                         showWarn=False)
                            if not tempStream:
                                raise StationException("Could not find station")
                            stream=tempStream
                            self.pytheas.stream=stream
                            self.pytheas.orCorrection=orCorrection
                        except:
                            logging.exception("Could not find %s for %s, while applying CA!" % (station,event))
                            with open(logfileCA,"a") as fid: 
                                state=" FAIL CHANNEL_NOT_FOUND\n"
                                fid.write(event.ljust(25)+station.rjust(6)+state)
                            continue
                        # get S pick
                        self.pytheas.sArr=self.pytheas.spickDict[event][station]
                        self.pytheas.sPick=self.pytheas.sArr-self.pytheas.stream[0].stats.starttime
                        #
                        self.ipairs+=1
                        # check if all timeseries have the same length
                        Wz=stream.select(component="Z")[0]
                        Wy=stream.select(component="N")[0]
                        Wx=stream.select(component="E")[0]
                        # First horizontals
                        Ax,Ay=tools.lengthcheck(Wx,Wy)
                        stream.select(component="E")[0].data=Ax 
                        stream.select(component="N")[0].data=Ay
                        # Now the vertical
                        Az,Ay=tools.lengthcheck(Wz,Wy)
                        stream.select(component="Z")[0].data=Az
                        stream.select(component="N")[0].data=Ay
                        # Horizontals second time 
                        Ax,Ay=tools.lengthcheck(Wx,Wy)
                        stream.select(component="E")[0].data=Ax 
                        stream.select(component="N")[0].data=Ay
                        ######################################
                        if abs(stream[1].stats.starttime - stream[2].stats.starttime) > stream[1].stats.delta:
                            logging.warning("Could not synchronize waveforms, skipping...")
                            with open(logfileCA,"a") as fid: 
                                state=" FAIL START_TIME_ERROR\n"
                                fid.write(event.ljust(25)+station.rjust(6)+state)
                            # add grade
                            self.pytheas.grade="F"            
                            self.pytheas.updateSplittingDict()
                            DB.addValues(self.dbConn,self.dbCur,self.pytheas.splittingDict)
                            continue
                        # TODO: make filter(s) defined in the TB config file
                        # get SNR
                        if not tools.isnone(self.pytheas.freqmin) or not tools.isnone(self.pytheas.freqmax):
                            stream.filter(type="bandpass",freqmin=self.pytheas.freqmin,freqmax=self.pytheas.freqmax,zerophase=True,corners=4)
                        SNR=self.pytheas.calcSNR(self.pytheas.sPick)
                        # deal with ain
                        if self.pytheas.actQT.isChecked():
                            logging.debug("Will rotate to LQT")
                            ain=stDict[station]['AN']
                        else:
                            logging.debug("Will rotate to ZRT")
                            ain=0.
                        self.pytheas.applyCA(self.meth[1:],stream,False,self.pytheas.freqmin,self.pytheas.freqmax,
                                     self.pytheas.sPick,stDict[station]['BAZ'],ain,False,True,
                                     eventCode=event)
                        while self.pytheas.CAthread.isRunning():
                            time.sleep(0.1) # stability hack
                            continue
                        if self.pytheas.CAthread.excFlag:
                            logging.error("Error while executing fully automatic clustering in pair: %s - %s" \
                                          % (event,station))
                            with open(logfileCA,"a") as fid:
                                try:
                                    _,_,exc_traceback=sys.exc_info()
                                    ss=traceback.extract_tb(exc_traceback)[0]
                                    state=" FAIL %s\n" % ss.line
                                except:
                                    state=" FAIL CA_THREAD_ERROR\n"
                                fid.write(event.ljust(25)+station.rjust(6)+state)
                            # add grade
                            self.pytheas.grade="F"
                            self.pytheas.updateSplittingDict()
                            DB.addValues(self.dbConn,self.dbCur,self.pytheas.splittingDict)
                            continue
                        phi=self.pytheas.CAthread.phi % 180
                        td=self.pytheas.CAthread.dt
                        optWindow=self.pytheas.CAthread.optWindow
                        filterRange=(self.pytheas.freqmin,self.pytheas.freqmax)
                        errors=(self.pytheas.CAthread.sphi,self.pytheas.CAthread.sdt)
                        _,_,cArray,pTest,dTest,\
                        pol,_,_,cEmax,\
                        CC_FS,CC_NE,Tvar,nContours=self.pytheas.CAthread.tempRes 
                        # grade the result
                        grade_score,grade=tools.autoGrading(
                            (phi,pol,self.pytheas.gradeCNF.polOff),
                            (errors[0],errors[1],SNR,CC_FS,CC_NE), # values
                            (self.pytheas.gradeCNF.error_bounds[0],
                            self.pytheas.gradeCNF.error_bounds[1],
                            self.pytheas.gradeCNF.snr_bound,
                            self.pytheas.gradeCNF.CC_FS_bound, 
                            self.pytheas.gradeCNF.CC_NE_bound),  # bounds
                            self.pytheas.gradeCNF.gradeDict
                                        )
                        # write results to splitting dict
                        self.pytheas.phi=phi; self.pytheas.td=td; self.pytheas.pol=pol
                        self.pytheas.nContours=nContours; self.pytheas.CC_FS=CC_FS; self.pytheas.CC_NE=CC_NE
                        self.pytheas.Tvar=Tvar; self.pytheas.errors=errors
                        self.pytheas.minPick,self.pytheas.maxPick=optWindow
                        self.pytheas.grade=grade; self.pytheas.grade_score=grade_score; self.pytheas.SNR=SNR
                        self.pytheas.pTest=pTest; self.pytheas.dTest=dTest
                        self.pytheas.cEmax=cEmax; self.pytheas.cArray=cArray
                        self.pytheas.initialClusters=self.pytheas.CAthread.initial; self.pytheas.initialClustersErrors=self.pytheas.CAthread.initial_err
                        self.pytheas.calinski=self.pytheas.CAthread.calinski
                        self.pytheas.clusters1=self.pytheas.CAthread.clusters1; self.pytheas.clusters2=self.pytheas.CAthread.clusters2                    
                        ## add values to database
                        self.pytheas.updateSplittingDict()
                        DB.addValues(self.dbConn,self.dbCur,self.pytheas.splittingDict)
                        self.ispairs+=1
                        with open(logfileCA,"a") as fid:
                            state=" DONE\n"
                            fid.write(event.ljust(25)+station.rjust(6)+state)
                    except:
                        logging.exception("Error in Cluster Analysis! Skipping %s - %s" % (event,station))
                        _,_,exc_traceback=sys.exc_info()
                        ss=traceback.extract_tb(exc_traceback)[0]
                        with open(logfileCA,"a") as fid:
                            state=" FAIL %s\n" % ss.line
                            fid.write(event.ljust(25)+station.rjust(6)+state)
                        # add grade
                        self.pytheas.grade="F"
                        self.pytheas.updateSplittingDict()
                        DB.addValues(self.dbConn,self.dbCur,self.pytheas.splittingDict)
                        continue
        except Exception as exc:
            self.excFlag=True
            self.iterFail.emit(exc)
        DB.close(self.dbConn)

## Subclassing QTableWidgetItem to enable sorting
class QTableWidgetNumItem(QtWidgets.QTableWidgetItem):
    """Set a different sorting method for number items."""
    def __lt__(self, other):
        return float(self.text()) < float(other.text())

class QTableWidgetDateItem(QtWidgets.QTableWidgetItem):
    """Set a different sorting method for date items."""
    def __lt__(self, other):
        return UTCDateTime(self.text()) < UTCDateTime(other.text())

## Custom Exceptions ##
class StationException(Exception):
    pass

class EventException(Exception):
    pass

class UserCancelException(Exception):
    pass

# grab silent exceptions
def my_excepthook(type,value,tback):
    sys.__excepthook__(type,value,tback)
sys.excepthook = my_excepthook

## LAUNCH
if __name__=='__main__':
    app=QtWidgets.QApplication(sys.argv)
    pytheasWin=Pytheas()
    pytheasWin.showMaximized()
    #pytheasWin.show() # don't show maximized. maybe add an option?
    sys.exit(app.exec_())