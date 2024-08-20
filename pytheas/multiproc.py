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

Authors: Spingos I. & Kaviris G. (c) 2019-2024
Special thanks to Millas C. for testing the software and providing valuable feedback from the 
very early stages of this endeavor!

For any issues, comments or suggestions please contact us at ispingos@geol.uoa.gr or through GitHub 
at https://www.github.com/ispingos/pytheas-splitting

"""
####################################################################
#                                                                  #
# This module offers parallelilzation support for                  #
# Catalogue Cluster Analysis.                                      #
#                                                                  #
####################################################################

import os, sys
from subprocess import call
import multiprocessing

core_count = int(sys.argv[1])

def call_cli(core_num):
	target = os.path.join(os.getcwd(), 'pytheas', 'cli.py')
	call(['python', target, ('%i' % core_num)])

if __name__ == '__main__':

	with multiprocessing.Pool() as pool:
		pool.map(call_cli, range(1, core_count + 1))
