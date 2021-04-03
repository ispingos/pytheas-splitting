#!/bin/sh

##############################
#
# version date: 01/12/2020
#
#############################

# upgrade pip
echo Upgrading pip, this might require su access..
sudo pip install --upgrade pip

# first upgrade older distr
echo Upgrading numpy...
pip install --upgrade numpy

echo Upgrading scipy...
pip install --upgrade scipy

echo Upgrading matplotlib...
pip install --upgrade matplotlib

# now install packages
echo Installing obspy...
pip install  --upgrade obspy

echo Installing scikit-learn...
pip install  --upgrade sklearn

echo Installing PyQt5...
pip install  --upgrade PyQt5

echo Installing configparser...
pip install  --upgrade configparser

 
echo Installation completed!