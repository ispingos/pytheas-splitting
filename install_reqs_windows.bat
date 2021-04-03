@echo off

rem ############################
rem version date: 01/12/2020
rem ############################

rem make sure pip is up to date
rem echo Upgrading pip...
rem pip install --upgrade pip

rem first upgrade older distr
echo Upgrading numpy...
pip install --upgrade numpy==1.19.3

echo Upgrading scipy...
pip install --upgrade scipy

echo Upgrading matplotlib...
pip install --upgrade matplotlib

rem now install packages
echo Installing obspy...
pip install --upgrade  obspy

echo Installing scikit-learn...
pip install --upgrade  sklearn

echo Installing PyQt5...
pip install --upgrade  PyQt5

echo Installing configparser...
pip install --upgrade  configparser

 
echo Installation completed!
pause