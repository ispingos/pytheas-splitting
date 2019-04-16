@echo off

:: first upgrade older distr
echo Upgrading numpy...
pip install --upgrade numpy

echo Upgrading scipy...
pip install --upgrade scipy

echo Upgrading matplotlib...
pip install --upgrade matplotlib

:: now install packages
echo Installing obspy...
pip install obspy

echo Installing scikit-learn...
pip install sklearn

echo Installing PyQt5...
pip install PyQt5

echo Installing configparser...
pip install configparser

 
echo Installation completed!
pause