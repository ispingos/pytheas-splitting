@echo off
if not "%1" == "max" start /MAX cmd /c %0 max & exit/b
python pytheas\pytheas.py
pause