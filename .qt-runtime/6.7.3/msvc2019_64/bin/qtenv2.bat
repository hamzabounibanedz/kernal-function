@echo off
echo Setting up environment for Qt usage...
set PATH=.qt-runtime\6.7.3\msvc2019_64\bin;%PATH%
cd /D .qt-runtime\6.7.3\msvc2019_64
echo Remember to call vcvarsall.bat to complete environment setup!
