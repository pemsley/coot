REM pyrogen version for windows

set RDKIT_DIR=C:\boost\RDKit_2014_03_1
set RDKIT_LIBDIR=%RDKIT_DIR%\build\lib
set RDKIT_BINDIR=%RDKIT_DIR%\build\bin
set RDBASE=C:\boost\RDKit_2014_03_1
set PYTHONPATH=%RDBASE%
set PATH=%PATH%;%RDBASE%\lib;C:\boost\boost_1_55_0\stage\lib

set prfx=%~dp0\..

echo prfx is %prfx%

set PYTHONHOME=%prfx%\python27

set PATH=%prfx%\python27;%prfx%\bin;%RDKIT_BINDIR%;%PATH%

python -m pyrogen %*
