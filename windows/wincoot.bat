@set LANG=en
@set LC_ALL=en
@set LC_NUMERIC=en
title WinCoot

REM First try to find a dictionary. Search in the registry for a
REM CCP4 installation and then setup CCP4.

if not exist "%CLIBD_MON%" (
   echo no $CLIBD_MON found trying to setup CCP4
   call :setup_ccp4
)

REM main line

set COOT_PREFIX=%~dp0
set COOT_GUILE_PREFIX=%COOT_PREFIX:\=/%

REM New COOT_HOME is in User directory now
set COOT_HOME=%USERPROFILE%\COOT
set COOT_BACKUP_DIR=%COOT_PREFIX%\coot-backup

set COOT_SHARE=%COOT_PREFIX%\share

if not exist "%CLIBD_MON%" (
  echo no $CLIBD_MON found setting COOT_REFMAC_LIB_DIR
  Call :setup_coot_dictionary
)

set COOT_SCHEME_DIR=%COOT_SHARE%/coot/scheme
set COOT_STANDARD_RESIDUES=%COOT_SHARE%\coot\standard-residues.pdb
set COOT_PIXMAPS_DIR=%COOT_SHARE%\coot\pixmaps
set COOT_RESOURCES_FILE=%COOT_SHARE%\coot\cootrc
set COOT_DATA_DIR=%COOT_SHARE%\coot
set COOT_REF_STRUCTS=%COOT_SHARE%\coot\reference-structures
set COOT_PYTHON_DIR=%COOT_PREFIX%\lib\python2.7\site-packages\coot
REM set PYTHONHOME=%COOT_PREFIX%\python27

set GUILE_LOAD_PATH=%COOT_GUILE_PREFIX%/share/guile/1.8;%COOT_GUILE_PREFIX%/share/guile;%COOT_GUILE_PREFIX%/share/guile/gtk-2.0;%COOT_GUILE_PREFIX%/share/guile/gui;%COOT_GUILE_PREFIX%/share/guile/www;%COOT_GUILE_PREFIX%/share/guile/site

set SYMINFO=%COOT_SHARE%\coot\syminfo.lib

set PATH=%COOT_PREFIX%\bin;%COOT_PREFIX%\lib;%PATH%;%COOT_PREFIX%\bin\extras

coot-bin.exe %*

Exit /B %ERRORLEVEL%


REM function to setup CCP4
:setup_ccp4
set /a CCP4_MAJOR=10

:loop_start
set /a CCP4_MAJOR=CCP4_MAJOR-1

call :find_ccp4 HKCU\SOFTWARE\CCP4-%CCP4_MAJOR%
if exist "%CLIBD_MON%" (goto :loop_end)

call :find_ccp4 HKLM\SOFTWARE\WOW6432Node\CCP4-%CCP4_MAJOR%
if exist "%CLIBD_MON%" (goto :loop_end)

if %CCP4_MAJOR% gtr 6 (goto :loop_start)
:loop_end

Exit /B 0


REM Function to find CCP4 in the registry and call ccp4 setup bat
:find_ccp4
REM Check in Registry for CCP4 InstallDir
for /f "skip=1 tokens=3" %%A in (
  'reg query %1 /v InstallDir 2^>Nul'
) do set CCP4DIR=%%~A
if not exist "%CCP4DIR%" (Exit /B 0)

REM find latest dir
for /f " tokens=*" %%i in (
  'dir "%CCP4DIR%\%CCP4_MAJOR%.*" /b /ad-h /t:c /od'
) do set SETUPFILE=%CCP4DIR%\%%i\ccp4.setup.bat
if not exist "%SETUPFILE%" (Exit /B 0)

REM Finally, setup CCP4
call "%SETUPFILE%"
Exit /B 0

:setup_coot_dictionary
set COOT_REFMAC_LIB_DIR=%COOT_SHARE%\coot\lib
set WARN=No valid dictionary found, no refinement possible!
if not exist "%COOT_REFMAC_LIB_DIR%" (
  REM check for Phenix place too
  if not exist "%COOT_MONOMER_LIB_DIR%" (
    msg "%USERNAME%" WinCoot WARNING:: %WARN%
  )
)
Exit /B 0
