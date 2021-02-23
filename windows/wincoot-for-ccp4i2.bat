@set LANG=en
title WinCoot


set COOT_PREFIX=yourWinCootdirectory
set COOT_GUILE_PREFIX=%COOT_PREFIX:\=/%

set COOT_HOME=%COOT_PREFIX%
set COOT_BACKUP_DIR=%COOT_PREFIX%\coot-backup

set COOT_SHARE=%COOT_PREFIX%\share

if not exist "%CLIBD_MON%" (
  echo no $CLIBD_MON found trying to setup CCP4
  Call :setup_ccp4
)

if not exist "%CLIBD_MON%" (
  echo no $CLIBD_MON found setting COOT_REFMAC_LIB_DIR
  set COOT_REFMAC_LIB_DIR=%COOT_SHARE%\coot\lib
)

set COOT_SCHEME_DIR=%COOT_SHARE%/coot/scheme
set COOT_STANDARD_RESIDUES=%COOT_SHARE%\coot\standard-residues.pdb
set COOT_PIXMAPS_DIR=%COOT_SHARE%\coot\pixmaps
set COOT_RESOURCES_FILE=%COOT_SHARE%\coot\cootrc
set COOT_DATA_DIR=%COOT_SHARE%\coot
set COOT_REF_STRUCTS=%COOT_SHARE%\coot\reference-structures
set COOT_PYTHON_DIR=%COOT_PREFIX%\lib\python2.7\site-packages\coot
REM set COOT_REF_SEC_STRUCTS=%COOT_SHARE%\coot\ss-reference-structures

REM set PYTHONHOME=%COOT_PREFIX%\python27

set GUILE_LOAD_PATH=%COOT_GUILE_PREFIX%/share/guile/1.8;%COOT_GUILE_PREFIX%/share/guile;%COOT_GUILE_PREFIX%/share/guile/gtk-2.0;%COOT_GUILE_PREFIX%/share/guile/gui;%COOT_GUILE_PREFIX%/share/guile/www;%COOT_GUILE_PREFIX%/share/guile/site

set SYMINFO=%COOT_SHARE%\coot\syminfo.lib

set PATH=%COOT_PREFIX%\bin;%COOT_PREFIX%\lib;%PATH%;%COOT_PREFIX%\bin\extras

coot-bin.exe %*

Exit /B %ERRORLEVEL%

:setup_ccp4
For /F "Skip=1 Tokens=2*" %%A In (
    'Reg Query "HKLM\SOFTWARE\WOW6432Node\CCP4-7"^
    /V "InstallDir" 2^>Nul'
) Do Set "CCP4Dir=%%~B"
REM find latest dir
FOR /F " tokens=*" %%i IN ('dir "%CCP4DIR%\7.*" /b /ad-h /t:c /od') DO SET vers=%%i
REM setup ccp4
call "%CCP4DIR%\%vers%\ccp4.setup.bat"
Exit /B 0
