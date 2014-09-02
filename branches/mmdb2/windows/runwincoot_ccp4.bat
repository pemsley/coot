@set LANG=en
title WinCoot


set COOT_PREFIX=%COOT_MASTER%
set COOT_GUILE_PREFIX=%COOT_PREFIX%

set COOT_HOME=%COOT_PREFIX%
set COOT_BACKUP_DIR=%COOT_PREFIX%\coot-backup

set COOT_SHARE=%COOT_PREFIX%\share
if not exist "%CLIBD_MON%" (
  echo no $CLIBD_MON found setting COOT_REFMAC_LIB_DIR
  set COOT_REFMAC_LIB_DIR=%COOT_SHARE%\coot\lib
)
set COOT_SCHEME_DIR=%COOT_GUILE_PREFIX%/share/coot/scheme
set COOT_STANDARD_RESIDUES=%COOT_SHARE%\coot\standard-residues.pdb
set COOT_PIXMAPS_DIR=%COOT_SHARE%\coot\pixmaps
set COOT_RESOURCES_FILE=%COOT_SHARE%\coot\cootrc
set COOT_DATA_DIR=%COOT_SHARE%\coot
set COOT_REF_STRUCTS=%COOT_SHARE%\coot\reference-structures
set COOT_PYTHON_DIR=%COOT_SHARE%\coot\python
set COOT_REF_SEC_STRUCTS=%COOT_SHARE%\coot\ss-reference-structures

if not "%PYTHONPATH%" == "" goto APPENDPATH
        set PYTHONPATH=%COOT_SHARE%\coot\python
        goto CONT
:APPENDPATH
        set PYTHONPATH=%COOT_SHARE%\coot\python;%PYTHONPATH%
:CONT
set PYTHONHOME=%COOT_PREFIX%\bin

set GUILE_LOAD_PATH=%COOT_GUILE_PREFIX%/share/guile/1.8;%COOT_GUILE_PREFIX%/share/guile;%COOT_GUILE_PREFIX%/share/guile/gtk-2.0;%COOT_GUILE_PREFIX%/share/guile/gui;%COOT_GUILE_PREFIX%/share/guile/www;%COOT_GUILE_PREFIX%/share/guile/site

set SYMINFO=%COOT_SHARE%\coot\syminfo.lib

set PATH=%COOT_PREFIX%\bin;%COOT_PREFIX%\lib;%PATH%

IF NOT EXIST %COOT_PREFIX%\etc\gtk-2.0\gdk-pixbuf.loaders (gdk-pixbuf-query-loaders.exe > %COOT_PREFIX%\etc\gtk-2.0\gdk-pixbuf.loaders)

coot-real.exe %*

