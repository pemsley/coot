REM this just runs the wincoot-autobuild-all bash script
REM we delay it by 3 min after start-up to make sure everything is running
sleep 180
REM C:\msys\bin\sh --login C:\msys\home\bernhard\make_wincoot_installer.bat
C:\msys\bin\sh --login C:\msys\home\bernhard\Projects\coot\windows\wincoot-autobuild-all
REM shutdown -i -s
shutdown -f -s
exit
