

# PATH=~/fed-4/bin:$PATH
# export PATH
# cd ~/Projects/coot/burn-up

pyextra=/home/emsley/autobuild/Linux-jackal-pre-release-gtk2-python-f10/bin
pyextra=/home/emsley/autobuild/f10/Linux-jackal-f10-pre-release-gtk2-python/bin
pyextra=$HOME/cycle-guile/bin

PATH=$pyextra:$PATH
export PATH

bash new-burn-point.sh
if Rscript burn-up-graph.r ; then

   # cp new-graph.png ~/public_html/coot/devel/burn-up.png
   # convert -antialias -scale 218x146 ~/public_html/coot/devel/burn-up.png ~/public_html/coot/devel/burn-up-icon.png

   # fink convert
   # convert -antialias -scale 218x146 burn-up.png burn-up-icon.png
   :

fi

if [ $(uname) = Darwin ] ; then
   open burn-up.png
else
   eog burn-up.png
fi
