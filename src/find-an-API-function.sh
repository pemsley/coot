
if [ -z "$1" ] ; then
   echo missing search string
   exit 2
fi

grep $1 *.h *.hh | grep -v ^mol | grep -v ^restraints-editor | grep -v widget | grep -v Widget | grep -v ^callbacks | grep -v ^graphics-info | grep -v rama_mouse | grep -v gtk-manual.h

