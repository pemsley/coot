# -*-shell-script-mode-*-

rev_no=$(svn info | awk '/^Revision:/{print $NF}')

if [ -z "$rev_no" ] ; then 
   # should never happen
   echo no svn status, using html
   wget --quiet -O svn-rev.html http://coot.googlecode.com/svn/
   rev_no=$(awk '/Revision/ {if ($1 == "<h2>coot") { sub(":", "", $4); print $4 } }' svn-rev.html)
fi

awk -v rev_no=$rev_no '
BEGIN {svn_revision_cc = "src/svn-revision.cc"
print "extern \"C\" {"         > svn_revision_cc 
print "int svn_revision() { "  > svn_revision_cc 
print "   return ", rev_no ";" > svn_revision_cc 
print "}"                      > svn_revision_cc
print "}"                      > svn_revision_cc
print ""                       > svn_revision_cc  
}'

if test ! -e src/svn-revision.cc ; then
   echo Missing src/svn-revision.cc
else 
   tail -4 src/svn-revision.cc | head -1 | awk '{split($2, arr, ";"); print arr[1]}'
fi

