
svn -qu status | tail -1 | awk '
BEGIN {svn_revision_cc = "src/svn-revision.cc"}
{
print "int svn_revision() { "  > svn_revision_cc 
print "   return ", $NF ";"    > svn_revision_cc 
print "}"                      > svn_revision_cc;
print ""                       > svn_revision_cc;  
}'
tail -3 src/svn-revision.cc | head -1 | awk '{split($2, arr, ";"); print arr[1]}'

