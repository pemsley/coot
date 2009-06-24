
svn -qu status | tail -1 | awk '
BEGIN {svn_revision_cc = "src/svn-revision.cc"}
{
# Adding this extern "C" is just magic.  Baah. 
# If it is not here then (in coot_version()):
#    c-interface.cc:174: undefined reference to "svn_revision"
#
print "extern \"C\" {"         > svn_revision_cc 
print "int svn_revision() { "  > svn_revision_cc 
print "   return ", $NF ";"    > svn_revision_cc 
print "}"                      > svn_revision_cc;
print "}"                      > svn_revision_cc;
print ""                       > svn_revision_cc;  
}'
tail -4 src/svn-revision.cc | head -1 | awk '{split($2, arr, ";"); print arr[1]}'

