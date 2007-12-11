
revision=`svn -qu status | tail -1 | awk '
BEGIN {svn_revision_cc = "src/svn-revision.cc"}
{
print "int svn_revision() { "  > svn_revision_cc 
print "   return ", $NF ";"    > svn_revision_cc 
print "}"                      > svn_revision_cc;
print ""                       > svn_revision_cc;  
}'`

