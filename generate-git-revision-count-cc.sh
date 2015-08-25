# -*-shell-script-mode-*-


rev_no=$(git rev-list --count HEAD)

awk -v rev_no=$rev_no '
BEGIN {git_revision_count_cc = "src/git-revision-count.cc"
print "extern \"C\" {"         > git_revision_count_cc 
print "   int git_revision_count() { " > git_revision_count_cc 
print "      return ", rev_no ";"      > git_revision_count_cc 
print "   }"                           > git_revision_count_cc
print "   int svn_revision() { "       > git_revision_count_cc 
print "      return git_revision_count();" > git_revision_count_cc 
print "   }"                           > git_revision_count_cc
print "}"                              > git_revision_count_cc
print ""                               > git_revision_count_cc  
}'

awk -v rev_no=$rev_no '
BEGIN {o = "pyrogen/coot_git.py"
   print "def revision_count():"    > o
   print "    return", rev_no       > o
}'		

if test ! -e src/git-revision-count.cc ; then
   echo Missing src/git-revision-count.cc
else 
   grep return src/git-revision-count.cc | head -1 | awk '{split($2, arr, ";"); print arr[1]}'
fi


