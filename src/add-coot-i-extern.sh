pre=coot_pre.i
post=coot.i
tmp1=tmp1.txt
tmp2=tmp2.txt

# first put all PyObject definitions in $tmp1
grep PyObject c-interface.h          > $tmp1
grep PyObject cc-interface.hh       >> $tmp1
#grep PyObject c-interface-mmdb.hh   >> $tmp1
grep PyObject c-interface-python.hh >> $tmp1

# now add the extern
sed -e 's/^/extern /' $tmp1 > $tmp2

# and put it in coot.i (together with USE_PYTHON)
cat $pre                    > $post
echo ""                     >> $post
echo "#ifdef USE_PYTHON"    >> $post
cat $tmp2                   >> $post
echo "#endif // USE_PYTHON" >> $post

rm $tmp1
rm $tmp2
