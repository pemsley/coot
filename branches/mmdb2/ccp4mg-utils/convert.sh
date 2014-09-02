
# for file in *.cc *.h
for file in mgtree.h mgtree.cc 
do
    # sed -e 's/[^_]*cartesian\.h/smn_cartesian.h/' $file > $file.tmp
    awk -f cartesian.awk $file > $file.tmp
    if test $? ; then
       mv $file.tmp $file
    fi
done

# currently we change mgtree.h cartesian.h to smn_cartesian.h by hand.

file=mgtree.h
sed -e 's/cartesian\.h/smn_cartesian.h/' $file > $file.tmp
if test $? ; then
   mv $file.tmp $file
fi
