
if [ $(uname -a) = Darwin ] ; then 
   libtoolize=glibtoolize
else 
   libtoolize=libtoolize
fi

echo aclocal -I m4 
aclocal -I m4 

echo $libtoolize --copy --no-warn
$libtoolize --copy --no-warn

echo autoconf
autoconf

echo automake --add-missing --copy
automake --add-missing --copy
