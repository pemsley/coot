
sys=`uname`

OS=`uname`
hostname=`hostname`
if [ -z "$HOST" ]; then
 HOST=$hostname
fi

if test -e ltmain.sh ; then
    # this should not be in the repo, I guess
    rm ltmain.sh
fi


# Bernie wants libtoolize
echo libtoolize --copy
libtoolize --copy

if test $sys = Darwin ; then
   echo aclocal -I macros -I /sw/share/aclocal
   aclocal -I macros -I /sw/share/aclocal
else

    # we can only add the directories if they exist, otherwise aclocal barfs (sigh)

    aclocal_extras=
    dir_list="$HOME/autobuild/$OS-$HOST-pre-release-gtk2-python/share/aclocal $HOME/autobuild/$OS-refinement-pre-release-gtk2-python/share/aclocal"


    for dir in $dir_list
    do
       if [ -e $dir ] ; then
	  aclocal_extras="$aclocal_extras -I $dir"
       fi
    done

    case "$sys" in
    MINGW32_NT-* )
	echo We have WIN
	aclocal_extras="-I /usr/local/share/aclocal -I /mingw/share/aclocal";;
    esac
    echo aclocal -I macros $aclocal_extras
    aclocal -I macros $aclocal_extras
fi

which git

echo autoconf
autoconf

which git

echo automake --add-missing --copy
automake --add-missing --copy


# configure needs to be newer that aclocal.m4, which on Paris
# checkout, it seems not to be.  Perhaps they are generated too close
# to distinguish.  So let's touch configure after automake, that
# should provide sufficient time distinction.

touch configure

# then in the build directory, configure.status needs to be newer that
# ../coot/configure or else it does a reconfigure.

# No longer can/should we do this... (all hail git)
# if [ ! -e src/svn-revision.cc ] ; then
#    bash generate-svn-revision-cc.sh
# fi


