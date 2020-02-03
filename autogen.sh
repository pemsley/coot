
sys=`uname`

OS=`uname`
hostname=`hostname`
if [ -z "$HOST" ]; then
 HOST=$hostname
fi


# Bernie wants libtoolize
echo libtoolize --copy --no-warn
libtoolize --copy --no-warn

if test $sys = Darwin ; then 
   echo aclocal -I macros -I /sw/share/aclocal
   aclocal -I macros -I /sw/share/aclocal
else
   # because on bubbles we have automake 1.9.6, we need 1.9.6 on kalypso, but system
   # is 1.9.5, so I build my own. But that needs m4 files in /usr/share/aclocal
   # It also needs libtool.m4 from a recent libtool
   #
   # but we can only add the directories if they exist, otherwise aclocal barfs (sigh)

    aclocal_extras=
    # dir_list="$HOME/libtool/share/aclocal $HOME/gettext/share/aclocal $HOME/glade/share/aclocal $HOME/autobuild/Linux/Coot-0.1/share/aclocal $HOME/autotools/share/aclocal-1.9 $HOME/automake/share/aclocal-1.9 $HOME/autobuild/$OS-$hostname/share/aclocal $HOME/autobuild/$OS-$hostname-pre-release/share/aclocal $HOME/autobuild/$OS-$hostname-pre-release-gtk2/share/aclocal /usr/share/aclocal $HOME/test/gtkglext/share/aclocal $HOME/build/share/aclocal $HOME/autobuild/Linux-$HOST/share/aclocal $HOME/gtk-1.2/share/aclocal $HOME/gtk-1/share/aclocal"
    dir_list="$HOME/autobuild/$OS-$HOST-pre-release-gtk2-python/share/aclocal"

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


