
sys=`uname`

OS=`uname`
hostname=`hostname`

if test $sys = Darwin ; then 
   aclocal -I macros -I /sw/share/aclocal
else
   # because on bubbles we have automake 1.9.6, we need 1.9.6 on kalypso, but system
   # is 1.9.5, so I build my own. But that needs m4 files in /usr/share/aclocal
   # It also needs libtool.m4 from a recent libtool
   #
   # but we can only add the directories if they exist, otherwise aclocal barfs (sigh)

   aclocal_extras=
   if [ -e $HOME/libtool/share/aclocal ] ; then
      aclocal_extras="-I $HOME/libtool/share/aclocal"
   fi
   if [ -d $HOME/autobuild/Linux/Coot-0.1/share/aclocal ] ; then
     aclocal_extras="$aclocal_extras -I $HOME/autobuild/Linux/Coot-0.1/share/aclocal"
   fi 
   if [ -d $HOME/autotools/share/aclocal-1.9 ] ; then
     aclocal_extras="$aclocal_extras -I $HOME/autotools/share/aclocal-1.9"
   fi
   if [ -e $HOME/automake/share/aclocal-1.9 ] ; then
      aclocal_extras="$aclocal_extras -I $HOME/automake/share/aclocal-1.9"
   fi
   if [ -e $HOME/autobuild/$OS-$hostname/share/aclocal ] ; then
      aclocal_extras="$aclocal_extras -I $HOME/autobuild/$OS-$hostname/share/aclocal"
   fi
   if [ -e $HOME/autobuild/$OS-$hostname-pre-release/share/aclocal ] ; then
      aclocal_extras="$aclocal_extras -I $HOME/autobuild/$OS-$hostname-pre-release/share/aclocal"
   fi
   if [ -e /usr/share/aclocal ] ; then
      aclocal_extras="$aclocal_extras -I /usr/share/aclocal"
   fi
   if [ -e $HOME/test/gtkglext ] ; then
      aclocal_extras="$aclocal_extras -I $HOME/test/gtkglext/share/aclocal"
   fi
   if [ -e $HOME/build/share/aclocal ] ; then
      aclocal_extras="$aclocal_extras -I $HOME/build/share/aclocal"
   fi
   if [ -e $HOME/autobuild/Linux-$HOST/share/aclocal ] ; then
      aclocal_extras="$aclocal_extras -I $HOME/autobuild/Linux-$HOST/share/aclocal"
   fi
   echo aclocal -I macros $aclocal_extras
   aclocal -I macros $aclocal_extras

fi

autoconf

automake --add-missing --copy


# configure needs to be newer that aclocal.m4, which on Paris
# checkout, it seems not to be.  Perhaps they are generated too close
# to distinguish.  So let's touch configure after automake, that
# should provide sufficient time distinction.

touch configure 

# then in the build directory, configure.status needs to be newer that
# ../coot/configure or else it does a reconfigure.


