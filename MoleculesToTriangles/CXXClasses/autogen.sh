
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


# we can only add the directories if they exist, otherwise aclocal barfs (sigh)

aclocal_extras=

# these need to correspond with the build directories
dir_list="$HOME/autobuild/$OS-$HOST-pre-release-gtk3/share/aclocal $HOME/autobuild/$OS-$HOST-gtk3/share/aclocal"

for dir in $dir_list
    do
       if [ -e $dir ] ; then
	  aclocal_extras="$aclocal_extras -I $dir"
       fi
    done

case "$sys" in
    MINGW32_NT-* )
	echo We have WIN
	aclocal_extras="-I /usr/local/share/aclocal -I /mingw/share/aclocal" ;;
esac
echo aclocal $aclocal_extras
aclocal $aclocal_extras

echo autoconf
autoconf

echo automake --add-missing --copy
automake --add-missing --copy

